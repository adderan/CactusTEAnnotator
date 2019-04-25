import subprocess
import argparse
import random
import os
import shutil
import networkx
import sys

from toil.job import Job
from toil.common import Toil

from sonLib.bioio import fastaRead, fastaWrite, catFiles, reverseComplement

def makeURL(path):
    return "file://%s" % path

import CactusTEAnnotator.treeBuilding as treeBuilding

def catFiles(fileList, target):
    with open(target, 'w') as targetWrite:
        for f in fileList:
            with open(f, 'r') as read:
                lines = read.readlines()
                targetWrite.write("".join(lines))
                targetWrite.write("\n")

class RepeatCandidate:
    def __init__(self, gffLine):
        chrom, annotationType, name, start, end, score, strand, a, group = gffLine.split()
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand
        self.group = group

def getGffLine(chrom, name, strand, start, end, family):
    return "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (chrom, "cactus_repeat_annotation", name, start, end, 0, strand, '.', family)


def getSequence(hal, chrom, strand, start, end, args):
    hal2fastaCmd = ["hal2fasta", hal, args.genome, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]
    fastaLines = subprocess.check_output(hal2fastaCmd).strip().split("\n")
    sequence = "\n".join(fastaLines[1:])
    if strand == "-":
        return reverseComplement(sequence)
    return sequence

def getFasta(hal, fastaFile, chrom, start, end, args):
    hal2fastaCmd = ["hal2fasta", hal, args.genome]
    if end and start and chrom:
        hal2fastaCmd.extend(["--sequence", chrom, "--start", str(start), "--length", str(end - start)])
    with open(fastaFile, 'w') as fh:
        subprocess.check_call(hal2fastaCmd, stdout = fh)

def getRootPath():
    import CactusTEAnnotator
    i = os.path.abspath(CactusTEAnnotator.__file__)
    return os.path.split(i)[0]

def highNFraction(seqs):
    for i in range(100):
        seq = random.choice(seqs)
        if seq == '':
            continue
        base = random.choice(seq)
        if base == 'N' or base == 'n':
            return True
    return False

def getInsertions(job, halID, args):
    """Extract candidate TE insertions from the Cactus alignment.
        """
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.getLocalTempFile()
    i = 0
    with open(gff, 'w') as gffWrite:
        for bedLine in subprocess.check_output(["halAlignedExtract", "--complement", hal, args.genome]).split("\n"):
            if len(bedLine.split()) != 3:
                continue
            chrom, start, end = bedLine.split()
            start = int(start)
            end = int(end)
            if end - start > args.maxInsertionSize:
                continue
            if end - start < args.minInsertionSize:
                continue
            sequence = getSequence(hal = hal, chrom = chrom, strand = "+", start = start, end = end, args = args)
            assert len(sequence) > 0

            if highNFraction(sequence):
                continue
            if args.maxInsertions > 0 and i >= args.maxInsertions:
                break
            gffWrite.write(getGffLine(chrom=chrom, start=start, end=end, name=str(i), strand="+", family=str(i)))
            gffWrite.write(getGffLine(chrom=chrom, start=start, end=end, name=str(i), strand="-", family=str(i)))
            i = i + 1

    print(sys.stderr, "Found %d insertions on branch" % i)

    return job.fileStore.writeGlobalFile(gff)

def runRepeatScout(job, halID, gffID, seqID, args):
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.readGlobalFile(gffID)
    seq = job.fileStore.readGlobalFile(seqID)

    seqs = job.fileStore.getLocalTempFile()
    with open(gff, 'r') as gffFile:
        with open(seqs, 'w') as seqsFile:
            for line in gffFile:
                repeatCandidate = RepeatCandidate(line)
                sequence = getSequence(hal = hal, chrom = repeatCandidate.chrom, strand = repeatCandidate.strand, start = repeatCandidate.start, end = repeatCandidate.end, args = args)
                seqsFile.write(">%s\n" % repeatCandidate.name)
                seqsFile.write(sequence)

    repeatScoutFreqs = job.fileStore.getLocalTempFile()
    subprocess.check_call(["build_lmer_table", "-sequence", seqs, "-freq", repeatScoutFreqs])

    repeatScoutLibrary = job.fileStore.getLocalTempFile()
    subprocess.check_call(["RepeatScout", "-sequence", seqs, "-output", repeatScoutLibrary, "-freq", repeatScoutFreqs])

    return job.fileStore.writeGlobalFile(repeatScoutLibrary)

def parseRepeatMaskerGFF(job, gffID, args):
    rmaskGffFile = job.fileStore.readGlobalFile(gffID)
    gffFile = job.fileStore.getLocalTempFile()
    with open(rmaskGffFile, 'r') as outputGFFRead:
        with open(gffFile, 'w') as gffWrite:
            for i, line in enumerate(outputGFFRead):
                try:
                    score, div, deletions, insertions, chrom, start, end, left, strand, family, repeatType, a, b, c, d = line.split()
                    start = int(start) + args.start
                    end = int(end) + args.start
                    gffLine = getGffLine(chrom = chrom, start = int(start), end = int(end), name = str(i), strand = strand, family = family)
                    gffWrite.write("%s" % gffLine)
                except ValueError:
                    continue
    return job.fileStore.writeGlobalFile(gffFile)


def runRepeatMasker(job, repeatLibraryID, seqID, args):
    repeatLibrary = job.fileStore.readGlobalFile(repeatLibraryID)
    seq = job.fileStore.readGlobalFile(seqID)
    repeatMaskerOutput = job.fileStore.getLocalTempDir()
    subprocess.check_call(["RepeatMasker", "-nolow", "-cutoff", "650", "-dir", repeatMaskerOutput, "-lib", repeatLibrary, seq])
    outputGff = "%s/%s.out" % (repeatMaskerOutput, os.path.basename(seq))
    job.fileStore.logToMaster("directory contents: %s" % os.listdir(repeatMaskerOutput))
    outputGffID = job.fileStore.writeGlobalFile(outputGff)

    return job.addChildJobFn(parseRepeatMaskerGFF, gffID = outputGffID, args = args).rv()

def minhashClustering(job, repeatCandidates, args):
    seqs = job.fileStore.getLocalTempFile()
    with open(seqs, "w") as seqsWrite:
        for seq in repeatCandidates:
            seqsWrite.write(">%s\n" % seq.name)
            seqsWrite.write("%s\n" % seq.seq)

    families = []
    nameToRepeatCandidate = {element.name:element for element in elements}
    for line in subprocess.check_output(["minhash", "--kmerLength", str(args.kmerLength), "--distanceThreshold", str(args.minhashDistanceThreshold), "--sequences", seqs]).split("\n"):
        if line == "":
            continue
        family = line.split()
        family = [nameToRepeatCandidate[name] for name in family]
        families.append(family)
    return families

####Partial order alignment###########

def writeSequencesToFile(gff, seqsFile):
    with open(seqsFile, "w") as seqsWrite:
        for repeatCandidate in repeatCandidates:
            seqsWrite.write(">%s\n" % repeatCandidate.name)
            seqsWrite.write("%s\n" % repeatCandidate.seq)

def runPoa(repeatCopies, heaviestBundle, familyNumber, args):
    seqs = os.path.join(args.outDir, "%d.fa" % int(familyNumber))
    with open(seqs, "w") as seqsWrite:
        for repeatCopy in repeatCopies:
            seqsWrite.write(">%s\n" % repeatCopy.name)
            seqsWrite.write("%s\n" % repeatCopy.seq)

    graph = getTempFile(rootDir = args.outDir)
    cmd = ["poa", "-read_fasta", seqs, "-po", graph, args.substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    subprocess.check_call(cmd)

    return graph

def clusterByHeaviestBundling(repeatCandidates, familyNumber, args):
    graph = runPoa(repeatCandidates, True, familyNumber, args)

    nameToRepeatCandidate = {repeatCandidate.name:repeatCandidate for repeatCandidate in repeatCopies}
    subfamilies = []
    for line in subprocess.check_output(["getHeaviestBundles", graph]).split("\n"):
        if line == "":
            continue
        seqsInBundle = [nameToRepeatCandidate[name] for name in line.split() if name in nameToRepeatCandidate]
        subfamilies.append(seqsInBundle)

    return subfamilies

def clusterByAlignmentDistances(repeatCandidates, familyNumber, args):
    graph = runPoa(repeatCopies, False, familyNumber, args)
    subfamilies = []
    nameToRepeatCandidate = {repeatCopy.name:repeatCandidate for repeatCandidate in repeatCandidates}
    for line in subprocess.check_output(["clusterByAlignmentDistances", graph, str(args.alignmentDistanceThreshold)]).split("\n"):
        line = line.strip().rstrip()
        if line == "":
            continue
        subfamilyNames = line.split()
        subfamilies.append([nameToRepeatCandidate[repeatCandidateName] for repeatCandidateName in subfamilyNames])

    return subfamilies

def joinSubfamiliesByMinhashDistance(repeatFamilies, args):
    graph = networkx.Graph()
    
    seqFileNames = []
    for i in range(len(repeatFamilies)):
        graph.add_node(i)
        seqFile = getTempFile(rootDir=args.outDir)
        writeSequencesToFile(repeatFamilies[i], seqFile)
        seqFileNames.append(seqFile)

        for j in range(i):
            if len(repeatFamilies[i]) <= 1 or len(repeatFamilies[j]) <= 1:
                continue
            dist = float(subprocess.check_output(["minhash", "--kmerLength", str(args.minhashClusterDistKmerLength), "--clusterA", seqFileNames[i], "--clusterB", seqFileNames[j]]))
            print("Found distance %f between clusters %d and %d of size %d and %d" % (dist, i, j, len(repeatFamilies[i]), len(repeatFamilies[j])))
            if dist < args.minhashClusterDistanceThreshold:
                assert graph.has_node(i)
                assert graph.has_node(j)
                graph.add_edge(i, j)

    newClusters = []
    for component in networkx.connected_components(graph):
        newCluster = []
        for node in component:
            newCluster.extend(repeatFamilies[node])
        newClusters.append(newCluster)
    return newClusters

def repeatScoutRepeatMasker(job, halID, args):
    hal = job.fileStore.readGlobalFile(halID)
    genomeFile = job.fileStore.getLocalTempFile()
    getFasta(hal = hal, fastaFile = genomeFile, chrom = args.chrom, start = args.start, end = args.end, args = args)
    seqID = job.fileStore.writeGlobalFile(genomeFile)
    getInsertionsJob = Job.wrapJobFn(getInsertions, halID = halID, args = args)

    repeatScoutJob = Job.wrapJobFn(runRepeatScout, halID = halID, gffID = getInsertionsJob.rv(), seqID = seqID, args = args)
    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID = repeatScoutJob.rv(), seqID = seqID, args = args)

    job.addChild(getInsertionsJob)
    getInsertionsJob.addFollowOn(repeatScoutJob)
    repeatScoutJob.addFollowOn(repeatMaskerJob)

    return {'finalGFF':repeatMaskerJob.rv(), 'library':repeatScoutJob.rv()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("genome", type=str)
    parser.add_argument("outGff", type=str)


    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)
    parser.add_argument("--maxInsertions", type=int, default=None)

    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)

    parser.add_argument("--repeatLibrary", type=str, default=None)

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    with Toil(args) as toil:
        halID = toil.importFile(makeURL(args.hal))

        repeatScoutRepeatMaskerJob = Job.wrapJobFn(repeatScoutRepeatMasker, halID = halID, args = args)
        results = toil.start(repeatScoutRepeatMaskerJob)
        toil.exportFile(results['finalGFF'], makeURL(args.outGff))
        if args.repeatLibrary:
            toil.exportFile(results['library'], makeURL(args.repeatLibrary))

if __name__ == "__main__":
    main()
