import argparse
import random
import os
import shutil
import networkx
import sys
import subprocess

from toil.job import Job
from toil.common import Toil

from sonLib.bioio import fastaRead, fastaWrite, catFiles, reverseComplement

def makeURL(path):
    return "file://%s" % path

import CactusTEAnnotator.treeBuilding as treeBuilding

dockerImage = "cactus-te-annotator:latest"

def runInContainer(parameters, streamfile=None):
    cmd = ["docker", "run", "-it", "--rm", "-v", "%s:/data" % os.getcwd(), dockerImage] + parameters
    if streamfile:
        subprocess.check_call(cmd, stdout=streamfile)
    else:
        output = subprocess.check_output(cmd)
        return output

def catFiles(fileList, target):
    with open(target, 'w') as targetWrite:
        for f in fileList:
            with open(f, 'r') as read:
                lines = read.readlines()
                targetWrite.write("".join(lines))
                targetWrite.write("\n")

class GtfInfo:
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


def getSequence(job, hal, genome, chrom, strand, start, end):
    hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]
    fastaLines = runInContainer(parameters=hal2fastaCmd).strip().split("\n")
    sequence = "\n".join(fastaLines[1:])
    if strand == "-":
        return reverseComplement(sequence)
    return sequence

def getFasta(job, hal, genome, chrom, start, end):
    hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome]
    if end and start and chrom:
        hal2fastaCmd.extend(["--sequence", chrom, "--start", str(start), "--length", str(end - start)])
    fastaFile = job.fileStore.getLocalTempFile()
    with open(fastaFile, 'w') as fh:
        runInContainer(parameters=hal2fastaCmd, streamfile=fh)
    return fastaFile

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

def getInsertionsOnBranch(job, halID, genome, args):
    """Extract candidate TE insertions from the Cactus alignment.
        """
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.getLocalTempFile()
    i = 0
    with open(gff, 'w') as gffWrite:
        for bedLine in runInContainer(parameters=["halAlignedExtract", "--complement", os.path.basename(hal), genome]).split("\n"):
            if len(bedLine.split()) != 3:
                continue
            chrom, start, end = bedLine.split()
            start = int(start)
            end = int(end)
            if end - start > args.maxInsertionSize:
                continue
            if end - start < args.minInsertionSize:
                continue
            sequence = getSequence(job=job, hal=hal, genome=genome, chrom=chrom, strand="+", start=start, end=end)
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

def runRepeatScout(job, genome, halID, gffID, seqID):
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.readGlobalFile(gffID)
    seq = job.fileStore.readGlobalFile(seqID)

    seqs = job.fileStore.getLocalTempFile()
    with open(gff, 'r') as gffFile:
        with open(seqs, 'w') as seqsFile:
            for line in gffFile:
                repeatCandidate = GtfInfo(line)
                sequence = getSequence(job=job, hal=hal, genome=genome, chrom=repeatCandidate.chrom, strand=repeatCandidate.strand, start=repeatCandidate.start, end=repeatCandidate.end)
                seqsFile.write(">%s\n" % repeatCandidate.name)
                seqsFile.write(sequence)

    repeatScoutFreqs = job.fileStore.getLocalTempFile()
    runInContainer(parameters=["build_lmer_table", "-sequence", os.path.basename(seqs), "-freq", os.path.basename(repeatScoutFreqs)])

    repeatScoutLibrary = job.fileStore.getLocalTempFile()
    runInContainer(parameters=["RepeatScout", "-sequence", os.path.basename(seqs), "-output", os.path.basename(repeatScoutLibrary), "-freq", os.path.basename(repeatScoutFreqs)])

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
    runInContainer(parameters=["RepeatMasker", "-nolow", "-cutoff", "650", "-dir", os.path.basename(repeatMaskerOutput), "-lib", os.path.basename(repeatLibrary), os.path.basename(seq)])
    outputGff = "%s/%s.out" % (repeatMaskerOutput, os.path.basename(seq))
    job.fileStore.logToMaster("directory contents: %s" % os.listdir(repeatMaskerOutput))
    outputGffID = job.fileStore.writeGlobalFile(outputGff)

    return job.addChildJobFn(parseRepeatMaskerGFF, gffID = outputGffID, args = args).rv()

def minhashClustering(job, repeatCandidates, args):
    graph = networkx.DiGraph()
    for repeatCandidate in repeatCandidates:
        graph.add_node(repeatCandidate.name)
    seqs = job.fileStore.getLocalTempFile()
    with open(seqs, "w") as seqsWrite:
        for seq in repeatCandidates:
            seqsWrite.write(">%s\n" % seq.name)
            seqsWrite.write("%s\n" % seq.seq)

    families = []
    nameToGtfInfo = {element.name:element for element in elements}
    for line in runInContainer(parameters=["minhash", "--kmerLength", str(args.kmerLength), "--distanceThreshold", str(args.minhashDistanceThreshold), "--sequences", os.path.basename(seqs)]).split("\n"):
        if line == "":
            continue
        family = line.split()
        family = [nameToGtfInfo[name] for name in family]
        families.append(family)
    return families


def writeSequencesToFile(gff, seqsFile):
    with open(seqsFile, "w") as seqsWrite:
        for repeatCandidate in repeatCandidates:
            seqsWrite.write(">%s\n" % repeatCandidate.name)
            seqsWrite.write("%s\n" % repeatCandidate.seq)

def splitFamilies(job, halID, genome, gffID, heaviestBundle, familyNumber, args):
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.readGlobalFile(gffID)
    families = {}
    with open(gff, "r") as gffRead:
        for line in gffRead:
            lineInfo = GtfInfo(line)
            families[lineInfo.family] = lineInfo
    for family in families:
        fasta = job.fileStore.getLocalTempFile()
        with open(fasta, "w") as fastaFile:
            for i in families[family]:
                fastaFile.write(">%s\n" % i.name)
                seq = getSequence(job, hal=hal, genome=genome, chrom=i.chrom, start=i.start, end=i.end)
                fastaFile.write("%s\n" % seq)
        fastaID = job.fileStore.writeGlobalFile(fasta)
        graphID = job.addChildJobFn(runPoa, fastaID=fastaID, args=args).rv()

def runPoa(job, fastaID, args):
    fasta = job.readGlobalFile(fastaID)
    graph = job.fileStore.getLocalTempFile()
    cmd = ["poa", "-read_fasta", os.path.basename(fasta), "-po", os.path.basename(graph), args.substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    runInContainer(parameters=cmd)
    return job.fileStore.writeGlobalFile(graph)

def clusterByHeaviestBundling(repeatCandidates, familyNumber, args):
    graph = runPoa(repeatCandidates, True, familyNumber, args)

    nameToGtfInfo = {repeatCandidate.name:repeatCandidate for repeatCandidate in repeatCopies}
    subfamilies = []
    for line in runInContainer(parameters=["getHeaviestBundles", graph]).split("\n"):
        if line == "":
            continue
        seqsInBundle = [nameToGtfInfo[name] for name in line.split() if name in nameToGtfInfo]
        subfamilies.append(seqsInBundle)

    return subfamilies

def clusterByAlignmentDistances(repeatCandidates, familyNumber, args):
    graph = runPoa(repeatCopies, False, familyNumber, args)
    subfamilies = []
    nameToGtfInfo = {repeatCopy.name:repeatCandidate for repeatCandidate in repeatCandidates}
    for line in runInContainer(parameters=["clusterByAlignmentDistances", os.path.basename(graph), str(args.alignmentDistanceThreshold)]).split("\n"):
        line = line.strip().rstrip()
        if line == "":
            continue
        subfamilyNames = line.split()
        subfamilies.append([nameToGtfInfo[repeatCandidateName] for repeatCandidateName in subfamilyNames])

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
            dist = float(runInContainer(parameters=["minhash", "--kmerLength", str(args.minhashClusterDistKmerLength), "--clusterA", seqFileNames[i], "--clusterB", seqFileNames[j]], stdout=True))
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

    genomeFile = getFasta(job=job, hal=hal, genome=args.genome, chrom=args.chrom, start=args.start, end=args.end)
    seqID = job.fileStore.writeGlobalFile(genomeFile)
    getInsertionsJob = Job.wrapJobFn(getInsertionsOnBranch, halID=halID, genome=args.genome, args=args)

    repeatScoutJob = Job.wrapJobFn(runRepeatScout, genome=args.genome, halID=halID, gffID = getInsertionsJob.rv(), seqID=seqID)
    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=repeatScoutJob.rv(), seqID=seqID, args=args)

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
