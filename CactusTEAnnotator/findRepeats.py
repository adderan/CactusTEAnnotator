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

def runCmd(parameters, args, streamfile=None):
    if args.localBinaries:
        cmd = parameters
    else:
        cmd = ["docker", "run", "-it", "--rm", "-v", "%s:/data" % os.getcwd(), dockerImage] + parameters
       
    if streamfile:
        subprocess.check_call(cmd, stdout=streamfile)
    else:
        output = subprocess.check_output(cmd)
        return output

#def catFiles(fileList, target):
#    with open(target, 'w') as targetWrite:
#        for f in fileList:
#            with open(f, 'r') as read:
#                lines = read.readlines()
#                targetWrite.write("".join(lines))
#                targetWrite.write("\n")

def catFilesJobFn(job, fileIDs):
    fileList = [job.fileStore.readGlobalFile(fileID) for fileID in fileIDs]
    combinedFile = job.fileStore.getLocalTempFile()
    catFiles(fileList, combinedFile)
    return job.fileStore.writeGlobalFile(combinedFile)

class GffInfo:
    def __init__(self, gffLine):
        chrom, annotationType, name, start, end, score, strand, a, family = gffLine.split()
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand
        self.family = family

    def printGff(self):
        return getGffLine(chrom=self.chrom, name=self.name, strand=self.strand, start=self.start, end=self.end, family=self.family)

def getGffLine(chrom, name, strand, start, end, family):
    return "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (chrom, "cactus_repeat_annotation", name, start, end, 0, strand, '.', family)


def getSequence(job, hal, genome, chrom, strand, start, end, args):
    hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]
    fastaLines = runCmd(parameters=hal2fastaCmd, args=args).strip().split("\n")
    sequence = "\n".join(fastaLines[1:])
    if strand == "-":
        sequence = reverseComplement(sequence)
    return sequence

def getFastaEntry(job, hal, genome, name, chrom, strand, start, end, args):
    sequence = getSequence(job, hal, genome, chrom, strand, start, end, args)
    fastaEntry = "%s\n%s\n" % (name, sequence)
    return fastaEntry

def getFasta(job, hal, genome, chrom, start, end, args, fastaFile=None):
    hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome]
    if end and start and chrom:
        hal2fastaCmd.extend(["--sequence", chrom, "--start", str(start), "--length", str(end - start)])
    if not fastaFile:
        fastaFile = job.fileStore.getLocalTempFile()
    with open(fastaFile, 'a') as fh:
        runCmd(parameters=hal2fastaCmd, streamfile=fh, args=args)
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
        for bedLine in runCmd(parameters=["halAlignedExtract", "--complement", os.path.basename(hal), genome], args=args).split("\n"):
            if len(bedLine.split()) != 3:
                continue
            chrom, start, end = bedLine.split()
            start = int(start)
            end = int(end)
            if end - start > args.maxInsertionSize:
                continue
            if end - start < args.minInsertionSize:
                continue
            sequence = getSequence(job=job, hal=hal, genome=genome, chrom=chrom, strand="+", start=start, end=end, args=args)
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
                repeatCandidate = GffInfo(line)
                fastaEntry = getFastaEntry(job=job, hal=hal, genome=genome, name=repeatCandidate.name, chrom=repeatCandidate.chrom, strand=repeatCandidate.strand, start=repeatCandidate.start, end=repeatCandidate.end, args=args)
                seqsFile.write(fastaEntry)

    repeatScoutFreqs = job.fileStore.getLocalTempFile()
    runCmd(parameters=["build_lmer_table", "-sequence", os.path.basename(seqs), "-freq", os.path.basename(repeatScoutFreqs)], args=args)

    repeatScoutLibrary = job.fileStore.getLocalTempFile()
    runCmd(parameters=["RepeatScout", "-sequence", os.path.basename(seqs), "-output", os.path.basename(repeatScoutLibrary), "-freq", os.path.basename(repeatScoutFreqs)], args=args)

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
    runCmd(parameters=["RepeatMasker", "-nolow", "-cutoff", "650", "-dir", os.path.basename(repeatMaskerOutput), "-lib", os.path.basename(repeatLibrary), os.path.basename(seq)], args=args)
    outputGff = "%s/%s.out" % (repeatMaskerOutput, os.path.basename(seq))
    job.fileStore.logToMaster("directory contents: %s" % os.listdir(repeatMaskerOutput))
    outputGffID = job.fileStore.writeGlobalFile(outputGff)

    return job.addChildJobFn(parseRepeatMaskerGFF, gffID = outputGffID, args = args).rv()

def minhashClustering(job, gffID, halID, genome, args):
    """Cluster a set of candidate repeat annotations
    by their pairwise Jaccard distances, estimated
    with minhash.
    """
    gff = job.fileStore.readGlobalFile(gffID)
    hal = job.fileStore.readGlobalFile(halID)
    fasta = job.fileStore.getLocalTempFile()

    repeatCandidates = {}
    with open(gff, "r") as gffRead:
        with open(fasta, "w") as fastaWrite:
            for line in gffRead:
                gffInfo = GffInfo(line)
                repeatCandidates[gffInfo.name] = gffInfo
                fastaEntry = getFastaEntry(job=job, hal=hal, genome=genome, name=gffInfo.name, chrom=gffInfo.chrom, strand=gffInfo.strand, start=gffInfo.start, end=gffInfo.end, args=args)
                fastaWrite.write(fastaEntry)
            

    graph = networkx.Graph()
    for i in repeatCandidates.keys():
        graph.add_node(i)
    for line in runCmd(parameters=["minhash", "--kmerLength", str(args.kmerLength), "--sequences", os.path.basename(fasta), "--distancesOnly"], args=args).split("\n"):
        if line == "":
            continue
        i, j, distance = line.split()
        assert graph.has_node(i)
        assert graph.has_node(j)
        if float(distance) < args.distanceThreshold:
            graph.add_edge(i, j)

    newGff = job.fileStore.getLocalTempFile()
    with open(newGff, "w") as gffWrite:
        for family_number, component in enumerate(networkx.connected_components(graph)):
            for i in component:
                repeatCandidates[i].family = family_number
                gffLine = repeatCandidates[i].printGff()
                gffWrite.write(gffLine)
    return job.fileStore.writeGlobalFile(newGff)


def writeSequencesToFile(gff, seqsFile):
    with open(seqsFile, "w") as seqsWrite:
        for repeatCandidate in repeatCandidates:
            seqsWrite.write(">%s\n" % repeatCandidate.name)
            seqsWrite.write("%s\n" % repeatCandidate.seq)

def buildLibrary_poa(job, halID, genome, gffID, args):
    """Build a poa graph for each family and extract
    repeat elements from each graph. Return a library of
    repeat elements in fasta format.
    """
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.readGlobalFile(gffID)
    repeatLibraryIDs = []
    families = {}
    with open(gff, "r") as gffRead:
        for line in gffRead:
            gffInfo = GffInfo(line)
            if gffInfo.family not in families:
                families[gffInfo.family] = []
            families[gffInfo.family].append(gffInfo)
    for family in families:
        fasta = job.fileStore.getLocalTempFile()
        if len(families[family]) < 2:
            continue
        with open(fasta, "w") as fastaFile:
            for i in families[family]:
                fastaEntry = getFastaEntry(job, hal=hal, genome=genome, name=i.name, chrom=i.chrom, strand=i.strand, start=i.start, end=i.end, args=args)
                fastaFile.write(fastaEntry)

        fastaID = job.fileStore.writeGlobalFile(fasta)
        poaJob = Job.wrapJobFn(runPoa, fastaID=fastaID, args=args)
        elementsJob = Job.wrapJobFn(getRepeatElementsFromGraph, graphID=poaJob.rv(), args=args)
        job.addChild(poaJob)
        poaJob.addFollowOn(elementsJob)
        repeatLibraryIDs.append(elementsJob.rv())

    return job.addFollowOnJobFn(catFilesJobFn, repeatLibraryIDs).rv()

def runPoa(job, fastaID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)
    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    cmd = ["poa", "-read_fasta", os.path.basename(fasta), "-po", os.path.basename(graph), substMatrix]
    if args.heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    runCmd(parameters=cmd, args=args)
    return job.fileStore.writeGlobalFile(graph)

def getRepeatElementsFromGraph(job, graphID, args):
    graph = job.fileStore.readGlobalFile(graphID)
    repeatLibrary = job.fileStore.getLocalTempFile()
    with open(repeatLibrary, "w") as fh:
        runCmd(parameters=["getHeaviestBundles", "--lpo", os.path.basename(graph), "--printElements"], streamfile=fh, args=args)

    return job.fileStore.writeGlobalFile(repeatLibrary)

def clusterByHeaviestBundling(repeatCandidates, args):
    graph = runPoa(repeatCandidates, args)

    nameToGffInfo = {repeatCandidate.name:repeatCandidate for repeatCandidate in repeatCopies}
    subfamilies = []
    for line in runCmd(parameters=["getHeaviestBundles", graph], args=args).split("\n"):
        if line == "":
            continue
        seqsInBundle = [nameToGffInfo[name] for name in line.split() if name in nameToGffInfo]
        subfamilies.append(seqsInBundle)

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
            dist = float(runCmd(parameters=["minhash", "--kmerLength", str(args.minhashClusterDistKmerLength), "--clusterA", seqFileNames[i], "--clusterB", seqFileNames[j]], args=args))
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

def poaPipeline(job, halID, genome, args):
    #Get candidate insertions from the cactus alignment
    getInsertionsJob = Job.wrapJobFn(getInsertionsOnBranch, halID=halID, genome=args.genome, args=args)

    #Cluster the sequences into large families with minhash
    #so that the graph construction step is computationally 
    #possible
    initialClusteringJob = Job.wrapJobFn(minhashClustering, gffID=getInsertionsJob.rv(), halID=halID, genome=genome, args=args)

    #Build a poa graph for each cluster and extract 
    #consensus repeat sequences from each graph
    buildLibraryJob = Job.wrapJobFn(buildLibrary_poa, halID=halID, genome=args.genome, gffID = initialClusteringJob.rv(), args=args)

    gff = job.fileStore.getLocalTempFile()
    gffID = job.fileStore.writeGlobalFile(gff)

    job.addChild(getInsertionsJob)
    getInsertionsJob.addFollowOn(initialClusteringJob)
    initialClusteringJob.addFollowOn(buildLibraryJob)

    return {'final.gff': gffID, 'library.fa': buildLibraryJob.rv(), 'clusters.gff': initialClusteringJob.rv()}


def repeatScoutPipeline(job, halID, args):
    hal = job.fileStore.readGlobalFile(halID)

    genomeFile = getFasta(job=job, hal=hal, genome=args.genome, chrom=args.chrom, start=args.start, end=args.end)
    seqID = job.fileStore.writeGlobalFile(genomeFile)
    getInsertionsJob = Job.wrapJobFn(getInsertionsOnBranch, halID=halID, genome=args.genome, args=args)

    repeatScoutJob = Job.wrapJobFn(runRepeatScout, genome=args.genome, halID=halID, gffID = getInsertionsJob.rv(), seqID=seqID)
    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=repeatScoutJob.rv(), seqID=seqID, args=args)

    job.addChild(getInsertionsJob)
    getInsertionsJob.addFollowOn(repeatScoutJob)
    repeatScoutJob.addFollowOn(repeatMaskerJob)

    return {'final.gff':repeatMaskerJob.rv(), 'library.fa':repeatScoutJob.rv()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("genome", type=str)
    parser.add_argument("outDir", type=str)


    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)
    parser.add_argument("--maxInsertions", type=int, default=None)

    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)

    parser.add_argument("--distanceThreshold", type=float, default=0.1)
    parser.add_argument("--kmerLength", type=int, default=5)
    parser.add_argument("--substMatrix", type=str, default="blosum80.mat")


    parser.add_argument("--usePoa", action="store_true", default=False, help="Use the POA pipeline")

    parser.add_argument("--localBinaries", action="store_true", default=False)


    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    if os.path.exists(args.outDir):
        print("Directory %s already exists" % args.outDir)
        exit()
    else:
        os.makedirs(args.outDir)

    with Toil(args) as toil:
        halID = toil.importFile(makeURL(args.hal))

        if args.usePoa:
            args.substMatrixID = toil.importFile(makeURL(os.path.join(os.path.dirname(__file__), args.substMatrix)))
            rootJob = Job.wrapJobFn(poaPipeline, halID=halID, genome=args.genome, args=args)
        else:
            rootJob = Job.wrapJobFn(repeatScoutPipeline, halID=halID, args=args)

        results = toil.start(rootJob)
        for filename in results:
            toil.exportFile(results[filename], makeURL(os.path.join(args.outDir, filename)))

if __name__ == "__main__":
    main()
