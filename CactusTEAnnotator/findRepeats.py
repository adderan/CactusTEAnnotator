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

def getFasta(job, hal, genome, chrom, start, end, args):
    hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome]
    if end and start and chrom:
        hal2fastaCmd.extend(["--sequence", chrom, "--start", str(start), "--length", str(end - start)])
    fastaFile = job.fileStore.getLocalTempFile()
    with open(fastaFile, 'a') as fh:
        runCmd(parameters=hal2fastaCmd, streamfile=fh, args=args)
    return fastaFile

def gffToFasta(job, hal, genome, gff, args):
    fasta = job.fileStore.getLocalTempFile()
    with open(fasta, "w") as fastaWrite:
        runCmd(parameters=["getSequencesFromHAL", hal, gff, genome], streamfile=fastaWrite, args=args)
    return fasta
   
def getRootPath():
    import CactusTEAnnotator
    i = os.path.abspath(CactusTEAnnotator.__file__)
    return os.path.split(i)[0]

def getTECandidatesOnBranch(job, halID, genome, args):
    """Use the HAL graph of the alignment to search for candidate TE insertions in this  \
    genome relative to its parent.
        """
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.getLocalTempFile()
    fasta = job.fileStore.getLocalTempFile()

    runCmd(parameters=["getTECandidates", os.path.basename(hal), genome, "--minLength", str(args.minTESize), "--maxLength", str(args.maxTESize), "--outGFF", os.path.basename(gff), "--outFasta", os.path.basename(fasta), "--maxSequences", str(args.maxInsertions)], args=args)

    return job.fileStore.writeGlobalFile(gff)

def runTRF(job, halID, gffID, genome, args):
    gff = job.fileStore.readGlobalFile(gffID)
    hal = job.fileStore.readGlobalFile(halID)
    fasta = gffToFasta(job=job, hal=hal, genome=genome, gff=gff, args=args)

    trfParameters = ["2", "5", "7", "80", "10", "50", "2000"]
    maskedFasta = fasta + ".".join(trfParameters) + ".mask"

    runCmd(parameters=["trf", fasta] + trfParameters + ["-m", "-h", "-ngs"], args=args)
    return job.fileStore.writeGlobalFile(maskedFasta)

def runRepeatScout(job, genome, halID, gffID, seqID):
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.readGlobalFile(gffID)
    seq = job.fileStore.readGlobalFile(seqID)

    seqs = gffToFasta(job=job, hal=hal, genome=genome, gff=gff, args=args)

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
    fasta = gffToFasta(job=job, hal=hal, genome=genome, gff=gff, args=args)

    with open(gff, "r") as gffRead:
        candidateRepeatsList = [GffInfo(line) for line in gffRead]
    candidateRepeats = {info.name:info for info in candidateRepeatsList}

    graph = networkx.Graph()
    for i in candidateRepeats.keys():
        graph.add_node(i)
    for line in runCmd(parameters=["minhash", "--kmerLength", str(args.kmerLength), "--sequences", os.path.basename(fasta), "--distancesOnly"], args=args).split("\n"):
        if line == "":
            continue
        i, j, distance = line.split()
        assert graph.has_node(i)
        assert graph.has_node(j)
        if float(distance) < args.distanceThreshold:
            graph.add_edge(i, j)

    newGffs = []
    for family_number, component in enumerate(networkx.connected_components(graph)):
        if len(component) < 2:
            continue
        newGff = job.fileStore.getLocalTempFile()
        job.fileStore.logToMaster("Found component %d size %d" % (family_number, len(component)))
        with open(newGff, "w") as gffWrite:
            for i in component:
                candidateRepeats[i].family = family_number
                gffLine = candidateRepeats[i].printGff()
                gffWrite.write(gffLine)
        newGffs.append(newGff)
    return [job.fileStore.writeGlobalFile(newGff) for newGff in newGffs]


def writeSequencesToFile(gff, seqsFile):
    with open(seqsFile, "w") as seqsWrite:
        for repeatCandidate in repeatCandidates:
            seqsWrite.write(">%s\n" % repeatCandidate.name)
            seqsWrite.write("%s\n" % repeatCandidate.seq)

def buildLibrary_poa(job, halID, genome, gffIDs, args):
    """Build a poa graph for each family and extract
    repeat elements from each graph. Return a library of
    repeat elements in fasta format.
    """
    hal = job.fileStore.readGlobalFile(halID)
    gffFiles = [job.fileStore.readGlobalFile(gffID) for gffID in gffIDs]
    elementsJobs = []
    fastaFiles = [gffToFasta(job=job, hal=hal, genome=genome, gff=gffFile, args=args) for gffFile in gffFiles]

    for clusterNum, fastaFile in enumerate(fastaFiles):
        fastaID = job.fileStore.writeGlobalFile(fastaFile)
        poaJob = Job.wrapJobFn(runPoa, fastaID=fastaID, args=args)
        elementsJob = Job.wrapJobFn(getRepeatElementsFromGraph, graphID=poaJob.rv(), clusterName=clusterNum, args=args)
        job.addChild(poaJob)
        poaJob.addFollowOn(elementsJob)
        elementsJobs.append(elementsJob)

    catFilesJob = Job.wrapJobFn(catFilesJobFn, fileIDs = [elementsJob.rv() for elementsJob in elementsJobs])
    for elementsJob in elementsJobs:
        elementsJob.addFollowOn(catFilesJob)
    job.addChild(catFilesJob)
    return catFilesJob.rv()

def runPoa(job, fastaID, args, heaviestBundle=True):
    fasta = job.fileStore.readGlobalFile(fastaID)
    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    cmd = ["poa", "-read_fasta", os.path.basename(fasta), "-po", os.path.basename(graph), substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    runCmd(parameters=cmd, args=args)
    return job.fileStore.writeGlobalFile(graph)

def getRepeatElementsFromGraph(job, graphID, clusterName, args):
    graph = job.fileStore.readGlobalFile(graphID)
    consensusSequences = job.fileStore.getLocalTempFile()
    with open(consensusSequences, "w") as fh:
        runCmd(parameters=["getHeaviestBundles", "--lpo", os.path.basename(graph)], streamfile=fh, args=args)

    #Give the sequences unique names
    repeatLibrary = job.fileStore.getLocalTempFile()
    with open(consensusSequences, "r") as seqRead:
        with open(repeatLibrary, "w") as repeatLibraryWrite:
            for name, sequence in fastaRead(seqRead):
                repeatLibraryWrite.write(">Cluster_%s_%s\n" % (str(clusterName), name))
                repeatLibraryWrite.write(sequence)
                repeatLibraryWrite.write("\n")


    return job.fileStore.writeGlobalFile(repeatLibrary)

def poaPipeline(job, halID, genome, args):
    #Get candidate TE insertions from the cactus alignment
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=args.genome, args=args)

    trfJob = Job.wrapJobFn(runTRF, halID=halID, gffID=getTECandidatesJob.rv(), genome=args.genome, args=args)

    #Cluster the sequences into large families with minhash
    #so that the graph construction step is computationally 
    #possible
    initialClusteringJob = Job.wrapJobFn(minhashClustering, gffID=trfJob.rv(), halID=halID, genome=genome, args=args)

    #Build a poa graph for each cluster and extract 
    #consensus repeat sequences from each graph
    buildLibraryJob = Job.wrapJobFn(buildLibrary_poa, halID=halID, genome=args.genome, gffIDs = initialClusteringJob.rv(), args=args)


    hal = job.fileStore.readGlobalFile(halID)
    genomeFile = getFasta(job=job, hal=hal, genome=args.genome, chrom=args.chrom, start=args.start, end=args.end, args=args)
    genomeID = job.fileStore.writeGlobalFile(genomeFile)

    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=buildLibraryJob.rv(), seqID=genomeID, args=args)


    job.addChild(getTECandidatesJob)
    getTECandidatesJob.addFollowOn(trfJob)
    trfJob.addFollowOn(initialClusteringJob)
    initialClusteringJob.addFollowOn(buildLibraryJob)
    buildLibraryJob.addFollowOn(repeatMaskerJob)

    #concatenate the gffs for each cluster and return them
    #for debugging
    catFilesJob = Job.wrapJobFn(catFilesJobFn, fileIDs=initialClusteringJob.rv())
    repeatMaskerJob.addFollowOn(catFilesJob)

    return {'clusters.gff': catFilesJob.rv(), 'masked_candidates.gff': trfJob.rv(), 'final.gff': repeatMaskerJob.rv(), 'library.fa': buildLibraryJob.rv()}


def repeatScoutPipeline(job, halID, args):
    hal = job.fileStore.readGlobalFile(halID)

    genomeFile = getFasta(job=job, hal=hal, genome=args.genome, chrom=args.chrom, start=args.start, end=args.end, args=args)
    genomeID = job.fileStore.writeGlobalFile(genomeFile)
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=args.genome, args=args)

    repeatScoutJob = Job.wrapJobFn(runRepeatScout, genome=args.genome, halID=halID, gffID = getTECandidatesJob.rv(), seqID=seqID)
    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=repeatScoutJob.rv(), seqID=genomeID, args=args)

    job.addChild(getTECandidatesJob)
    getTECandidatesJob.addFollowOn(repeatScoutJob)
    repeatScoutJob.addFollowOn(repeatMaskerJob)

    return {'final.gff':repeatMaskerJob.rv(), 'library.fa':repeatScoutJob.rv()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("genome", type=str)
    parser.add_argument("outDir", type=str)


    parser.add_argument("--minTESize", type=int, default=100)
    parser.add_argument("--maxTESize", type=int, default=10000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)
    parser.add_argument("--maxInsertions", type=int, default=None)

    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)

    parser.add_argument("--distanceThreshold", type=float, default=0.1)

    parser.add_argument("--heaviestBundlingThreshold", type=float, default=0.8)
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
