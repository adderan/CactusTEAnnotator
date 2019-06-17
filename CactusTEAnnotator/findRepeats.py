import argparse
import random
import os
import shutil
import sys
import subprocess

from toil.job import Job
from toil.common import Toil

from sonLib.bioio import fastaRead, fastaWrite, catFiles, reverseComplement

def makeURL(path):
    return "file://%s" % path

dockerImage = "cactus-te-annotator:latest"

def runCmd(parameters, args, outfile=None):
    if args.localBinaries:
        cmd = parameters
    else:
        cmd = ["docker", "run", "-it", "--rm", "-v", "%s:/data" % os.getcwd(), dockerImage] + parameters
       
    if outfile:
        with open(outfile, "w") as outfileWrite:
            subprocess.check_call(cmd, stdout=outfileWrite)
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

def gffToFasta(job, hal, genome, gff, args):
    fasta = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getSequencesFromHAL", hal, gff, genome], outfile=fasta, args=args)
    return fasta
   
def getRootPath():
    import CactusTEAnnotator
    i = os.path.abspath(CactusTEAnnotator.__file__)
    return os.path.split(i)[0]

def getTECandidatesOnBranch(job, halID, genome, args, includeReverse=True):
    """Use the HAL graph of the alignment to search for candidate TE insertions in this  \
    genome relative to its parent.
        """
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.getLocalTempFile()
    fasta = job.fileStore.getLocalTempFile()

    cmd = ["getTECandidates", os.path.basename(hal), genome, "--minLength", str(args.minTESize), "--maxLength", str(args.maxTESize), "--outGFF", os.path.basename(gff), "--outFasta", os.path.basename(fasta), "--maxSequences", str(args.maxInsertions)]
    if not includeReverse:
        cmd.extend(["--ignoreReverse"])
    runCmd(parameters=cmd, args=args)

    return {'fasta': job.fileStore.writeGlobalFile(fasta), 'gff': job.fileStore.writeGlobalFile(gff)}

def runTRF(job, fastaID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)

    trfParameters = ["2", "5", "7", "80", "10", "50", "2000"]
    maskedFasta = os.path.basename(fasta) + "." + ".".join(trfParameters) + ".mask"

    runCmd(parameters=["trf", os.path.basename(fasta)] + trfParameters + ["-m", "-h", "-ngs"], args=args)

    #TRF masks the low-complexity sequences with Ns, so discard
    #sequences with too many
    NFilteredFasta = job.fileStore.getLocalTempFile()
    runCmd(parameters=["filterNs", maskedFasta, str(args.maxNFraction)], outfile=NFilteredFasta, args=args)
    return job.fileStore.writeGlobalFile(NFilteredFasta)

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

def runLastz(job, fastaID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)
    alignments = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--notrivial", "--format=cigar", "%s[multiple]" % os.path.basename(fasta), os.path.basename(fasta)], outfile=alignments, args=args)

    return job.fileStore.writeGlobalFile(alignments)

def sampleLastzAlignments(job, fastaID, args):
    levels = [0.01, 0.1, 1.0]
    fasta = job.fileStore.readGlobalFile(fastaID)

    returnValues = {}
    samplingRatesID = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    alignmentsID = None
    nLevels = len(levels)
    alignmentJobs = []
    for i in range(nLevels):
        alignmentJobs.append(Job.wrapJobFn(runLastzAndGetCoveredSeeds, fastaID=fastaID, samplingRatesID = samplingRatesID, baseSamplingRate = levels[i], args=args))

        alignmentsID = alignmentJobs[i].rv('alignments')
        samplingRatesID = alignmentJobs[i].rv('samplingRates')

        returnValues["alignments_%f.cigar" % levels[i]] = alignmentsID
        returnValues["samplingRates_%f.txt" % levels[i]] = samplingRatesID

        if i > 0:
            alignmentJobs[i-1].addFollowOn(alignmentJobs[i])
    job.addChild(alignmentJobs[0])

    returnValues["alignments"] = alignmentsID
    return returnValues

def runLastzAndGetCoveredSeeds(job, fastaID, samplingRatesID, baseSamplingRate, args):
    fasta = job.fileStore.readGlobalFile(fastaID)

    samplingRatesFile = job.fileStore.readGlobalFile(samplingRatesID)

    alignments = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--format=cigar", "--notrivial", "--samplingRates=%s" % os.path.basename(samplingRatesFile), "--baseSamplingRate=%f" % baseSamplingRate, "%s[multiple]" % os.path.basename(fasta), os.path.basename(fasta)], outfile=alignments, args=args)

    alignmentsID = job.fileStore.writeGlobalFile(alignments)

    samplingRates = {}
    with open(samplingRatesFile, 'r') as f:
        for line in f:
            seed, rate = line.split()
            samplingRates[seed] = float(rate)

    seedCountsFile = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--tableonly=count", "%s[multiple]" % os.path.basename(fasta)], outfile=seedCountsFile, args=args)

    for line in runCmd(parameters=["getCoveredSeeds", "--seeds", os.path.basename(seedCountsFile), "--alignments", os.path.basename(alignments), "--sequences", os.path.basename(fasta)], args=args).split("\n"):
        if len(line.split()) != 3:
            continue
        seed, count, length = line.split()
        count = int(count)
        length = int(length)
        updatedSamplingRate = baseSamplingRate * (1.0/(count*length))
        if not seed in samplingRates or samplingRates[seed] > updatedSamplingRate:
            samplingRates[seed] = updatedSamplingRate

    newSamplingRatesFile = job.fileStore.getLocalTempFile()
    with open(newSamplingRatesFile, "w") as fh:
        for seed in samplingRates:
            fh.write("%s %lf\n" % (seed, samplingRates[seed]))
    samplingRatesID = job.fileStore.writeGlobalFile(newSamplingRatesFile)

    return {"alignments": alignmentsID, "samplingRates": samplingRatesID}

def runRepeatMasker(job, repeatLibraryID, seqID, args):
    repeatLibrary = job.fileStore.readGlobalFile(repeatLibraryID)
    seq = job.fileStore.readGlobalFile(seqID)
    repeatMaskerOutput = job.fileStore.getLocalTempDir()
    runCmd(parameters=["RepeatMasker", "-nolow", "-cutoff", "650", "-dir", os.path.basename(repeatMaskerOutput), "-lib", os.path.basename(repeatLibrary), os.path.basename(seq)], args=args)
    outputGff = "%s/%s.out" % (repeatMaskerOutput, os.path.basename(seq))
    job.fileStore.logToMaster("directory contents: %s" % os.listdir(repeatMaskerOutput))
    outputGffID = job.fileStore.writeGlobalFile(outputGff)

    return outputGffID

def minhashClustering(job, fastaID, args):
    """Cluster a set of candidate repeat annotations
    by their pairwise Jaccard distances, estimated
    with minhash.
    """
    fasta = job.fileStore.readGlobalFile(fastaID)
    distances = job.fileStore.getLocalTempFile()

    runCmd(parameters=["minhash", "--kmerLength", str(args.kmerLength), "--sequences", os.path.basename(fasta)], outfile=distances, args=args)

    clusters = job.fileStore.getLocalTempFile()
    runCmd(parameters=["build_clusters", os.path.basename(distances), "--distanceThreshold", str(args.distanceThreshold)], outfile=clusters, args=args)

    return job.fileStore.writeGlobalFile(clusters)

def buildLibrary_poa(job, fastaID, clustersID, args):
    """Build a poa graph for each family and extract
    repeat elements from each graph. Return a library of
    repeat elements in fasta format.
    """
    fasta = job.fileStore.readGlobalFile(fastaID)

    clustersFile = job.fileStore.readGlobalFile(clustersID)
    with open(clustersFile, "r") as clustersRead:
        clusters = [line.split() for line in clustersRead]

    elementsJobs = []
    for i, seqList in enumerate(clusters):
        cluster_i_fasta = job.fileStore.getLocalTempFile()
        runCmd(parameters=["samtools", "faidx", os.path.basename(fasta)] + seqList, outfile=cluster_i_fasta, args=args)

        cluster_i_fastaID = job.fileStore.writeGlobalFile(cluster_i_fasta)
        poaJob = Job.wrapJobFn(runPoa, fastaID=cluster_i_fastaID, args=args)
        elementsJob = Job.wrapJobFn(getRepeatElementsFromGraph, graphID=poaJob.rv(), clusterName=i, args=args)
        job.addChild(poaJob)
        poaJob.addFollowOn(elementsJob)
        elementsJobs.append(elementsJob)

    catFilesJob = Job.wrapJobFn(catFilesJobFn, fileIDs=[elementsJob.rv() for elementsJob in elementsJobs])
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
    runCmd(parameters=["getHeaviestBundles", "--lpo", os.path.basename(graph)], outfile=consensusSequences, args=args)

    #Give the sequences unique names
    repeatLibrary = job.fileStore.getLocalTempFile()
    with open(consensusSequences, "r") as seqRead:
        with open(repeatLibrary, "w") as repeatLibraryWrite:
            for name, sequence in fastaRead(seqRead):
                repeatLibraryWrite.write(">Cluster_%s_%s\n" % (str(clusterName), name))
                repeatLibraryWrite.write(sequence)
                repeatLibraryWrite.write("\n")


    return job.fileStore.writeGlobalFile(repeatLibrary)

def lastzPipeline(job, halID, genome, args):
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=genome, includeReverse=False, args=args)
    job.addChild(getTECandidatesJob)

    trfJob = Job.wrapJobFn(runTRF, fastaID=getTECandidatesJob.rv('fasta'), args=args)
    getTECandidatesJob.addFollowOn(trfJob)

    lastzJob = Job.wrapJobFn(sampleLastzAlignments, fastaID=trfJob.rv(), args=args)
    trfJob.addFollowOn(lastzJob)

    return {"masked_candidates.fa": trfJob.rv(), "alignments.cigar": lastzJob.rv('alignments'), "alignmentFiles": lastzJob.rv()}

def poaPipeline(job, halID, genome, args):
    #Get candidate TE insertions from the cactus alignment
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=genome, args=args)
    job.addChild(getTECandidatesJob)

    #Mask low complexity regions and simple repeats with TRF
    trfJob = Job.wrapJobFn(runTRF, fastaID=getTECandidatesJob.rv('fasta'), args=args)
    getTECandidatesJob.addFollowOn(trfJob)

    #Cluster the sequences into large families with minhash
    #so that the graph construction step is computationally 
    #possible
    initialClusteringJob = Job.wrapJobFn(minhashClustering, fastaID=trfJob.rv(), args=args)
    trfJob.addFollowOn(initialClusteringJob)

    #Build a poa graph for each cluster and extract 
    #consensus repeat sequences from each graph
    buildLibraryJob = Job.wrapJobFn(buildLibrary_poa, clustersID=initialClusteringJob.rv(), fastaID=trfJob.rv(), args=args)
    initialClusteringJob.addFollowOn(buildLibraryJob)

    if args.skipRepeatMasker:
        finalGffID = None
    else:
        #Download a section of the genome to annotate
        #with RepeatMasker, using the library we constructed
        hal = job.fileStore.readGlobalFile(halID)
        genomeFile = job.fileStore.getLocalTempFile()
        hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome]
        if args.chrom and args.start and args.end:
            hal2fastaCmd.extend(["--sequence", chrom, "--start", str(start), "--length", str(end - start)])
        runCmd(parameters=hal2fastaCmd, outfile=genomeFile, args=args)

        genomeID = job.fileStore.writeGlobalFile(genomeFile)
        repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=buildLibraryJob.rv(), seqID=genomeID, args=args)
        buildLibraryJob.addFollowOn(repeatMaskerJob)
        finalGffID = repeatMaskerJob.rv()

    return {'candidates.gff': getTECandidatesJob.rv('gff'), 'masked_candidates.fa': trfJob.rv(),'clusters.txt': initialClusteringJob.rv(), 'final.gff': finalGffID, 'library.fa': buildLibraryJob.rv()}


def repeatScoutPipeline(job, halID, args):
    hal = job.fileStore.readGlobalFile(halID)

    #genomeFile = getFasta(job=job, hal=hal, genome=args.genome, chrom=args.chrom, start=args.start, end=args.end, args=args)
    genomeID = job.fileStore.writeGlobalFile(genomeFile)
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=args.genome, args=args)

    repeatScoutJob = Job.wrapJobFn(runRepeatScout, genome=args.genome, halID=halID, gffID = getTECandidatesJob.rv(), seqID=seqID)
    repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=repeatScoutJob.rv(), seqID=genomeID, args=args)

    job.addChild(getTECandidatesJob)
    getTECandidatesJob.addFollowOn(repeatScoutJob)
    repeatScoutJob.addFollowOn(repeatMaskerJob)

    return {'final.gff':repeatMaskerJob.rv(), 'library.fa':repeatScoutJob.rv()}

def exportResultsFiles(toil, results, outDir):
    for item in results:
        if results[item] is None:
            continue
        if isinstance(results[item], dict):
            #Store any extra files returned from a job in their own directory
            subDir = os.path.join(outDir, item)
            os.makedirs(subDir)
            exportResultsFiles(toil=toil, results=results[item], outDir=subDir)
        else:
            toil.exportFile(results[item], makeURL(os.path.join(outDir, item)))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("genome", type=str)
    parser.add_argument("outDir", type=str)


    parser.add_argument("--minTESize", type=int, default=100)
    parser.add_argument("--maxTESize", type=int, default=10000)
    parser.add_argument("--maxNFraction", type=float, default=0.2)
    parser.add_argument("--maxInsertions", type=int, default=None)

    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)

    parser.add_argument("--distanceThreshold", type=float, default=0.01)

    parser.add_argument("--heaviestBundlingThreshold", type=float, default=0.8)
    parser.add_argument("--kmerLength", type=int, default=5)
    parser.add_argument("--substMatrix", type=str, default="blosum80.mat")


    parser.add_argument("--usePoa", action="store_true", default=False, help="Use the POA pipeline")
    parser.add_argument("--useLastz", action="store_true", default=False, help="")
    parser.add_argument("--skipRepeatMasker", action="store_true", default=False)

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
        elif args.useLastz:
            rootJob = Job.wrapJobFn(lastzPipeline, halID=halID, genome=args.genome, args=args)
        else:
            rootJob = Job.wrapJobFn(repeatScoutPipeline, halID=halID, args=args)

        results = toil.start(rootJob)
        exportResultsFiles(toil=toil, results=results, outDir=args.outDir)

if __name__ == "__main__":
    main()
