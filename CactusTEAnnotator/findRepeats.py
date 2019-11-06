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

    return {"alignments": job.fileStore.writeGlobalFile(alignments)}

def sortAlignments(job, alignmentsID, args):
    """Sort a set of alignments in cigar format by highest to lowest
    score.
    """
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sortedAlignments = job.fileStore.getLocalTempFile()
    runCmd(parameters=["sort", os.path.basename(alignments), "-k10", "-n", "-r"], outfile=sortedAlignments, args=args)
    return job.fileStore.writeGlobalFile(sortedAlignments)

def sampleLastzAlignments(job, fastaID, args):
    returnValues = {}

    #calculate how shallow we should begin sampling the alignments based on
    #an estimate of the largest repeat family size
    fasta = job.fileStore.readGlobalFile(fastaID)
    largestFamilySize = 0.0
    for line in runCmd(parameters=["lastz", "--tableonly=count", "%s[multiple]" % os.path.basename(fasta)], args=args).split("\n"):
        if line.split() != 2:
            continue
        other, count = line.split()
        count = float(count)
        if count > largestFamilySize:
            largestFamilySize = count

    #assume an average family length (in number of seeds covered)
    averageFamilyLength = 5000.0
    lowestSamplingLevel = 1.0/((largestFamilySize*largestFamilySize) * averageFamilyLength)
    
    #Fill out the other levels, incrementing by an order of magnitude
    levels = []
    level = lowestSamplingLevel
    while (level < 1.0):
        levels.append(level)
        level *= 10.0

    #Always sample at full depth on the last iteration
    levels.append(1.0)

    #Store the levels for debugging
    levelsFile = job.fileStore.getLocalTempFile()
    with open(levelsFile, "w") as f:
        f.write("Largest family size = %d\n" % largestFamilySize)
        f.write(" ".join(levels) + "\n")
    returnValues["levels.txt"] = job.fileStore.writeGlobalFile(levelsFile)

    #Start with every seed treated equally
    seedWeightsID = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    alignmentsID = None
    nLevels = len(levels)
    alignmentJobs = []
    for i in range(nLevels):
        alignmentJobs.append(Job.wrapJobFn(runLastzAndGetCoveredSeeds, fastaID=fastaID, seedWeightsID = seedWeightsID, baseSamplingRate = levels[i], args=args))

        alignmentsID = alignmentJobs[i].rv('alignments')
        seedWeightsID = alignmentJobs[i].rv('seedWeightsRates')

        returnValues["alignments_%f.cigar" % levels[i]] = alignmentsID
        returnValues["seedWeightsRates_%f.txt" % levels[i]] = seedWeightsID

        if i > 0:
            alignmentJobs[i-1].addFollowOn(alignmentJobs[i])
    job.addChild(alignmentJobs[0])

    returnValues["alignments"] = alignmentsID
    return returnValues

def runLastzAndGetCoveredSeeds(job, fastaID, seedWeightsID, baseSamplingRate, args):
    fasta = job.fileStore.readGlobalFile(fastaID)

    seedWeightsFile = job.fileStore.readGlobalFile(seedWeightsID)

    alignments = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--format=cigar", "--notrivial", "--seedWeights=%s" % os.path.basename(seedWeightsFile), "--baseSamplingRate=%f" % baseSamplingRate, "%s[multiple]" % os.path.basename(fasta), os.path.basename(fasta)], outfile=alignments, args=args)

    alignmentsID = job.fileStore.writeGlobalFile(alignments)

    seedWeights = {}
    with open(seedWeightsFile, 'r') as f:
        for line in f:
            seed, rate = line.split()
            seedWeights[seed] = float(rate)

    seedCountsFile = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--tableonly=count", "%s[multiple]" % os.path.basename(fasta)], outfile=seedCountsFile, args=args)

    for line in runCmd(parameters=["getCoveredSeeds", "--seeds", os.path.basename(seedCountsFile), "--alignments", os.path.basename(alignments), "--sequences", os.path.basename(fasta)], args=args).split("\n"):
        if len(line.split()) != 3:
            continue
        seed, count, length = line.split()
        count = float(count)
        length = float(length)

        #Percentage of the n*(n-1)/2 pairwise alignments between a 
        #family of n repeat copies that will be sampled
        #in order to build a transitively complete set of 
        #alignments representing the family
        alignmentSampleFraction = 0.1

        #Adjust for the fact that lastz will explore the entire alignment
        #if any of the seeds contained within it is matched.
        updatedSeedWeight = alignmentSampleFraction * (1.0/(count*length))

        if not seed in seedWeights or seedWeights[seed] > updatedSeedWeight:
            seedWeights[seed] = updatedSeedWeight

    newSeedWeightsFile = job.fileStore.getLocalTempFile()
    with open(newSeedWeightsFile, "w") as fh:
        for seed in seedWeights:
            fh.write("%s %lf\n" % (seed, seedWeights[seed]))
    seedWeightsID = job.fileStore.writeGlobalFile(newSeedWeightsFile)

    return {"alignments": alignmentsID, "seedWeights": seedWeightsID}

def repeatLibraryFromPinchGraph(job, alignmentsID, sequencesID, args):
    """Construct a pinch graph from the set of pairwise alignments
    representing the repeat family. Then use the graph to define repeat
    element boundaries and extract the consensus sequence for each 
    defined element.
    """
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sequences = job.fileStore.readGlobalFile(sequencesID)

    repeatLibrary = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getElementsFromPinchGraph", os.path.basename(alignments), os.path.basename(sequences)], outfile=repeatLibrary, args=args)
    return job.fileStore.writeGlobalFile(repeatLibrary)

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

    if args.lastzExact:
        lastzJob = Job.wrapJobFn(runLastz, fastaID=trfJob.rv(), args=args)
    else:
        lastzJob = Job.wrapJobFn(sampleLastzAlignments, fastaID=trfJob.rv(), args=args)
    trfJob.addFollowOn(lastzJob)
    alignmentsID = lastzJob.rv('alignments')


    #Sort the alignments so that high-scoring homologies are prioritized
    #when building the pinch graph
    sortAlignmentsJob = Job.wrapJobFn(sortAlignments, alignmentsID=alignmentsID, args=args)
    lastzJob.addFollowOn(sortAlignmentsJob)
    alignmentsID = sortAlignmentsJob.rv()

    return {"masked_candidates.fa": trfJob.rv(), "alignments.cigar": alignmentsID, "alignmentFiles": lastzJob.rv()}

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
    parser.add_argument("--lastzExact", action="store_true", default=False)
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
