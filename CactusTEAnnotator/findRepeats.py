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
    return "file://%s" % os.path.abspath(path)

dockerImage = "cactus-te-annotator:latest"

def runCmd(parameters, args, outfile=None, mode="w"):
    if args.localBinaries:
        cmd = parameters
    else:
        cmd = ["docker", "run", "-it", "--rm", "-v", "%s:/data" % os.getcwd(), dockerImage] + parameters
       
    if outfile:
        with open(outfile, mode=mode) as outfileWrite:
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
    returnValues = {}
    if args.precomputedFiles:
        assert "candidate_TEs.gff" in args.precomputedFileIDs
        assert "candidate_TEs.fa" in args.precomputedFileIDs
        assert "masked_candidate_TEs.fa" in args.precomputedFileIDs
        job.fileStore.logToMaster("Using precomputed masked TEs file")
        returnValues["gff"] = args.precomputedFileIDs["candidate_TEs.gff"]
        returnValues["fasta"] = args.precomputedFileIDs["candidate_TEs.fa"]
        returnValues["masked_fasta"] = args.precomputedFileIDs["masked_candidate_TEs.fa"]
        return returnValues
    hal = job.fileStore.readGlobalFile(halID)
    gff = job.fileStore.getLocalTempFile()
    fasta = job.fileStore.getLocalTempFile()

    cmd = ["getTECandidates", os.path.basename(hal), genome, "--minLength", str(args.minTESize), "--maxLength", str(args.maxTESize), "--outGFF", os.path.basename(gff), "--outFasta", os.path.basename(fasta), "--maxSequences", str(args.maxTECandidatesToProcess)]
    if not includeReverse:
        cmd.extend(["--ignoreReverse"])
    runCmd(parameters=cmd, args=args)

    fastaID = job.fileStore.writeGlobalFile(fasta)
    gffID = job.fileStore.writeGlobalFile(gff)
    returnValues["gff"] = gffID
    returnValues["fasta"] = fastaID
    returnValues["masked_fasta"] = job.addChildJobFn(runTRF, fastaID=fastaID, args=args).rv()
    return returnValues

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

def runLastz(job, fastaID, args, querydepth=None):
    if args.precomputedFiles and "alignments.cigar" in args.precomputedFileIDs:
        job.fileStore.logToMaster("Using precomputed lastz alignments file.")
        return args.precomputedFileIDs["alignments.cigar"]
    fasta = job.fileStore.readGlobalFile(fastaID)
    alignments = job.fileStore.getLocalTempFile()

    parameters = ["lastz", "--notrivial", "--format=cigar", "%s[multiple]" % os.path.basename(fasta), os.path.basename(fasta)]
    if querydepth:
        parameters.extend(["--querydepth=keep,nowarn:%i" % querydepth])
    runCmd(parameters=parameters, outfile=alignments, args=args)

    return job.fileStore.writeGlobalFile(alignments)

def sortAlignments(job, alignmentsID, args):
    """Sort a set of alignments in cigar format by highest to lowest
    score.
    """
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sortedAlignments = job.fileStore.getLocalTempFile()
    runCmd(parameters=["sort", os.path.basename(alignments), "-k10", "-n", "-r"], outfile=sortedAlignments, args=args)
    return job.fileStore.writeGlobalFile(sortedAlignments)

def seedSampleProb(N, n, L):
    #N : number of alignments sampled per family
    #p : probability of starting an alignment at any given 
    #    seed hit
    #n : size of repeat family
    #N = n*(n-1) * (1 - (1 - p)^L)
    #p = 1 - (1 - N/(n*(n-1))) ^ (1/L)
    p = 1.0 - (1.0 - N/(n*(n-1.0)))**(1.0/L)
    return p

def sampleLastzAlignments(job, fastaID, args, files=None, nIterations=0):
    fasta = job.fileStore.readGlobalFile(fastaID)

    if files:
        ignoredSeeds = job.fileStore.readGlobalFile(files["ignoredSeeds"])
        alignments = job.fileStore.readGlobalFile(files["alignments"])
        levels = job.fileStore.readGlobalFile(files["levels"])
    else:
        files = {}
        ignoredSeeds = job.fileStore.getLocalTempFile()
        alignments = job.fileStore.getLocalTempFile()
        levels = job.fileStore.getLocalTempFile()

    ignoredSeedsSet = set()
    with open(ignoredSeeds, "r") as ignoredSeedsFile:
        for line in ignoredSeedsFile:
            seed = line.strip().rstrip()
            ignoredSeedsSet.add(str(seed))

    seedCounts = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--tableonly=count", "%s[multiple]" % os.path.basename(fasta)], outfile=seedCounts, args=args)
    highestMultiplicity = 0
    with open(seedCounts, "r") as seedCountsFile:
        for line in seedCountsFile:
            info, count = line.split()
            seed, unpacked = info.split("/")
            if not seed in ignoredSeedsSet and count > highestMultiplicity:
                highestMultiplicity = count
    
    p = seedSampleProb(N = highestMultiplicity, n = highestMultiplicity, L = args.assumedTELength)
    with open(levels, "a") as levelsFile:
        levelsFile.write("%s\n", str(p))

    runCmd(parameters=["lastz", "--format=cigar", "--notrivial", "--ignoredSeeds=%s" % os.path.basename(ignoredSeeds), "--baseSamplingRate=%E" % p, "%s[multiple]" % os.path.basename(fasta), os.path.basename(fasta)], outfile=alignments, mode="a", args=args)


    seedCountsFile = job.fileStore.getLocalTempFile()
    runCmd(parameters=["lastz", "--tableonly=count", "%s[multiple]" % fasta], outfile=seedCountsFile, args=args)
    
    runCmd(parameters=["getCoveredSeeds", "--seeds", os.path.basename(seedCountsFile), "--alignments", os.path.basename(alignments), "--sequences", os.path.basename(fasta)], outfile=ignoredSeeds, mode="a", args=args)

    uniqueSeeds = job.fileStore.getLocalTempFile()
    runCmd(parameters=["uniq", os.path.basename(ignoredSeeds)], outfile=uniqueSeeds, args=args)
    
    files["ignoredSeeds"] = job.fileStore.writeGlobalFile(uniqueSeeds)
    files["alignments"] = job.fileStore.writeGlobalFile(alignments)
    files["levels"] = job.fileStore.writeGlobalFile(levels)

    if nIterations > 5 or p > 0.5:
        return files
    else:
        return job.addChildJobFn(sampleLastzAlignments, fastaID=fastaID, files=files, nIterations=nIterations+1).rv()

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

def minhashDistances(job, fastaID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)

    distances = job.fileStore.getLocalTempFile()
    runCmd(parameters=["minhash", "--kmerLength", str(args.kmerLength), "--sequences", os.path.basename(fasta)], outfile=distances, args=args)
    return job.fileStore.writeGlobalFile(distances)

def alignmentDistances(job, sequencesID, alignmentsID, args):
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sequences = job.fileStore.readGlobalFile(sequencesID)

    distances = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getAlignmentDistances", "--sequences", sequences, "--alignments", alignments], outfile=distances, args=args)
    return job.fileStore.writeGlobalFile(distances)

def clusterByDistance(job, distancesID, args):
    distances = job.fileStore.readGlobalFile(distancesID)
    clusters = job.fileStore.getLocalTempFile()
    runCmd(parameters=["buildClusters", "--distances", os.path.basename(distances), "--distanceThreshold", str(args.distanceThreshold)], outfile=clusters, args=args)

    return job.fileStore.writeGlobalFile(clusters)

def makeFamilySequenceFiles(job, clustersID, fastaID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)

    clustersFile = job.fileStore.readGlobalFile(clustersID)
    with open(clustersFile, "r") as clustersRead:
        clusters = [line.split() for line in clustersRead]

    fastaIDs = {}
    for i, seqList in enumerate(clusters):
        if len(seqList) < args.minClusterSize:
            continue
        if args.maxTEsPerFamily and len(seqList) > args.maxTEsPerFamily:
            seqList = random.sample(seqList, args.maxTEsPerFamily)
        cluster_i_fasta = job.fileStore.getLocalTempFile()
        runCmd(parameters=["samtools", "faidx", os.path.basename(fasta)] + seqList, outfile=cluster_i_fasta, args=args)

        cluster_i_fastaID = job.fileStore.writeGlobalFile(cluster_i_fasta)
        clusterName = "cluster_%d" % i
        fastaIDs[clusterName] = cluster_i_fastaID
    return fastaIDs

def chooseRandomConsensusEachCluster(job, fastaID, clustersID, args):
    fasta = job.fileStore.readGlobalFile(fastaID)
    clusters = job.fileStore.readGlobalFile(clustersID)

    library = job.fileStore.getLocalTempFile()
    with open(clusters, "r") as clustersFile:
        for line in clustersFile:
            seqList = line.split()
            seqName = random.choice(seqList)
            runCmd(parameters=["samtools", "faidx", os.path.basename(fasta), seqName], outfile=library, mode="a", args=args)
    
    return job.fileStore.writeGlobalFile(library)

def buildLibraryFromClusters(job, clusterSeqFileIDs, args):
    """Build a poa graph for each family and extract
    repeat elements from each graph. Return a library of
    repeat elements in fasta format.
    """
    returnValues = {}
    consensusJobs = []
    for clusterName in clusterSeqFileIDs:
        fastaID = clusterSeqFileIDs[clusterName]

        if args.consensusMethod == "poa-heaviest-path":

            poaJob = Job.wrapJobFn(runPoa, fastaID=fastaID, heaviestBundle=True, args=args)
            consensusJob = Job.wrapJobFn(getConsensusPOAHeaviestPath, graphID=poaJob.rv(), clusterName=clusterName, args=args)
        elif args.consensusMethod == "poa-dense-path":
            poaJob = Job.wrapJobFn(runPoa, fastaID=fastaID, heaviestBundle=False, args=args)
            consensusJob = Job.wrapJobFn(getConsensusPOADensePath, graphID=poaJob.rv(), clusterName=clusterName, args=args)

        job.addChild(poaJob)
        poaJob.addFollowOn(consensusJob)
        consensusJobs.append(consensusJob)
        returnValues["%s.po" % clusterName] = poaJob.rv()
        returnValues["%s_consensus.fa" % clusterName] = consensusJob.rv()

    catFilesJob = Job.wrapJobFn(catFilesJobFn, fileIDs=[consensusJob.rv() for consensusJob in consensusJobs])
    for consensusJob in consensusJobs:
        consensusJob.addFollowOn(catFilesJob)
    job.addChild(catFilesJob)
    returnValues["library.fa"] = catFilesJob.rv()
    return returnValues

def runPoa(job, fastaID, args, heaviestBundle=False):
    fasta = job.fileStore.readGlobalFile(fastaID)
    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    cmd = ["poa", "-read_fasta", os.path.basename(fasta), "-po", os.path.basename(graph), substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    runCmd(parameters=cmd, args=args)
    return job.fileStore.writeGlobalFile(graph)

def getConsensusPOAHeaviestPath(job, graphID, clusterName, args):
    graph = job.fileStore.readGlobalFile(graphID)
    consensusSequences = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getHeaviestBundles", "--namePrefix", clusterName, "--lpo", os.path.basename(graph)], outfile=consensusSequences, args=args)
    return job.fileStore.writeGlobalFile(consensusSequences)

def getConsensusPOADensePath(job, graphID, clusterName, args):
    graph = job.fileStore.readGlobalFile(graphID)
    consensusSequences = job.fileStore.getLocalTempFile()
    job.fileStore.logToMaster("Using dense path consensus method")
    runCmd(parameters=["getConsensusPOA", "--lpo", os.path.basename(graph), "--minConsensusScore", str(args.minConsensusScore), "--namePrefix", clusterName], outfile=consensusSequences, args=args)
    return job.fileStore.writeGlobalFile(consensusSequences)

def workflow(job, halID, genome, args):
    returnValues = {}

    parameters = job.fileStore.getLocalTempFile()
    with open(parameters, "w") as parametersFile:
        parametersFile.write("parameters for this run\n" \
            + "alignment distance threshold for clustering: %f\n" % args.distanceThreshold \
            + "initial clustering method: %s\n" % args.initialClusteringMethod \
            + "consensus method: %s\n" % args.consensusMethod)
    returnValues["parameters.txt"] = job.fileStore.writeGlobalFile(parameters)
    
    #Get candidate TE insertions from the cactus alignment
    getTECandidatesJob = Job.wrapJobFn(getTECandidatesOnBranch, halID=halID, genome=genome, includeReverse=False, args=args)
    job.addChild(getTECandidatesJob)
    returnValues["candidate_TEs.gff"] = getTECandidatesJob.rv('gff')
    returnValues["candidate_TEs.fa"] = getTECandidatesJob.rv('fasta')
    returnValues["masked_candidate_TEs.fa"] = getTECandidatesJob.rv('masked_fasta')
    fastaID = getTECandidatesJob.rv('masked_fasta')
    if args.getCandidateTEsOnly:
        return returnValues

    if args.initialClusteringMethod == 'lastz':
        alignmentsJob = Job.wrapJobFn(runLastz, fastaID=fastaID, args=args, querydepth=2)
        alignmentsID = alignmentsJob.rv()
        getTECandidatesJob.addFollowOn(alignmentsJob)
        
        returnValues["alignments.cigar"] = alignmentsID
        distancesJob = Job.wrapJobFn(alignmentDistances, sequencesID=fastaID, alignmentsID=alignmentsID, args=args)
        alignmentsJob.addFollowOn(distancesJob)
        initialClusteringJob = Job.wrapJobFn(clusterByDistance, distancesID=distancesJob.rv(), args=args)
        distancesJob.addFollowOn(initialClusteringJob)
        returnValues["distances.txt"] = distancesJob.rv()

    elif args.initialClusteringMethod == 'minhash':
        distancesJob = Job.wrapJobFn(minhashDistances, fastaID=fastaID, args=args)
        getTECandidatesJob.addFollowOn(distancesJob)
        initialClusteringJob = Job.wrapJobFn(clusterByDistance, distancesID=distancesJob.rv(), args=args)
        distancesJob.addFollowOn(initialClusteringJob)
        returnValues["distances.txt"] = distancesJob.rv()

    returnValues["clusters.txt"] = initialClusteringJob.rv()
    makeFamilySequenceFilesJob = Job.wrapJobFn(makeFamilySequenceFiles, clustersID=initialClusteringJob.rv(), fastaID=fastaID, args=args)
    initialClusteringJob.addFollowOn(makeFamilySequenceFilesJob)
    returnValues["families"] = makeFamilySequenceFilesJob.rv()

    if args.skipConsensus:
        return returnValues

    if args.consensusMethod.startswith("poa"):
        #Build a poa graph for each cluster and extract 
        #consensus repeat sequences from each graph
        buildLibraryJob = Job.wrapJobFn(buildLibraryFromClusters, clusterSeqFileIDs=makeFamilySequenceFilesJob.rv(), args=args)
        makeFamilySequenceFilesJob.addFollowOn(buildLibraryJob)
        returnValues["library.fa"] = buildLibraryJob.rv("library.fa")
        returnValues["consensus"] = buildLibraryJob.rv()
    elif args.consensusMethod == "random":
        buildLibraryJob = Job.wrapJobFn(chooseRandomConsensusEachCluster, fastaID=fastaID, clustersID=initialClusteringJob.rv(), args=args)
        makeFamilySequenceFilesJob.addFollowOn(buildLibraryJob)
        returnValues["library.fa"] = buildLibraryJob.rv()
    if args.skipRepeatMasker:
        finalGffID = None
    else:
        #Download a section of the genome to annotate
        #with RepeatMasker, using the library we constructed
        hal = job.fileStore.readGlobalFile(halID)
        genomeFile = job.fileStore.getLocalTempFile()
        hal2fastaCmd = ["hal2fasta", os.path.basename(hal), genome]
        if args.chrom and args.start and args.end:
            hal2fastaCmd.extend(["--sequence", args.chrom, "--start", str(start), "--length", str(end - start)])
        runCmd(parameters=hal2fastaCmd, outfile=genomeFile, args=args)

        genomeID = job.fileStore.writeGlobalFile(genomeFile)
        repeatMaskerJob = Job.wrapJobFn(runRepeatMasker, repeatLibraryID=buildLibraryJob.rv(), seqID=genomeID, args=args)
        buildLibraryJob.addFollowOn(repeatMaskerJob)
        finalGffID = repeatMaskerJob.rv()
    returnValues["final_annotations.gff"] = finalGffID

    return returnValues

def importPrecomputedFiles(toil, args):
    args.precomputedFileIDs = {}
    files = os.listdir(args.precomputedFiles)
    for file in files:
        if os.path.isdir(os.path.join(args.precomputedFiles, file)):
            continue
        args.precomputedFileIDs[file] = toil.importFile(makeURL(os.path.join(args.precomputedFiles, file)))


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
    parser.add_argument("--maxTECandidatesToProcess", type=int, default=None)

    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)

    parser.add_argument("--distanceThreshold", type=float, default=0.005)

    parser.add_argument("--heaviestBundlingThreshold", type=float, default=0.8)
    parser.add_argument("--kmerLength", type=int, default=5)
    parser.add_argument("--substMatrix", type=str, default="blosum80.mat")


    parser.add_argument("--getCandidateTEsOnly", action="store_true", default=False, help="")
    parser.add_argument("--lastzExact", action="store_true", default=False)
    parser.add_argument("--skipRepeatMasker", action="store_true", default=False)
    parser.add_argument("--skipConsensus", action="store_true", default=False)

    parser.add_argument("--localBinaries", action="store_true", default=False)

    parser.add_argument("--precomputedFiles", type=str, default=None)

    parser.add_argument("--initialClusteringMethod", type=str, default="lastz")
    parser.add_argument("--minClusterSize", type=int, default=3)
    parser.add_argument("--consensusMethod", type=str, default="poa-dense-path")
    parser.add_argument("--maxTEsPerFamily", type=int, default=100)
    parser.add_argument("--minConsensusScore", type=int, default=200)
    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    if os.path.exists(args.outDir):
        print("Directory %s already exists" % args.outDir)
        exit()
    else:
        os.makedirs(args.outDir)

    with Toil(args) as toil:
        halID = toil.importFile(makeURL(args.hal))
        args.substMatrixID = toil.importFile(makeURL(os.path.join(os.path.dirname(__file__), args.substMatrix)))

        if args.precomputedFiles:
            importPrecomputedFiles(toil, args)

        rootJob = Job.wrapJobFn(workflow, halID=halID, genome=args.genome, args=args)

        results = toil.start(rootJob)
        exportResultsFiles(toil=toil, results=results, outDir=args.outDir)

if __name__ == "__main__":
    main()
