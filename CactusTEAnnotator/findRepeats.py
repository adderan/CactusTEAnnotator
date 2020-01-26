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

def buildRepeatLibrary(job, alignmentsID, fastaID, args):
    """Construct a pinch graph from the set of pairwise alignments
    representing the repeat family. Then use the graph to define repeat
    element boundaries and extract the consensus sequence for each 
    defined family.
    """
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sequences = job.fileStore.readGlobalFile(fastaID)

    repeatLibrary = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getConsensusTEs", "--alignments", os.path.basename(alignments), "--sequences", os.path.basename(sequences)], outfile=repeatLibrary, args=args)
    return job.fileStore.writeGlobalFile(repeatLibrary)

def runRepeatMasker(job, repeatLibraryID, seqID, args):
    """Use RepeatMasker to annotate the genome using a
    custom repeat library.
    """
    repeatLibrary = job.fileStore.readGlobalFile(repeatLibraryID)
    seq = job.fileStore.readGlobalFile(seqID)
    repeatMaskerOutput = job.fileStore.getLocalTempDir()
    runCmd(parameters=["RepeatMasker", "-nolow", "-cutoff", "650", "-dir", os.path.basename(repeatMaskerOutput), "-lib", os.path.basename(repeatLibrary), os.path.basename(seq)], args=args)
    outputGff = "%s/%s.out" % (repeatMaskerOutput, os.path.basename(seq))
    job.fileStore.logToMaster("directory contents: %s" % os.listdir(repeatMaskerOutput))
    outputGffID = job.fileStore.writeGlobalFile(outputGff)

    return outputGffID

def alignmentDistances(job, sequencesID, alignmentsID, args):
    alignments = job.fileStore.readGlobalFile(alignmentsID)
    sequences = job.fileStore.readGlobalFile(sequencesID)

    distances = job.fileStore.getLocalTempFile()
    runCmd(parameters=["getAlignmentDistances", "--sequences", sequences, "--alignments", alignments], outfile=distances, args=args)
    return job.fileStore.writeGlobalFile(distances)

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

    buildRepeatLibraryJob = Job.wrapJobFn(buildRepeatLibrary, alignmentsID=alignmentsID, fastaID=fastaID, args=args)
    getTECandidatesJob.addFollowOn(buildRepeatLibraryJob)
    returnValues["library.fa"] = buildRepeatLibraryJob.rv()
    
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
        buildRepeatLibraryJob.addFollowOn(repeatMaskerJob)
        finalGffID = repeatMaskerJob.rv()
        returnValues["annotations.gff"] = finalGffID

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
    parser.add_argument("--skipRepeatMasker", action="store_true", default=False)
    parser.add_argument("--localBinaries", action="store_true", default=False)
    parser.add_argument("--precomputedFiles", type=str, default=None)

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
