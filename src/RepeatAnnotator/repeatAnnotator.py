from toil.common import Toil
from toil.job import Job
import subprocess
import argparse
import random
import os

from toil.lib.bioio import logger

def makeURL(path):
    return "file://" + path

def getRootPath():
    import RepeatAnnotator
    i = os.path.abspath(RepeatAnnotator.__file__)
    return os.path.split(i)[0]

def highNFraction(seqs):
    for i in range(10):
        seq = random.choice(seqs)
        base = random.choice(seq)
        if base == 'N' or base == 'n':
            return True
    return False

def annotateInsertions(job, halID, args):
    getInsertionsJob = Job.wrapJobFn(getInsertions, halID = halID, args = args)

    extractFastaSequencesJob = Job.wrapJobFn(extractFastaSequences, halID=halID, insertionsBedID = getInsertionsJob.rv(), args = args)
    getDistancesJob = Job.wrapJobFn(getDistances, extractFastaSequencesJob.rv())

    getInsertionsJob.addFollowOn(extractFastaSequencesJob)
    extractFastaSequencesJob.addFollowOn(getDistancesJob)
    job.addChild(getInsertionsJob)

    return getDistancesJob.rv()
def getInsertions(job, halID, args):
    hal = job.fileStore.readGlobalFile(halID)
    insertions = job.fileStore.getLocalTempFile()
    with open(insertions, 'w') as insertionsWrite:
        subprocess.check_call(["halAlignedExtract", "--complement", hal, args.reference], stdout=insertionsWrite)
    return job.fileStore.writeGlobalFile(insertions)

def extractFastaSequences(job, halID, insertionsBedID, args):
    insertionsBed = job.fileStore.readGlobalFile(insertionsBedID)
    fasta = job.fileStore.getLocalTempFile()
    hal = job.fileStore.readGlobalFile(halID)
    i = 0
    with open(fasta, 'w') as fastaWrite:
        with open(insertionsBed, 'r') as bedRead:
            for bedLine in bedRead:
                chrom, start, end = bedLine.split()
                fastaLines = subprocess.check_output(["hal2fasta", hal, args.reference, "--sequence", chrom, "--start", start, "--length", str(int(end) - int(start))])

                seqLen = sum([len(line) for line in fastaLines[1:]])
                if seqLen > args.maxInsertionSize:
                    continue
                if seqLen < args.minInsertionSize:
                    continue
                if highNFraction(fastaLines[1:]):
                    continue
                for fastaLine in fastaLines:
                    fastaWrite.write(fastaLine)
                i = i + 1
                if i > args.maxInsertions:
                    break
    return job.fileStore.writeGlobalFile(fasta)

class TEInsertion:
    def __init__(self, chrom, start, end, group):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.group = group

def buildClusters(job, distancesID, insertionsBedID, args):
    distancesFile = job.fileStore.readGlobalFile(distancesID)
    insertionsFile = job.fileStore.readGlobalFile(insertionsBedID)
   
    insertions = []
    with open(insertionsFile, 'r') as insertionsRead:
        for line, i in enumerate(insertionsRead):
            chrom, start, end = line.split()
            insertions.append(TEInsertion(chrom = chrom, start = start, end = end, group = i))

    clusterToSeq = {i: [insertions[i]] for i in range(len(insertions))}
    with open(distancesFile, 'r') as distancesRead:
        for line in distancesRead:
            i, j, dist = line.split()
            assert i > j
            clusterToSeq[i].extend(clusterToSeq[j])

def getDistances(job, insertionsFastaID):
    insertionsFasta = job.fileStore.readGlobalFile(insertionsFastaID)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["pairwise_distances", "--sequences", insertionsFasta], stdout=distancesWrite)
    assert False
    return job.fileStore.writeGlobalFile(distances)

def buildClusters(job, distancesID, insertionsFastaID):
    clusters = {}
    insertionsFasta = job.fileStore.readGlobalFile(insertionsFastaID)

def extractConsensus(job, seqsID, args):
    seqs = job.fileStore.readGlobalFile(seqsID)
    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    subprocess.check_call(["poa", seqs, "-read_fasta", seqs, "-po", graph, substMatrix])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignmentFile", type=str)
    parser.add_argument("gffFile", type=str)
    parser.add_argument("reference", type=str)
    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--reference", type=str)
    parser.add_argument("--maxNFraction", type=float, default=0.1)

    parser.add_argument("--kmerLength", type=int, default=10)

    parser.add_argument("--maxInsertions", type=int, default=1000)

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()
    with Toil(args) as toil:
        args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))
        halID = toil.importFile(args.alignmentFile)
        rootJob = Job.wrapJobFn(annotateInsertions, halID=halID, args=args)
        gffID = toil.start(rootJob)
        toil.exportFile(gffID, args.gffFile)


if __name__ == "__main__":
    main()
