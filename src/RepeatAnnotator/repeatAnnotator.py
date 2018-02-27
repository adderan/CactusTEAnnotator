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
    getDistancesJob = Job.wrapJobFn(getDistances, insertionsFastaID=extractFastaSequencesJob.rv(0))
    buildClustersJob = Job.wrapJobFn(buildClusters, distancesID=getDistancesJob.rv(), insertionsBedID=extractFastaSequencesJob.rv(1), args=args)
    writeGFFJob = Job.wrapJobFn(writeGFF, elements=buildClustersJob.rv())

    getInsertionsJob.addFollowOn(extractFastaSequencesJob)
    extractFastaSequencesJob.addFollowOn(getDistancesJob)
    getDistancesJob.addFollowOn(buildClustersJob)
    buildClustersJob.addFollowOn(writeGFFJob)
    job.addChild(getInsertionsJob)

    return writeGFFJob.rv()

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
    filteredBed = job.fileStore.getLocalTempFile()
    i = 0
    with open(filteredBed, 'w') as filteredBedWrite:
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
                    filteredBedWrite.write(bedLine)
                    i = i + 1
                    if args.maxInsertions > 0 and i > args.maxInsertions:
                        break
    return (job.fileStore.writeGlobalFile(fasta), job.fileStore.writeGlobalFile(filteredBed))

def getDistances(job, insertionsFastaID):
    insertionsFasta = job.fileStore.readGlobalFile(insertionsFastaID)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["pairwise_distances", "--sequences", insertionsFasta], stdout=distancesWrite)
    return job.fileStore.writeGlobalFile(distances)

class TEInsertion:
    def __init__(self, chrom, start, end, group):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.group = group

def buildClusters(job, distancesID, insertionsBedID, args):
    distances = job.fileStore.readGlobalFile(distancesID)
    insertionsBed = job.fileStore.readGlobalFile(insertionsBedID)
   
    insertions = []
    with open(insertionsBed, 'r') as insertionsRead:
        for i, line in enumerate(insertionsRead):
            chrom, start, end = line.split()
            insertions.append(TEInsertion(chrom = chrom, start = start, end = end, group = i))

    clusterToSeq = {i: [insertions[i]] for i in range(len(insertions))}
    with open(distances, 'r') as distancesRead:
        for line in distancesRead:
            i, j, dist = line.split()
            assert i > j
            if dist > args.distanceThreshold:
                clusterToSeq[i].extend(clusterToSeq[j])
                insertions[j].group = i
    return insertions

def writeGFF(job, elements):
    gff = job.fileStore.getLocalTempFile()
    with open(gff, 'w') as gffWrite:
        for element in elements:
            gffWrite.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\n" % 
                    (element.chrom, "cactus_repeat_annotation",
                    "repeat_copy", start, end, 0, '+', '.', group))

    return job.fileStore.writeGlobalFile(writeGFF)

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
    parser.add_argument("--distanceThreshold", type=float, default=0.2)

    parser.add_argument("--maxInsertions", type=int, default=0)

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()
    with Toil(args) as toil:
        args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))
        halID = toil.importFile(args.alignmentFile)
        rootJob = Job.wrapJobFn(annotateInsertions, halID=halID, args=args)
        if args.restart:
            gffID = toil.restart()
        else:
            gffID = toil.start(rootJob)
        toil.exportFile(gffID, args.gffFile)


if __name__ == "__main__":
    main()
