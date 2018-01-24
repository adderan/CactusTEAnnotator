from toil.common import Toil
from toil.job import Job
import subprocess
import argparse
import random
import os

from sonLib.bioio import fastaRead, fastaWrite
from toil.lib.bioio import logger

class Element:
    def __init__(self, chrom, start, end, group, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.group = group
        self.name = name

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

def extractFastaSequence(hal, reference, chrom, start, end):
    fastaLines = subprocess.check_output(["hal2fasta", hal, reference, "--sequence", chrom, "--start", str(start), "--length", str(int(end) - int(start))]).split("\n")
    return fastaLines

def annotateInsertions(job, halID, args):
    """Produce a set of TE annotations for one branch of a Cactus alignment.
        """
    getInsertionsJob = Job.wrapJobFn(getInsertions, halID=halID, args = args)
    job.addChild(getInsertionsJob)
    insertions = getInsertionsJob.rv()

    getDistancesJob = Job.wrapJobFn(getDistances, halID=halID, elements=insertions, args=args)
    getInsertionsJob.addFollowOn(getDistancesJob)

    joinByDistanceJob = Job.wrapJobFn(joinByDistance, distancesID=getDistancesJob.rv(), elements=insertions, args=args)
    getDistancesJob.addFollowOn(joinByDistanceJob)

    writeGFFJob = Job.wrapJobFn(writeGFF, elements=joinByDistanceJob.rv())
    joinByDistanceJob.addFollowOn(writeGFFJob)

    return writeGFFJob.rv()

def getInsertions(job, halID, args):
    """Extract candidate TE insertions from the Cactus alignment.
        """
    hal = job.fileStore.readGlobalFile(halID)

    bed = job.fileStore.getLocalTempFile()
    with open(bed, 'w') as bedWrite:
        subprocess.check_call(["halAlignedExtract", "--complement", hal, args.reference], stdout=bedWrite)
    insertions = []
    i = 0
    with open(bed, 'r') as bedRead:
        for bedLine in bedRead:
            if len(bedLine.split()) != 3:
                continue
            chrom, start, end = bedLine.split()
            start = int(start)
            end = int(end)
            if end - start > args.maxInsertionSize:
                continue
            if end - start < args.minInsertionSize:
                continue
            fastaLines = subprocess.check_output(["hal2fasta", hal, args.reference, "--sequence", chrom, "--start", str(start), "--length", str(end - start)])
            assert len(fastaLines) > 0
            assert len(fastaLines[0]) > 0

            if highNFraction(fastaLines[1:]):
                continue
            if args.maxInsertions > 0 and i > args.maxInsertions:
                break
            insertions.append(Element(chrom=chrom, start=start, end=end, group=i, name="cte%d" % i))
            i = i + 1
    return insertions

def getDistances(job, halID, elements, args):
    hal = job.fileStore.readGlobalFile(halID)
    seqs = job.fileStore.getLocalTempFile()
    getFastaSequences(elements=elements, hal=hal, seqs=seqs, args=args)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["pairwise_distances", "--sequences", seqs], stdout=distancesWrite)
    return job.fileStore.writeGlobalFile(distances)


def joinByDistance(job, distancesID, elements, args):
    """Join candidate TE insertions into coarse families based on Jaccard distance. 
    TE insertions transitively joined into the same family if their distance is 
    below the configured threshold.
    """
    distances = job.fileStore.readGlobalFile(distancesID)

    clusterToElement = {i: [elements[i]] for i in range(len(elements))}
    with open(distances, 'r') as distancesRead:
        for line in distancesRead:
            i, j, dist = line.split()
            i = int(i)
            j = int(j)
            dist = float(dist)
            assert i > j
            if dist > args.distanceThreshold:
                group_i = elements[i].group
                group_j = elements[j].group
                if group_i == group_j:
                    continue
                group_a = min(group_i, group_j)
                group_b = max(group_i, group_j)

                clusterToElement[group_a].extend(clusterToElement[group_b])
                for element in clusterToElement[group_b]:
                    element.group = group_a
                del clusterToElement[group_b]

    filteredElements = []
    for cluster in clusterToElement:
        if len(clusterToElement[cluster]) >= args.minClusterSize:
            filteredElements.extend(clusterToElement[cluster])
    return filteredElements

def writeGFF(job, elements):
    """Write a set of TE annotations in GFF format.
    """
    gff = job.fileStore.getLocalTempFile()
    with open(gff, 'w') as gffWrite:
        for element in elements:
            gffWrite.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (element.chrom, "cactus_repeat_annotation", element.name, element.start, element.end, 0, '+', '.', element.group))
    return job.fileStore.writeGlobalFile(gff)

def readGFF(job, gffID):
    elements = []
    gff = job.fileStore.readGlobalFile(gffID)
    with open(gff, 'r') as gffRead:
        for line in gffRead:
            chrom, source, name, start, end, score, strand, a, group = line.split()
            elements.append(Element(chrom=chrom, name=name, start=int(start), end=int(end), group=group))
    return elements

####Partial order alignment###########

def runPoa(job, elements, halID, args):
    """Run partial order alignment with heaviest bundling to assign each TE sequence
    to a consensus.
    """
    logger.info("Aligning group with %d elements" % len(elements))
    logger.info("Elements: %s" % [element.name for element in elements])
    hal = job.fileStore.readGlobalFile(halID)
    seqs = job.fileStore.getLocalTempFile()
    getFastaSequences(elements=elements, hal=hal, seqs=seqs, args=args)
    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    subprocess.check_call(["poa", "-hb", "-read_fasta", seqs, "-po", graph, substMatrix])
    seqToBundle = parsePartialOrder(graph)
    logger.info("Parsed families: %s" % seqToBundle)
    elementsAfterBundling = []
    for element in elements:
        if element.name in seqToBundle:
            element.group = "%s_%s" % (element.group, seqToBundle[element.name])
        elementsAfterBundling.append(element)

    return elementsAfterBundling

def parsePartialOrder(poFile):
    sequences = {}
    name = ""
    with open(poFile, 'r') as poRead:
        for line in poRead:
            if line.startswith("SOURCENAME"):
                name = line.split("=")[1].rstrip()
            elif line.startswith("SOURCEINFO"):
                info = line.split("=")[1].rstrip()
                bundle = info[3]
                if bundle == -1:
                    continue
                #Ignore the consensus sequences
                if name.startswith("CONSENS"):
                    continue
                sequences[name] = bundle
    return sequences

def buildSubfamilies(job, elements, halID, args):
    families = {}
    for element in elements:
        if element.group not in families:
            families[element.group] = []
        families[element.group].append(element)

    updatedElements = []
    for family in families:
        elementsInFamily = families[family]
        updatedElements.append(job.addChildJobFn(runPoa, elements=elementsInFamily, halID=halID, args=args).rv())

    return job.addFollowOnJobFn(buildSubfamilies2, updatedElements).rv()

def buildSubfamilies2(job, updatedElements):
    elements = []
    logger.info("Processing updated elements for %d families" % len(updatedElements))
    for elementsList in updatedElements:
        logger.info("Processing %d elements in family" % len(elementsList))
        elements.extend(elementsList)
    return elements

def getFastaSequences(elements, hal, seqs, args):
    logger.info("Extracting sequences for %d elements" % len(elements))
    with open(seqs, 'w') as seqsWrite:
        for element in elements:
            fastaLines = extractFastaSequence(hal=hal, reference=args.reference, chrom=element.chrom, start=element.start, end=element.end)
            #Write header
            seqsWrite.write(">%s\n" % element.name)
            #Write sequence
            seqsWrite.write(" ".join(fastaLines[1:]))

def addRepeatAnnotatorOptions(parser):
    parser.add_argument("reference", type=str)
    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)

    parser.add_argument("--kmerLength", type=int, default=10)
    parser.add_argument("--distanceThreshold", type=float, default=0.2)
    parser.add_argument("--minClusterSize", type=int, default=2)

    parser.add_argument("--maxInsertions", type=int, default=0)
    parser.add_argument("--buildSubfamilies", action="store_true", default=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignmentFile", type=str)
    parser.add_argument("outGFF", type=str)
    addRepeatAnnotatorOptions(parser)

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()
    with Toil(args) as toil:
        args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))
        halID = toil.importFile(makeURL(args.alignmentFile))
        rootJob = Job.wrapJobFn(annotateInsertions, halID=halID, args=args)
        gffID = toil.start(rootJob)

        toil.exportFile(gffID, makeURL(args.outGFF))


if __name__ == "__main__":
    main()
