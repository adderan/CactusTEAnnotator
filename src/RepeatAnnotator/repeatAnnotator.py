from toil.common import Toil
from toil.job import Job
import subprocess
import argparse
import random
import os
import shutil

from sonLib.bioio import fastaRead, fastaWrite, catFiles, reverseComplement
from toil.lib.bioio import logger
import networkx

import RepeatAnnotator.treeBuilding as treeBuilding
from sonLib.nxnewick import NXNewick

def catFiles(fileList, target):
    with open(target, 'w') as targetWrite:
        for f in fileList:
            with open(f, 'r') as read:
                lines = read.readlines()
                targetWrite.write("".join(lines))
                targetWrite.write("\n")

class Element:
    def __init__(self, chrom, start, end, family, subfamily, name, seqID, revCompSeqID, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.family = family
        self.subfamily = None
        self.name = name
        self.seqID = seqID
        self.revCompSeqID = revCompSeqID
        self.strand = strand

def makeURL(path):
    return "file://" + path

def getFastaSequence(hal, chrom, start, end, args):
    fastaLines = subprocess.check_output(["hal2fasta", hal, args.reference, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]).strip().split("\n")
    return fastaLines

def getRootPath():
    import RepeatAnnotator
    i = os.path.abspath(RepeatAnnotator.__file__)
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

def annotateInsertions(job, halID, args):
    """Produce a set of TE annotations for one branch of a Cactus alignment.
        """
    getInsertionsJob = Job.wrapJobFn(getInsertions, halID=halID, args = args)
    job.addChild(getInsertionsJob)
    insertions = getInsertionsJob.rv()

    getDistancesJob = Job.wrapJobFn(getDistances, elements=insertions, args=args)
    getInsertionsJob.addFollowOn(getDistancesJob)

    joinByDistanceJob = Job.wrapJobFn(joinByDistance, distancesID=getDistancesJob.rv(), elements=insertions, threshold=args.distanceThreshold, args=args)
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
            fastaLines = getFastaSequence(hal, chrom, start, end, args)
            assert len(fastaLines) > 0
            assert len(fastaLines[0]) > 0

            if highNFraction(fastaLines[1:]):
                continue
            if args.maxInsertions > 0 and i > args.maxInsertions:
                break

            seqFasta = job.fileStore.getLocalTempFile()
            with open(seqFasta, 'w') as seqFastaWrite:
                seqFastaWrite.write(">%d\n" % i)
                seqFastaWrite.write("\n".join(fastaLines[1:]))
                seqFastaWrite.write("\n")


            insertions.append(Element(chrom=chrom, start=start, end=end, family=None, subfamily=None, name=str(i), seqID=job.fileStore.writeGlobalFile(seqFasta), revCompSeqID=None, strand="+"))

            i = i + 1
    job.fileStore.logToMaster("Found %d insertions on branch" % len(insertions))
    return insertions

def getDistances(job, elements, args):
    seqs = job.fileStore.getLocalTempFile()
    seqFiles = [job.fileStore.readGlobalFile(element.seqID) for element in elements]
    catFiles(seqFiles, seqs)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["pairwise_distances", "--kmerLength", "5", "--sequences", seqs], stdout=distancesWrite)
    return job.fileStore.writeGlobalFile(distances)

def joinByDistance(job, distancesID, elements, threshold, args, subfamily=False):
    """Join candidate TE insertions into coarse families based on 
    Jaccard distance. TE insertions are transitively joined into 
    the same family if their distance is below the configured threshold.
    """

    def complementName(name):
        return name + "-"

    distances = job.fileStore.readGlobalFile(distancesID)
    nameToElement = {element.name:element for element in elements}

    graph = networkx.Graph()
    for element in elements:
        graph.add_node(element.name)
        graph.add_node(complementName(element.name))
    with open(distances, 'r') as distancesRead:
        for line in distancesRead:
            i, j, dist, strand = line.split()
            if strand == "-":
                j = complementName(j)
            dist = float(dist)
            assert graph.has_node(i)
            assert graph.has_node(j)
            if dist > threshold:
                graph.add_edge(i, j, weight=dist)
    i = 0
    newElements = []
    for component in networkx.connected_components(graph):
        if len(component) < args.minClusterSize:
            continue
        for elementName in component:
            if elementName.endswith("-"):
                elementName = elementName[:-1]
                strand = "-"
            else:
                strand = "+"

            oldElement = nameToElement[elementName]
            newElements.append(Element(chrom=oldElement.chrom, start=oldElement.start, end=oldElement.end, family=i, subfamily=oldElement.subfamily, name=oldElement.name, seqID=oldElement.seqID, revCompSeqID=oldElement.revCompSeqID, strand=strand))
        i = i + 1

    return newElements

def writeGFF(job, elements):
    """Write a set of TE annotations in GFF format.
    """
    gff = job.fileStore.getLocalTempFile()
    with open(gff, 'w') as gffWrite:
        for element in elements:
            if element.subfamily:
                family = "%s_%s" % (element.family, element.subfamily)
            else:
                family = element.family

            gffWrite.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (element.chrom, "cactus_repeat_annotation", element.name, element.start, element.end, 0, element.strand, '.', element.family))
    return job.fileStore.writeGlobalFile(gff)

def readGFF(job, halID, gffID, args):
    elements = []
    gff = job.fileStore.readGlobalFile(gffID)
    hal = job.fileStore.readGlobalFile(halID)
    with open(gff, 'r') as gffRead:
        for line in gffRead:
            chrom, source, name, start, end, score, strand, a, family = line.split()

            fastaLines = getFastaSequence(hal=hal, chrom=chrom, start=int(start), end=int(end), args=args)


            #Replace chromosome with name
            fastaLines[0] = ">%s" % name
            seqFile = job.fileStore.getLocalTempFile()
            with open(seqFile, 'w') as seqFileWrite:
                seqFileWrite.write("\n".join(fastaLines))

            revCompSeqFile = job.fileStore.getLocalTempFile()
            with open(revCompSeqFile, 'w') as seqFileWrite:
                seqFileWrite.write(">%s\n" % name)
                seqFileWrite.write(reverseComplement("\n".join(fastaLines)))
                seqFileWrite.write("\n")

            elements.append(Element(chrom=chrom, name=name, start=int(start), end=int(end), family=family, subfamily=None, seqID=job.fileStore.writeGlobalFile(seqFile), revCompSeqID=job.fileStore.writeGlobalFile(revCompSeqFile), strand=strand))
    return elements

####Partial order alignment###########

def runPoa(job, elements, heaviestBundle, familyNumber, args):
    logger.info("Aligning family with %d elements" % len(elements))
    logger.info("Elements: %s" % [element.name for element in elements])
    seqs = job.fileStore.getLocalTempFile()
    seqFiles = []
    for element in elements:
        if element.strand == "+":
            seqFiles.append(job.fileStore.readGlobalFile(element.seqID))
        elif element.strand == "-":
            seqFiles.append(job.fileStore.readGlobalFile(element.revCompSeqID))

    catFiles(seqFiles, seqs)

    graph = job.fileStore.getLocalTempFile()
    substMatrix = job.fileStore.readGlobalFile(args.substMatrixID)
    cmd = ["poa", "-read_fasta", seqs, "-po", graph, substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    subprocess.check_call(cmd)

    #Store the graph for debugging
    if args.graphsDir:
        if not os.path.isdir(args.graphsDir):
            os.mkdir(args.graphsDir)
        shutil.copyfile(graph, os.path.join(args.graphsDir, "%d.po" % int(familyNumber)))
        shutil.copyfile(seqs, os.path.join(args.graphsDir, "%d.fa" % int(familyNumber)))

    return job.fileStore.writeGlobalFile(graph)

def updateElements(job, elements, partitioning):
    nameToElement = {element.name:element for element in elements}
    updatedElements = []
    logger.info("Partitioning = %s" % partitioning)
    logger.info("element names = %s" % nameToElement.keys())
    assert(sum([len(partition) for partition in partitioning]) == len(nameToElement))
    for i, partition in enumerate(partitioning):
        for elementName in partition:
            element = nameToElement[elementName]
            element.family = element.family + "_%d" % i
            updatedElements.append(element)
    return updatedElements

def flatten(job, listOfLists):
    flattenedList = []
    for l in listOfLists:
        flattenedList.extend(l)
    return flattenedList

def runTreeBuilding(job, graphID, args):
    """Runs the partitioning tree building method on a sequence graph in 
    PO format.
    """
    graphFile = job.fileStore.readGlobalFile(graphID)
    
    graph = treeBuilding.POGraph(graphFile)
    partitions = graph.getPartitions()
    job.fileStore.logToMaster("Parsed partitions %s" % partitions)
    job.fileStore.logToMaster("tmp file: %s" % os.path.dirname(graphFile))
    tree = treeBuilding.buildTree(graph.threads, partitions)
    partitioning = treeBuilding.getLeafPartitioning(tree)

    return partitioning
    
def getAlignmentDistances(job, graphID, args):
    graph = job.fileStore.readGlobalFile(graphID)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["getAlignmentDistances", graph], stdout=distancesWrite)
    return job.fileStore.writeGlobalFile(distances)


def runNeighborJoining(job, elements, graphID, args):
    if len(elements) < 3:
        return frozenset([frozenset([element.name for element in elements])])

    seqs = job.fileStore.getLocalTempFile()
    seqFiles = [job.fileStore.readGlobalFile(element.seqID) for element in elements]
    graph = job.fileStore.readGlobalFile(graphID)
    distances = job.fileStore.getLocalTempFile()
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["getAlignmentDistances", graph], stdout=distancesWrite)

    partitioning = {}
    for line in subprocess.check_output(["neighborJoining", distances, str(len(elements))]).split("\n"):
        if line == "":
            continue
        nodeNum, family = line.split()
        family = int(family)
        nodeNum = int(nodeNum)
        if not family in partitioning:
            partitioning[family] = set()
        partitioning[family].add(elements[nodeNum].name)


    partitioning = frozenset([frozenset(partition) for partition in partitioning.values()])
    job.fileStore.logToMaster("Partitioning = %s" % partitioning)
    job.fileStore.logToMaster("work dir = %s" % os.path.dirname(seqs))

    return partitioning

def buildSubfamilies(job, elements, args):
    families = {}
    for element in elements:
        if element.family not in families:
            families[element.family] = []
        families[element.family].append(element)

    updatedElements = []
    for family in families:
        elementsInFamily = families[family]
        poaJob = Job.wrapJobFn(runPoa, elements=elementsInFamily, heaviestBundle=False, familyNumber=family, args=args)

        alignmentDistancesJob = Job.wrapJobFn(getAlignmentDistances, graphID=poaJob.rv(), args=args)
        joinByDistanceJob = Job.wrapJobFn(joinByDistance, distancesID=alignmentDistancesJob.rv(), elements=elementsInFamily, threshold=args.alignmentDistanceThreshold, args=args)

        job.addChild(poaJob)
        poaJob.addFollowOn(alignmentDistancesJob)
        alignmentDistancesJob.addFollowOn(joinByDistanceJob)
        updatedElements.append(joinByDistanceJob.rv())

    return job.addFollowOnJobFn(flatten, updatedElements).rv()

def addRepeatAnnotatorOptions(parser):
    parser.add_argument("reference", type=str)
    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)

    parser.add_argument("--kmerLength", type=int, default=10)
    parser.add_argument("--distanceThreshold", type=float, default=0.7)
    parser.add_argument("--minClusterSize", type=int, default=1)

    parser.add_argument("--maxInsertions", type=int, default=0)

    #Methods for splitting families of insertions
    parser.add_argument("--heaviestBundling", action="store_true")
    parser.add_argument("--partialOrderTreeBuilding", action="store_true")
    parser.add_argument("--neighborJoining", action="store_true")
    parser.add_argument("--alignmentDistance", action="store_true")
    #POA default is 0.9, which leaves most sequences un-bundled
    parser.add_argument("--heaviestBundlingThreshold", type=float, default=0.5)

    parser.add_argument("--alignmentDistanceThreshold", type=float, default=0.7)

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
