import subprocess
import argparse
import random
import os
import shutil

from sonLib.bioio import fastaRead, fastaWrite, catFiles, reverseComplement,getTempFile
from toil.lib.bioio import logger

import CactusTEAnnotator.treeBuilding as treeBuilding

def catFiles(fileList, target):
    with open(target, 'w') as targetWrite:
        for f in fileList:
            with open(f, 'r') as read:
                lines = read.readlines()
                targetWrite.write("".join(lines))
                targetWrite.write("\n")

class Element:
    def __init__(self, chrom, start, end, family, name, seq, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.family = family
        self.name = name
        self.seq = seq
        self.strand = strand

def makeURL(path):
    return "file://" + path

def getFastaSequence(hal, chrom, start, end, args):
    fastaLines = subprocess.check_output(["hal2fasta", hal, args.reference, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]).strip().split("\n")
    return fastaLines

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

def getInsertions(hal, args):
    """Extract candidate TE insertions from the Cactus alignment.
        """
    insertions = []
    i = 0
    for bedLine in subprocess.check_output(["halAlignedExtract", "--complement", hal, args.reference]).split("\n"):
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

        seq = "\n".join(fastaLines[1:])
        forward = Element(chrom=chrom, start=start, end=end, family=None, name=str(i), seq=seq, strand="+")
        backward = Element(chrom=chrom, start=start, end=end, family=None, name="%d_comp" % i, seq=reverseComplement(seq), strand="-")

        forward.reverseName = backward.name
        backward.reverseName = forward.name
        insertions.append(forward)
        insertions.append(backward)

        i = i + 1
    logger.info("Found %d insertions on branch" % len(insertions))
    return insertions

def minhashClustering(elements, args):
    seqs = getTempFile(rootDir = args.workDir)
    with open(seqs, "w") as seqsWrite:
        for element in elements:
            seqsWrite.write(">%s\n" % element.name)
            seqsWrite.write("%s\n" % element.seq)


    partitioning = set()
    for line in subprocess.check_output(["pairwise_distances", "--kmerLength", "5", "--confidenceLevel", str(args.minhashConfidenceLevel), "--sequences", seqs]).split("\n"):
        cluster = frozenset(line.split())
        partitioning.add(cluster)
    return frozenset(partitioning)

def writeGFF(elements, args):
    """Write a set of TE annotations in GFF format.
    """
    gff = getTempFile(rootDir = args.workDir)
    with open(gff, 'w') as gffWrite:
        for element in elements:
            gffWrite.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (element.chrom, "cactus_repeat_annotation", element.name, element.start, element.end, 0, element.strand, '.', element.family))
    return gff

####Partial order alignment###########

def runPoa(elements, heaviestBundle, familyNumber, args):
    logger.info("Aligning family with %d elements" % len(elements))
    logger.info("Elements: %s" % [element.name for element in elements])
    seqs = os.path.join(args.workDir, "%d.fa" % int(familyNumber))
    with open(seqs, "w") as seqsWrite:
        for element in elements:
            seqsWrite.write(">%s\n" % element.name)
            seqsWrite.write("%s\n" % element.seq)

    graph = os.path.join(args.workDir, "%d.po" % int(familyNumber))
    cmd = ["poa", "-read_fasta", seqs, "-po", graph, args.substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    subprocess.check_call(cmd)

    return graph

def updateElements(elements, partitioning):
    nameToElement = {element.name:element for element in elements}
    updatedElements = []
    logger.info("Partitioning = %s" % partitioning)
    logger.info("element names = %s" % nameToElement.keys())
    assert(sum([len(partition) for partition in partitioning]) == len(nameToElement))
    for i, partition in enumerate(partitioning):
        for elementName in partition:
            element = nameToElement[elementName]
            if element.family is None:
                element.family = str(i)
            else:
                element.family = element.family + "_%d" % i
            updatedElements.append(element)
    return updatedElements

def flatten(listOfLists):
    flattenedList = []
    for l in listOfLists:
        flattenedList.extend(l)
    return flattenedList

def runTreeBuilding(graphFile, args):
    """Runs the partitioning tree building method on a sequence graph in 
    PO format.
    """
    graph = treeBuilding.POGraph(graphFile)
    partitions = graph.getPartitions()
    logger.info("Parsed partitions %s" % partitions)
    logger.info("tmp file: %s" % os.path.dirname(graphFile))
    tree = treeBuilding.buildTree(graph.threads, partitions)
    partitioning = treeBuilding.getLeafPartitioning(tree)

    return partitioning
    
def getAlignmentDistances(graph, args):
    distances = getTempFile(rootDir = args.workDir)
    with open(distances, 'w') as distancesWrite:
        subprocess.check_call(["getAlignmentDistances", graph], stdout=distancesWrite)
    return distances


def runNeighborJoining(elements, graph, args):
    if len(elements) < 3:
        return frozenset([frozenset([element.name for element in elements])])

    seqs = job.fileStore.getLocalTempFile()
    seqFiles = [job.fileStore.readGlobalFile(element.seqID) for element in elements]
    graph = getTempFile(rootDir = args.workDir)
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

def buildSubfamilies(elements, args):
    families = {}
    for element in elements:
        if element.family not in families:
            families[element.family] = []
        families[element.family].append(element)

    updatedElements = []
    for family in families:
        elementsInFamily = families[family]
        graph = runPoa(elements=elementsInFamily, heaviestBundle=False, familyNumber=family, args=args)

        alignmentDistances= getAlignmentDistances(graph=graph, args=args)
        partitioning = joinByDistance(distances=alignmentDistances, elements=elementsInFamily, threshold=args.alignmentDistanceThreshold, args=args)
        elementsInFamily = updateElements(elements = elementsInFamily, partitioning = partitioning)


        updatedElements.extend(elementsInFamily)

    return updatedElements

def addRepeatAnnotatorOptions(parser):
    parser.add_argument("reference", type=str)
    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)

    parser.add_argument("--kmerLength", type=int, default=10)
    parser.add_argument("--minhashConfidenceLevel", type=float, default=0.05)
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
    parser.add_argument("--substMatrix", type=str, default = os.path.join(getRootPath(), "blosum80.mat"))

def main():
    parser = argparse.ArgumentParser()
    addRepeatAnnotatorOptions(parser)
    parser.add_argument("hal", type=str)
    parser.add_argument("outGFF", type=str)
    parser.add_argument("workDir", type=str)

    args = parser.parse_args()

    if os.path.exists(args.workDir):
        print("Work dir already exists, exiting")
        exit()
    else:
        os.mkdir(args.workDir)

    insertions = getInsertions(hal=args.hal, args=args)
    partitioning = minhashClustering(elements = insertions, args = args)
    elements = updateElements(elements = insertions, partitioning = partitioning)

    familiesGFF = writeGFF(elements = elements, args = args)
    shutil.copyfile(familiesGFF, os.path.join(args.workDir, "clusters.gff"))

    #elements = buildSubfamilies(elements = elements, args = args)

    gff = writeGFF(elements = elements, args = args)
    shutil.copyfile(gff, args.outGFF)

if __name__ == "__main__":
    main()
