import subprocess
import argparse
import random
import os
import shutil
import networkx

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

class RepeatInstance:
    def __init__(self, chrom, start, end, family, name, seq, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.seq = seq
        self.strand = strand

def getSequence(hal, chrom, start, end, args):
    fastaLines = subprocess.check_output(["hal2fasta", hal, args.reference, "--sequence", chrom, "--start", str(start), "--length", str(end - start)]).strip().split("\n")
    sequence = fastaLines[1:]
    return "\n".join(sequence)

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
        sequence = getSequence(hal, chrom, start, end, args)
        assert len(sequence) > 0

        if highNFraction(sequence):
            continue
        if args.maxInsertions > 0 and i > args.maxInsertions:
            break
        forward = RepeatInstance(chrom=chrom, start=start, end=end, family=None, name=str(i), seq=sequence, strand="+")
        backward = RepeatInstance(chrom=chrom, start=start, end=end, family=None, name="%d_comp" % i, seq=reverseComplement(sequence), strand="-")

        forward.reverseName = backward.name
        backward.reverseName = forward.name
        insertions.append(forward)
        insertions.append(backward)

        i = i + 1
    logger.info("Found %d insertions on branch" % len(insertions))
    return insertions

def minhashClustering(elements, args):
    seqs = getTempFile(rootDir = args.outDir)
    with open(seqs, "w") as seqsWrite:
        for element in elements:
            seqsWrite.write(">%s\n" % element.name)
            seqsWrite.write("%s\n" % element.seq)

    families = []
    nameToElement = {element.name:element for element in elements}
    for line in subprocess.check_output(["minhash", "--kmerLength", str(args.kmerLength), "--distanceThreshold", str(args.minhashDistanceThreshold), "--sequences", seqs]).split("\n"):
        if line == "":
            continue
        family = line.split()
        family = [nameToElement[name] for name in family]
        families.append(family)
    return families

def readGFF(gff, hal, args):
    elements = []
    with open(gff, 'r') as gffRead:
        for line in gffRead:
            chrom, annotationType, name, start, end, score, strand, a, group = line.split()
            sequence = getSequence(hal, chrom, int(start), int(end), args)
            elements.append(RepeatInstance(chrom = chrom, start = int(start), end = int(end), family = group, name = name, seq = sequence, strand = strand))
    return elements

def writeGFF(families, args):
    """Write a set of TE annotations in GFF format.
    """
    gff = getTempFile(rootDir = args.outDir)
    with open(gff, 'w') as gffWrite:
        for i, family in enumerate(families):
            for element in family:
                gffWrite.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (element.chrom, "cactus_repeat_annotation", element.name, element.start, element.end, 0, element.strand, '.', i))
    return gff

####Partial order alignment###########

def writeSequencesToFile(repeatInstances, seqsFile):
    with open(seqsFile, "w") as seqsWrite:
        for repeatCopy in repeatInstances:
            seqsWrite.write(">%s\n" % repeatCopy.name)
            seqsWrite.write("%s\n" % repeatCopy.seq)

def runPoa(repeatCopies, heaviestBundle, familyNumber, args):
    seqs = os.path.join(args.outDir, "%d.fa" % int(familyNumber))
    with open(seqs, "w") as seqsWrite:
        for repeatCopy in repeatCopies:
            seqsWrite.write(">%s\n" % repeatCopy.name)
            seqsWrite.write("%s\n" % repeatCopy.seq)

    graph = getTempFile(rootDir = args.outDir)
    cmd = ["poa", "-read_fasta", seqs, "-po", graph, args.substMatrix]
    if heaviestBundle:
        cmd.extend(["-hb", "-hbmin", str(args.heaviestBundlingThreshold)])
    subprocess.check_call(cmd)

    return graph

def clusterByHeaviestBundling(repeatCopies, familyNumber, args):
    graph = runPoa(repeatCopies, True, familyNumber, args)

    nameToRepeatInstance = {repeatInstance.name:repeatInstance for repeatInstance in repeatCopies}
    subfamilies = []
    for line in subprocess.check_output(["getHeaviestBundles", graph]).split("\n"):
        if line == "":
            continue
        seqsInBundle = [nameToRepeatInstance[name] for name in line.split() if name in nameToRepeatInstance]
        subfamilies.append(seqsInBundle)

    return subfamilies

def clusterByAlignmentDistances(repeatCopies, familyNumber, args):
    graph = runPoa(repeatCopies, False, familyNumber, args)
    subfamilies = []
    nameToRepeatCopy = {repeatCopy.name:repeatCopy for repeatCopy in repeatCopies}
    for line in subprocess.check_output(["clusterByAlignmentDistances", graph, str(args.alignmentDistanceThreshold)]).split("\n"):
        line = line.strip().rstrip()
        if line == "":
            continue
        subfamilyNames = line.split()
        subfamilies.append([nameToRepeatCopy[repeatCopyName] for repeatCopyName in subfamilyNames])
        
    return subfamilies

def joinSubfamiliesByMinhashDistance(repeatFamilies, args):
    graph = networkx.Graph()
    
    seqFileNames = []
    for i in range(len(repeatFamilies)):
        graph.add_node(i)
        seqFile = getTempFile(rootDir=args.outDir)
        writeSequencesToFile(repeatFamilies[i], seqFile)
        seqFileNames.append(seqFile)

        for j in range(i):
            if len(repeatFamilies[i]) <= 1 or len(repeatFamilies[j]) <= 1:
                continue
            dist = float(subprocess.check_output(["minhash", "--kmerLength", str(args.minhashClusterDistKmerLength), "--clusterA", seqFileNames[i], "--clusterB", seqFileNames[j]]))
            print("Found distance %f between clusters %d and %d of size %d and %d" % (dist, i, j, len(repeatFamilies[i]), len(repeatFamilies[j])))
            if dist < args.minhashClusterDistanceThreshold:
                assert graph.has_node(i)
                assert graph.has_node(j)
                graph.add_edge(i, j)

    newClusters = []
    for component in networkx.connected_components(graph):
        newCluster = []
        for node in component:
            newCluster.extend(repeatFamilies[node])
        newClusters.append(newCluster)
    return newClusters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("reference", type=str)
    parser.add_argument("outDir", type=str)

    parser.add_argument("--minInsertionSize", type=int, default=100)
    parser.add_argument("--maxInsertionSize", type=int, default=50000)
    parser.add_argument("--maxNFraction", type=float, default=0.1)

    parser.add_argument("--kmerLength", type=int, default=5)
    parser.add_argument("--minhashClusterDistKmerLength", type=int, default=15)
    parser.add_argument("--minhashDistanceThreshold", type=float, default=0.3)
    parser.add_argument("--minClusterSize", type=int, default=1)

    parser.add_argument("--maxInsertions", type=int, default=0)

    #Methods for splitting families of insertions
    #POA default is 0.9, which leaves most sequences un-bundled
    parser.add_argument("--subfamilySplitMethod", type=str, default="alignmentDistance")
    parser.add_argument("--heaviestBundlingThreshold", type=float, default=0.9)

    parser.add_argument("--alignmentDistanceThreshold", type=float, default=0.1)
    parser.add_argument("--substMatrix", type=str, default = os.path.join(getRootPath(), "blosum80.mat"))
    parser.add_argument("--inGFF", type=str, default = None)
    parser.add_argument("--skipSubfamilies", action = "store_true")
    parser.add_argument("--minhashClusterDistanceThreshold", type=float, default=0.1)

    args = parser.parse_args()

    if os.path.exists(args.outDir):
        print("Work dir already exists, exiting")
        exit()
    else:
        os.mkdir(args.outDir)


    subfamilySplitFn = None
    if args.subfamilySplitMethod == "alignmentDistance":
        subfamilySplitFn = clusterByAlignmentDistance
    elif args.subfamilySplitMethod == "heaviestBundling":
        subfamilySplitFn = clusterByHeaviestBundling

    insertions = getInsertions(hal=args.hal, args=args)
    coarseFamilies = minhashClustering(elements = insertions, args = args)


    #Store an initial GFF from the families created by minhash clustering
    #coarseFamiliesDict = {str(i):family for i, family in enumerate(coarseFamilies)}
    familiesGff = writeGFF(families = coarseFamilies, args = args)
    shutil.copyfile(familiesGff, os.path.join(args.outDir, "families.gff"))

    families = []
    for i, family in enumerate(coarseFamilies):
        subfamilies = subfamilySplitFn(repeatCopies = family, familyNumber = i, args = args)
        families.extend(subfamilies)
    subfamiliesGFF = writeGFF(families = families, args = args)
    shutil.copyfile(subfamiliesGFF, os.path.join(args.outDir, "subfamilies.gff"))

    families = joinSubfamiliesByMinhashDistance(families, args)

    gff = writeGFF(families = families, args = args)
    shutil.copyfile(gff, os.path.join(args.outDir, "final.gff"))

if __name__ == "__main__":
    main()
