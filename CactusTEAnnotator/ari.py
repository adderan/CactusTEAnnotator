import sys
import argparse

"""Compute the Adjusted Rand Index between two sets of features.
"""

class Feature:
    def __init__(self, chrom, start, end, name, strand, family):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
        self.family = family

class Element:
    def __init__(self, name, refName, family, refFamily):
        self.name = name
        self.refName = refName
        self.family = family
        self.refFamily = refFamily

def readGff(gff):
    features = {}
    for line in gff:
        info = line.split()
        if len(info) != 9:
            continue
        chrom, source, name, start, end, score, strand, a, family = info
        if start > end:
            continue          
        if not chrom in features:
            features[chrom] = []
        features[chrom].append(Feature(chrom=chrom, start=int(start), end=int(end), name=name, strand=strand, family=family))
    return(features)

def readRmaskGff(gff):
    features = {}
    for line in gff:
        info = line.split()
        if len(info) != 16:
            continue
        chrom, source, annotationType, start, end, score, strand, a, b, geneID, c, transcriptID, e, familyID, g, h = info
        geneID = geneID[1:len(geneID) - 2]
        familyID = familyID[1:len(familyID) - 2]
        transcriptID = transcriptID[1:len(transcriptID) - 2]
        if start > end:
            continue
        if not chrom in features:
            features[chrom] = []

        #print("family = %s" % family)
        features[chrom].append(Feature(chrom=chrom, start=int(start), end=int(end), name=transcriptID, strand=strand, family=geneID))
    return(features)


def findIntersectionEfficient(features, refFeatures, args):
    elements = []
    for chrom in features:
        if chrom not in refFeatures:
            #print("Skipping chromosome %s" % chrom)
            continue
        #sort by start position
        features[chrom].sort(key=lambda x: x.start)
        refFeatures[chrom].sort(key=lambda x: x.start)

        i = 0
        j = 0
        print("Finding intersection")
        while i < len(features[chrom]) and j < len(refFeatures[chrom]):
            f1 = features[chrom][i]
            f2 = refFeatures[chrom][j]

            #move to beginning of first overlap
            while f1.end < f2.start and (i + 1 < len(features[chrom])):
                i = i + 1
                print("Skipping feature %s" % f1.name)
                f1 = features[chrom][i]
            while f2.end < f1.start and j + 1 < (len(refFeatures[chrom])):
                j = j + 1
                print("Skipping feature %s" % f2.name)
                f2 = refFeatures[chrom][j]

            if (f1.start <= f2.end) and (f2.start <= f1.end) and f1.strand == f2.strand:
                #overlap
                overlap_fraction_f2 = (min(f1.end, f2.end) - max(f1.start, f2.start))/float(f2.end - f2.start)
                overlap_fraction_f1 = (min(f1.end, f2.end) - max(f1.start, f2.start))/float(f1.end - f1.start)
                if overlap_fraction_f1 > args.minOverlap and overlap_fraction_f2 > args.minOverlap:
                    print("Found overlapping features %s %s" % (f1.name, f2.name))
                    elements.append(Element(name=f1.name, refName=f2.name, family=f1.family, refFamily=f2.family))
                    i = i + 1
                    j = j + 1
                    continue
            if f1.start < f2.start:
                i = i + 1
            else:
                j = j + 1

    return(elements)

def findIntersectionSlow(features, refFeatures):
    elements = []
    for chrom in features:
        if chrom not in refFeatures:
            continue
        for i in range(len(features[chrom])):
            for j in range(len(refFeatures[chrom])):
                if overlap(features[chrom][i], refFeatures[chrom][j]) == 0:
                    elements.append(Element(name=features[chrom][i].name, refName=refFeatures[chrom][j].name, family=features[chrom][i].family, refFamily=refFeatures[chrom][i].family))
                    #print("Found overlapping feature %d %d and %s" % (len(elements), i, refFeatures[chrom][j].name))
                    #print("%d %d" % (features[chrom][i].start, features[chrom][i].end))
                    #print("%d %d" % (refFeatures[chrom][j].start, refFeatures[chrom][j].end))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--rmask", type=file)
    parser.add_argument("--features", type=file)
    parser.add_argument("--minOverlap", type=int, default=0.6)
    args = parser.parse_args()

    
    features = readGff(args.features)

    refFeatures = readRmaskGff(args.rmask)
    
    #print("Found %d chromosomes" % len(features))
    #print("Found %d chromosomes in reference" % len(refFeatures))
    
    elements = findIntersectionEfficient(features, refFeatures, args)

    for element in elements:
        print("%s\t%s\t%s\t%s\n" % (element.family, element.refFamily, element.refName, element.name))

    #print("Found %d elements" % len(elements))
    #Calculate rand index
    a = 0.0
    b = 0.0
    c = 0.0
    d = 0.0

    k = 0
    for i in range(len(elements)):
        for j in range(i):
            if k % 1000 != 0:
                continue
            if (elements[i].family == elements[j].family) and (elements[i].refFamily == elements[j].refFamily):
                a += 1.0
            elif (elements[i].family != elements[j].family) and (elements[i].refFamily != elements[j].refFamily):
                b += 1.0
            elif (elements[i].family == elements[j].family) and (elements[i].refFamily != elements[j].refFamily):
                c += 1.0
            elif (elements[i].family != elements[j].family) and (elements[i].refFamily == elements[j].refFamily):
                d += 1.0

    #print("Found %d elements" % len(elements))

    #print("a = %f" % a)
    #print("b = %f" % b)
    #print("c = %f" % c)
    #print("d = %f" % d)

    print("Found %d overlapping elements" % len(elements))
    r = (a + b)/(a + b + c + d)
    overCollapsed = float(c)/(a + b + c + d)
    underCollapsed = float(d)/(a + b + c + d)

    #clusters = set([feature.name for feature in features.values()])
    #print("Number of clusters = %d" % len(clusters))
    print("Rand index = %f" % r)
    print("Overcollapsed = %f" % overCollapsed)
    print("Undercollapsed = %f" % underCollapsed)


if __name__ == "__main__":
    main()
