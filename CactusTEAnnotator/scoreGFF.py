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

def readRmaskGff(gffFilename):
    features = []
    with open(gffFilename, 'r') as gffFile:
        for line in gffFile:
            info = line.split()
            if len(info) != 16:
                continue
            chrom, source, annotationType, start, end, score, strand, a, b, geneID, c, transcriptID, e, familyID, g, h = info
            geneID = geneID[1:len(geneID) - 2]
            familyID = familyID[1:len(familyID) - 2]
            transcriptID = transcriptID[1:len(transcriptID) - 2]
            if start > end:
                continue

            #print("family = %s" % family)
            features.append(Feature(chrom=chrom, start=int(start), end=int(end), name=transcriptID, strand=strand, family=geneID))
    return(features)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--query", type=str)
    parser.add_argument("--target", type=file)
    parser.add_argument("--minOverlap", type=int, default=0.6)
    args = parser.parse_args()

    
    queryAnnotations = readRmaskGff(args.query)
    targetAnnotations = readRmaskGff(args.target)
    

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
    print("Fraction of true matchings found: %f" % float(a/(a + d)))
    print("Fraction false matchings: %f" % float(c/(a + c)))
    print("Total matchings: %f" % (a + c))
    print("Total matchings found by reference: %f" % (a + d))


if __name__ == "__main__":
    main()
