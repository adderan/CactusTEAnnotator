import sys
import argparse

"""Compute the Adjusted Rand Index between two sets of features.
"""

class Feature:
    def __init__(self, chrom, start, end, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name

class Element:
    def __init__(self, name, refName):
        self.name = name
        self.refName = refName

def readGff(gff):
    features = {}
    for line in gff:
        info = line.split()
        if len(info) != 9:
            continue
        chrom, source, annotationType, start, end, score, strand, a, name = info
        if start > end:
            continue          
        if not chrom in features:
            features[chrom] = []
        features[chrom].append(Feature(chrom=chrom, start=int(start), end=int(end), name=name))
    return(features)

def readRmaskGff(gff):
    features = {}
    for line in gff:
        info = line.split()
        if len(info) != 16:
            continue
        chrom, source, annotationType, start, end, score, strand, a, b, name, c, d, e, f, g, h = info
        if start > end:
            continue
        if not chrom in features:
            features[chrom] = []
        features[chrom].append(Feature(chrom, int(start), int(end), name))
    return(features)

def overlap(f1, f2):
    assert f1.start < f1.end
    assert f2.start < f2.end
    assert f1.chrom == f2.chrom

    if (f1.start <= f2.end) and (f2.start <= f1.end):
        #Overlap
        return 0
    elif f1.start < f2.start:
        #No overlap and f1 occurs before f2
        return -1
    else:
        #No overlap and f1 occurs after f2
        return 1

def findIntersectionEfficient(features, refFeatures):
    elements = []
    for chrom in features:
        if chrom not in refFeatures:
            #print("Skipping chromosome %s" % chrom)
            continue
        #sort by start position
        features[chrom].sort(key=lambda x: x.start)
        refFeatures[chrom].sort(key=lambda x: x.start)

        i = 0
        for feature in features[chrom]:
            if i >= len(refFeatures[chrom]):
                break

            while overlap(feature, refFeatures[chrom][i]) > 0:
                i = i + 1
                if i >= len(refFeatures[chrom]):
                    break
            if i >= len(refFeatures[chrom]):
                break

            if overlap(feature, refFeatures[chrom][i]) == 0:
                #print "Found overlap betwen " + feature.name + " and " + refFeatures[chrom][i].name
                elements.append(Element(name=feature.name, refName=refFeatures[chrom][i].name))
                #print("Found overlapping feature %d %d and %s" % (len(elements), i, refFeatures[chrom][i].name))
                #print("%d %d" % (feature.start, feature.end))
                #print("%d %d" % (refFeatures[chrom][i].start, refFeatures[chrom][i].end))

    return(elements)

def findIntersectionSlow(features, refFeatures):
    elements = []
    for chrom in features:
        if chrom not in refFeatures:
            continue
        for i in range(len(features[chrom])):
            for j in range(len(refFeatures[chrom])):
                if overlap(features[chrom][i], refFeatures[chrom][j]) == 0:
                    elements.append(Element(name=features[chrom][i].name, refName=refFeatures[chrom][j].name))
                    #print("Found overlapping feature %d %d and %s" % (len(elements), i, refFeatures[chrom][j].name))
                    #print("%d %d" % (features[chrom][i].start, features[chrom][i].end))
                    #print("%d %d" % (refFeatures[chrom][j].start, refFeatures[chrom][j].end))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--rmask", type=file)
    parser.add_argument("--features", type=file)
    args = parser.parse_args()

    
    features = readGff(args.features)

    refFeatures = readRmaskGff(args.rmask)
    
    #print("Found %d chromosomes" % len(features))
    #print("Found %d chromosomes in reference" % len(refFeatures))
    
    elements = findIntersectionEfficient(features, refFeatures)

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
            if (elements[i].name == elements[j].name) and (elements[i].refName == elements[j].refName):
                a += 1.0
            elif (elements[i].name != elements[j].name) and (elements[i].refName != elements[j].refName):
                b += 1.0
            elif (elements[i].name == elements[j].name) and (elements[i].refName != elements[j].refName):
                c += 1.0
            elif (elements[i].name != elements[j].name) and (elements[i].refName == elements[j].refName):
                d += 1.0

    #print("Found %d elements" % len(elements))

    #print("a = %f" % a)
    #print("b = %f" % b)
    #print("c = %f" % c)
    #print("d = %f" % d)

    r = (a + b)/(a + b + c + d)
    overCollapsed = float(c)/(a + b + c + d)
    underCollapsed = float(d)/(a + b + c + d)

    clusters = set([feature.name for feature in features.values()])
    print("Number of clusters = %d" % len(clusters))
    print("Rand index = %f" % r)
    print("Overcollapsed = %f" % overCollapsed)
    print("Undercollapsed = %f" % underCollapsed)


if __name__ == "__main__":
    main()
