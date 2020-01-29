from __future__ import print_function
import argparse
import sys
import os
import matplotlib.pyplot as pyplot


class GffLine:
    def __init__(self, chrom, start, end, name, family):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.family = family
        self.match = None

def readGff(filename):
    annotations = []
    with open(filename, "r") as file:
        for line in file:
            if len(line.split()) != 15:
                continue
            a, b, c, d, chrom, start, end, score, strand, family, e, f, g, h, id = line.split()
            gffLine = GffLine(chrom=chrom, start=start, end=end, name=id, family=family)
            annotations.append(gffLine)
    return annotations


def overlap(reference, query):
    assert reference.start < reference.end
    assert query.start < query.end
    if reference.start > query.end or query.start > reference.end:
        return 0.0
    overlapStart = reference.start if reference.start > query.start else query.start
    overlapEnd = reference.end if reference.end < query.end else query.end

    overlappingBases = overlapEnd - overlapStart + 1
    overlapFraction = float(overlappingBases)/float(reference.end - reference.start)
    return overlapFraction
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--queryGff", type=str)
    parser.add_argument("--referenceGff", type=str)
    parser.add_argument("--out", type=str)
    parser.add_argument("--skipConcordance", action="store_true", default=False)

    args = parser.parse_args()

    if os.path.exists(args.out):
        print("Output folder %s already exists." % args.out, file=sys.stderr)
        exit(1)
    args.out = os.path.abspath(args.out)
    os.mkdir(args.out)

    queryAnnotations = readGff(args.queryGff)
    referenceAnnotations = readGff(args.referenceGff)

    queryAnnotations.sort(key = lambda x: (x.start, x.end))
    referenceAnnotations.sort(key = lambda x: (x.start, x.end))

    matchingAnnotations = []
    referenceIter = iter(referenceAnnotations)
    queryIter = iter(queryAnnotations)
    reference = next(referenceIter)
    query = next(queryIter)
    while reference and query:
        if overlap(query, reference) > 0.8 and overlap(reference, query) > 0.8:
            matchingAnnotations.append((query, reference))
            query.match = reference
            reference.match = query
            reference = next(referenceIter, None)
            query = next(queryIter, None)
        elif reference.start < query.start:
            reference = next(referenceIter, None)
        else:
            query = next(queryIter, None)

    print("Found %d overlapping annotations" % len(matchingAnnotations), file=sys.stderr)

    if not args.skipConcordance:
        nConcordant = 0
        total = 0
        for i in range(len(matchingAnnotations)):
            if i % 50 != 0:
                continue
            reference_i, query_i = matchingAnnotations[i]
            for j in range(i):
                reference_j, query_j = matchingAnnotations[j]
                a = reference_i.family == reference_j.family
                b = query_i.family == query_j.family
                if ((a and b) or ((not a) and (not b))):
                    nConcordant = nConcordant + 1
                total = total + 1
        concordance = float(nConcordant)/float(total)
        print("Concordance = %f" % concordance, file=sys.stderr)


    referenceFamilies = {}
    for reference in referenceAnnotations:
        if reference.family not in referenceFamilies:
            referenceFamilies[reference.family] = []
        referenceFamilies[reference.family].append(reference)
    familyStats = []
    for family in referenceFamilies:
        coveredRepeatCopies = 0
        for repeatCopy in referenceFamilies[family]:
            if repeatCopy.match is not None:
                coveredRepeatCopies = coveredRepeatCopies + 1
        familyCoverage = float(coveredRepeatCopies)/float(len(referenceFamilies[family]))

        averageLength = float(sum([copy.end - copy.start for copy in referenceFamilies[family]]))/float(len(referenceFamilies[family]))
        concordance = 0
        
        familyStats.append((family, len(referenceFamilies[family]), coveredRepeatCopies, familyCoverage, averageLength))
    familyStats.sort(key = lambda x: x[1], reverse=True)
    familyCoverageStats = os.path.join(args.out, "family_coverage.txt")
    with open(familyCoverageStats, "w") as coverageStatsFile:
        for familyInfo in familyStats:
            print("\t".join([str(item) for item in familyInfo]), file=coverageStatsFile)
    
    #Plot of number of repeat copies per family and how
    #many are covered by the query annotations
    nFamilies = 20
    familyNames = [stats[0] for stats in familyStats]
    familySize = [stats[1] for stats in familyStats]
    coveredRepeatCopies = [stats[2] for stats in familyStats]
    pyplot.bar(familyNames[:nFamilies], familySize[:nFamilies], label="RepeatMasker")
    pyplot.bar(familyNames[:nFamilies], coveredRepeatCopies[:nFamilies], label="TE Annotations from Cactus")
    pyplot.ylabel("Number of copies")
    pyplot.xlabel('RepBase Family')
    pyplot.title("Coverage of Cactus TE annotations on RepeatMasker families")
    pyplot.xticks(rotation=45)
    pyplot.legend()
    pyplot.savefig(os.path.join(args.out, "family_coverage.png"))
    pyplot.clf()


if __name__ == "__main__":
    main()