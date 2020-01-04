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


def overlap(target, query):
    assert target.start < target.end
    assert query.start < query.end
    if target.start > query.end or query.start > target.end:
        return 0.0
    overlapStart = target.start if target.start > query.start else query.start
    overlapEnd = target.end if target.end < query.end else query.end

    overlappingBases = overlapEnd - overlapStart + 1
    overlapFraction = float(overlappingBases)/float(target.end - target.start)
    return overlapFraction

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--queryGff", type=str)
    parser.add_argument("--targetGff", type=str)
    parser.add_argument("--out", type=str)

    args = parser.parse_args()

    if os.path.exists(args.out):
        print("Output folder %s already exists.", file=sys.stderr)

    queryAnnotations = readGff(args.queryGff)
    targetAnnotations = readGff(args.targetGff)

    queryAnnotations.sort(key = lambda x: (x.start, x.end))
    targetAnnotations.sort(key = lambda x: (x.start, x.end))

    matchingAnnotations = []
    targetIter = iter(targetAnnotations)
    queryIter = iter(queryAnnotations)
    target = next(targetIter)
    query = next(queryIter)
    while target and query:
        if overlap(query, target) > 0.5:
            matchingAnnotations.append((query, target))
            query.match = target
            target.match = query
            target = next(targetIter, None)
            query = next(queryIter, None)
        elif target.start < query.start:
            target = next(targetIter, None)
        else:
            query = next(queryIter, None)

    print("Found %d overlapping annotations" % len(matchingAnnotations), file=sys.stderr)

    nConcordant = 0
    total = 0
    for i in range(len(matchingAnnotations)):
        if i % 50 != 0:
            continue
        target_i, query_i = matchingAnnotations[i]
        for j in range(i):
            target_j, query_j = matchingAnnotations[j]
            a = target_i.family == target_j.family
            b = query_i.family == query_j.family
            if ((a and b) or ((not a) and (not b))):
                nConcordant = nConcordant + 1
            total = total + 1
    concordance = float(nConcordant)/float(total)
    print("Concordance = %f" % concordance, file=sys.stderr)


    if args.familyCoverageStats:
        targetFamilies = {}
        for target in targetAnnotations:
            if target.family not in targetFamilies:
                targetFamilies[target.family] = []
            targetFamilies[target.family].append(target)
        familyStats = []
        for family in targetFamilies:
            coveredRepeatCopies = 0
            for repeatCopy in targetFamilies[family]:
                if repeatCopy.match is not None:
                    coveredRepeatCopies = coveredRepeatCopies + 1
            familyCoverage = float(coveredRepeatCopies)/float(len(targetFamilies[family]))

            averageLength = float(sum([copy.end - copy.start for copy in targetFamilies[family]]))/float(len(targetFamilies[family]))
            concordance = 0
            
            familyStats.append((family, len(targetFamilies[family]), familyCoverage, averageLength))
        familyStats.sort(key = lambda x: x[1], reverse=True)
        with open(args.familyCoverageStats, "w") as coverageStatsFile:
            for familyInfo in familyStats:
                print("\t".join([str(item) for item in familyInfo]), file=coverageStatsFile)
if __name__ == "__main__":
    main()