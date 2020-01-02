import argparse


class GffLine:
    def __init__(self, chrom, start, end, name, family):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.family = family

def readGff(filename):
    with open(filename, "r") as file:
        for line in file:
            if len(line.split()) != 15:
                continue
            a, b, c, d, chrom, start, end, score, strand, family, e, f, g, h, id = line.split()
            gffLine = GffLine(chrom=chrom, start=start, end=end, name=id, family=family)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--queryGff", type=str)
    parser.add_argument("--targetGff", type=str)

    args = parser.parse_args()

    queryAnnotations = readGff(args.queryGff)
    targetAnnotations = readGff(args.targetGff)

    sort(queryAnnotations, key = lambda x: x.start)
    sort(targetAnnotations, key = lambda x: x.start)



if __name__ == "__main__":
    main()