import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rmaskGff", type=str)

    args = parser.parse_args()

    with open(args.rmaskGff, "r") as rmaskGffFile:
        for line in rmaskGffFile:
            chrom = line[4]
            start = line[5]
            end = line[6]
            print "%s\t%s\t%s\n" % (chrom, start, end)

if __name__ == "__main__":
    main()