import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rmaskGff", type=str)

    args = parser.parse_args()

    with open(args.rmaskGff, "r") as rmaskGffFile:
        for line in rmaskGffFile:
            info = line.split()
            chrom = info[4]
            start = info[5]
            end = info[6]
            print("%s\t%s\t%s\n" % (chrom, start, end))

if __name__ == "__main__":
    main()