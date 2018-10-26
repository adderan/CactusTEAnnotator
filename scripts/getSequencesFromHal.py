import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hal", type=str)
    parser.add_argument("gff", type=str)
    parser.add_argument("genome", type=str)
    args = parser.parse_args()

    with open(args.gff, "r") as gffRead:
        for line in gffRead:
            info = line.split()
            start = info[3]
            end = info[4]
            chrom = info[0]
            name = info[11]
            name = name[1:]
            name = name[:-2]

            seq = subprocess.check_output(["hal2fasta", args.hal, args.genome, "--sequence", chrom, "--start", start, "--length", str(int(end) - int(start))])
            seq = seq.split("\n")
            seq = seq[1:]
            seq = "\n".join(seq)
            print(">%s" % name)
            print(seq)

if __name__ == "__main__":
    main()
