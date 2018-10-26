import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff")
    parser.add_argument("--fasta")
    
    args = parser.parse_args()

    with open(args.gff, 'r') as gffRead:
        for line in gffRead:
            info = line.split()

            chrom = info[0]
            start = info[3]
            end = info[4]
            name = info[11]
            name = name[1:]
            name = name[:-2]

            seq = subprocess.check_output(["samtools", "faidx", args.fasta, "%s:%d-%d" % (chrom, int(start), int(end))])
            seq = seq.split("\n")
            seq = seq[1:]
            seq = "\n".join(seq)
            print(">%s" % name)
            print(seq)
if __name__ == "__main__":
    main()

