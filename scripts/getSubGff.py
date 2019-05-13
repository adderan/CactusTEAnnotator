import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff")
    parser.add_argument("chrom")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)

    args = parser.parse_args()
    with open(args.gff, "r") as gffFile:
        for line in gffFile:
            chrom, source, annotationType, start, end, score, strand, a, gene_id_, gene_id, transcript_id_, transcript_id, family_id_, family_id, class_id_, class_id = line.split()

            if chrom == args.chrom and (int(start) >= args.start) and (int(end) <= args.end):
                print line

if __name__ == "__main__":
    main()
