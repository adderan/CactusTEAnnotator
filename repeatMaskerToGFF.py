import argparse
import sys


def main():
    parser = argparse.ArgumentParser()

    for line in sys.stdin:
        if not len(line.split()) == 15:
            continue
        a, b, c, d, chromosome, start, end, left, strand, repeatFamily, repeatClass, e, f, g, h= line.split()
        if strand == 'C':
            strand = '-'
        print "%s repeatMasker\t%s\t%s\t%s\t0\t%s\t.\t%s" % (chromosome, repeatClass, start, end, strand, repeatFamily)
        
if __name__ == "__main__":
    main()
