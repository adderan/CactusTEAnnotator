import argparse
import random

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--numSeqs", type=int, default=10)
    parser.add_argument("--seqLength", type=int, default=100)
    args = parser.parse_args()

    bases = ['a', 'c', 't', 'g', 'A', 'C', 'T', 'G']
    for i in range(args.numSeqs):
        seq = [random.choice(bases) for j in range(args.seqLength)]
        print(">seq_%d" % i)
        print("".join(seq))

if __name__ == "__main__":
    main()
