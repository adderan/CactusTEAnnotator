#include <stdint.h>


uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, int numHashes);
uint32_t hashKmer(char *seq, int length);
double minhash_jaccard(int a, int b, uint32_t **minhashValues, int numSeeds);
