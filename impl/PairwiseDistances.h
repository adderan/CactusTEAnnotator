#include <stdint.h>


uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, uint32_t *seeds, int numSeeds);

double minhash_jaccard(int a, int b, uint32_t **minhashValues, int numSeeds);
