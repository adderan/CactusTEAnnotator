#include <stdint.h>


uint32_t **precompute_minhash(char **seqs, int numSeqs, int kmerLength, int numHashes);
uint32_t hashKmer(char *seq, int length);

double **getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes);
double **getDistancesExact(char **seqs, int numSeqs, int kmerLength);
double exactJaccardDistance(char *a, char*b, int kmerLength);

