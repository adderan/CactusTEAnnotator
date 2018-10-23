#include <stdint.h>

double **getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes);
double **getDistancesExact(char **seqs, int numSeqs, int kmerLength);
double exactJaccardDistance(char *a, char*b, int kmerLength);

