#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <getopt.h>

#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
#include "MurmurHash3.h"


using namespace std;
//using namespace boost::numeric::ublas;

uint32_t hashKmer(char *seq, int length, int seed) {
  for (int i = 0; i < length; i++) {
    if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G') return -1;
  }
  char data[8];
  MurmurHash3_x86_32(seq, length, seed, data);
  return *(uint32_t*)data;
}

uint32_t minhash(char *seq, int kmerLength, int seed) {
  uint32_t min = INT_MAX;
  //cerr << "seq = " << string(seq) << endl;
  for (int i = 0; i < strlen(seq) - kmerLength; i++) {
    uint32_t hash = hashKmer(seq + i, kmerLength, seed);
    //cerr << "hash = " << hash << endl;
    if (hash < min) {
      min = hash;
    }
  }
  return min;
}

uint32_t *precompute_minhash(struct List *seqs, int kmerLength, int seed) {
  uint32_t *hashes = (uint32_t*)malloc(sizeof(uint32_t)* seqs->length); 
  for (int i = 0; i < seqs->length; i++) {
    hashes[i] = minhash((char*)seqs->list[i], kmerLength, seed);
  }
  return hashes;
}


double minhash_jaccard(int a, int b, uint32_t **minhashValues, int numSeeds) {
  int matches = 0;
  for (int i = 0; i < numSeeds; i++) {
    if (minhashValues[i][a] == minhashValues[i][b]) {
      matches++;
    }
  }
  return (double) matches/(double) numSeeds;
}

int main(int argc, char **argv) {
  char *sequencesFilename = NULL;
  int kmerLength = 10;

  int seeds[10] = {0, 10, 50, 89, 10000, 45, 8, 20, 56, 22};
  int numSeeds = 10;

  /*
   * Parse the options.
   */
  int i;
  while (1) {
    static struct option long_options[] = { 
      { "sequences", required_argument, 0, 'a' }, 
      { "kmerLength", required_argument, 0, 'b' }, 
      { 0, 0, 0, 0 } };

    int option_index = 0;

    int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

    if (key == -1) {
      break;
    }

    switch (key) {
      case 'a':
        sequencesFilename = stString_copy(optarg);
        break;
      case 'b':
        i = sscanf(optarg, "%d", &kmerLength);
        assert(i == 1);
        break;
      default:
        return 1;
    }
  }

  struct List *sequences = constructEmptyList(0, NULL);
  struct List *seqNames = constructEmptyList(0, free);
  struct List *seqLengths = constructEmptyList(0, free);
  FILE *sequencesFile = fopen(sequencesFilename, "r");
  fastaRead(sequencesFile, sequences, seqLengths, seqNames);
  cerr << "Read " << seqLengths->length << " sequences" << endl;

  //convert to uppercase
  for (int i = 0; i < seqLengths->length; i++) {
    char *seq = (char*)sequences->list[i];
    for (int k = 0; k < strlen(seq); k++) {
      seq[k] = (char)toupper(seq[k]);
    }
  }

  uint32_t **minhashValues = (uint32_t**)malloc(sizeof(char*)*numSeeds);
  //precompute hash values
  for (int i = 0; i < numSeeds; i++) {
    minhashValues[i] = precompute_minhash(sequences, kmerLength, seeds[i]);
  }

  for (int i = 0; i < seqLengths->length; i++) {
    for (int j = 0; j < i; j++) {
      //cerr << "seq = " << string((char*)sequences->list[i]) << endl;
      double dist = minhash_jaccard(i, j, minhashValues, numSeeds);
      cout << i << " " << j << " " << dist << endl;
    }
  }

}
