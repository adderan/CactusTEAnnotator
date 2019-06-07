#include <stdint.h>
#include <vector>

using namespace std;

vector<tuple<int, int, double> > getDistances(char **seqs, int numSeqs, int kmerLength, int sketchSize);
vector<tuple<int, int, double> > getDistancesExact(char **seqs, int numSeqs, int kmerLength);
double exactJaccardDistance(char *a, char*b, int kmerLength);
vector<uint32_t> buildSketch(char *seq, int kmerLength);
double minhashJaccard(vector<uint32_t> &a, vector<uint32_t> &b, int length_a, int length_b, int kmerLength, int sketchSize);
double getDistanceBetweenFamilies(char **seqs1, char **seqs2, int numSeqs1, int numSeqs2, int kmerLength, int sketchSize);


