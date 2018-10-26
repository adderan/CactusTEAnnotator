#include <stdint.h>
#include <vector>
#include <set>

using namespace std;

vector<tuple<int, int, double> > getDistances(char **seqs, int numSeqs, int kmerLength, int numHashes);
vector<tuple<int, int, double> > getDistancesExact(char **seqs, int numSeqs, int kmerLength);
double exactJaccardDistance(char *a, char*b, int kmerLength);
set<set<long> > buildClusters(vector<tuple<int, int, double> > &pValues, int nSeqs, double confidenceLevel);
