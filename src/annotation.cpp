#include <vector>
#include <string>
#include "insertions.h"
#include "annotation.h"
#include "clustering.h"


using namespace std;

void buildClusters(AlignmentConstPtr alignment, RepeatAnnotatorOpts opts) {
  vector<Insertion*> insertions;
  if (opts.referenceName != "") {
    const Genome *reference = alignment->openGenome(opts.referenceName);
    insertions = getInsertionsOnBranch(reference, opts);
  }
}


set<string> getSeeds(string seq, RepeatAnnotatorOpts opts) {
  set<string> seeds;
  if (seq.length() < opts.seedLength) return seeds;
  for (unsigned int i = 0; i < (seq.length() - opts.seedLength); i++) {
    seeds.insert(seq.substr(i, i + opts.seedLength));
  }
  return seeds;
}


double getDistance(string seq1, string seq2, RepeatAnnotatorOpts opts) {
  set<string> seq1Seeds = getSeeds(seq1, opts);
  set<string> seq2Seeds = getSeeds(seq2, opts);
  set<string> intersectionSeeds;
  set<string> unionSeeds;
  set_intersection(seq1Seeds.begin(), seq1Seeds.end(), seq2Seeds.begin(), seq2Seeds.end(), inserter(intersectionSeeds, intersectionSeeds.begin()));
  set_union(seq1Seeds.begin(), seq1Seeds.end(), seq2Seeds.begin(), seq2Seeds.end(), inserter(unionSeeds, unionSeeds.begin()));
  int intersectionSize = distance(intersectionSeeds.begin(), intersectionSeeds.end());
  int unionSize = distance(unionSeeds.begin(), unionSeeds.end());
  return (double)(intersectionSize)/unionSize;
}
  
