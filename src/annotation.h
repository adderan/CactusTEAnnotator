#ifndef _ANNOTATION_H
#define _ANNOTATION_H
#include "insertions.h"


class GenomeIterator {
private:
  hal::AlignmentConstPtr alignment;
  const hal::Genome *root;
  std::stack<const hal::Genome *> visited;
public:
  GenomeIterator(hal::AlignmentConstPtr _alignment);
  const hal::Genome *next();
};

double kmerDistance(string a, string b);
double **buildDistanceMatrix(vector<string> seqs, int kmerLength);
vector<Sequence*> annotateRepeatsOnBranch(const hal::Genome *reference, InsertionIterator &insertionIt);

#endif
  
