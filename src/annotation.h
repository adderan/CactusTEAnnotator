#ifndef _ANNOTATION_H
#define _ANNOTATION_H
#include "insertions.h"


class GenomeIterator {
private:
  AlignmentConstPtr alignment;
  const Genome *root;
  std::stack<const Genome *> visited;
public:
  GenomeIterator(AlignmentConstPtr _alignment);
  const Genome *next();
};

double insertionDistance(Insertion *a, Insertion *seq2);
double kmerDistance(string a, string b);
double **buildDistanceMatrix(vector<string> seqs, int kmerLength);
vector<Insertion*> annotateInsertionsOnBranch(const Genome *reference, InsertionIterator &insertionIt);

#endif
  
