#ifndef _ANNOTATION_H
#define _ANNOTATION_H
#include "insertions.h"

double insertionDistance(Insertion *a, Insertion *seq2);
double kmerDistance(string a, string b);
double **buildDistanceMatrix(vector<string> seqs, int kmerLength);

#endif
  
