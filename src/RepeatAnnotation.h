#ifndef _ANNOTATION_H
#define _ANNOTATION_H

#include <stack>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "hal.h"

using namespace hal;
using namespace std;


class Seq {
public:
  hal_size_t start;
  hal_size_t end;
  string seqName;
  char *seq;
  string repeatFamily;
  int group;
  char strand;
  int score;

  Seq() {};
  //Create annotation string
  void toGFF(ostream* gffStream);

  //Metric for similarity between insertions
  double distance(Seq *other);
};



class InsertionIterator {
private:
  double maxNFraction;
  hal_size_t insertionJoinDistance;
  hal_size_t minInsertionSize;
  hal_size_t maxInsertionSize;
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  bool filter(char *seq);
public:
  InsertionIterator() {};
  string toGFF();
  void goToGenome(const Genome *genome);
  InsertionIterator(double _maxNFraction,
		  hal_size_t _insertionJoinDistance,
		  hal_size_t _minInsertionSize, hal_size_t _maxInsertionSize): maxNFraction(_maxNFraction), insertionJoinDistance(_insertionJoinDistance), minInsertionSize(_minInsertionSize), maxInsertionSize(_maxInsertionSize) {}
  Seq* next();
  Seq* nextGappedInsertion();
};

class GenomeIterator {
private:
  AlignmentConstPtr alignment;
  const Genome *root;
  stack<const Genome *> visited;
public:
  GenomeIterator(AlignmentConstPtr _alignment);
  const Genome *next();
};


void getInsertionLengthsOnBranch(const Genome* genome, InsertionIterator &insertionIt);
vector<Seq*> annotateRepeatsOnBranch(const Genome *genome, InsertionIterator &insertionIter, hal_size_t maxInsertions);

boost::numeric::ublas::mapped_matrix<double> buildDistanceMatrix(vector<Seq*> &seqs, int kmerLength);
map<Seq*, vector<Seq*> > buildGroups(vector<Seq*> &seqs, boost::numeric::ublas::mapped_matrix<double> &distanceMatrix, double similarityThreshold);

#endif
