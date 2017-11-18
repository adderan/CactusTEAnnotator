#ifndef _REPEATS_H
#define _REPEATS_H

#include <stack>
#include <vector>
#include "hal.h"

using namespace hal;
using namespace std;


class Insertion {
public:
  hal_size_t start;
  hal_size_t end;
  string seqName;
  string seq;
  unsigned int group;

  Insertion() {};

  //Create annotation string
  string toGFF();

  //Metric for similarity between insertions
  double distance(Insertion *other);
};

class GenomeIterator {
private:
  AlignmentConstPtr alignment;
  const Genome *root;
  std::stack<const Genome *> visited;
public:
  GenomeIterator(AlignmentConstPtr _alignment);
  const Genome *next();
};

class InsertionIterator {
private:
  double maxNFraction;
  hal_size_t insertionJoinDistance;
  hal_size_t minInsertionSize;
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  bool filter(string seq);
public:
  InsertionIterator() {};
  string toGFF();
  void goToGenome(const Genome *genome);
  InsertionIterator(double _maxNFraction,
				       hal_size_t _insertionJoinDistance,
				       hal_size_t _minInsertionSize): maxNFraction(_maxNFraction), insertionJoinDistance(_insertionJoinDistance), minInsertionSize(_minInsertionSize) {}
  Insertion* next();
  Insertion* nextGappedInsertion();
};

#endif
