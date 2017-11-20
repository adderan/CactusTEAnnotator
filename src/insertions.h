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
  string repeatFamily;
  int group;
  char strand;
  int score;

  Insertion() {};
  
  //Create annotation string
  void toGFF(ostream* gffStream);

  //Metric for similarity between insertions
  double distance(Insertion *other);
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
  bool filter(string seq);
public:
  InsertionIterator() {};
  string toGFF();
  void goToGenome(const Genome *genome);
  InsertionIterator(double _maxNFraction,
				       hal_size_t _insertionJoinDistance,
		    hal_size_t _minInsertionSize, hal_size_t _maxInsertionSize): maxNFraction(_maxNFraction), insertionJoinDistance(_insertionJoinDistance), minInsertionSize(_minInsertionSize), maxInsertionSize(_maxInsertionSize) {}
  Insertion* next();
  Insertion* nextGappedInsertion();
};

#endif
