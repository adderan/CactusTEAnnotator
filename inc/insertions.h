#ifndef _REPEATS_H
#define _REPEATS_H

#include <stack>
#include "hal.h"

using namespace hal;
using namespace std;



typedef struct RepeatAnnotatorOpts {
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
  CLParserPtr optionsParser;
  bool joinNeighborInsertions;
} RepeatAnnotatorOpts;

void getInsertions(AlignmentConstPtr alignment, RepeatAnnotatorOpts &opts);

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
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  RepeatAnnotatorOpts opts;
public:
  InsertionIterator() {};
  InsertionIterator(const Genome *_genome, RepeatAnnotatorOpts &_opts);
  virtual string next();
};

class InsertionIteratorJoinNeighbors: public InsertionIterator {
private:
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  RepeatAnnotatorOpts opts;
public:
  InsertionIteratorJoinNeighbors(const Genome *_genome, RepeatAnnotatorOpts &_opts);
  string next();
};
#endif
