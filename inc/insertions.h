#ifndef _REPEATS_H
#define _REPEATS_H

#include <stack>
#include "hal.h"

using namespace hal;

typedef struct RepeatAnnotatorOpts {
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
} RepeatAnnotatorOpts;

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
  InsertionIterator(const Genome *_genome, RepeatAnnotatorOpts _opts);
  std::string next();
};

class SmoothedInsertionIterator {
private:
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  RepeatAnnotatorOpts opts;
public:
  SmoothedInsertionIterator(const Genome *_genome, RepeatAnnotatorOpts _opts);
  std::string next();
};
#endif
