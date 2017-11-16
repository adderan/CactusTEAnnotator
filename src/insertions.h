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
  string toGFF();
};

typedef struct RepeatAnnotatorOpts {
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
  CLParserPtr optionsParser;
  string referenceName;
  double maxNFraction;
  hal_size_t seedLength;
} RepeatAnnotatorOpts;

void getInsertions(AlignmentConstPtr alignment, RepeatAnnotatorOpts opts);
vector<Insertion*> getInsertionsOnBranch(const Genome *genome, RepeatAnnotatorOpts opts);

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
  bool filter(string seq);
public:
  InsertionIterator() {};
  InsertionIterator(const Genome *_genome, RepeatAnnotatorOpts &_opts);
  virtual Insertion* next();
};

class InsertionIteratorJoinNeighbors: public InsertionIterator {
private:
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  RepeatAnnotatorOpts opts;
public:
  InsertionIteratorJoinNeighbors(const Genome *_genome, RepeatAnnotatorOpts &_opts);
  //Insertion* next();
};
#endif
