#include <stdio.h>
#include <stack>
#include "hal.h"
#include "omp.h"
#include "insertions.h"

using namespace hal;
using namespace std;


GenomeIterator::GenomeIterator(AlignmentConstPtr _alignment) {
  alignment = _alignment;
  root = alignment->openGenome(alignment->getRootName());
  visited.push(root);
}

const Genome * GenomeIterator::next() {
  if (visited.empty()) return NULL;

  root = visited.top();
  visited.pop();

  for (hal_size_t childIndex = 0; childIndex < root->getNumChildren(); childIndex++) {
    const Genome* child = root->getChild(childIndex);
    visited.push(child);
  }
  return root;
  
}

InsertionIterator::InsertionIterator(const Genome *_genome, RepeatAnnotatorOpts _opts) {

  genome = _genome;

  topSeg = genome->getTopSegmentIterator();
  endSeg = genome->getTopSegmentEndIterator();
  opts = _opts;
}


std::string InsertionIterator::next() {
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent() && topSeg->getLength() > opts.minInsertionSize) {
      string seq;
      topSeg->getString(seq);
      topSeg->toRight();
      return seq;
    }
    topSeg->toRight();
  }
  return "";
}
      



SmoothedInsertionIterator::SmoothedInsertionIterator(const Genome *_genome, RepeatAnnotatorOpts _opts) {

  genome = _genome;

  topSeg = genome->getTopSegmentIterator();
  endSeg = genome->getTopSegmentEndIterator();
  opts = _opts;
}


std::string SmoothedInsertionIterator::next()

{

  vector<string> insertion;
  hal_size_t gapLength = 0;

  while (topSeg->equals(endSeg) == false) {
    string seq;
    topSeg->getString(seq);

    if (!topSeg->hasParent()) {
      insertion.push_back(seq);
      gapLength = 0;
    }
    else if (seq.length() + gapLength < opts.insertionJoinDistance) {
      gapLength += seq.length();
      insertion.push_back(seq);
    }

    topSeg->toRight();

    if (insertion.size() > 0) {
      if (insertion.size() > 1) cerr << "Joining " << insertion.size()
				       << " segments to form an insertion." << endl;
      string insertionSeq;
      for (int i = 0; i < insertion.size(); i++) {
	insertionSeq += insertion.at(i);
      }
      if (insertionSeq.length() > opts.minInsertionSize) {
	return insertionSeq;
      }
    }
    insertion.clear();
    gapLength = 0;

      
  }
  return "";
  
}

