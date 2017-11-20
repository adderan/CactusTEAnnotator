#include <stdio.h>
#include <stack>
#include "hal.h"
#include "insertions.h"

using namespace hal;
using namespace std;


double fractionN(string seq) {
  double numN = 0.0;
  for (unsigned int i = 0; i < seq.length(); i++) {
    if (seq[i] == 'N') {
      numN += 1.0;
    }
  }
  return numN/seq.length();
}

void Insertion::toGFF(ostream* gffStream) {
  *gffStream << seqName << " cactus_repeat_annotation repeat_copy " << start << " " << end <<  " " << score << " " << strand << " " << "." << " " << group << endl;
}

void InsertionIterator::goToGenome(const Genome * _genome) {
  genome = _genome;
  topSeg = genome->getTopSegmentIterator();
  endSeg = genome->getTopSegmentEndIterator();
}

bool InsertionIterator::filter(string seq) {
  if (seq.length() < minInsertionSize || seq.length() > maxInsertionSize) return false;
  if (fractionN(seq) > maxNFraction) {
    return false;
  }
  return true;
}


Insertion* InsertionIterator::next() {
  if (insertionJoinDistance > 0) {
    return NULL;
  }
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      string seq;
      topSeg->getString(seq);
      topSeg->toRight();
      if (filter(seq)) {
	Insertion *insertion = new Insertion;
	insertion->seq = seq;
	insertion->seqName = topSeg->getSequence()->getName();
	insertion->start = topSeg->getStartPosition();
	insertion->end = topSeg->getEndPosition();
	insertion->strand = '+';
	insertion->score = 0;
	return insertion;
      }
    }
    topSeg->toRight();
  }
  return NULL;
}
/*
Insertion *InsertionIterator::nextGappedInsertion()

{

  vector<TopSegmentIteratorConstPtr> insertion;
  hal_size_t gapLength = 0;
  
  while (topSeg->equals(endSeg) == false) {
    string seq;
    topSeg->getString(seq);

    if (!topSeg->hasParent()) {
      insertion.push_back(seq);
      gapLength = 0;
    }
    else if ((seq.length() + gapLength < insertionJoinDistance) && (insertion.size() > 0)) {
      gapLength += seq.length();
      insertion.push_back(seq);
      cerr << "Appending aligned sequence of length " << seq.length() << endl;
    }

    topSeg->toRight();

    if (insertion.size() > 0) {
      string insertionSeq;
      for (int i = 0; i < insertion.size(); i++) {
	insertionSeq += insertion.at(i);
      }
      if (filter(insertionSeq)) {
	Insertion *insertion = new Insertion;
	insertion->seq = insertionSeq;
	insertion->end = topSeg->getEndPosition();
	return insertion;
      }
    }
    insertion.clear();
    gapLength = 0;
      
  }
  return NULL;
  
}

*/
