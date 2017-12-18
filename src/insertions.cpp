#include <stdio.h>
#include <stack>
#include "hal.h"
#include "insertions.h"

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

void Sequence::toGFF(ostream* gffStream) {
  *gffStream << seqName << "\tcactus_repeat_annotation\trepeat_copy" << "\t" << start << "\t" << end <<  "\t" << score << "\t" << strand << "\t" << "." << "\t" << group << endl;
}

void InsertionIterator::goToGenome(const hal::Genome * _genome) {
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


Sequence* InsertionIterator::next() {
  if (insertionJoinDistance > 0) {
    return NULL;
  }
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      string seq;
      seq.clear();
      topSeg->getString(seq);
      if (filter(seq)) {
	Sequence *insertion = new Sequence;
	insertion->seq = seq;
	insertion->seqName = topSeg->getSequence()->getName();
	hal_size_t seqStart = topSeg->getSequence()->getStartPosition();
	insertion->start = topSeg->getStartPosition() - seqStart;
	insertion->end = topSeg->getEndPosition() - seqStart;
	insertion->strand = '+';
	insertion->score = 0;
	topSeg->toRight();
	return insertion;
      }
    }
    topSeg->toRight();
  }
  return NULL;
}

/*
Sequence *SequenceIterator::nextGappedSequence()

{
  
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      //Found an insertion. Scan forward for another insertion within a certain distance
      string insertion;
      topSeg->getString(insertion);
      while(topSeg->equals(endSeg) == false) {
	topSeg->toRight();
	if (!topSeg->hasParent) {
	  string seq 
  }
  
}
*/
