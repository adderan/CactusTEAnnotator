#include <stdio.h>
#include <stack>
#include "hal.h"
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

double fractionN(string seq) {
  double numN = 0.0;
  for (unsigned int i = 0; i < seq.length(); i++) {
    if (seq[i] == 'N') {
      numN += 1.0;
    }
  }
  return numN/seq.length();
}

string Insertion::toGFF() {
  return "";
}

void InsertionIterator::goToGenome(const Genome * _genome) {
  genome = _genome;
  topSeg = genome->getTopSegmentIterator();
  endSeg = genome->getTopSegmentEndIterator();
}

bool InsertionIterator::filter(string seq) {
  if (seq.length() < minInsertionSize) return false;
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
