#include <stdio.h>
#include "hal.h"
#include "halRepeats.h"

using namespace hal;
using namespace std;

HalRepeats::HalRepeats()
{
}

HalRepeats::~HalRepeats()
{
}


void HalRepeats::analyzeBranch(AlignmentConstPtr alignment,
			       hal_size_t gapThreshold,
			       double nThreshold)
{
  
  const Genome *reference = alignment->openGenome("Human");
  assert(reference != NULL);
  if (reference->getParent() == NULL)
  {
    throw hal_exception("Reference genome must have parent");
  }

  TopSegmentIteratorConstPtr top  = reference->getTopSegmentIterator();
  top->toSite(0);


  hal_index_t end = reference->getSequenceLength();
  RearrangementPtr rearrangement = reference->getRearrangement(top->getArrayIndex(),
                                               gapThreshold, nThreshold);

  const Sequence *sequence;
  do {
    sequence = reference->getSequenceBySite( rearrangement->getLeftBreakpoint()->getStartPosition());
    if (rearrangement->getID() == Rearrangement::Insertion) {
      cout << rearrangement->getLeftBreakpoint()->getStartPosition() << " " <<
	rearrangement->getRightBreakpoint()->getEndPosition() + 1 << endl;

      TopSegmentIteratorConstPtr segment = rearrangement->getLeftBreakpoint()->copy();
      //segment->toSite(rearrangement->getLeftBreakpoint()->getStartPosition());
      while (true) {
	string seq;
	segment->getString(seq);
	cout << seq << endl;
	if (segment == rearrangement->getRightBreakpoint()) {
	  break;
	}
	segment->toRight();
      }
    }
  }
  while (rearrangement->identifyNext() == true &&
	 rearrangement->getLeftBreakpoint()->getStartPosition() <= end);
}

