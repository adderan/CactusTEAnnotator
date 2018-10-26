#ifndef _REPEAT_ANNOTATOR_H
#define _REPEAT_ANNOTATOR_H

extern "C" {
	#include "lpo.h"
	#include "seq_util.h"
}

double **getAlignmentDistances(LPOSequence_T *graph);

#endif
