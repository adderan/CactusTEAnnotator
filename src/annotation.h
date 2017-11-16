#ifndef _ANNOTATION_H
#define _ANNOTATION_H

void buildClusters(AlignmentConstPtr alignment, RepeatAnnotatorOpts opts);
set<string> getSeeds(string seq, RepeatAnnotatorOpts opts);
double getDistance(string seq1, string seq2, RepeatAnnotatorOpts opts);

#endif
  
