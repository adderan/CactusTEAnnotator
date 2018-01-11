#ifndef _ANNOTATION_H
#define _ANNOTATION_H

#include <stack>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "hal.h"

using namespace hal;
using namespace std;


class CRASequence {
public:
  hal_size_t start;
  hal_size_t end;
  string seqName;
  char *seq;
  string repeatFamily;
  int group;
  char strand;
  int score;

  CRASequence() {};
  //Create annotation string
  void toGFF(ostream* gffStream);

  //Metric for similarity between insertions
  double distance(CRASequence *other);
};



class InsertionIterator {
private:
  double maxNFraction;
  hal_size_t insertionJoinDistance;
  hal_size_t minInsertionSize;
  hal_size_t maxInsertionSize;
  TopSegmentIteratorConstPtr topSeg;
  TopSegmentIteratorConstPtr endSeg;
  const Genome *genome;
  bool filter(char *seq);
public:
  InsertionIterator() {};
  string toGFF();
  void goToGenome(const Genome *genome);
  InsertionIterator(double _maxNFraction,
		  hal_size_t _insertionJoinDistance,
		  hal_size_t _minInsertionSize, hal_size_t _maxInsertionSize): maxNFraction(_maxNFraction), insertionJoinDistance(_insertionJoinDistance), minInsertionSize(_minInsertionSize), maxInsertionSize(_maxInsertionSize) {}
  CRASequence* next();
  CRASequence* nextGappedInsertion();
};

class GenomeIterator {
private:
  AlignmentConstPtr alignment;
  const Genome *root;
  stack<const Genome *> visited;
public:
  GenomeIterator(AlignmentConstPtr _alignment);
  const Genome *next();
};


template <typename Object> map<Object*, vector<Object*> > buildTransitiveClusters(vector<Object*> objects, boost::numeric::ublas::mapped_matrix<double> similarityMatrix, double similarityThreshold) {

  map<Object *,vector<Object *> > clusterToObj;
  map<Object *,Object *> objToCluster;

  for (uint i = 0; i < objects.size(); i++) {
    objToCluster[objects[i]] = objects[i];
    clusterToObj[objects[i]].push_back(objects[i]);
  }
  for (uint i = 0; i < objects.size(); i++) {
    for (uint j = 0; j < i; j++) {
      Object *a = objects[i];
      Object *b = objects[j];
      double similarity = similarityMatrix (i, j);
      if (similarity > similarityThreshold) {
        //Combine the clusters
        Object* cluster_a = objToCluster[a];
        Object* cluster_b = objToCluster[b];
        if (cluster_a != cluster_b) {
          Object* new_cluster = (cluster_a > cluster_b) ? cluster_a : cluster_b;
          Object* cluster_to_delete = (cluster_a > cluster_b) ? cluster_b : cluster_a;
          vector<Object*> objectsInCluster = clusterToObj[cluster_to_delete];
          typename vector<Object*>::iterator it;
          for (it = objectsInCluster.begin(); it != objectsInCluster.end(); it++) {
            objToCluster[*it] = new_cluster;
            clusterToObj[new_cluster].push_back(*it);
          }
          clusterToObj.erase(cluster_to_delete);
        }
      }
    }
  }
  return clusterToObj;
};


void getInsertionLengthsOnBranch(const Genome* genome, InsertionIterator &insertionIt);
vector<CRASequence*> annotateRepeatsOnBranch(const Genome *genome, InsertionIterator &insertionIter, hal_size_t maxInsertions);

boost::numeric::ublas::mapped_matrix<double> buildDistanceMatrix(vector<char*> seqs, int kmerLength);

#endif
