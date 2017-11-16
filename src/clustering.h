#ifndef _CLUSTERING_H
#define _CLUSTERING_H

#include <map>
#include <iterator>
#include <vector>
#include <stdlib.h>

#include "omp.h"

using namespace std;

template <typename Object> map<Object*, vector<Object*> > buildTransitiveClusters(vector<Object*> objects, double distanceThreshold) {

  map<Object *,vector<Object *> > clusterToObj;
  map<Object *,Object *> objToCluster;

  for (uint i = 0; i < objects.size(); i++) {
    objToCluster[objects[i]] = objects[i];
    clusterToObj[objects[i]].push_back(objects[i]);
  }
  
#pragma omp parallel shared(objects, objToCluster, clusterToObj) num_threads(8)
  {
    #pragma omp parallel for
    for (uint i = 0; i < objects.size(); i++) {
      for (uint j = 0; j < i; j++) {
	Object *a = objects[i];
	Object *b = objects[j];
	double distance = a->distance(b);
	if (distance < distanceThreshold) {
	  //Combine the clusters
#pragma omp critical
	  {
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
    }
  }
  return clusterToObj;
}

#endif
