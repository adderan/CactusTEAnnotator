#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <iostream>
#include <climits>

#include "PairwiseDistances.h"

using namespace std;

class Point {
    public:
        double x,y;
        Point(double _x, double _y) : x(_x), y(_y) {}

};


static double point_distance(Point *a, Point *b) {
    double dx = (a->x - b->x);
    double dy = (a->y - b->y);
    return sqrt(dx*dx + dy*dy);
}


/*
   void clusterBuilderTest() {
   cerr << "Running clustering tests" << endl;

   vector<Point*> points;
   points.push_back(new Point(0,0));
   points.push_back(new Point(1,0));
   points.push_back(new Point(3,3));
   points.push_back(new Point(3,4));
   points.push_back(new Point(10,10));
   points.push_back(new Point(9,10));
   points.push_back(new Point(3,5));
   mapped_matrix *distanceMatrix = new double*[points.size()];
   for (uint i = 0; i < points.size(); i++) {
   distanceMatrix[i] = (double*) calloc(i, sizeof(double));
   for (uint j = 0; j < i; j++) {
   distanceMatrix[i][j] = point_distance(points[i], points[j]);
   }
   }
   map<Point*, vector<Point*> > clusters = buildTransitiveClusters<Point>(points, distanceMatrix, 1.1);
   cerr << "Built " << clusters.size() << " clusters" << endl;
   assert(clusters.size() == 3);
   }
   */

void insertionIteratorTest() {
}

double exactKmerDistance(char *a, char*b, int kmerLength) {
    int nMatches = 0;
    for (int i = 0; i < strlen(a) - kmerLength; i++) {
        for (int j = 0; j < strlen(b) - kmerLength; j++) {
            if (strncmp(a + i, b + j, kmerLength) == 0) {
                nMatches++;
                break;
            }
        }
    }
    int nKmers = strlen(a) + strlen(b) - 2*kmerLength;
    return 2*(float)nMatches/nKmers;
}

void kmerDistanceTest() {
    int len = 1000;
    int kmerLength = 4;
    char *a = (char*)malloc(len*sizeof(char));
    char *b = (char*)malloc(len*sizeof(char));

    strcpy(a, "AGAGCATGTGACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");
    strcpy(b, "AGAGCATGCGGCATATCGATCGTAGCACTATGTGTATGTCGATGTGATGCTGATGCTAGCTAGCTAGCTGATCGATGC\0");


    char **seqs = (char**)malloc(2*sizeof(char*));
    seqs[0] = a;
    seqs[1] = b;

    int numSeeds = 100000;
    uint32_t *seeds = (uint32_t*) malloc(numSeeds*sizeof(uint32_t));
    for (int i = 0; i < numSeeds; i++) {
        seeds[i] = (uint32_t) (INT_MAX*rand());
    }

    uint32_t **minhashValues = precompute_minhash(seqs, 2, kmerLength, seeds, numSeeds);

    double dist = minhash_jaccard(0, 1, minhashValues, numSeeds);

    double exactDist = exactKmerDistance(a, b, kmerLength);
    cerr << "Length = " << strlen(a) << endl;
    cerr << "dist = " << dist << endl;
    cerr << "exact dist = " << exactDist << endl;

}


int main(int argc, char **argv) {
    kmerDistanceTest();
}
