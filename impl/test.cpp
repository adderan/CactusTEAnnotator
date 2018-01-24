#include <math.h>
#include <vector>
#include <stdlib.h>
#include "RepeatAnnotation.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>

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

/*
void testDistanceMatrix() {
  char *kmer1 = "AAGTCAGTACATAGGGACAGTCAG";
  char *kmer2 = "GTACGTACAATGAGAAGGGAATCA";
  char *a = (char*)malloc(sizeof(char)*(strlen(kmer1)*3 + strlen(kmer2)*3));
  char *b = (char*)malloc(sizeof(char)*(strlen(kmer2)*3));
  char *c = (char*)malloc(sizeof(char)*(strlen(kmer1)*4 + strlen(kmer2)*1));
  sprintf(a, "%s%s%s%s%s%s", kmer1, kmer2, kmer2, kmer1, kmer1, kmer2);
  sprintf(b, "%s%s%s", kmer2, kmer2,  kmer2);
  sprintf(c, "%s%s%s%s%s", kmer1, kmer1, kmer1, kmer1, kmer2);

  vector<char*> seqs;
  seqs.push_back(a);
  seqs.push_back(b);
  seqs.push_back(c);
  boost::numeric::ublas::mapped_matrix<double> dist = buildDistanceMatrix(seqs, 5);
  cerr << "Distance = " << dist(1,0) << endl;
  assert(dist(1,0) > 0);
  assert(dist(1,0) < dist(2,0));
  assert(dist(2,1) > 0);
}

*/
/*

void kmerDistanceTest() {
  string kmer1 = "AAGTCAGTACATAGGGACAGTCAG";
  string kmer2 = "GTACGTACAATGAGAAGGGAATCA";
  string a = kmer1 + kmer2 + kmer2 + kmer1 + kmer1 + kmer2;
  string b = kmer2 + kmer2 + kmer2;
  assert(kmerDistance(a, b) > 0.0);
  assert(kmerDistance(a, b) < 1.0);
  assert(kmerDistance(a, b) == kmerDistance(b, a));
}
*/


int main(int argc, char **argv) {
  //clusterBuilderTest();
  insertionIteratorTest();
  //testDistanceMatrix();
}
