#include <math.h>
#include <vector>
#include "clustering.h"
#include "insertions.h"
#include "annotation.h"

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
  map<Point*, vector<Point*> > clusters = buildTransitiveClusters<Point>(points, &point_distance, 1.1);
  cerr << "Built " << clusters.size() << " clusters" << endl;
  assert(clusters.size() == 3);
}

void insertionIteratorTest() {
}

void testDistanceMatrix() {
  string kmer1 = "AAGTCAGTACATAGGGACAGTCAG";
  string kmer2 = "GTACGTACAATGAGAAGGGAATCA";
  string a = kmer1 + kmer2 + kmer2 + kmer1 + kmer1 + kmer2;
  string b = kmer2 + kmer2 + kmer2;
  string c = kmer1 + kmer1 + kmer1 + kmer1 + kmer2;

  vector<string> seqs;
  seqs.push_back(a);
  seqs.push_back(b);
  seqs.push_back(c);
  double **dist = buildDistanceMatrix(seqs, 5);
  cerr << "Distance = " << dist[1][0] << endl;
  assert(dist[1][0] > 0);
  assert(dist[1][0] < 1.0);
  assert(dist[2][1] > 0);

}

void kmerDistanceTest() {
  string kmer1 = "AAGTCAGTACATAGGGACAGTCAG";
  string kmer2 = "GTACGTACAATGAGAAGGGAATCA";
  string a = kmer1 + kmer2 + kmer2 + kmer1 + kmer1 + kmer2;
  string b = kmer2 + kmer2 + kmer2;
  assert(kmerDistance(a, b) > 0.0);
  assert(kmerDistance(a, b) < 1.0);
  assert(kmerDistance(a, b) == kmerDistance(b, a));
}
  


int main(int argc, char **argv) {
  clusterBuilderTest();
  insertionIteratorTest();
  kmerDistanceTest();
  testDistanceMatrix();
}
