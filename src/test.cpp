#include <math.h>
#include <vector>
#include "clustering.h"
#include "insertions.h"

using namespace std;

class Point {
public:
  double x,y;
  Point(double _x, double _y) : x(_x), y(_y) {}
  double distance(Point *other) {
    double dx = (x - other->x);
    double dy = (y - other->y);
    return sqrt(dx*dx + dy*dy);
  }
};

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
  map<Point*, vector<Point*> > clusters = buildTransitiveClusters<Point>(points, 1.1);
  cerr << "Built " << clusters.size() << " clusters" << endl;
  assert(clusters.size() == 3);
}

void insertionIteratorTest() {
}


int main(int argc, char **argv) {
  clusterBuilderTest();
}
