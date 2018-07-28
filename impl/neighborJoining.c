#include <stdlib.h>
#include <stdio.h>
#include "sonLib.h"

int main(int argc, char **argv) {
  FILE *distancesFile = fopen(argv[1], "r");
  int size = atoi(argv[2]);

  stMatrix *distanceMatrix = stMatrix_construct(size, size);

  char *line = NULL;
  int i;
  int j;
  double dist;
  size_t len;

  while ((getline(&line, &len, distancesFile)) != -1) {
    sscanf(line, "%d %d %lf\n", &i, &j, &dist);
    *stMatrix_getCell(distanceMatrix, i, j) = dist;
  }

  stTree *tree = stPhylogeny_neighborJoin(distanceMatrix, NULL);

  char nodeName[4];
  for (int i = 0; i < size; i++) {
    sprintf(nodeName, "%d", i);
    stTree *leaf_i = stTree_findChild(tree, nodeName);
    stTree *parent = stTree_getParent(leaf_i);
    if (parent != NULL) {
      printf("%d %s\n", i, stTree_getLabel(parent));
    }
  }

}
