#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"


int main(int argc, char **argv) {
	FILE *alignmentsFile = fopen(argv[1], "r");

	struct PairwiseAlignment *pairwiseAlignment = NULL;

	stUnionFind *families = stUnionFind_construct();
	while((pairwiseAlignment = cigarRead(alignmentsFile)) != NULL) {
		printf("start coordinate = %ld\n", pairwiseAlignment->start1);
		printf("end coordinate = %ld\n", pairwiseAlignment->end1);
		printf("contig name = %s\n", pairwiseAlignment->contig1);
		if (!stUnionFind_find(families, (void*) pairwiseAlignment->contig1)) {
			stUnionFind_add(families, (void*) pairwiseAlignment->contig1);
		}
		if (!stUnionFind_find(families, (void*) pairwiseAlignment->contig2)) {
			stUnionFind_add(families, (void*) pairwiseAlignment->contig2);
		}
		stUnionFind_union(families, (void*) pairwiseAlignment->contig1, (void*) pairwiseAlignment->contig2);
	}
	fclose(alignmentsFile);

	stUnionFindIt *it = stUnionFind_getIterator(families);
	stSet *family;
	while ((family = stUnionFindIt_getNext(it)) != NULL) {

		stSetIterator *setIt = stSet_getIterator(family);
		char *item;
		while ((item = stSet_getNext(setIt)) != NULL) {
			printf("%s ", item);
		}
		printf("\n");
		stSet_destructIterator(setIt);
	}
	stUnionFind_destructIterator(it);
}
