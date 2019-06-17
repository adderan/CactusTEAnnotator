#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"

int main(int argc, char **argv) {
	FILE *alignmentsFile = fopen(argv[1], "r");

	struct PairwiseAlignment *pairwiseAlignment = NULL;

	stUnionFind *families = stUnionFind_construct();
	stSet *alignments = stSet_construct();
	stHash *seqNames = stHash_construct();

	//stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	while((pairwiseAlignment = cigarRead(alignmentsFile)) != NULL) {
		stSet_insert(alignments, pairwiseAlignment);
		stUnionFind_add(families, (void*) stHash_stringKey(pairwiseAlignment->contig1));
		stUnionFind_add(families, (void*) stHash_stringKey(pairwiseAlignment->contig2));

		stHash_insert(seqNames, (void*) stHash_stringKey(pairwiseAlignment->contig1), (void*) pairwiseAlignment->contig1);
		stHash_insert(seqNames, (void*) stHash_stringKey(pairwiseAlignment->contig2), (void*) pairwiseAlignment->contig2);


	}

	stSetIterator *alignmentsIt = stSet_getIterator(alignments);
	while((pairwiseAlignment = stSet_getNext(alignmentsIt)) != NULL) {
		stUnionFind_union(families, (void*) stHash_stringKey(pairwiseAlignment->contig1), (void*) stHash_stringKey(pairwiseAlignment->contig2));
	}
	fclose(alignmentsFile);

	stUnionFindIt *it = stUnionFind_getIterator(families);
	stSet *family;
	while ((family = stUnionFindIt_getNext(it)) != NULL) {

		stSetIterator *setIt = stSet_getIterator(family);
		void *item;
		while ((item = stSet_getNext(setIt)) != NULL) {
			printf("%s ", (char*) stHash_search(seqNames, item));
		}
		printf("\n");
		stSet_destructIterator(setIt);
	}
	stUnionFind_destructIterator(it);
}
