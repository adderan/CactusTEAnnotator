#include <stdio.h>

#include "lpo.h"
#include "seq_util.h"

#include "sonLib.h"

int h(int x) {
	return (x * 2654435761) % 2 << 32;
}

int hash_partition(stSet *partition) {
	int partitionHash = 0;
	stSetIterator *it = stSet_getIterator(partition);
	stSet *threadSet;
	while((threadSet = stSet_getNext(it)) != NULL) {
		int threadSetHash = 0;
		stSetIterator *threadIt = stSet_getIterator(threadSet);
		int thread;
		while((thread = stSet_getNext(threadSet)) != NULL) {
			threadSetHash += h(thread);
		}
		stSetIterator_destruct(threadIt);
		partitionHash += h(threadSetHash);
	}
	stSetIterator_destruct(it);
	return partitionHash;
}


int main(int argc, char **argv) {
	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;

	stHash *partitionMultiplicity = stHash_construct();
	stHash *partitions = stHash_construct();

	int *containsPosition = (int*)calloc(sizeof(int), graph->nsource_seq);
	for (int i = 0; i < graph->length; i++) {
		LPOSource_t *source = &seq[i].source;

		memset(containsPosition, 0, nsource_seq*sizeof(int));
		do {
			containsPosition[source->iseq] = source->ipos+1;
		}
		while (source = source->more);


		stSet *partition = stSet_construct2(stSet_destruct);

		//Whether each sequence is present at this position in the group

		for (right = &seq[i].right; right && right->ipos>=0; right=right->more) {
			//Find sequences shared in node i and right adjacent node
			stSet *threadSet = stSet_construct();
			for(LPOSource_T *rightSource = &right.source; rightSource != NULL; rightSource = rightSource->more) {
				if (containsPosition[rightSource->iseq] == rightSource->ipos) {
					stSet_insert(threadSet, rightSource->iseq);
				}
			}
			if (stSet_size(threadSet) > 0) {
				stSet_insert(partition, threadSet);
			}
			else {
				stSet_destruct(threadSet);
			}
		}

		int partitionHash = hash_partition(partition);
		if (!stHash_search(partitionMultiplicity, (void*)partitionHash)) {
			int *multiplicity = (int*)malloc(sizeof(int));
			multiplicity[0] = 0;
			stHash_insert(partitionMultiplicity, (void*)partitionHash, multiplicity);
			stHash_insert(partitions, partition);
		}
		int *multiplicity = stHash_search(partitionMultiplicity, (void*)partitionHash);
		*multiplicity++;

	}
	//convert to a sorted set with tuples of (partition_hash, multiplicity)
	stHashIterator *it = stHash_getIterator(multiplicity);



}
