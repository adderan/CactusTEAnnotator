#include <stdio.h>
#include <stdlib.h>
#include "lpo.h"


int is_consensus(char *seqName) {
	return !strncmp(seqName, "CONSENS", 6);
}

int main(int argc, char **argv) {

	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;

	for (int i = 0; i < graph->nsource_seq; i++) {
		if (!is_consensus(graph->source_seq[i].name)) {
			printf("%s %d\n", graph->source_seq[i].name, graph->source_seq[i].bundle_id);
		}
		else {
			printf("%s %d\n", graph->source_seq[i].name, graph->source_seq[i].length);
		}
	}
}
