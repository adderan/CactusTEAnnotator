#include <stdio.h>

#include "lpo.h"
#include "seq_util.h"


int main(int argc, char **argv) {
	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;


	int *containsPosition = (int*)calloc(sizeof(int), graph->nsource_seq);
	for (int i = 0; i < graph->length; i++) {
		LPOLetterSource_T *source = &seq[i].source;

		memset(containsPosition, 0, graph->nsource_seq*sizeof(int));
		do {
			containsPosition[source->iseq] = source->ipos+1;
		}
		while (source = source->more);


		//Whether each sequence is present at this position in the group

		for (LPOLetterLink_T *right = &seq[i].right; right && right->ipos>=0; right=right->more) {
			//Find sequences shared in node i and right adjacent node
			for(LPOLetterSource_T *rightSource = &seq[right->ipos].source; rightSource != NULL; rightSource = rightSource->more) {
				if (containsPosition[rightSource->iseq] == rightSource->ipos) {
					printf("%d ", rightSource->iseq);
				}
			}
			printf(", ");

		}
		printf("\n");


	}

}
