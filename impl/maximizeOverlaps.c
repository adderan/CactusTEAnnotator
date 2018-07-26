#include <stdio.h>

#include "lpo.h"
#include "seq_util.h"


int main(int argc, char **argv) {

    FILE *lpoFile = fopen(argv[1], "r");
    LPOSequence_T *graph = read_lpo(lpoFile);
    fclose(lpoFile);

    fprintf(stderr, "Length = %d\n", graph->length);

    LPOLetter_T *seq = graph->letter;

    //bool *seenThreads = (bool*)calloc(sizeof(bool), seq->nsource_seq);

    LPOLetterSource_T *source = NULL;
    for (int i = 0; i < graph->length; i++) {
	fprintf(stderr, "Base = %s\n", &seq[i].letter);
	source = &seq[i].source;
	int degree = 0;
	while(source != NULL) {
	    degree++;
	    fprintf(stderr, "Contains sequence %d\n", source->iseq);
	    source = source->more;
	}
	fprintf(stderr, "Node degree %d\n", degree);
    }
}
