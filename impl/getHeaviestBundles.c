#include <stdio.h>
#include <getopt.h>
#include "lpo.h"


//Print sequence i in the graph
void printSequence(LPOSequence_T *graph, int seqNum, FILE *output) {
	LPOLetter_T *seq = graph->letter;
	int pos = 0;
	for (int i = 0; i < graph->length; i++) {
		LPOLetterSource_T *source = &seq[i].source;
		do {
			if (source->iseq == seqNum && source->ipos == pos) {
				fprintf(output, "%c", seq[i].letter);
				pos++;
				break;
			}
		}
		while((source = source->more));
	}

}

int is_consensus(char *seqName) {
	return !strncmp(seqName, "CONSENS", 6);
}

int main(int argc, char **argv) {
	char *lpoFilename = NULL;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                lpoFilename = strdup(optarg);
                break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	for (int i = 0; i < graph->nsource_seq; i++) {
		if (is_consensus(graph->source_seq[i].name)) {
			printf(">%s\n", graph->source_seq[i].name);
			printSequence(graph, i, stdout);
			printf("\n");
		}
	}

}
