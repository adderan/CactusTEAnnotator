#include "lpo.h"
#include "seq_util.h"
#include <getopt.h>


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

	printf("digraph {\n");

	for (int i = 0; i < graph->length; i++) {
		for (LPOLetterLink_T *right = &graph->letter[i].right; right && right->ipos > 0; right = right->more) {
			printf("\t%s_%d -> %s_%d\n", &graph->letter[i].letter, i, &graph->letter[right->ipos].letter, right->ipos);
		}

	}
	
	printf("}\n");
	

}
