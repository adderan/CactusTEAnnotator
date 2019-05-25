#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <iostream>
#include <getopt.h>

extern "C" {
#include "lpo.h"
}

using namespace std;

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
		while(source = source->more);
	}

}

int is_consensus(char *seqName) {
	return !strncmp(seqName, "CONSENS", 6);
}

int main(int argc, char **argv) {
	char *lpoFilename = NULL;
	bool printElements = false;
    int i;
    while (1) {
        static struct option long_options[] = {
            { "lpo", required_argument, 0, 'a' }, 
			{ "printElements", no_argument, 0, 'b'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                lpoFilename = strdup(optarg);
                break;
            case 'b':
				printElements = true;
                break;
            default:
                return 1;
        }
    }

	FILE *lpoFile = fopen(lpoFilename, "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;

	if (printElements) {
		//just print the consensus sequences
		for (int i = 0; i < graph->nsource_seq; i++) {
			if (is_consensus(graph->source_seq[i].name)) {
				printf(">%s\n", graph->source_seq[i].name);
				printSequence(graph, i, stdout);
				printf("\n");
			}
		}
		exit(0);
	}

	map<int, vector<string> > families;
	for (int i = 0; i < graph->nsource_seq; i++) {
		if (!is_consensus(graph->source_seq[i].name)) {
			families[graph->source_seq[i].bundle_id].push_back(graph->source_seq[i].name);
		}
	}

	for (map<int, vector<string> >::iterator it = families.begin(); it != families.end(); ++it) {
		for (auto &seq: it->second) {
			cout << seq << " ";
		}
		cout << endl;
	}
}
