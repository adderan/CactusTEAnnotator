#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <iostream>

extern "C" {
#include "lpo.h"
}

using namespace std;

int is_consensus(char *seqName) {
	return !strncmp(seqName, "CONSENS", 6);
}

int main(int argc, char **argv) {

	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

	LPOLetter_T *seq = graph->letter;

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
