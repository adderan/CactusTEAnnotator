#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include "sonLib.h"
#include "bioioC.h"

stUnionFind *clusters;
stHash *seqs;

void addSeq(char *seqName) {
    if (stHash_search(seqs, (void*) stHash_stringKey(seqName))) return;
    char *seqNameTmp = malloc(sizeof(char) *(strlen(seqName) + 1));
    strcpy(seqNameTmp, seqName);
    stHash_insert(seqs, (void*) stHash_stringKey(seqName), (void*) seqNameTmp);
    stUnionFind_add(clusters, (void*) stHash_stringKey(seqName));
}


int main(int argc, char **argv) {
    char *distancesFilename = NULL;
    double distanceThreshold = 0.1;
	while (1) {
		static struct option long_options[] = {
			{ "distances", required_argument, 0, 'a' }, 
            { "distanceThreshold", required_argument, 0, 'b'},
			{ 0, 0, 0, 0 } };

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

		if (key == -1) {
			break;
		}

		switch (key) {
			case 'a':
				distancesFilename = strdup(optarg);
				break;
            case 'b':
                sscanf(optarg, "%lf", &distanceThreshold);
                break;
			default:
				return 1;
		}
	}

    clusters = stUnionFind_construct();
    seqs = stHash_construct2(NULL, free);

    FILE *distancesFile = fopen(distancesFilename, "r");

    char *line;
    char seq1[50];
    char seq2[50];
    double dist;
    size_t nBytes;
    while(getline(&line, &nBytes, distancesFile) != EOF) {
        sscanf(line, "%s %s %lf", seq1, seq2, &dist);
        addSeq(seq1);
        addSeq(seq2);
        if (fabs(dist) < distanceThreshold) {
            stUnionFind_union(clusters, (void*) stHash_stringKey(seq1), (void*) stHash_stringKey(seq2));
        }
    }
    fclose(distancesFile);

    stUnionFindIt *it = stUnionFind_getIterator(clusters);
    stSet *cluster;
    while((cluster = stUnionFindIt_getNext(it)) != NULL) {
        stSetIterator *clusterIt = stSet_getIterator(cluster);
        void *seqID;
        while((seqID = stSet_getNext(clusterIt)) != NULL) {
            char *seqName = stHash_search(seqs, seqID);
            fprintf(stdout, "%s ", seqName);
        }
        fprintf(stdout, "\n");
        stSet_destructIterator(clusterIt);
        
    }
    stUnionFind_destructIterator(it);
    stUnionFind_destruct(clusters);
    stHash_destruct(seqs);
}