#include <stdio.h>

#include "lpo.h"
#include "seq_util.h"

int main(int argc, char **argv) {
    FILE *lpoFile = fopen(argv[1], "r");
    LPOSequence_T *graph = read_lpo(lpoFile);
    fclose(lpoFile);


    LPOLetter_T *seq = graph->letter;
    
    int **nAligned = (int**)calloc(sizeof(int*), graph->nsource_seq);
    for (int i = 0; i < graph->nsource_seq; i++) {
	nAligned[i] = (int*)calloc(sizeof(int), graph->nsource_seq);
    }

    //count all aligned positions
    for (int i = 0; i < graph->length; i++) {
	for (LPOLetterSource_T *seq_a = &seq[i].source; seq_a != NULL; seq_a = seq_a->more) {
	    for (LPOLetterSource_T *seq_b = &seq[i].source; seq_b != seq_a; seq_b = seq_b->more) {
		int a = seq_a->iseq;
		int b = seq_b->iseq;
		if (a < b) {
		    int temp = a;
		    a = b;
		    b = temp;
		}
		nAligned[a][b]++;
	    }
	}
    }

    //output all pairwise jukes cantor distances
    for (int i = 0; i < graph->nsource_seq; i++) {
	for (int j = 0; j < i; j++) {
	    int N = nAligned[i][j];
	    int len_i = graph->source_seq[i].length;
	    int len_j = graph->source_seq[j].length;
	    double p = 2*N/(double)(len_i + len_j);

	    printf("%s %s %f\n", graph->source_seq[i].name, graph->source_seq[j].name, p);
	}
    }
}
	    
	    
	
	
    
