#include <stdio.h>
#include <stdlib.h>
#include "RepeatAnnotator.h"


int main(int argc, char **argv) {
	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	printf("length = %d\n", graph->length);
	return(0);
}



