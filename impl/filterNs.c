#include <stdio.h>
#include "sonLib.h"
#include "bioioC.h"

double nContentThreshold;

static void checkNContent(const char *name, const char *seq, int64_t seqLength) {
	int numN = 0;
	for (int64_t i = 0; i < seqLength; i++) {
		if (seq[i] == 'N') numN++;
	}
	double nFraction = (double)numN / (double) seqLength;
	
	if (nFraction < nContentThreshold) {
		printf(">%s\n", name);
		printf("%s\n", seq);
	}
}

int main(int argc, char **argv) {

	FILE *fasta = fopen(argv[1], "r");

	sscanf(argv[2], "%lf", &nContentThreshold);

	fastaReadToFunction(fasta, checkNContent);
	fclose(fasta);
}
