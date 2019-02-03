#include <iostream>
#include <getopt.h>
#include <vector>
#include <stdlib.h>

extern "C" {
#include "sonLib.h"
#include "commonC.h"
#include "bioioC.h"
}


using namespace std;

int main(int argc, char **argv) {
	char *rmaskGffFilename;
	char *featuresGffFilename;
    int i;
    while (1) {
        static struct option long_options[] = { 
            { "--rmask", required_argument, 0, 'a' }, 
			{ "--features", required_argument, 0, 'b'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                rmaskGffFilename = stString_copy(optarg);
                break;
			case 'b':
				featuresGffFilename = stString_copy(optarg);
				break;
            default:
                return 1;
        }
    }

	FILE *featuresGffFile = fopen(featuresGffFilename, "r");
	char *line;
	char chrom[50];
	char annotationType[100];
	char name[100];
	int start, end;
	int score;
	char strand;
	char b;
	int group;
	while ((line = getline(featuresGffFile)) != NULL) {
		sscanf(line, "%s %s %s %d %d %d %c %c %d\n", chrom, annotationType, name, &start, &end, &score, &strand, &b, &group);
		printf("Name = %s\n", name);
	}
}

