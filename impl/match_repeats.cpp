#include <getopt.h>
#include <string>
#include "hal.h"


using namespace hal;

int main(int argc, char **argv) {
	char *halFilename = NULL;
	char *genomeName = NULL;
	int maxLength = 10000;
	int minLength = 100;
	char *fastaFilename = NULL;
	char *gffFilename = NULL;
    while (1) {
        static struct option long_options[] = {
            { "hal", required_argument, 0, 'a' }, 
			{ "genome", required_argument, 0, 'b'},
			{ "maxLength", required_argument, 0, 'c'},
			{ "minLength", required_argument, 0, 'd'},
			{ "outFasta", required_argument, 0, 'e'},
			{ "outGFF", required_argument, 0, 'f'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:f:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
               	halFilename = strdup(optarg);
                break;
			case 'b':
				genomeName = strdup(optarg);
				break;
			case 'c':
				sscanf(optarg, "%d", &maxLength);
				break;
			case 'd':
				sscanf(optarg, "%d", &minLength);
				break;
			case 'e':
				fastaFilename = strdup(optarg);
				break;
			case 'f':
				gffFilename = strdup(optarg);
				break;
            default:
                return 1;
        }
    }

	CLParserPtr parser = hdf5CLParserInstance(false);
	AlignmentConstPtr alignment = openHalAlignmentReadOnly(halFilename, parser);

	const Genome *genome = alignment->openGenome(genomeName);
	if (genome == NULL) {
		fprintf(stderr, "Can't open genome.");
		exit(1);
	}
	if (genome->getParent() == NULL) {
		fprintf(stderr, "Specified genome has no parent. Can't find candidate transposon insertions.");
		exit(1);
	}

	FILE *fastaFile = NULL;
	FILE *gffFile = stdout;
	if (fastaFilename) {
		fastaFile = fopen(fastaFilename, "w");
	}

	if (gffFilename) {
		gffFile = fopen(gffFilename, "w");
	}


	TopSegmentIteratorConstPtr topSeg = genome->getTopSegmentIterator();
	TopSegmentIteratorConstPtr endSeg = genome->getTopSegmentEndIterator();

	int i = 0;
	for (; !topSeg->equals(endSeg); topSeg->toRight()) {
		if (topSeg->hasParent()) continue;
		if (topSeg->getLength() > maxLength) continue;
		if (topSeg->getLength() < minLength) continue;

		const Sequence *sequence = topSeg->getSequence();
		int start = topSeg->getStartPosition() - sequence->getStartPosition();
		int end = topSeg->getEndPosition() - sequence->getStartPosition();

		//Print the GFF line
		fprintf(gffFile, "%s\tcandidate_transposon\t%d\t%d\t%d\t0\t+\t.\t%d\n", sequence->getName().c_str(), i, start, end, i);


		if (fastaFile) {
			std::string seqBuffer;
			sequence->getSubString(seqBuffer, start, topSeg->getLength());
			fprintf(fastaFile, ">%s\n", i);
			fprintf(fastaFile, "%s\n", seqBuffer.c_str());
		}

		i++;
	}

}
