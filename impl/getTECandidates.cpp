#include <getopt.h>
#include <stdlib.h>
#include <string>
#include "hal.h"


using namespace hal;
using namespace std;

double getNFraction(DNAIteratorConstPtr it, hal_size_t start, hal_size_t length) {
	int numN = 0;
	hal_size_t genomeLength = it->getGenome()->getSequenceLength();

	it->jumpTo(start);
	for (hal_size_t i = 0; i < length; i++) {
		if (it->getArrayIndex() >= genomeLength) break;
		if (it->getChar() == 'N') {
			numN++;
		}
		it->toRight();
	}

	return (double)numN/(double)length;
}

int main(int argc, char **argv) {
	string halPath;
	string genomeName;
	hal_size_t minLength = 100;
	hal_size_t maxLength = 10000;
	string fastaFilename;
	string gffFilename;
	int maxSequences = 0;

	CLParserPtr parser = hdf5CLParserInstance(false);
	parser->addArgument("halPath", "Input hal file");
	parser->addArgument("genome", "Name of genome to scan");
	parser->addOption("outGFF", "Output GFF representing the candidate TEs", "");
	parser->addOption("outFasta", "Ouput fasta file containing the candidate TEs.", "");
	parser->addOption("minLength", "Minimum length of candidate TEs", 100);
	parser->addOption("maxLength", "Maximum length of candidate TEs", 10000);
	parser->addOption("maxSequences", "Maximium number of candidate TEs to output", 0);

	try {
		parser->parseOptions(argc, argv);
		halPath = parser->getArgument<string>("halPath");
		genomeName = parser->getArgument<string>("genome");
		gffFilename = parser->getOption<string>("outGFF");
		fastaFilename = parser->getOption<string>("outFasta");
		minLength = parser->getOption<hal_size_t>("minLength");
		maxLength = parser->getOption<hal_size_t>("maxLength");
		maxSequences = parser->getOption<int>("maxSequences");
	}
	catch (exception& e) {
		cerr << e.what() << endl;
	}

	AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, parser);

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
	if (fastaFilename != "") {
		fastaFile = fopen(fastaFilename.c_str(), "w");
	}

	if (gffFilename != "") {
		gffFile = fopen(gffFilename.c_str(), "w");
	}


	TopSegmentIteratorConstPtr topSeg = genome->getTopSegmentIterator();
	TopSegmentIteratorConstPtr endSeg = genome->getTopSegmentEndIterator();

	DNAIteratorConstPtr dnaIt = genome->getDNAIterator();
	int i = 0;
	for (; topSeg->equals(endSeg) == false; topSeg->toRight()) {
		if (topSeg->hasParent()) continue;
		if (topSeg->getLength() > maxLength) continue;
		if (topSeg->getLength() < minLength) continue;

		double nFraction = getNFraction(dnaIt, topSeg->getStartPosition(), topSeg->getLength());

		if (nFraction > 0.1) continue;

		const Sequence *sequence = topSeg->getSequence();
		int start = topSeg->getStartPosition() - sequence->getStartPosition();
		int end = topSeg->getEndPosition() - sequence->getStartPosition();

		//Print the GFF line
		fprintf(gffFile, "%s\tcandidate_transposon\tcte_%d\t%d\t%d\t0\t+\t.\t%d\n", sequence->getName().c_str(), i, start, end, i);


		if (fastaFile) {
			string seqBuffer;
			sequence->getSubString(seqBuffer, start, topSeg->getLength());
			fprintf(fastaFile, ">cte_%d\n", i);
			fprintf(fastaFile, "%s\n", seqBuffer.c_str());
		}

		i++;
		if (maxSequences > 0 && i >= maxSequences) break;
	}

	if (fastaFile) fclose(fastaFile);
	if (gffFile) fclose(gffFile);
}
