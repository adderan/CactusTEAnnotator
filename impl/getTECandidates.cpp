#include <getopt.h>
#include <stdlib.h>
#include <string>
#include "hal.h"
#include <fstream>

extern "C" {
#include "sonLib.h"
}


using namespace hal;
using namespace std;


void printGffEntry(ostream *output, const vector<TopSegmentIteratorConstPtr> &segments, int seqNum, bool forward) {

	if (forward) {
		for (int i = 0; i < segments.size(); i++) {
			const Sequence *sequence = segments[i]->getSequence();
			hal_size_t start = sequence->getStartPosition() + segments[i]->getStartPosition();
			hal_size_t end = sequence->getStartPosition() + segments[i]->getEndPosition();

			char strand = '+';
			*output <<
				sequence->getName() << "\t"
				"candidate_transposon" << "\t" <<
				"cte" << seqNum << "_" << i << "\t" <<
				start << "\t" <<
				end << "\t" <<
				"0" << "\t" <<
				strand << "\t" <<
				"." << "\t" <<
				seqNum << endl;
		}
	}
}

void printFastaEntry(ostream *output, const vector<TopSegmentIteratorConstPtr> &segments, int seqNum, bool forward) {

	*output << ">" << "cte" << seqNum << endl;
	for (int i = 0; i < segments.size(); i++) {
		const Sequence *sequence = segments[i]->getSequence();
		string seq;
		hal_size_t start = sequence->getStartPosition() +segments[i]->getStartPosition();
		hal_size_t length = segments[i]->getLength();
		sequence->getSubString(seq, start, length);
		*output << seq;
	}
	*output << endl;
}

int main(int argc, char **argv) {
	string halPath;
	string genomeName;
	hal_size_t minLength;
	hal_size_t maxLength;
	string fastaFilename;
	string gffFilename;
	hal_size_t maxDiscardedSequence;
	int maxSequences = 0;

	CLParserPtr parser = hdf5CLParserInstance(false);
	parser->addArgument("halPath", "Input hal file");
	parser->addArgument("genome", "Name of genome to scan");
	parser->addOption("outGFF", "Output GFF representing the candidate TEs", "");
	parser->addOption("outFasta", "Ouput fasta file containing the candidate TEs.", "");
	parser->addOption("minLength", "Minimum length of candidate TEs", 100);
	parser->addOption("maxLength", "Maximum length of candidate TEs", 10000);
	parser->addOption("maxSequences", "Maximium number of candidate TEs to output", 0);

	parser->addOption("maxDiscardedSequence", "Maximum distance to join multiple HAL insertions as one candidate TE", 100);

	try {
		parser->parseOptions(argc, argv);
		halPath = parser->getArgument<string>("halPath");
		genomeName = parser->getArgument<string>("genome");
		gffFilename = parser->getOption<string>("outGFF");
		fastaFilename = parser->getOption<string>("outFasta");
		minLength = parser->getOption<hal_size_t>("minLength");
		maxLength = parser->getOption<hal_size_t>("maxLength");
		maxSequences = parser->getOption<int>("maxSequences");
		maxDiscardedSequence = parser->getOption<hal_size_t>("maxDiscardedSequence");
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

	ostream *gffFile = &cout;
	ofstream gffStream;
	ostream *fastaFile = NULL;
	ofstream fastaStream;
	if (fastaFilename != "") {
		fastaStream.open(fastaFilename);
		fastaFile = &fastaStream;
	}

	if (gffFilename != "") {
		gffStream.open(gffFilename);
		gffFile = &gffStream;
	}

	TopSegmentIteratorConstPtr topSeg = genome->getTopSegmentIterator();
	TopSegmentIteratorConstPtr endSeg = genome->getTopSegmentEndIterator();

	int i = 0;
	vector<TopSegmentIteratorConstPtr> segments;
	while(true) {
		hal_size_t discardedSequence = 0;
		hal_size_t seqLength = 0;
		//scan to the next insert
		while(!topSeg->equals(endSeg) && topSeg->hasParent()) {
			topSeg->toRight();
		}

		while (true) {
			if (topSeg->hasParent()) {
				discardedSequence += topSeg->getLength();
				if (discardedSequence > maxDiscardedSequence) break;
			}
			else {
				segments.push_back(topSeg->copy());
				seqLength += topSeg->getLength();
			}
			if (topSeg->equals(endSeg)) break;
			topSeg->toRight();
		}

		if (seqLength < minLength) continue;
		if (seqLength > maxLength) continue;

		printGffEntry(gffFile, segments, i, true);

		if (fastaFile) {
			//forward strand
			printFastaEntry(fastaFile, segments, i, true);
		}

		segments.clear();
		i++;
		if (maxSequences > 0 && i >= maxSequences) break;
		if (topSeg->equals(endSeg)) break;
	}

	if (fastaFile) fastaStream.close();
	if (gffFile) gffStream.close();
}
