#include <stdio.h>
#include <fstream>
#include "hal.h"

using namespace hal;
using namespace std;

int main(int argc, char **argv) {
	
	CLParserPtr parser = hdf5CLParserInstance(false);
	AlignmentConstPtr alignment = openHalAlignmentReadOnly(argv[1], parser);
	
	ifstream gffFile(argv[2]);

	if (argv[3] == NULL) {
		cerr << "Must specify genome" << endl;
		exit(1);
	}

	const Genome *genome = alignment->openGenome(argv[3]);
	DNAIteratorConstPtr it = genome->getDNAIterator(0);

	string chrom, annotationType, name, strand, score, a, group;

	hal_size_t start, end;
	const Sequence *sequence;
	while (gffFile >> chrom >> annotationType >> name >> start >> end >> score >> strand >> a >> group) {
		sequence = genome->getSequence(chrom);

		cout << ">" << name << endl;
		it->jumpTo(sequence->getStartPosition() + start);
		while(it->getArrayIndex() <= sequence->getStartPosition() + end) {
			cout << it->getChar();
			it->toRight();
		}
		cout << endl;
	}

}
