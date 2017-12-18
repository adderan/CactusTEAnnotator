#include <stdio.h>
#include <stdlib.h>
#include "insertions.h"
#include "annotation.h"
#include "clustering.h"

using namespace std;


static hal::CLParserPtr initParser()
{
  hal::CLParserPtr optionsParser = hal::hdf5CLParserInstance();
  optionsParser->addOption("alignment", "input hal file", "");
  optionsParser->addOption("insertionsFasta", "Input fasta file containing insertions to annotate with repeat classes", "");
                           
  optionsParser->setDescription("Identify mutations on branch between given "
                                "genome and its parent.");

  optionsParser->addOption("minInsertionSize", 
                           "minimum insertion length to consider.",
                           100);
  optionsParser->addOption("maxInsertionSize",
			   "maximum insertion length to consider.",
			   50000);
  optionsParser->addOption("reference",
			   "Genome to get insertions for",
			   "");
  optionsParser->addOption("insertionJoinDistance", "Maximum length of gaps in insertions.", 0);
  optionsParser->addOptionFlag("getLengths", "Get insertion lengths", true);
  optionsParser->addOption("maxNFraction", "Maximum fraction of Ns in insertion", 0.1);
  
  optionsParser->addOption("seedLength", "Length of seeds to use for insertion clustering", 20);

  optionsParser->addOptionFlag("getInsertionLengths", "", false);
  optionsParser->addOptionFlag("getInsertions", "", false);
  optionsParser->addOptionFlag("annotateInsertions", "", false);
  return optionsParser;
}

int main(int argc, char** argv)
{
  hal::CLParserPtr optionsParser = initParser();

  string halPath;
  string insertionsFasta;
  
  hal_size_t minInsertionSize;
  hal_size_t maxInsertionSize;
  hal_size_t insertionJoinDistance;
  string referenceName;
  double maxNFraction;

  
  hal_size_t seedLength;

  bool getInsertionLengths;
  bool annotateInsertions;
  bool getInsertions;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getOption<string>("alignment");
    insertionsFasta = optionsParser->getOption<string>("insertionsFasta");

    //Insertions
    referenceName = optionsParser->getOption<string>("reference");
    minInsertionSize = optionsParser->getOption<hal_size_t>("minInsertionSize");
    maxInsertionSize = optionsParser->getOption<hal_size_t>("maxInsertionSize");
    insertionJoinDistance = optionsParser->getOption<hal_size_t>("insertionJoinDistance");
    maxNFraction = optionsParser->getOption<double>("maxNFraction");

    //Clustering
    seedLength = optionsParser->getOption<hal_size_t>("seedLength");

    //Execution modes
    getInsertionLengths = optionsParser->getFlag("getInsertionLengths");
    annotateInsertions = optionsParser->getFlag("annotateInsertions");
    getInsertions = optionsParser->getFlag("getInsertions");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  try
  {
    hal::AlignmentConstPtr alignment = hal::openHalAlignmentReadOnly(halPath,
							   optionsParser);

    InsertionIterator insertionIt = InsertionIterator(maxNFraction, insertionJoinDistance, minInsertionSize, maxInsertionSize);
    if (getInsertionLengths) {
      if (referenceName != "") {
	const hal::Genome *reference = alignment->openGenome(referenceName);
	getInsertionLengthsOnBranch(reference, insertionIt);
      }
      else {
	GenomeIterator genomeIt(alignment);
	const hal::Genome *genome;
	while ((genome = genomeIt.next()) != NULL) {
	  getInsertionLengthsOnBranch(genome, insertionIt);
	}
      }
    }

    else if (annotateInsertions) {
      if (referenceName != "") {
	const hal::Genome *reference = alignment->openGenome(referenceName);
	vector<Sequence*> repeats = annotateRepeatsOnBranch(reference, insertionIt);
	for (uint i = 0; i < repeats.size(); i++) {
	  repeats[i]->toGFF(&cout);
	}
      }
    }

  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  
  return 0;
}

