#include <stdio.h>
#include <stdlib.h>
#include "insertions.h"
#include "annotation.h"
#include "clustering.h"

using namespace std;
using namespace hal;


static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
                           
  optionsParser->setDescription("Identify mutations on branch between given "
                                "genome and its parent.");

  optionsParser->addOption("minInsertionSize", 
                           "minimum insertion length to consider.",
                           50);
  optionsParser->addOption("reference",
			   "Genome to get insertions for",
			   "");
  optionsParser->addOption("insertionJoinDistance", "Maximum length of gaps in insertions.", 0);
  optionsParser->addOptionFlag("getLengths", "Get insertion lengths", true);
  optionsParser->addOption("maxNFraction", "Maximum fraction of Ns in insertion", 0.1);
  
  optionsParser->addOption("seedLength", "Length of seeds to use for insertion clustering", 20);

  optionsParser->addOptionFlag("getInsertionLengths", "", false);
  optionsParser->addOptionFlag("buildClusters", "", false);
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
  string referenceName;
  double maxNFraction;

  
  hal_size_t seedLength;

  bool getInsertionLengths;
  bool buildClusters;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");

    //Insertions
    referenceName = optionsParser->getOption<string>("reference");
    minInsertionSize = optionsParser->getOption<hal_size_t>("minInsertionSize");
    insertionJoinDistance = optionsParser->getOption<hal_size_t>("insertionJoinDistance");
    maxNFraction = optionsParser->getOption<double>("maxNFraction");

    //Clustering
    seedLength = optionsParser->getOption<hal_size_t>("seedLength");

    //Execution modes
    getInsertionLengths = optionsParser->getFlag("getInsertionLengths");
    buildClusters = optionsParser->getFlag("buildClusters");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath,
							   optionsParser);

    InsertionIterator insertionIt = InsertionIterator(maxNFraction, insertionJoinDistance, minInsertionSize);
    if (getInsertionLengths) {
      if (referenceName != "") {
	const Genome *reference = alignment->openGenome(referenceName);
	insertionIt.goToGenome(reference);
	Insertion *insertion;
	while((insertion = insertionIt.next())) {
	  cerr << (insertion->seq).length() << endl;
	}
	delete insertion;
      }
    }

    else if (buildClusters) {
      const Genome *reference = alignment->openGenome(referenceName);
      insertionIt.goToGenome(reference);
      Insertion *insertion;
      vector<Insertion*> insertions;
      int nInsertions = 0;
      int maxInsertions = 100;
      while((insertion = insertionIt.next())) {
	if (nInsertions == maxInsertions) break;
	insertions.push_back(insertion);
	nInsertions++;
      }
      map<Insertion*, vector<Insertion*> > clusters = buildTransitiveClusters<Insertion>(insertions, &insertionDistance, 0.5);
      cerr << "Built " << clusters.size() << " clusters" << endl;
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

