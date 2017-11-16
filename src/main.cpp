#include <stdio.h>
#include <stdlib.h>
#include "insertions.h"
#include "annotation.h"

using namespace std;
using namespace hal;


static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("executionMode", "");
                           
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
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string executionMode;
  
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
  string referenceName;
  double maxNFraction;
  hal_size_t seedLength;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    executionMode = optionsParser->getArgument<string>("executionMode");
    
    referenceName = optionsParser->getOption<string>("reference");
    minInsertionSize = optionsParser->getOption<hal_size_t>("minInsertionSize");
    insertionJoinDistance = optionsParser->getOption<hal_size_t>("insertionJoinDistance");
    maxNFraction = optionsParser->getOption<double>("maxNFraction");
    
    seedLength = optionsParser->getOption<hal_size_t>("seedLength");
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
    RepeatAnnotatorOpts opts = {minInsertionSize,
				insertionJoinDistance,
				optionsParser,
				referenceName,
				maxNFraction,
				seedLength};
    if (executionMode == "getLengths") {
      getInsertions(alignment, opts);
    }
    else if (executionMode == "buildClusters") {
      buildClusters(alignment, opts);
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

