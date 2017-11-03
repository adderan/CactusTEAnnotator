#include <stdio.h>
#include <stdlib.h>
#include "insertions.h"

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
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  hal_size_t minInsertionSize;
  hal_size_t insertionJoinDistance;
  string referenceName;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    referenceName = optionsParser->getOption<string>("reference");
    minInsertionSize = optionsParser->getOption<hal_size_t>("minInsertionSize");
    insertionJoinDistance = optionsParser->getOption<hal_size_t>("insertionJoinDistance");
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
				optionsParser};
    if (referenceName == "") {
      getInsertions(alignment, opts);

    }
    else {
      const Genome *reference = alignment->openGenome(referenceName);
      InsertionIterator insertionIterator(reference, opts);
      std::string insertionSeq;
      while ((insertionSeq = insertionIterator.next()) != "") {
	cout << reference->getName() << " " << insertionSeq.length() << endl;
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

