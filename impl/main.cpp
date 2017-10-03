#include <stdio.h>
#include <stdlib.h>
#include "halRepeats.h"

using namespace std;
using namespace hal;


static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
                           
  optionsParser->setDescription("Identify mutations on branch between given "
                                "genome and its parent.");

  optionsParser->addOption("start",
                           "coordinate within reference genome (or sequence"
                           " if specified) to start at",
                           0);
  optionsParser->addOption("length",
                           "length of the reference genome (or sequence"
                           " if specified) to convert.  If set to 0,"
                           " the entire thing is converted",
                           0);
  optionsParser->addOption("maxGap", 
                           "maximum indel length to be considered a gap.  Gaps "
                           " can be nested within other rearrangements.", 
                           20);
  optionsParser->addOption("maxNFraction",
                           "maximum fraction of Ns in a rearranged segment "
                           "for it to not be ignored as missing data.",
                           1.0);
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  hal_index_t start;
  hal_size_t length;
  hal_size_t maxGap;
  double nThreshold;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    maxGap = optionsParser->getOption<hal_size_t>("maxGap");
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

    HalRepeats repeats;
    repeats.analyzeBranch(alignment, maxGap, nThreshold);
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

