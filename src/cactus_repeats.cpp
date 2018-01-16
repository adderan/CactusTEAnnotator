#include <stdio.h>
#include <vector>
#include <map>
#include "hal.h"
#include "RepeatAnnotation.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace std;
using namespace hal;
using namespace boost::numeric::ublas;


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

  optionsParser->addOption("kmerSize", "Length of seeds to use for insertion clustering", 20);
  optionsParser->addOption("chunkSize", "Number of insertions to cluster in one chunk", 1000);

  optionsParser->addOption("similarityThreshold", "Similarity level for joining sequences into a cluster", 0.4);
  optionsParser->addOptionFlag("getInsertionLengths", "", false);
  optionsParser->addOptionFlag("getInsertions", "", false);
  optionsParser->addOptionFlag("annotateInsertions", "", false);
  optionsParser->addOption("maxInsertions", "Maximum number of insertions to process", 0);
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
  int maxInsertions;
  string referenceName;
  double maxNFraction;


  int kmerSize;
  int chunkSize;
  double similarityThreshold;

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
    maxInsertions = optionsParser->getOption<int>("maxInsertions");

    //Clustering
    kmerSize = optionsParser->getOption<int>("kmerSize");
    chunkSize = optionsParser->getOption<int>("chunkSize");
    similarityThreshold = optionsParser->getOption<double>("similarityThreshold");

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
    const hal::Genome *reference = alignment->openGenome(referenceName);
    insertionIt.goToGenome(reference);
    Seq *insertion;
    std::vector<std::vector<Seq*> *> chunks;
    int n_insertions = 0;
    std::vector<Seq*> *chunk;
    while((insertion = insertionIt.next())) {
      if (n_insertions%chunkSize == 0) {
        cerr << "Read " << n_insertions << " insertions" << endl;
        chunk = new std::vector<Seq*>();
        chunks.push_back(chunk);
      }
      n_insertions++;
      if (maxInsertions != 0 && n_insertions >= maxInsertions) break;
      chunk->push_back(insertion);
    }
    cerr << "Found " << n_insertions << " candidate insertions on branch " << reference->getName() << endl;

    map<Seq*, std::vector<Seq*> > groups;
    for (int i = 0; i < chunks.size(); i++) {
      cerr << "Processing chunk of size " << chunks[i]->size() << endl;
      mapped_matrix<double> similarityMatrix = buildDistanceMatrix(*chunks[i], kmerSize);
      map<Seq*, std::vector<Seq*> > groups_chunk_i = buildGroups(*chunks[i], similarityMatrix, 
          similarityThreshold);

      for (map<Seq*, std::vector<Seq*> >::iterator it = groups_chunk_i.begin(); it != groups_chunk_i.end(); it++) {
        groups[it->first] = it->second;
      }
    }


    //produce final output
    std::vector<Seq*> repeats;
    map<Seq*, std::vector<Seq*> >::iterator clusterIter;
    int familyNumber = 0;
    for (clusterIter = groups.begin(); clusterIter != groups.end(); clusterIter++) {
      std::vector<Seq*> insertionsInCluster = clusterIter->second;
      if (insertionsInCluster.size() > 1) {
        for(uint i = 0; i < insertionsInCluster.size(); i++) {
          Seq* insertion = insertionsInCluster[i];
          insertion->repeatFamily = "cactus";
          insertion->group = familyNumber;
          repeats.push_back(insertion);
        }
        familyNumber++;
      }
    }

    for (uint i = 0; i < repeats.size(); i++) {
      repeats[i]->toGFF(&cout);
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
