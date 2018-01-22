#include <stdio.h>
#include <vector>
#include <map>
#include "hal.h"
#include "omp.h"
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

  optionsParser->addOption("numThreads", "Maximum number of threads to use", 8);

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

  optionsParser->addOption("kmerLength", "Length of seeds to use for insertion clustering", 20);
  optionsParser->addOption("chunkSize", "Number of insertions to cluster in one chunk", 1000);

  optionsParser->addOption("similarityThreshold", "Similarity level for joining sequences into a cluster", 0.4);
  optionsParser->addOption("groupJoinThreshold", "Jaccard distance threshold for joining clusters", 0.1);
  optionsParser->addOption("minGroupSize", "Minimium size of repeat family to output", 5);
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
  int numThreads;

  hal_size_t minInsertionSize;
  hal_size_t maxInsertionSize;
  hal_size_t insertionJoinDistance;
  int maxInsertions;
  string referenceName;
  double maxNFraction;


  int kmerLength;
  int chunkSize;
  double similarityThreshold;
  double groupJoinThreshold;
  int minGroupSize;

  bool getInsertionLengths;
  bool annotateInsertions;
  bool getInsertions;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getOption<string>("alignment");
    insertionsFasta = optionsParser->getOption<string>("insertionsFasta");
    numThreads = optionsParser->getOption<int>("numThreads");

    //Insertions
    referenceName = optionsParser->getOption<string>("reference");
    minInsertionSize = optionsParser->getOption<hal_size_t>("minInsertionSize");
    maxInsertionSize = optionsParser->getOption<hal_size_t>("maxInsertionSize");
    insertionJoinDistance = optionsParser->getOption<hal_size_t>("insertionJoinDistance");
    maxNFraction = optionsParser->getOption<double>("maxNFraction");
    maxInsertions = optionsParser->getOption<int>("maxInsertions");

    //Clustering
    kmerLength = optionsParser->getOption<int>("kmerLength");
    chunkSize = optionsParser->getOption<int>("chunkSize");
    similarityThreshold = optionsParser->getOption<double>("similarityThreshold");
    groupJoinThreshold = optionsParser->getOption<double>("groupJoinThreshold");
    minGroupSize = optionsParser->getOption<int>("minGroupSize");

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

    int nTotalInsertions = 0;
    std::vector<std::vector<Seq*> > groups;
    #pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < chunks.size(); i++) {

      #pragma omp critical
      {
        cerr << "Processing chunk # " << i << " of " << chunks.size() << ", of size " << chunks[i]->size() << endl;
      }
      mapped_matrix<double> similarityMatrix = buildDistanceMatrix(*(chunks[i]), kmerLength);
      map<Seq*, std::vector<Seq*> > groups_chunk_i = buildGroups(*(chunks[i]), similarityMatrix, 
          similarityThreshold);

      #pragma omp critical
      {
        for (map<Seq*, std::vector<Seq*> >::iterator it = groups_chunk_i.begin(); 
            it != groups_chunk_i.end(); it++) {
          if (it->second.size() > minGroupSize) {
            groups.push_back(it->second);
            nTotalInsertions += it->second.size();
          }
        }
      }
    }
    cerr << "Found " << nTotalInsertions << " total insertions before joining clusters." << endl;

    //join groups based on kmer distance
    joinGroups(groups, kmerLength, groupJoinThreshold);


    //produce final output
    std::vector<Seq*> repeats;
    std::vector<std::vector<Seq*> >::iterator clusterIter;
    int familyNumber = 0;
    for (clusterIter = groups.begin(); clusterIter != groups.end(); clusterIter++) {
      std::vector<Seq*> insertionsInCluster = *clusterIter;
      for(uint i = 0; i < insertionsInCluster.size(); i++) {
        Seq* insertion = insertionsInCluster[i];
        insertion->repeatFamily = "cactus";
        insertion->group = familyNumber;
        repeats.push_back(insertion);
      }
      familyNumber++;
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
