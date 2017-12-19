#include <vector>
#include <string>
#include <map>
#include <stack>

#include "lpo.h"
#include "align_score.h"

#include "RepeatAnnotation.h"
#include "hal.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>


using namespace std;
using namespace hal;
using namespace boost::numeric::ublas;



double fractionN(string seq) {
  double numN = 0.0;
  for (unsigned int i = 0; i < 10; i++) {
    if (seq[(i*seq.length()/10] == 'N') {
      numN += 1.0;
    }
  }
  return numN/10.0;
}

void CRASequence::toGFF(ostream* gffStream) {
  *gffStream << seqName << "\tcactus_repeat_annotation\trepeat_copy" << "\t" << start << "\t" << end <<  "\t" << score << "\t" << strand << "\t" << "." << "\t" << group << endl;
}

void InsertionIterator::goToGenome(const hal::Genome * _genome) {
  genome = _genome;
  topSeg = genome->getTopSegmentIterator();
  endSeg = genome->getTopSegmentEndIterator();
}

bool InsertionIterator::filter(string seq) {
  if (seq.length() < minInsertionSize || seq.length() > maxInsertionSize) return false;
  if (fractionN(seq) > maxNFraction) {
    return false;
  }
  return true;
}


CRASequence* InsertionIterator::next() {
  if (insertionJoinDistance > 0) {
    return NULL;
  }
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      string seq;
      seq.clear();
      topSeg->getString(seq);
      if (filter(seq)) {
	CRASequence *insertion = new CRASequence;
	insertion->seq = seq;
	insertion->seqName = topSeg->getSequence()->getName();
	hal_size_t seqStart = topSeg->getSequence()->getStartPosition();
	insertion->start = topSeg->getStartPosition() - seqStart;
	insertion->end = topSeg->getEndPosition() - seqStart;
	insertion->strand = '+';
	insertion->score = 0;
	topSeg->toRight();
	return insertion;
      }
    }
    topSeg->toRight();
  }
  return NULL;
}

/*
CRASequence *CRASequenceIterator::nextGappedCRASequence()

{
  
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      //Found an insertion. Scan forward for another insertion within a certain distance
      string insertion;
      topSeg->getString(insertion);
      while(topSeg->equals(endSeg) == false) {
	topSeg->toRight();
	if (!topSeg->hasParent) {
	  string seq 
  }
  
}
*/


GenomeIterator::GenomeIterator(hal::AlignmentConstPtr _alignment) {
  alignment = _alignment;
  root = alignment->openGenome(alignment->getRootName());
  visited.push(root);
}

const hal::Genome * GenomeIterator::next() {
  //dfs on the hal species tree
  if (visited.empty()) return NULL;

  root = visited.top();
  visited.pop();

  for (hal_size_t childIndex = 0; childIndex < root->getNumChildren(); childIndex++) {
    const hal::Genome* child = root->getChild(childIndex);
    visited.push(child);
  }
  return root;
  
}

double **buildDistanceMatrix(vector<string> seqs, int kmerLength) {
  map<string, vector<int> > kmerIndex;
  //Build index from kmers to sequences containing that kmer
  for (uint i = 0; i < seqs.size(); i++) {
    string seq = seqs[i];
    if (seq.length() < kmerLength) continue;
    for (int j = 0; j < (seq.length() - kmerLength); j++) {
      kmerIndex[seq.substr(j, j + kmerLength)].push_back(i);
    }
  }
  cerr << "Finished building kmer index " << endl;

  //Allocate only above the diagonal of the distance matrix
  int N = seqs.size();
  double **dist = new double*[N];
  for (uint i = 0; i < N; i++) {
    dist[i] = (double*) calloc(i, sizeof(double));
  }
  
  for (map<string, vector<int> >::iterator it = kmerIndex.begin(); it != kmerIndex.end(); it++) {
    vector<int> seqsWithKmer = it->second;
    for (uint i = 0; i < seqsWithKmer.size(); i++) {
      for (uint j = 0; j < i; j++) {
	if (seqsWithKmer[i] == seqsWithKmer[j]) continue;
	int a = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[i] : seqsWithKmer[j];
	int b = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[j] : seqsWithKmer[i];
	dist[a][b] += 1.0;
      }
    }
  }

  cerr << "Finished computing pairwise distances" << endl;

  //Divide by the number of kmers in each sequence
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < i; j++) {
      string a = seqs[i];
      string b = seqs[j];
      if (a.length() < kmerLength || b.length() < kmerLength) {
	dist[i][j] = 1.0;
      }
      double nKmers = (double)(a.length())/kmerLength + (double)(b.length())/kmerLength;
      dist[i][j] = 1.0 - dist[i][j]/nKmers;
    }
  }
  return dist;
}

vector<CRASequence*> annotateRepeatsOnBranch(const hal::Genome *genome, InsertionIterator &insertionIter) {
  insertionIter.goToGenome(genome);
  CRASequence *insertion;
  vector<CRASequence*> insertions;
  while((insertion = insertionIter.next())) {
    insertions.push_back(insertion);
  }
  cerr << "Found " << insertions.size() << " candidate insertions on branch " << genome->getName() << endl;

  vector<string> seqs;
  for (uint i = 0; i < insertions.size(); i++) {
    seqs.push_back(insertions[i]->seq);
  }
  cerr << "Built distance matrix of size " << seqs.size() << endl;
  double **distanceMatrix = buildDistanceMatrix(seqs, 5);
  cerr << "Finished building distance matrix" << endl;
  map<CRASequence*, vector<CRASequence*> > clusters = buildTransitiveClusters<CRASequence>(insertions, distanceMatrix, 0.8);
  cerr << "Built " << clusters.size() << " clusters from " << insertions.size() << " insertions " << endl;
  
  vector<CRASequence*> repeats;
  map<CRASequence*, vector<CRASequence*> >::iterator clusterIter;
  int familyNumber = 0;
  for (clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++) {
    vector<CRASequence*> insertionsInCluster = clusterIter->second;
    if (insertionsInCluster.size() > 1) {
      for(uint i = 0; i < insertionsInCluster.size(); i++) {
	CRASequence* insertion = insertionsInCluster[i];
	insertion->repeatFamily = "cactus";
	insertion->group = familyNumber;
	repeats.push_back(insertion);
      }
      familyNumber++;
    }
  }
  return repeats;
  
}

void getInsertionLengthsOnBranch(const hal::Genome* genome, InsertionIterator &insertionIt) {
  insertionIt.goToGenome(genome);
  CRASequence *insertion;
  while((insertion = insertionIt.next())) {
    cout << genome->getName() << " " << (insertion->seq).length() << endl;
  }
  delete insertion;
}

vector<CRASequence*> liftoverRepeatAnnotations(vector<CRASequence*> repeats, const hal::Genome *source, const hal::Genome *target) {
  
}





set<string> getSeeds(string seq, int seedLength) {
  assert(seq.length() < seedLength);
  set<string> seeds;
  int i = 0;
  while (i < (seq.length() - seedLength)) {
    seeds.insert(seq.substr(i, i + seedLength));
    i = i + seedLength;
  }
  return seeds;
}

double kmerDistance(string seq1, string seq2) {
  int seedLength = 5;
  if (seq1.length() < seedLength || seq2.length() < seedLength) return 1.0;
  set<string> seq1Seeds = getSeeds(seq1, seedLength);
  set<string> seq2Seeds = getSeeds(seq2, seedLength);
  set<string> intersectionSeeds;
  set<string> unionSeeds;
  set_intersection(seq1Seeds.begin(), seq1Seeds.end(), seq2Seeds.begin(), seq2Seeds.end(), inserter(intersectionSeeds, intersectionSeeds.begin()));
  set_union(seq1Seeds.begin(), seq1Seeds.end(), seq2Seeds.begin(), seq2Seeds.end(), inserter(unionSeeds, unionSeeds.begin()));
  int intersectionSize = distance(intersectionSeeds.begin(), intersectionSeeds.end());
  int unionSize = distance(unionSeeds.begin(), unionSeeds.end());
  double dist = 1.0 - (double)(intersectionSize)/unionSize;
  cerr << "Computed distance " << dist << endl;
  return dist;
}
  
/*
LPOCRASequence_T *buildSequenceGraph(vector<Sequence*> &sequences, char *matrixFilename) {
  ResidueScoreMatrix_T scoreMatrix;
  
  read_score_matrix(matrixFilename, &scoreMatrix);
  
  LPOSequence_T **inputSeqs =  (LPOSequence_T**)calloc(sequences.size(), sizeof(LPOSequence_T*));
  for (uint i = 0; i < sequences.size(); i++) {
    char *seqName = new char[(sequences[i]->seqName).length() + 1];
    strcpy(seqName, (sequences[i]->seqName).c_str());
    char *seq = new char[(sequences[i]->seq).length() + 1];
    strcpy(seq, (sequences[i]->seq).c_str());
    create_seq(i, inputSeqs, seqName, seqName, seq, false)
  }

  LPOSequence_T *align = buildup_progressive_lpo(sequences.size(), inputSeqs, &scoreMatrix,
						 false, true, NULL, matrix_scoring_function, false, true);
  return align;
}
  
*/
