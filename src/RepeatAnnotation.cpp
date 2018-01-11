#include <vector>
#include <string>
#include <stack>

#include "lpo.h"
#include "align_score.h"

#include "RepeatAnnotation.h"
#include "hal.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

#include "MurmurHash3.h"


using namespace std;
using namespace hal;

char *getSequenceFromHal(const Genome *genome, hal_size_t start, hal_size_t end) {
  DNAIteratorConstPtr dnaIt = genome->getDNAIterator(start);
  char *seq = new char[end - start + 1];
  for (hal_size_t i = 0; i < end - start; i++) {
    seq[i] = toupper(dnaIt->getChar());
    dnaIt->toRight();
  }
  seq[end - start] = '\0';
  return seq;
}

uint32_t hashKmer(char *seq, int length) {
  for (int i = 0; i < length; i++) {
    if (seq[i] == 'N') return -1;
  }
  char data[8];
  MurmurHash3_x86_32(seq, length, 0, data);
  return *(uint32_t*)data;
}


double fractionN(char *seq) {
  double numN = 0.0;
  int length = strlen(seq);
  for (uint i = 0; i < 10; i++) {
    if (seq[(i*length)/10] == 'N') {
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

bool InsertionIterator::filter(char *seq) {
  if (strlen(seq) < minInsertionSize || strlen(seq) > maxInsertionSize) return false;
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
      char *seq = getSequenceFromHal(topSeg->getGenome(), topSeg->getStartPosition(), topSeg->getEndPosition());
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
      else {
        delete seq;
      }
    }
    topSeg->toRight();
  }
  return NULL;
}

/*

CRASequence *InsertionIterator::nextGappedInsertion() {
  hal_size_t start = topSeg->getStartPosition();
  hal_size_t gapLength = 0;
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      hal_size_t length = topSeg->getEndPosition() - start;
      if (length >= minInsertionSize) {
        char *seq = getSequenceFromHal(topSeg->getGenome(), start, topSeg->getEndPosition());
        if (filter(seq)) {
          CRASequence *insertion = new CRASequence;
          insertion->seq = seq;
          insertion->seqName = topSeg->getSequence()->getName();
          hal_size_t seqStart = topSeg->getSequence()->getStartPosition();
          insertion->start = start - seqStart;
          insertion->end = topSeg->getEndPosition() - seqStart;
          insertion->strand = '+';
          insertion->score = 0;
          topSeg->toRight();
          return insertion;
        }
        else {
          delete seq;
          topSeg->toRight();
          start = topSeg->getStartPosition();
          gapLength = 0;
        }
      }
    }
    else {
      //has has parent
      hal_size_t length = topSeg->getEndPosition() - topSeg->getStartPosition();
      if (length < insertionJoinDistance) {

      }
    }
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

boost::numeric::ublas::mapped_matrix<double> buildDistanceMatrix(vector<char*> seqs, int kmerLength) {
  typedef boost::unordered_map<uint32_t, vector<int> > KmerIndex;
  KmerIndex index;
  //Build index from kmers to sequences containing that kmer
  for (uint i = 0; i < seqs.size(); i++) {
    char *seq = seqs[i];
    if (i%1000 == 0) {
      cerr << "Indexed " << i << " sequences" << endl;
    }
    if (strlen(seq) < kmerLength) continue;
    for (int j = 0; j < (strlen(seq) - kmerLength); j++) {
      uint32_t kmerHash = hashKmer(seq + j, kmerLength);
      if (kmerHash == -1) continue;
      index[kmerHash].push_back(i);
    }
  }
  cerr << "Finished building kmer index " << endl;

  int N = seqs.size();
  boost::numeric::ublas::mapped_matrix<double> dist(N, N, N);

  int npairs = 0;
  BOOST_FOREACH(KmerIndex::value_type kv, index) {
    vector<int> seqsWithKmer = kv.second;
    for (uint i = 0; i < seqsWithKmer.size(); i++) {
      for (uint j = 0; j < i; j++) {
        if (seqsWithKmer[i] == seqsWithKmer[j]) continue;
        npairs++;
        if (npairs % 1000000 == 0) cerr << "Filled " << npairs << " matrix cells" << endl;
        int a = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[i] : seqsWithKmer[j];
        int b = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[j] : seqsWithKmer[i];
        dist (a, b) += 1.0;
      }
    }
  }
  cerr << "Computed " << npairs << " pairwise distances" << endl;

  cerr << "Finished computing pairwise distances" << endl;

  //Divide by the number of kmers in each sequence
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < i; j++) {
      if (dist(i, j) == 0.0) {
        continue;
      }
      string a = seqs[i];
      string b = seqs[j];
      if (a.length() < kmerLength || b.length() < kmerLength) {
        dist (i, j) = 1.0;
      }
      double nKmers = (double)(a.length())/kmerLength + (double)(b.length())/kmerLength;
      dist (i, j) = dist (i, j) / nKmers;
    }
  }
  cerr << "Finished normalizing distances" << endl;
  return dist;
}

vector<CRASequence*> annotateRepeatsOnBranch(const hal::Genome *genome, InsertionIterator &insertionIter, hal_size_t maxInsertions) {
  insertionIter.goToGenome(genome);
  CRASequence *insertion;
  vector<CRASequence*> insertions;
  int i = 0;
  while((insertion = insertionIter.next())) {
    if (i%1000 == 0) cerr << "Read " << i << " insertions" << endl;
    i++;
    if (maxInsertions != 0 && i >= maxInsertions) break;
    insertions.push_back(insertion);
  }
  cerr << "Found " << insertions.size() << " candidate insertions on branch " << genome->getName() << endl;

  vector<char*> seqs;
  for (uint i = 0; i < insertions.size(); i++) {
    seqs.push_back(insertions[i]->seq);
  }
  cerr << "Built distance matrix of size " << seqs.size() << endl;
  boost::numeric::ublas::mapped_matrix<double> distanceMatrix = buildDistanceMatrix(seqs, 30);
  cerr << "Finished building distance matrix" << endl;
  map<CRASequence*, vector<CRASequence*> > clusters = buildTransitiveClusters<CRASequence>(insertions, distanceMatrix, 0.3);
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
    cout << genome->getName() << " " << insertion->end - insertion->start << endl;
  }
  delete insertion;
}

vector<CRASequence*> liftoverRepeatAnnotations(vector<CRASequence*> repeats, const hal::Genome *source, const hal::Genome *target) {

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
