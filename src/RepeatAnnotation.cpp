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
//using namespace boost::numeric::ublas;

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
    if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G') return -1;
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

void Seq::toGFF(ostream* gffStream) {
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

map<Seq*, vector<Seq*> > buildGroups(vector<Seq*> &seqs, boost::numeric::ublas::mapped_matrix<double> &similarityMatrix, double similarityThreshold) {

  map<Seq *,vector<Seq *> > clusterToSeq;
  map<Seq *,Seq *> seqToCluster;

  for (uint i = 0; i < seqs.size(); i++) {
    seqToCluster[seqs[i]] = seqs[i];
    clusterToSeq[seqs[i]].push_back(seqs[i]);
  }
  for (uint i = 0; i < seqs.size(); i++) {
    for (uint j = 0; j < i; j++) {
      Seq *a = seqs[i];
      Seq *b = seqs[j];
      double similarity = similarityMatrix (i, j);
      if (similarity > similarityThreshold) {
        //Combine the clusters
        Seq* cluster_a = seqToCluster[a];
        Seq* cluster_b = seqToCluster[b];
        if (cluster_a != cluster_b) {
          Seq* new_cluster = (cluster_a > cluster_b) ? cluster_a : cluster_b;
          Seq* cluster_to_delete = (cluster_a > cluster_b) ? cluster_b : cluster_a;
          vector<Seq*> seqsInCluster = clusterToSeq[cluster_to_delete];
          vector<Seq*>::iterator it;
          for (it = seqsInCluster.begin(); it != seqsInCluster.end(); it++) {
            seqToCluster[*it] = new_cluster;
            clusterToSeq[new_cluster].push_back(*it);
          }
          clusterToSeq.erase(cluster_to_delete);
        }
      }
    }
  }
  return clusterToSeq;
};


Seq* InsertionIterator::next() {
  if (insertionJoinDistance > 0) {
    return NULL;
  }
  while (topSeg->equals(endSeg) == false) {
    if (!topSeg->hasParent()) {
      char *seq = getSequenceFromHal(topSeg->getGenome(), topSeg->getStartPosition(), topSeg->getEndPosition());
      if (filter(seq)) {
        Seq *insertion = new Seq;
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

   Seq *InsertionIterator::nextGappedInsertion() {
   hal_size_t start = topSeg->getStartPosition();
   hal_size_t gapLength = 0;
   while (topSeg->equals(endSeg) == false) {
   if (!topSeg->hasParent()) {
   hal_size_t length = topSeg->getEndPosition() - start;
   if (length >= minInsertionSize) {
   char *seq = getSequenceFromHal(topSeg->getGenome(), start, topSeg->getEndPosition());
   if (filter(seq)) {
   Seq *insertion = new Seq;
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

boost::numeric::ublas::mapped_matrix<double> buildDistanceMatrix(vector<Seq*> &seqs, int kmerLength) {
  typedef boost::unordered_map<uint32_t, vector<int> > KmerIndex;
  KmerIndex index;
  //Build index from kmers to sequences containing that kmer
  for (uint i = 0; i < seqs.size(); i++) {
    char *seq = seqs[i]->seq;
    if (strlen(seq) < kmerLength) continue;
    for (int j = 0; j < (strlen(seq) - kmerLength); j++) {
      uint32_t kmerHash = hashKmer(seq + j, kmerLength);
      if (kmerHash == -1) continue;
      index[kmerHash].push_back(i);
    }
  }
  int N = seqs.size();
  boost::numeric::ublas::mapped_matrix<double> dist(N, N, N);

  int npairs = 0;
  int nrows = 0;
  BOOST_FOREACH(KmerIndex::value_type kv, index) {
    vector<int> seqsWithKmer = kv.second;
    nrows++;
    for (uint i = 0; i < seqsWithKmer.size(); i++) {
      for (uint j = 0; j < i; j++) {
        if (seqsWithKmer[i] == seqsWithKmer[j]) continue;

        int a = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[i] : seqsWithKmer[j];
        int b = (seqsWithKmer[i] > seqsWithKmer[j]) ? seqsWithKmer[j] : seqsWithKmer[i];
        dist (a, b) += 1.0;
      }
    }
  }
  //Divide by the number of kmers in each sequence
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < i; j++) {
      if (dist(i, j) == 0.0) {
        continue;
      }
      Seq *a = seqs[i];
      Seq *b = seqs[j];
      if (strlen(a->seq) < kmerLength || strlen(b->seq) < kmerLength) {
        dist (i, j) = 1.0;
      }
      double nKmers = ((double)strlen(a->seq) - (double)kmerLength) + ((double)strlen(b->seq) - (double)kmerLength);

      //avoid double-counting shared kmers
      double unionSize = nKmers - dist(i, j);
      dist (i, j) = dist (i, j) / unionSize;
    }
  }
  cerr << "Finished normalizing distances" << endl;
  return dist;
}


set<uint32_t> getKmers(vector<Seq*> &seqs, int kmerLength, int maxKmers) {

  set<uint32_t> kmers;
  for (int i = 0; i < maxKmers; i++) {
    int seqnum = (int) (((double)seqs.size()/RAND_MAX) * rand());
    int seqLen = strlen(seqs[seqnum]->seq);
    if (seqLen < kmerLength) continue;
    int pos = (int)(((double)(seqLen - kmerLength - 1)/RAND_MAX) * rand());
    uint32_t kmer = hashKmer(seqs[seqnum]->seq + pos, kmerLength);
    if (kmer != -1) {
      kmers.insert(kmer);
    }
  }
  return kmers;

}

double kmerDistance(vector<Seq*> &a, vector<Seq*> &b, int kmerLength) {
  set<uint32_t> a_kmers = getKmers(a, kmerLength, 10000);
  set<uint32_t> b_kmers = getKmers(b, kmerLength, 10000);
  set<uint32_t> a_u_b;
  set<uint32_t> a_int_b;

  set_union(a_kmers.begin(), a_kmers.end(), b_kmers.begin(), b_kmers.end(), std::inserter(a_u_b, a_u_b.begin()));
  set_intersection(a_kmers.begin(),a_kmers.end(),b_kmers.begin(), b_kmers.end(), std::inserter(a_int_b, a_int_b.begin()));
  return (double) (a_int_b.size()) / (double) (a_u_b.size());
}

void joinGroups(vector<vector<Seq*> > &groups, int kmerLength, double similarityThreshold) {
  for (int i = 0; i < groups.size(); i++) {
    for (int j = 0; j < i; j++) {
      if (groups[i].size() <= 1 || groups[j].size() <= 1) continue;
      double similarity = kmerDistance(groups[i], groups[j], kmerLength);
      cerr << "Similarity between groups " << i << " and " << j << " = " << similarity << endl;
      if (similarity > similarityThreshold) {
        cerr << "Joining group " << i << " to group " << j << endl;
        for (vector<Seq*>::iterator it = groups[i].begin(); it != groups[i].end(); it++) {
          groups[j].push_back(*it);
        } 
        groups[i] = vector<Seq*>();
      }
    }

  }
}

vector<Seq*> liftoverRepeatAnnotations(vector<Seq*> repeats, const hal::Genome *source, const hal::Genome *target) {

}

/*
   LPOSeq_T *buildSequenceGraph(vector<Sequence*> &sequences, char *matrixFilename) {
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
