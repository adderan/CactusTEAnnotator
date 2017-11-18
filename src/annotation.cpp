#include <vector>
#include <string>
#include <map>
#include "insertions.h"
#include "annotation.h"
#include "clustering.h"


using namespace std;


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

double insertionDistance(Insertion *a, Insertion *b) {
  return kmerDistance(a->seq, b->seq);
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
  
