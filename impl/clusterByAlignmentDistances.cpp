#include <vector>
#include <tuple>
#include <set>
#include <stdlib.h>
#include <iostream>
#include "PairwiseDistances.h"

extern "C" {
#include "sonLib.h"
#include "bioioC.h"
#include "commonC.h"
#include "lpo.h"
#include "seq_util.h"
}

using namespace std;

int main(int argc, char **argv) {
	FILE *lpoFile = fopen(argv[1], "r");
	LPOSequence_T *graph = read_lpo(lpoFile);
	fclose(lpoFile);

    double distanceThreshold;
    sscanf(argv[2], "%lf", &distanceThreshold);


	LPOLetter_T *seq = graph->letter;

	int **nAligned = (int**)calloc(sizeof(int*), graph->nsource_seq);
	for (int i = 0; i < graph->nsource_seq; i++) {
		nAligned[i] = (int*)calloc(sizeof(int), graph->nsource_seq);
	}

	//count all aligned positions
	for (int i = 0; i < graph->length; i++) {
		for (LPOLetterSource_T *seq_a = &seq[i].source; seq_a != NULL; seq_a = seq_a->more) {
			for (LPOLetterSource_T *seq_b = &seq[i].source; seq_b != seq_a; seq_b = seq_b->more) {
				int a = seq_a->iseq;
				int b = seq_b->iseq;
				if (a < b) {
					int temp = a;
					a = b;
					b = temp;
				}
				nAligned[a][b]++;
			}
		}
	}

    vector<tuple<int, int, double> > distances;
	for (int i = 0; i < graph->nsource_seq; i++) {
		for (int j = 0; j < i; j++) {
			int N = nAligned[i][j];
			int len_i = graph->source_seq[i].length;
			int len_j = graph->source_seq[j].length;
			double p = 2*N/(double)(len_i + len_j);

            //jukes-cantor distance
            cerr << "p = " << p << endl;
            double d;
            if (p <= 0.25) {
                d = 1000.0;
            }
            else {
                d = (-3.0/4.0)*log(1.0 - (4.0/3.0)*(1.0 - p));
            }
            cerr << "Found jukes cantor distance " << d << endl;

            distances.push_back(make_tuple(i, j, d));
		}
	}
    set<set<long> > clusters = buildClusters(distances, graph->nsource_seq, distanceThreshold);

    for (auto &t: clusters) {
        for (auto &seqNum: t) {
            printf("%s ", graph->source_seq[seqNum].name);
        }
        printf("\n\n");
    }
}
