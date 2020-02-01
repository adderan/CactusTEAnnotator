#include "repeatGraphs.h"
#include "sonLib.h"
#include "bioioC.h"
#include <getopt.h>

#include "stPinchIterator.h"

int main(int argc, char **argv) {
	char *sequencesFilename = NULL;
	char *alignmentsFilename = NULL;
	char *namePrefix = "";
	int64_t minConsensusLength= 30;
	int64_t gapPenalty = 1;
	char *gvizDebugFilename = NULL;
	double minConsensusDegree = 10.0;
	while (1) {
        static struct option long_options[] = {
            { "sequences", required_argument, 0, 'a' }, 
			{ "alignments", required_argument, 0, 'b'},
			{ "namePrefix", required_argument, 0, 'c'},
			{ "minConsensusLength", required_argument, 0, 'd'},
			{ "gapPenalty", required_argument, 0, 'e'},
			{ "gvizDebugFilename", required_argument, 0, 'f'},
			{ "minConsensusDegree", required_argument, 0, 'g'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:f:g:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                sequencesFilename = strdup(optarg);
                break;
			case 'b':
				alignmentsFilename = strdup(optarg);
				break;
			case 'c':
				namePrefix = strdup(optarg);
				break;
			case 'd':
				sscanf(optarg, "%ld", &minConsensusLength);
				break;
			case 'e':
				sscanf(optarg, "%ld", &gapPenalty);
				break;
			case 'f':
				gvizDebugFilename = strdup(optarg);
				break;
			case 'g':
				sscanf(optarg, "%lf", &minConsensusDegree);
				break;
            default:
                return 1;
        }
    }

	struct List *seqs = constructEmptyList(0, NULL);
	struct List *seqLengths = constructEmptyList(0, free);
	struct List *headers = constructEmptyList(0, free);
	FILE *sequencesFile = fopen(sequencesFilename, "r");
	fastaRead(sequencesFile, seqs, seqLengths, headers);
	fclose(sequencesFile);


	stPinchThreadSet *graph = stPinchThreadSet_construct();

	stHash *sequences = stHash_construct();
	
	for (int64_t i = 0; i < seqs->length; i++) {
		char *threadName = headers->list[i];
		int64_t threadID;
		sscanf(threadName, "%ld", &threadID);
		int64_t seqLength = strlen(seqs->list[i]);
		stHash_insert(sequences, (void*) threadID, seqs->list[i]);
		stPinchThreadSet_addThread(graph, threadID, 0, seqLength);
	}
	stPinchIterator *pinchIterator = 
		stPinchIterator_constructFromFile(alignmentsFilename);

	stPinch *pinch = NULL;
	while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		assert(stHash_search(sequences, (void*)pinch->name1));
		assert(stHash_search(sequences, (void*)pinch->name2));
		stPinchThread *thread1 = stPinchThreadSet_getThread(graph, pinch->name1);
		stPinchThread *thread2 = stPinchThreadSet_getThread(graph, pinch->name2);
		assert(stPinchThread_getName(thread1) == pinch->name1);
		assert(stPinchThread_getName(thread2) == pinch->name2);
		//if (pinch->length < 10) continue;
		if (pinch->strand == 0) continue;
		stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, 
			pinch->length, pinch->strand, singleCopyFilterFn);
	}


	if (gvizDebugFilename) {
		printBiedgedGraph(graph, gvizDebugFilename);
	}

	fprintf(stderr, "Graph has %ld blocks\n", stPinchThreadSet_getTotalBlockNumber(graph));
	
	stList *blockOrdering = getBlockOrdering(graph);
	assert(stList_length(blockOrdering) == stPinchThreadSet_getTotalBlockNumber(graph));

	int64_t pathScore;
	int64_t consensusNum = 0;

	int64_t N = stList_length(blockOrdering);
	int64_t *scores = calloc(N, sizeof(int64_t));
	int64_t *directions = calloc(N, sizeof(int64_t));
	getHeaviestPathScores(blockOrdering, 1, scores, directions);

	while (true) {

		stList *path = tracebackHeaviestPath(blockOrdering, scores, directions, &pathScore);
		if (!path) break;


		char *consensusSeq = getConsensusSequence(path, sequences);

		stList_destruct(path);

		double consensusDegree = (double)pathScore/(double)strlen(consensusSeq);

		if ((strlen(consensusSeq) < minConsensusLength) || (consensusDegree < minConsensusDegree)) {
			free(consensusSeq);
			continue;
		}

		fprintf(stdout, ">%s_consensus_%ld length=%ld score=%ld\n", namePrefix, consensusNum, strlen(consensusSeq), pathScore);
		fprintf(stdout, "%s\n", consensusSeq);
		free(consensusSeq);

		consensusNum++;
	}
	destructList(seqs);
	destructList(headers);
	destructList(seqLengths);

	stList_destruct(blockOrdering);
	stPinchThreadSet_destruct(graph);
}
