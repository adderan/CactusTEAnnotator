#include <stdlib.h>
#include "sonLib.h"
#include "sonLibGlobalsTest.h"
#include "stPinchGraphs.h"
#include "stPinchPhylogeny.h"
#include "stPhylogeny.h"

#include "Consensus.h"


stPinchThreadSet *graph;
stPinchThread *thread1;
stPinchThread *thread2;

void setup() {
	graph = stPinchThreadSet_construct();
	thread1 = stPinchThreadSet_addThread(graph, 1, 0, 100);
	thread2 = stPinchThreadSet_addThread(graph, 2, 0, 100);
}

void teardown() {
	stPinchThreadSet_destruct(graph);
}

void testDirectedWalk(CuTest *testCase) {
	setup();

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);

	stPinchSegment *seg1 = stPinchThread_getSegment(thread1, 25);
	stPinchSegment *seg2 = stPinchThread_getSegment(thread2, 25);

	CuAssertTrue(testCase, directedWalk(seg1, seg2, _3PRIME));
	CuAssertTrue(testCase, directedWalk(seg1, seg2, _5PRIME));

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	CuAssertTrue(testCase, directedWalk(seg1, seg2, _5PRIME));
	CuAssertTrue(testCase, directedWalk(seg1, seg2, _3PRIME));


	stPinchThread_pinch(thread1, thread2, 80, 10, 10, 1);

	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	CuAssertTrue(testCase, directedWalk(seg1, seg2, _3PRIME));
	CuAssertTrue(testCase, directedWalk(seg1, seg2, _5PRIME));

	teardown();
}

static void testAcyclic(CuTest *testCase) {
	setup();

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);
	stPinchThread_pinch(thread1, thread2, 40, 40, 10, 1);

	CuAssertTrue(testCase, graphIsAcyclic(graph));

	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 0);
	CuAssertTrue(testCase, !graphIsAcyclic(graph));

	teardown();
}

static void testConnectingThreads(CuTest *testCase) {
	setup();

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 1);
	stPinchThread_pinch(thread1, thread2, 90, 90, 5, 1);

	stPinchBlock *block1 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50));
	stPinchBlock *block2 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 70));
	stPinchBlock *block3 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 90));
	
	stPinchEnd *end1 = stPinchEnd_construct(block1, _3PRIME);
	stPinchEnd *end2 = stPinchEnd_construct(block2, _5PRIME);
	stPinchEnd *end3 = stPinchEnd_construct(block3, _5PRIME);

	CuAssertTrue(testCase, stSortedSet_size(getConnectingThreads(end1, end2)) == 2);
	CuAssertTrue(testCase, stSortedSet_size(getConnectingThreads(end1, end3)) == 0);

	CuAssertTrue(testCase, stPinchEnd_getNumberOfConnectedPinchEnds(end1) == 1);
	CuAssertTrue(testCase, stList_length(stPinchEnd_getSubSequenceLengthsConnectingEnds(end1, end2)) == 2);
	
	//why?
	CuAssertTrue(testCase, stList_length(stPinchEnd_getSubSequenceLengthsConnectingEnds(end1, end3)) == 2);

	teardown();
}

static void testBlockOrdering(CuTest *testCase) {
	setup();
	stPinchThread *thread3 = stPinchThreadSet_addThread(graph, 3, 0, 100);
	stPinchThread *thread4 = stPinchThreadSet_addThread(graph, 4, 0, 100);

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 1);
	stPinchThread_pinch(thread3, thread4, 50, 50, 10, 1);
	stPinchThread_pinch(thread3, thread4, 70, 70, 10, 1);


	stList *ordering = getBlockOrdering(graph);
	assert(stList_length(ordering) == 4);
	stList_destruct(ordering);

	teardown();
}

char *generateSeq(int64_t length) {
	char *seq = malloc(sizeof(char)*length);
	char bases[4] = {'A', 'C', 'T', 'G'};
	for (int64_t i = 0; i < length; i++) {
		int64_t baseIndex = (int64_t) (((double)rand())/RAND_MAX * 4);
		seq[i] = bases[baseIndex];
	}
	return seq;
}

static void testGetConsensus(CuTest *testCase) {
	setup();
	char *seq1 = generateSeq(100);
	char *seq2 = strdup(seq1);
	stHash *sequences = stHash_construct();
	stHash_insert(sequences, (void*)1, seq1);
	stHash_insert(sequences, (void*)2, seq2);

	stPinchThread_pinch(thread1, thread2, 20, 20, 10, 1);
	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);

	stPinchSegment *seg1 = stPinchThread_getSegment(thread1, 25);
	stPinchBlock *block1 = stPinchSegment_getBlock(seg1);
	stPinchEnd *end1 = stPinchEnd_construct(block1, _5PRIME);

	stPinchSegment *seg2 = stPinchThread_getSegment(thread1, 55);
	stPinchBlock *block2 = stPinchSegment_getBlock(seg2);
	stPinchEnd *end2 = stPinchEnd_construct(block2, _5PRIME);

	stList *path = stList_construct();
	stList_append(path, end1);
	stList_append(path, end2);

	char *consensusSeq = getConsensusSequence(path, sequences);

	assert(strlen(consensusSeq) == 40);
	assert(strncmp(consensusSeq, seq1 + 20, 10) == 0);

	free(consensusSeq);
	free(seq1);
	free(seq2);
	teardown();
}

static void testTreeBuilding(CuTest *testCase) {
	setup();
	stPinchThread *thread3 = stPinchThreadSet_addThread(graph, 3, 0, 100);
	char *seq1 = generateSeq(100);
	char *seq2 = generateSeq(100);
	char *seq3 = generateSeq(100);

	stPinchThread_pinch(thread1, thread2, 20, 20, 10, 1);
	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread2, thread3, 50, 50, 10, 1);
	//stPinchBlock *block1 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 25));
	//stPinchBlock *block2 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 55));

	
	stList *chain = stList_construct();

	stHash *pinchThreadsToStrings = stHash_construct();
	stHash_insert(pinchThreadsToStrings, thread1, seq1);
	stHash_insert(pinchThreadsToStrings, thread2, seq2);
	stHash_insert(pinchThreadsToStrings, thread3, seq3);
	stList *featureBlocks = stFeatureBlock_getContextualFeatureBlocksForChainedBlocks(chain, 10, 10, true, true, pinchThreadsToStrings);

	printf("Blocks: %ld\n", stList_length(featureBlocks));
	stList_destruct(featureBlocks);

	teardown();
}

int main(int argc, char **argv) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testConnectingThreads);
	SUITE_ADD_TEST(suite, testDirectedWalk);
	SUITE_ADD_TEST(suite, testAcyclic);
	SUITE_ADD_TEST(suite, testBlockOrdering);
	SUITE_ADD_TEST(suite, testGetConsensus);
	SUITE_ADD_TEST(suite, testTreeBuilding);
	CuSuiteRun(suite);
}
