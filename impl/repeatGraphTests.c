#include "repeatGraphs.h"
#include "sonLibGlobalsTest.h"

void testDirectedWalk(CuTest *testCase) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 0, 0, 100);
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);


	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 0);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 1);

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
}

static void testAcyclic(CuTest *testCase) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);
	stPinchThread_pinch(thread1, thread2, 40, 40, 10, 1);

	CuAssertTrue(testCase, graphIsAcyclic(threadSet));

	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 0);
	CuAssertTrue(testCase, !graphIsAcyclic(threadSet));
}

static void testPoGraph(CuTest *testCase) {
    stPinchThreadSet *graph = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(graph, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(graph, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 1);

    stPinchBlock *block1 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50));
	stPinchBlock *block2 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 70));

	stList *poGraph = getPartialOrderGraph(graph);

    CuAssertTrue(testCase, stList_length(poGraph) == 2);
	PartialOrderNode *node1 = stList_get(poGraph, 0);
	PartialOrderNode *node2 = stList_get(poGraph, 1);
	CuAssertTrue(testCase, node1->data == block1);
	CuAssertTrue(testCase, node2->data == block2);
	CuAssertTrue(testCase, stList_length(node2->incomingNodes) == 1);

}

static void testHeaviestPath(CuTest *testCase) {
	stPinchThreadSet *graph = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(graph, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(graph, 2, 0, 100);
	//stPinchThread *thread3 = stPinchThreadSet_addThread(graph, 3, 0, 100);

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 1);

	//stPinchThread_pinch(thread1, thread3, 90, 90, 5, 1);

    stPinchBlock *block1 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50));
	stPinchBlock *block2 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 70));

	//stPinchEnd *end1 = stPinchEnd_construct(block1, _3PRIME);
	//stPinchEnd *end2 = stPinchEnd_construct(block2, _5PRIME);

	stList *poGraph = getPartialOrderGraph(graph);

	stList *heaviestPath = getHeaviestPath(poGraph, 1.0);
	CuAssertTrue(testCase, stList_length(heaviestPath) == 2);

	PartialOrderNode *node1 = stList_get(heaviestPath, 0);
	PartialOrderNode *node2 = stList_get(heaviestPath, 1);

	CuAssertTrue(testCase, node1->data == block1);
	CuAssertTrue(testCase, node2->data == block2);



}

static void testConnectingThreads(CuTest *testCase) {
	stPinchThreadSet *graph = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(graph, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(graph, 2, 0, 100);

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
}

int main(int argc, char **argv) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testConnectingThreads);
	SUITE_ADD_TEST(suite, testDirectedWalk);
	SUITE_ADD_TEST(suite, testAcyclic);
    SUITE_ADD_TEST(suite, testPoGraph);
	SUITE_ADD_TEST(suite, testHeaviestPath);
	CuSuiteRun(suite);
}
