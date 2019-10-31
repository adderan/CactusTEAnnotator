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

static void testHeaviestPath(CuTest *testCase) {
	stPinchThreadSet *graph = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(graph, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(graph, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 1);
	stPinchBlock *block1 = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50));
	stPinchBlock *block2 = stPinchSegment_getBlock(stPinchThread_getSegment(thread2, 70));
	stPinchEnd *end1 = stPinchEnd_construct(block1, 0);
	stPinchEnd *end2 = stPinchEnd_construct(block2, 1);

	stList *subsequences = stPinchEnd_getSubSequenceLengthsConnectingEnds(end1, end2);

	fprintf(stderr, "Number of connecting adjacencies: %ld\n", stList_length(subsequences));
	stIntTuple *adj = stList_pop(subsequences);
	fprintf(stderr, "Tuple length: %ld\n", stIntTuple_length(adj));
	fprintf(stderr, "tuple: %ld\n", stIntTuple_get(adj, 0));

}

int main(int argc, char **argv) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testDirectedWalk);
	SUITE_ADD_TEST(suite, testAcyclic);
	SUITE_ADD_TEST(suite, testHeaviestPath);
	CuSuiteRun(suite);
}
