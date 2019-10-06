#include "repeatGraphs.h"
#include "sonLibGlobalsTest.h"

void testDirectedWalk(CuTest *testCase) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);


	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);

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

static void testOrdering(CuTest *testCase) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);
	stPinchThread_pinch(thread1, thread2, 40, 40, 10, 1);

	CuAssertTrue(testCase, getOrdering(threadSet) != NULL);

	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 0);
	CuAssertTrue(testCase, getOrdering(threadSet) == NULL);
}

int main(int argc, char **argv) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testDirectedWalk);
	SUITE_ADD_TEST(suite, testOrdering);
	CuSuiteRun(suite);
}
