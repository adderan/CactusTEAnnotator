#include "repeatGraphs.h"

void testDirectedWalk() {
	fprintf(stderr, "Testing directedWalk\n");
	assert((true ^ true) == false);
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);


	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);

	stPinchSegment *seg1 = stPinchThread_getSegment(thread1, 25);
	stPinchSegment *seg2 = stPinchThread_getSegment(thread2, 25);

	assert(pinchCreatesCycle(seg1, seg2, 0));
	assert(!pinchCreatesCycle(seg1, seg2, 1));

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);
	stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 10));
	stPinchEnd blockALeftEnd = stPinchEnd_constructStatic(blockA, _5PRIME);
	assert(stPinchEnd_getNumberOfConnectedPinchEnds(&blockALeftEnd) == 0);

	assert(pinchCreatesCycle(seg1, seg2, 0));
	assert(!pinchCreatesCycle(seg1, seg2, 1));


	stPinchThreadSet_destruct(threadSet);
	threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	thread1 = stPinchThreadSet_getThread(threadSet, 1);
	thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 0);

	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	assert(pinchCreatesCycle(seg1, seg2, 0));
	assert(pinchCreatesCycle(seg1, seg2, 1));
	fprintf(stderr, "Passed\n");
}

void testOrdering() {
	fprintf(stderr, "Testing graph ordering\n");
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);
	stPinchThread_pinch(thread1, thread2, 40, 40, 10, 1);

	assert(getOrdering(threadSet, NULL));

	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 0);
	assert(!getOrdering(threadSet, NULL));

	fprintf(stderr, "Passed\n");
}

int main(int argc, char **argv) {
	testDirectedWalk();
	testOrdering();
}
