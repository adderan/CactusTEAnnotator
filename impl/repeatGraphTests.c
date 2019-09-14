#include "repeatGraphs.h"

void testDirectedWalk() {
	fprintf(stderr, "Testing directedWalk\n");
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);


	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);

	stPinchSegment *seg1 = stPinchThread_getSegment(thread1, 25);
	stPinchSegment *seg2 = stPinchThread_getSegment(thread2, 25);

	//moving toward 5' end
	assert(!directedWalk(seg1, seg2, _3PRIME));
	//moving toward 3' end
	assert(!directedWalk(seg1, seg2, _5PRIME));

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);
	stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 10));
	stPinchEnd blockALeftEnd = stPinchEnd_constructStatic(blockA, _5PRIME);
	assert(stPinchEnd_getNumberOfConnectedPinchEnds(&blockALeftEnd) == 0);

	assert(!directedWalk(seg1, seg2, _3PRIME));
	assert(!directedWalk(seg1, seg2, _5PRIME));


	stPinchThreadSet_destruct(threadSet);
	threadSet = stPinchThreadSet_construct();
	stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	thread1 = stPinchThreadSet_getThread(threadSet, 1);
	thread2 = stPinchThreadSet_getThread(threadSet, 2);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 0);

	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	//should be able to go backwards on thread 1, traverse
	//the reverse block, and then forwards on thread 2
	assert(directedWalk(seg1, seg2, _5PRIME));
	assert(!directedWalk(seg1, seg2, _3PRIME));
	assert(pinchCreatesCycle(seg1, seg2, 0));
	assert(!pinchCreatesCycle(seg1, seg2, 1));
	fprintf(stderr, "Passed\n");
}

void testOrdering() {
	fprintf(stderr, "Testing graph ordering\n");
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stPinchThread *thread1 = stPinchThreadSet_addThread(threadSet, 1, 0, 100);
	stPinchThread *thread2 = stPinchThreadSet_addThread(threadSet, 2, 0, 100);

	stPinchThread_pinch(thread1, thread2, 10, 10, 10, 1);
	stPinchThread_pinch(thread1, thread2, 40, 40, 10, 0);

	assert(getOrdering(threadSet));

	stPinchThread_pinch(thread1, thread2, 70, 70, 10, 0);
	assert(!getOrdering(threadSet));

	fprintf(stderr, "Passed\n");
}

int main(int argc, char **argv) {
	testDirectedWalk();
	testOrdering();
}
