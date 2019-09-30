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

	assert(!directedWalk(seg1, seg2, _3PRIME));
	assert(!directedWalk(seg1, seg2, _5PRIME));

	stPinchThread_pinch(thread1, thread2, 50, 50, 10, 1);
	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	assert(!directedWalk(seg1, seg2, _5PRIME));
	assert(!directedWalk(seg1, seg2, _3PRIME));


	stPinchThread_pinch(thread1, thread2, 80, 10, 10, 1);

	seg1 = stPinchThread_getSegment(thread1, 25);
	seg2 = stPinchThread_getSegment(thread2, 25);

	assert(directedWalk(seg1, seg2, _3PRIME));
	assert(directedWalk(seg1, seg2, _5PRIME));
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
