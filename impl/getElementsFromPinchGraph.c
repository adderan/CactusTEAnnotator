#include "repeatGraphs.h"

int main(int argc, char **argv) {
	char *sequencesFilename = argv[1];
	char *alignmentsFilename = argv[2];
	char *gvizDebugFilename = argv[3];

	stPinchThreadSet *threadSet = buildRepeatGraph(sequencesFilename, alignmentsFilename, gvizDebugFilename);

	stPinchThreadSet_destruct(threadSet);
}
