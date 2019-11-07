#include <stdio.h>
#include "bioioC.h"
#include "sonLibTypes.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"
#include "stPinchIterator.h"

#include "repeatGraphs.h"

stSortedSet *getThreads(stPinchSegment *segment) {
	stSortedSet *threads = stSortedSet_construct();
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (block) {
		stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
		stPinchSegment *otherSeg = NULL;
		while((otherSeg = stPinchBlockIt_getNext(&blockIt)) != NULL) {
			stSortedSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(otherSeg)));
		}
	}
	else {
		stSortedSet_insert(threads, (void*) stPinchThread_getName(stPinchSegment_getThread(segment)));
	}
	return threads;
}

bool singleCopyFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (stPinchSegment_getThread(seg1) == stPinchSegment_getThread(seg2)) return true;
	bool filter = false;

	stSortedSet *threads1 = getThreads(seg1);
	stSortedSet *threads2 = getThreads(seg2);

	stSortedSet *intersect = stSortedSet_getIntersection(threads1, threads2);
	if (stSortedSet_size(intersect) > 0) filter = true;
	stSortedSet_destruct(threads1);
	stSortedSet_destruct(threads2);
	stSortedSet_destruct(intersect);
	return filter;
}

/*
void printBiedgedGraph(stPinchThreadSet *graph, char *gvizFilename) {
	stList *ordering = getOrdering(graph);
	stHash *blockToEnd = stHash_construct();
	stListIterator *orderingIt = stList_getIterator(ordering);
	stPinchEnd *end;
	while((end = stList_getNext(orderingIt))) {
		stHash_insert(blockToEnd, stPinchEnd_getBlock(end), end);
	}

	FILE *gvizFile = fopen(gvizFilename, "w");
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	stPinchBlock *block;

	fprintf(gvizFile, "strict graph {\n");
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt))) {
		stPinchEnd *blackEnd = stHash_search(blockToEnd, block);
		stPinchBlock *block = stPinchEnd_getBlock(blackEnd);
		stPinchEnd *redEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(blackEnd));

		//end nodes
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=grey];\n", stPinchEnd_hashFn(blackEnd));
		fprintf(gvizFile, "\t%ld[style=filled fillcolor=red];\n", stPinchEnd_hashFn(redEnd));


		//block edge
		fprintf(gvizFile, "\t%ld -- %ld[dir=none color=blue penwidth=10];\n", stPinchEnd_hashFn(blackEnd), stPinchEnd_hashFn(redEnd));

		//fprintf(gvizFile, "{rank = same; %ld; %ld;}\n", stPinchEnd_hashFn(end1), stPinchEnd_hashFn(end2));

		//near adjacencies
		stSet *nearAdjacent = stPinchEnd_getConnectedPinchEnds(blackEnd);
		stSetIterator *nearAdjacentIt = stSet_getIterator(nearAdjacent);
		stPinchEnd *adjEnd;
		while((adjEnd = stSet_getNext(nearAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(blackEnd), stPinchEnd_hashFn(adjEnd));

		}
		stSet_destructIterator(nearAdjacentIt);

		stSet *farAdjacent = stPinchEnd_getConnectedPinchEnds(redEnd);
		stSetIterator *farAdjacentIt = stSet_getIterator(farAdjacent);
		while((adjEnd = stSet_getNext(farAdjacentIt))) {
			fprintf(gvizFile, "\t%ld -- %ld;\n", stPinchEnd_hashFn(redEnd), stPinchEnd_hashFn(adjEnd));
		}
	}
	fprintf(gvizFile, "}\n");

	stList_destruct(ordering);
	fclose(gvizFile);
}
*/

stList *getPartialOrderGraph(stPinchThreadSet *graph) {
	stList *stack = stList_construct();

	stSet *seen = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void*) stPinchEnd_destruct);

	stHash *blockIndex = stHash_construct();

	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	stPinchBlock *startBlock;
	int64_t nodeID = 0;

	stList *poGraph = stList_construct();

	while((startBlock = stPinchThreadSetBlockIt_getNext(&blockIt))) {

		stPinchEnd *startEnd = stPinchEnd_construct(startBlock, _5PRIME);

		if (stSet_search(seen, startEnd)) continue;

		stList_append(stack, startEnd);
		stPinchBlock *block;
		stPinchEnd *end;
		while(stList_length(stack) > 0) {
			end = stList_pop(stack);
			block = stPinchEnd_getBlock(end);

			//this only has to be true because negative pinches
			//are discarded
			assert(stPinchEnd_getOrientation(end) == _5PRIME);

			if (stSet_search(seen, end)) continue;

			stSet *incomingAdjacencies = stPinchEnd_getConnectedPinchEnds(end);
			stSetIterator *incomingEndIt = stSet_getIterator(incomingAdjacencies);
			stPinchEnd *incomingRedEnd;
			stList *unseenPredecessors = stList_construct();
			while((incomingRedEnd = stSet_getNext(incomingEndIt)) != NULL) {
				stPinchEnd *incomingBlackEnd = stPinchEnd_construct(stPinchEnd_getBlock(incomingRedEnd), !stPinchEnd_getOrientation(incomingRedEnd));
				if (!stSet_search(seen, incomingBlackEnd)) 
					stList_append(unseenPredecessors, incomingBlackEnd);
			}
			stSet_destructIterator(incomingEndIt);

			if (stList_length(unseenPredecessors) > 0) {
				//re-visit this block after visiting all predecessors
				stList_append(stack, end);
				stList_appendAll(stack, unseenPredecessors);
			}

			else {
				stSet_insert(seen, end);
				//seen all the predecessors, so add this block
				//to the partial order graph
				PONode *node = malloc(sizeof(PONode));
				node->nodeID = nodeID;
				node->weight = stPinchBlock_getDegree(block) * stPinchBlock_getLength(block);
				node->incomingNodes = malloc(sizeof(int64_t)*stPinchEnd_getNumberOfConnectedPinchEnds(end));
				node->nIncomingNodes = stPinchEnd_getNumberOfConnectedPinchEnds(end);
				stList *incomingEndList = stSet_getList(incomingAdjacencies);
				for (int k = 0; k < stList_length(incomingEndList); k++) {
					stPinchBlock *incomingBlock = stPinchEnd_getBlock(stList_get(incomingEndList, k));
					//We should have already encountered this block
					node->incomingNodes[k] = (int64_t) stHash_search(blockIndex, incomingBlock);
				}
				node->data = end;
				stHash_insert(blockIndex, block, (void*)nodeID);
				stList_append(poGraph, node);
				nodeID++;

				//visit next blocks
				stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
				stSet *outgoingAdjacencies = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
				stSetIterator *outgoingEndIt = stSet_getIterator(outgoingAdjacencies);
				stPinchEnd *outgoingEnd;
				while((outgoingEnd = stSet_getNext(outgoingEndIt))) {
					if (!stSet_search(seen, outgoingEnd)) stList_append(stack, outgoingEnd);
				}
				stPinchEnd_destruct(oppositeEnd);
			}
			stList_destruct(unseenPredecessors);
		}
	}
	return poGraph;
}

/*If the graph is not acyclic, returns NULL. If the graph is acyclic, 
 * return a list of edges in the DAG corresponding to this pinch graph,
 * in ordering implied by traversing the graph. */
bool graphIsAcyclic(stPinchThreadSet *graph) {
	stSet *seen = stSet_construct();

	stHash *coloring = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void*)) stPinchEnd_destruct, NULL);

	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	stPinchBlock *startBlock;
	while((startBlock = stPinchThreadSetBlockIt_getNext(&blockIt))) {
		if (stSet_search(seen, startBlock)) continue;
		stList *stack = stList_construct();

		stPinchEnd *leftEnd = stPinchEnd_construct(startBlock, _5PRIME);

		stList_append(stack, leftEnd);
		stHash_insert(coloring, leftEnd, BLACK);

		stPinchBlock *block;
		stPinchEnd *end;
		while (stList_length(stack) > 0) {
			end = stList_pop(stack);
			block = stPinchEnd_getBlock(end);

			stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));

			assert(stHash_search(coloring, end));
			if (stHash_search(coloring, end) ==
					stHash_search(coloring, oppositeEnd)) {
				//encountered cycle
				stList_destruct(stack);
				return false;
			}
			stHash_insert(coloring, oppositeEnd, OP(stHash_search(coloring, end)));

			if (!stSet_search(seen, block)) {
				stSet_insert(seen, block);
				//fprintf(stderr, "Found %ld connected ends on near side\n", stSet_size(nearAdjacentEnds));
				stSet *incomingEnds = stPinchEnd_getConnectedPinchEnds(end);
				stSetIterator *incomingEndIt = stSet_getIterator(incomingEnds);
				stPinchEnd *adjEnd;
				while((adjEnd = stSet_getNext(incomingEndIt))) {
					void *colorToTagEnd = OP(stHash_search(coloring, end));
					if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
						fprintf(stderr, "About to encounter cycle\n");
						return NULL;
					}
					stHash_insert(coloring, adjEnd, colorToTagEnd);
					stList_append(stack, adjEnd);
				}

				stSet_destructIterator(incomingEndIt);
				stSet *outgoingEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
				stSetIterator *outgoingEndIt = stSet_getIterator(outgoingEnds);
				outgoingEndIt = stSet_getIterator(outgoingEnds);
				while((adjEnd = stSet_getNext(outgoingEndIt))) {
					void *colorToTagEnd = OP(stHash_search(coloring, oppositeEnd));
					if (stHash_search(coloring, adjEnd) == OP(colorToTagEnd)) {
						fprintf(stderr, "About to encounter cycle\n");
						return NULL;
					}
					stHash_insert(coloring, adjEnd, colorToTagEnd);
					stList_append(stack, adjEnd);
				}
				stSet_destructIterator(outgoingEndIt);
			}
		}
		stList_destruct(stack);
	}
	stSet_destruct(seen);
	return coloring;
}

stPinchBlock *getFirstBlock(stPinchThread *thread) {
	stPinchSegment *segment = stPinchThread_getFirst(thread);
	while(segment && (stPinchSegment_getBlock(segment) == NULL)) {
		segment = stPinchSegment_get3Prime(segment);
	}
	if (!segment) {
		segment = stPinchThread_getFirst(thread);
		while(segment && (stPinchSegment_getBlock(segment) == NULL)) {
			segment = stPinchSegment_get3Prime(segment);
		}
	}
	if (!segment) return NULL;
	return stPinchSegment_getBlock(segment);
}

stPinchEnd *getAdjacentEnd(stPinchSegment *segment, bool direction) {
	if (stPinchSegment_getBlock(segment)) {
		return stPinchEnd_construct(stPinchSegment_getBlock(segment), direction);
	}
	while (segment && (stPinchSegment_getBlock(segment) == NULL)) {
		if (direction == _3PRIME) {
			segment = stPinchSegment_get3Prime(segment);
		}
		else {
			assert(direction == _5PRIME);
			segment = stPinchSegment_get5Prime(segment);
		}
	}
	if (!segment) return NULL;
	stPinchBlock *block = stPinchSegment_getBlock(segment);
	if (!block) return NULL;

	bool orientation = direction ^ stPinchSegment_getBlockOrientation(segment);
	return stPinchEnd_construct(block, orientation);
}


//inefficient function for checking whether there is a directed
//walk (a path of alternating block and adjacency edges) between
//two segments.
bool directedWalk(stPinchSegment *seg1, stPinchSegment *seg2, bool startDirection) {
	stPinchBlock *block;
	stSet *seenBlocks = stSet_construct();
	stList *stack = stList_construct();
	stList *directions = stList_construct();
	stPinchSegment *segment = seg1;
	stList_append(stack, segment);
	stList_append(directions, (void*)startDirection);
	bool curDirection = startDirection;
	bool walk = false;
	while(segment || stList_length(stack) > 0) {
		if (!segment) {
			segment = stList_pop(stack);
			curDirection = stList_pop(directions);
		}
		if (segment == seg2) {
			walk = true;
			goto out;
		}
		block = stPinchSegment_getBlock(segment);
		if (block && !stSet_search(seenBlocks, block)) {
			stSet_insert(seenBlocks, block);
			stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
			stPinchSegment *nextSeg;
			while((nextSeg = stPinchBlockIt_getNext(&blockIt))) {
				if (nextSeg != segment) {
					stList_append(stack, nextSeg);

					//switch directions if the next segment 
					//or the current segment is reversed relative 
					//to the block
					bool nextSegDirection = curDirection 
						^ (stPinchSegment_getBlockOrientation(nextSeg)==0)
						^ (stPinchSegment_getBlockOrientation(segment)==0);
					stList_append(directions, (void*)nextSegDirection);
				}
			}

		}
		segment = (curDirection == _5PRIME) ? 
			stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);

	}
out:
	stList_destruct(stack);
	stList_destruct(directions);
	stSet_destruct(seenBlocks);
	return walk;
}

bool acyclicFilterFn(stPinchSegment *seg1, stPinchSegment *seg2) {
	if (singleCopyFilterFn(seg1, seg2)) return true;
	return (directedWalk(seg1, seg2, _3PRIME) || directedWalk(seg1, seg2, _5PRIME));
}


stPinchThreadSet *initializeGraph(stHash *sequences) {
	stPinchThreadSet *threadSet = stPinchThreadSet_construct();
	stHashIterator *sequencesIt = stHash_getIterator(sequences);
	void *threadName;
	while((threadName = stHash_getNext(sequencesIt))) {
		char *sequence = stHash_search(sequences, threadName);
		stPinchThreadSet_addThread(threadSet, (int64_t) threadName, 0, strlen(sequence));
	}
	return threadSet;
}

void applyPinch(stPinchThreadSet *threadSet, stPinch *pinch) {
	if (pinch->strand == 0) return;
	stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
	stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
	stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, acyclicFilterFn);
}

stPinchThreadSet *buildRepeatGraph(stHash *sequences, char *alignmentsFilename) {

	stPinchThreadSet *threadSet = initializeGraph(sequences);
#ifndef NDEBUG
	//maintain another copy of the graph to keep track of what it looked
	//like before applying a bad pinch
	stPinchThreadSet *threadSet_debug = initializeGraph(sequences);
#endif

	stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(alignmentsFilename);
	stPinch *pinch = NULL;

	while((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
		applyPinch(threadSet, pinch);

#ifndef NDEBUG
		/*
		   stList *ordering = getOrdering(threadSet);
		   if (!ordering) {
		   fprintf(stderr, "Pinch created cycle.\n");
		   printBiedgedGraph(threadSet, "after_bad_pinch.gvz");
		   printBiedgedGraph(threadSet_debug, "before_bad_pinch.gvz");
		   exit(1);
		   }
		   */
		applyPinch(threadSet_debug, pinch);
#endif
	}
	stPinchIterator_destruct(pinchIterator);
	return threadSet;
}

stPinchSegment *getSegment(stPinchBlock *block, int64_t threadName) {
	stPinchSegment *segment;
	stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
	while((segment = stPinchBlockIt_getNext(&blockIt))) {
		if (stPinchThread_getName(stPinchSegment_getThread(segment)) == threadName) break;
	}
	return segment;
}

stList *traversePath(stPinchThreadSet *threadSet, stList *path, stHash *sequences) {
	stList *pathSeqList = stList_construct();

	int64_t pathLength = 0;
	for (int64_t i = 0; i < stList_length(path); i++) {
		if (i > 0) {
			//Fill in the adjacency connecting the previous block
			//to this one
			PONode *prevNode = stList_get(path, i - 1);
			stPinchEnd *prevEnd = prevNode->data;
			stPinchEnd *end1 = stPinchEnd_construct(stPinchEnd_getBlock(prevEnd), !stPinchEnd_getOrientation(prevEnd));
			PONode *node = stList_get(path, i);
			stPinchEnd *end2 = node->data;
			assert(end1);
			assert(end2);

			stPinchBlock *block1 = stPinchEnd_getBlock(end1);
			stPinchBlock *block2 = stPinchEnd_getBlock(end2);
			stSortedSet *threads1 = 
				getThreads(stPinchBlock_getFirst(block1));

			stSortedSet *threads2 = 
				getThreads(stPinchBlock_getFirst(block2));

			stSortedSet *sharedThreadsSet = stSortedSet_getIntersection(threads1, threads2);
			stList *sharedThreads = stSortedSet_getList(sharedThreadsSet);

			stList *candidateAdjacencySeqs = stList_construct();
			for (int j = 0; j < stList_length(sharedThreads); j++) {
				int64_t threadName = (int64_t) stList_get(sharedThreads, j);

				stPinchSegment *seg1 = getSegment(block1, threadName);
				stPinchSegment *seg2 = getSegment(block2, threadName);

				int64_t adjacencyLength = stPinchSegment_getStart(seg2) - stPinchSegment_getStart(seg1) - stPinchSegment_getLength(seg1);
				assert(adjacencyLength >= 0);

				char *adj = malloc(sizeof(char)* (adjacencyLength + 1));
				char *seq = stHash_search(sequences, (void*)threadName);
				strncpy(seq, adj, adjacencyLength);
				adj[adjacencyLength] = '\0';
				stList_append(candidateAdjacencySeqs, adj);

			}

			//pick the shortest one
			int64_t minLength = INT_MAX;
			char *adjacencySequence = NULL;
			for (int j = 0; j < stList_length(candidateAdjacencySeqs); j++) {
				char *adj = stList_get(candidateAdjacencySeqs, j);
				if (strlen(adj) < minLength) {
					minLength = strlen(adj);
					free(adjacencySequence);
					adjacencySequence = adj;
				}
				else {
					free(adj);
				}
			}
			assert(adjacencySequence);

			pathLength += strlen(adjacencySequence);
			stList_append(pathSeqList, adjacencySequence);

		}

		PONode *node = stList_get(path, i);
		stPinchEnd *end = node->data;
		assert(end);
		stPinchBlock *block = stPinchEnd_getBlock(end);
		stPinchSegment *segment = stPinchBlock_getFirst(block);
		int64_t threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
		char *blockSequence = malloc(sizeof(char) * (stPinchSegment_getLength(segment) + 1));
		strncpy(stHash_search(sequences, (void*)threadName), blockSequence, stPinchSegment_getLength(segment));
		blockSequence[stPinchSegment_getLength(segment) + 1] = '\0';

		pathLength += strlen(blockSequence);
		stList_append(pathSeqList, blockSequence);
	}

	/*
	char *pathSequence = malloc(sizeof(char)*(pathLength + 1));
	char *pos = pathSequence;
	for (int i = 0; i < stList_length(pathSeqList); i++) {
		char *pathSubsequence = stList_get(pathSeqList, i);
		strncpy(pos, pathSubsequence, strlen(pathSubsequence));
		pos += strlen(pathSubsequence);
		free(pathSubsequence);
	}
	pathSequence[pathLength] = '\0';
	assert(pathLength == strlen(pathSequence));
	*/
	return pathSeqList;
}

//basic dp method for getting the consensus sequence
//that maximizes the total weight of blocks visited
stList *heaviestPath(stList *poGraph) {
	int N = stList_length(poGraph);
	int64_t *paths = (int64_t*) calloc(N, sizeof(int64_t));
	int64_t *score = (int64_t*) calloc(N, sizeof(int64_t));

	int64_t bestScore = 0;
	int64_t bestEndpoint = -1;

	for (int64_t i = 0; i < N; i++) {
		PONode *node = stList_get(poGraph, i);

		//fprintf(stderr, "Node ID: %ld, incoming nodes: %ld node weight: %ld\n", node->nodeID, node->nIncomingNodes, node->weight);

		int64_t bestLeft = -1;
		int64_t maxWeight = 0;

		for (int64_t k = 0; k < node->nIncomingNodes; k++) {
			PONode *incomingNode = stList_get(poGraph, node->incomingNodes[k]);
			if (incomingNode->weight > maxWeight) {
				bestLeft = incomingNode->nodeID;
				maxWeight = incomingNode->weight;
			}
		}

		paths[i] = bestLeft;
		score[i] = maxWeight + score[bestLeft];

		if (score[i] > bestScore) {
			bestEndpoint = i;
			bestScore = score[i];
		}
		fprintf(stderr, "Score: %ld\n", score[i]);
	}

	//traceback
	fprintf(stderr, "Tracing back path\n");
	stList *path = stList_construct();
	for (int64_t i = bestEndpoint; i < N; i = paths[i]) {
		stList_append(path, stList_get(poGraph, i));
	}

	return path;
}
