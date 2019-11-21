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

stSortedSet *getConnectingThreads(stPinchEnd *end1, stPinchEnd *end2) {
	stSortedSet *connectingThreads = stSortedSet_construct2((void*)stIntTuple_destruct);
	stPinchBlock *block1 = stPinchEnd_getBlock(end1);
	stPinchBlock *block2 = stPinchEnd_getBlock(end2);

	stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block1);
	stPinchSegment *seg_block1;
	
	while((seg_block1 = stPinchBlockIt_getNext(&blockIt))) {
		bool traverse5Prime = stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(end1), seg_block1);
	
		//because no negative pinches
		assert(traverse5Prime == false);
		stPinchSegment *seg_block2 = seg_block1;
		do {
			seg_block2 = (traverse5Prime) ? 
				stPinchSegment_get5Prime(seg_block2) : stPinchSegment_get3Prime(seg_block2);
		}
		while(seg_block2 && (stPinchSegment_getBlock(seg_block2) == NULL));
		if (!seg_block2) continue;
		if (stPinchSegment_getBlock(seg_block2) == block2) {
			int64_t threadName = stPinchThread_getName(stPinchSegment_getThread(seg_block2));
			int64_t adjStart = stPinchSegment_getStart(seg_block1) 
				+ stPinchSegment_getLength(seg_block1);
			int64_t adjLength = stPinchSegment_getStart(seg_block2) - adjStart;
			stSortedSet_insert(connectingThreads, stIntTuple_construct3(threadName, adjStart, adjLength));
		}
	}

	//Check
	//stList *subsequences = stPinchEnd_getSubSequenceLengthsConnectingEnds(end1, end2);
	//assert(stList_length(subsequences) == stSortedSet_size(connectingThreads));
	return connectingThreads;
}


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
				PartialOrderNode *node = malloc(sizeof(PartialOrderNode));
				node->nodeID = nodeID;
				node->degree = stPinchBlock_getDegree(block);
				node->length =  stPinchBlock_getLength(block);
				stList *incomingEndList = stSet_getList(incomingAdjacencies);

				node->nIncomingNodes = stList_length(incomingEndList);
				node->incomingNodes = calloc(node->nIncomingNodes, sizeof(int64_t));
				node->incomingEdgeWeights = calloc(node->nIncomingNodes, sizeof(int64_t));
				for (int k = 0; k < node->nIncomingNodes; k++) {
					stPinchEnd *incomingEnd = stList_get(incomingEndList, k);
					stPinchBlock *incomingBlock = stPinchEnd_getBlock(incomingEnd);
					//We should have already encountered this block
					node->incomingNodes[k] = (int64_t) stHash_search(blockIndex, incomingBlock);
					stSortedSet *connectingThreads = getConnectingThreads(incomingEnd, end);
					node->incomingEdgeWeights[k] = stSortedSet_size(connectingThreads);
					stSortedSet_destruct(connectingThreads);
				}
				node->data = block;
				node->orientation = stPinchEnd_getOrientation(end);
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

stList *traversePath(stPinchThreadSet *graph, stList *path, stHash *sequences) {
	stList *pathSeqList = stList_construct();

	for (int64_t i = 0; i < stList_length(path); i++) {
		if (i > 0) {
			//Fill in the adjacency connecting the previous block
			//to this one
			PartialOrderNode *prevNode = stList_get(path, i - 1);
			stPinchBlock *block1 = prevNode->data;

			PartialOrderNode *node = stList_get(path, i);
			stPinchBlock *block2 = node->data;

			//can only assume the orientations because there 
			//are no negative pinches
			stPinchEnd *end1 = stPinchEnd_construct(block1, _3PRIME);
			stPinchEnd *end2 = stPinchEnd_construct(block2, _5PRIME);

			stSortedSet *adjacencies = getConnectingThreads(end1, end2);
			stSet *connectedEnds = stPinchEnd_getConnectedPinchEnds(end1);
			assert(stSet_size(connectedEnds) > 0);
			assert(stSortedSet_size(adjacencies) > 0);

			stSortedSetIterator *adjIt = stSortedSet_getIterator(adjacencies);
			stIntTuple *adjTuple;
			stIntTuple *bestAdj = NULL;
			int64_t bestAdjLength = INT_MAX;
			while ((adjTuple = stSortedSet_getNext(adjIt)) != NULL) {
				int64_t adjLength = stIntTuple_get(adjTuple, 2);

				if (adjLength < bestAdjLength) {
					bestAdj = adjTuple;
					bestAdjLength = adjLength;
				}
			}

			int64_t threadName = stIntTuple_get(bestAdj, 0);
			int64_t adjStart = stIntTuple_get(bestAdj, 1);
			int64_t adjLength = stIntTuple_get(bestAdj, 2);

			if (adjLength > 0) {
				char *adjacencySequence = calloc(adjLength + 1, sizeof(char));
				char *sequence = stHash_search(sequences, (void*) threadName);
				assert(sequence);
				assert(strlen(sequence) >= adjStart + adjLength);

				strncpy(adjacencySequence, sequence + adjStart, adjLength);
				assert(adjacencySequence);

				stList_append(pathSeqList, adjacencySequence);
			}
			stSortedSet_destruct(adjacencies);

		}

		PartialOrderNode *node = stList_get(path, i);
		stPinchBlock *block = (stPinchBlock*) node->data;
		assert(block);
		stPinchSegment *segment = stPinchBlock_getFirst(block);
		int64_t threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
		char *threadSequence = stHash_search(sequences, (void*) threadName);
		int64_t blockLength = stPinchSegment_getLength(segment);
		int64_t blockStart = stPinchSegment_getStart(segment);
		char *blockSequence = calloc(blockLength + 1, sizeof(char));
		strncpy(blockSequence, threadSequence + blockStart, blockLength);
		stList_append(pathSeqList, blockSequence);
	}
	return pathSeqList;
}

stList *getHeaviestPath(stList *graph) {
	int64_t N = stList_length(graph);
	int64_t *score = calloc(N, sizeof(int64_t));
	int64_t *paths = calloc(N, sizeof(int64_t));

	int64_t bestScore = 0;
	int64_t endpoint = -1;
	for (int64_t i = 0; i < N; i++) {
		PartialOrderNode *node = stList_get(graph, i);
		int64_t blockWeight = node->degree * node->length;

		int64_t maxEdgeWeight = 0;
		paths[i] = -1;
		for (int64_t k = 0; k < node->nIncomingNodes; k++) {
			PartialOrderNode *incomingNode = 
				stList_get(graph, node->incomingNodes[k]);
			
			int64_t edgeWeight = node->incomingEdgeWeights[k];

			if (edgeWeight > maxEdgeWeight) {
				paths[i] = k;
				assert(paths[i] < N);
				maxEdgeWeight = edgeWeight;
				score[i] = score[incomingNode->nodeID] + blockWeight;
			}
		}

		if (score[i] > bestScore) {
			bestScore = score[i];
			endpoint = i;
		}
	}

	if (endpoint == -1) return NULL;

	//traceback and adjust weights
	stList *heaviestPath = stList_construct();
	int64_t i = endpoint;
	while(true) {
		PartialOrderNode *node = stList_get(graph, i);
		assert(node);
		stList_append(heaviestPath, node);

		int64_t k = paths[i];
		if (k == -1) break;

		//zero out the edge chosen in this traversal
		node->incomingEdgeWeights[k] = 0;

		i = node->incomingNodes[k];
	}
	free(paths);
	free(score);
	stList_reverse(heaviestPath);
	return heaviestPath;
}