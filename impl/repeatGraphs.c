#include <stdio.h>
#include "bioioC.h"
#include "sonLibTypes.h"
#include "pairwiseAlignment.h"

#include "stPinchGraphs.h"

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

	stPinchBlock *block1 = stPinchSegment_getBlock(seg1);
	stPinchBlock *block2 = stPinchSegment_getBlock(seg2);
	if (block1 && stPinchBlock_getDegree(block1) > MAX_BLOCK_DEGREE) return true;
	if (block2 && stPinchBlock_getDegree(block2) > MAX_BLOCK_DEGREE) return true;

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

void printBiedgedGraph(stPinchThreadSet *graph, char *gvizFilename) {
	stList *ordering = getBlockOrdering(graph);
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

stPinchSegment *getSegment(stPinchBlock *block, int64_t threadName) {
	stPinchSegment *segment;
	stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
	while((segment = stPinchBlockIt_getNext(&blockIt))) {
		if (stPinchThread_getName(stPinchSegment_getThread(segment)) == threadName) break;
	}
	return segment;
}

char *getConsensusSequence(stList *path, stHash *sequences) {
	int64_t consensusSeqLength = 0;
	char *consensusSeq = malloc(sizeof(char));
	int64_t pos = 0;

	for (int64_t i = 0; i < stList_length(path); i++) {
		if (i > 0) {
			//Fill in the adjacency connecting the previous block
			//to this one
			stPinchEnd *oppositeEnd1 = stList_get(path, i - 1);
			stPinchBlock *block1 = stPinchEnd_getBlock(oppositeEnd1);
			stPinchEnd *end1 = stPinchEnd_construct(block1, !stPinchEnd_getOrientation(oppositeEnd1));
			stPinchEnd *end2 = stList_get(path, i);

			stSortedSet *adjacencies = getConnectingThreads(end1, end2);
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
				consensusSeqLength += adjLength;
				consensusSeq = realloc(consensusSeq, (consensusSeqLength + 1)*sizeof(char));
				char *sequence = stHash_search(sequences, (void*) threadName);
				assert(sequence);
				assert(strlen(sequence) >= adjStart + adjLength);

				strncpy(consensusSeq + pos, sequence + adjStart, adjLength);
				pos += adjLength;
			}
			stSortedSet_destructIterator(adjIt);
			stSortedSet_destruct(adjacencies);
			stPinchEnd_destruct(end1);

		}

		stPinchEnd *end = stList_get(path, i);
		stPinchBlock *block = stPinchEnd_getBlock(end);
		assert(block);
		stPinchSegment *segment = stPinchBlock_getFirst(block);
		assert(segment);
		int64_t threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
		char *threadSequence = (char*) stHash_search(sequences, (void*) threadName);
		assert(threadSequence);
		int64_t blockLength = stPinchSegment_getLength(segment);
		int64_t blockStart = stPinchSegment_getStart(segment);
		consensusSeqLength += blockLength;
		consensusSeq = realloc(consensusSeq, (consensusSeqLength + 1)*sizeof(char));
		strncpy(consensusSeq + pos, threadSequence + blockStart, blockLength);
		pos += blockLength;
	}
	consensusSeq[consensusSeqLength] = '\0';
	return consensusSeq;
}

stList *getBlockOrdering3(stPinchEnd *startEnd, stSet *seen) {
	stList *ordering = stList_construct();
	stList *stack = stList_construct();
	stList_append(stack, startEnd);
	while (stList_length(stack) > 0) {
		stPinchEnd *end = stList_pop(stack);
		stPinchBlock *block = stPinchEnd_getBlock(end);

		if (stSet_search(seen, block)) {
			stPinchEnd_destruct(end);
			continue;
		}
		stSet_insert(seen, block);

		stList_append(ordering, end);

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
		stSet *connectedEnds = stPinchEnd_getConnectedPinchEnds(oppositeEnd);
		stPinchEnd *connectedEnd;
		stSetIterator *connectedEndsIt = stSet_getIterator(connectedEnds);
		while ((connectedEnd = stSet_getNext(connectedEndsIt)) != NULL) {
			stPinchBlock *connectedBlock = stPinchEnd_getBlock(connectedEnd);

			//duplicate the end because this copy will be destructed
			//along with the connected ends set
			stPinchEnd *connectedEndCopy = stPinchEnd_construct(connectedBlock, stPinchEnd_getOrientation(connectedEnd));
			stList_append(stack, connectedEndCopy);
		}
		stPinchEnd_destruct(oppositeEnd);
		stSet_destructIterator(connectedEndsIt);
		stSet_destruct(connectedEnds);

	}
	stList_destruct(stack);
	return ordering;
}

stList *getBlockOrdering2(stPinchThreadSet *graph, stSet *seen) {
	//start DFS from the highest-weight block
	stPinchBlock *startBlock = NULL;
	int64_t highestWeight = 0;
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	stPinchBlock *block;
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
		if (stSet_search(seen, block)) continue;
		int64_t blockWeight = stPinchBlock_getDegree(block) * stPinchBlock_getLength(block);
		if (blockWeight > highestWeight) {
			highestWeight = blockWeight;
			startBlock = block;
		}
	}
	if (startBlock == NULL) return NULL;
	
	stPinchEnd *rightStart = stPinchEnd_construct(startBlock, _5PRIME);
	stPinchEnd *leftStart = stPinchEnd_construct(startBlock, _3PRIME);

	stList *orderingToRight = getBlockOrdering3(rightStart, seen);
	stSet_remove(seen, startBlock);

	stList *orderingToLeft = getBlockOrdering3(leftStart, seen);

	stList *ordering = stList_construct();
	//Reverse the direction and trim off first block to avoid double counting
	for (int64_t i = stList_length(orderingToLeft) - 1; i >= 1; i--) {
		stPinchEnd *end = stList_get(orderingToLeft, i);
		stPinchBlock *block = stPinchEnd_getBlock(end);
		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
		stList_append(ordering, oppositeEnd);
		stPinchEnd_destruct(end);
	}
	stList_appendAll(ordering, orderingToRight);
	stList_destruct(orderingToLeft);
	stList_destruct(orderingToRight);
	return ordering;
}

//get a partial ordering of the blocks implied by
//the order they are seen in a DFS from the highest-weight
//block
stList *getBlockOrdering(stPinchThreadSet *graph) {
	stSet *seen = stSet_construct();
	stList *ordering = stList_construct3(0, (void (*)(void *))stPinchEnd_destruct);

	stList *componentOrdering = NULL;
	while((componentOrdering = getBlockOrdering2(graph, seen)) != NULL) {
		stList_appendAll(ordering, componentOrdering);
	}
	return ordering;
}

int64_t getMinAdjacencyLength(stPinchEnd *end1, stPinchEnd *end2) {
	stSortedSet *connectingThreads = getConnectingThreads(end1, end2);
	int64_t minAdjacencyLength = INT_MAX;
	stIntTuple *connectingThreadInfo;
	stSortedSetIterator *connectingThreadsIt = stSortedSet_getIterator(connectingThreads);
	while((connectingThreadInfo = stSortedSet_getNext(connectingThreadsIt))) {
		int64_t adjLength = stIntTuple_get(connectingThreadInfo, 2);
		if (adjLength < minAdjacencyLength) {
			minAdjacencyLength = adjLength;
		}
	}
	stSortedSet_destruct(connectingThreads);
	stSortedSet_destructIterator(connectingThreadsIt);
	return minAdjacencyLength;
}

stList *getHeaviestPath(stList *blockOrdering, int64_t gapPenalty, stSortedSet *pathThreads, stSet *ignoredBlocks, int64_t *pathScore) {
	int64_t N = stList_length(blockOrdering);
	int64_t *scores = calloc(N, sizeof(int64_t));
	int64_t *directions = calloc(N, sizeof(int64_t));

	stHash *blockIndex = stHash_construct();

	int64_t bestPathStart = -1;
	int64_t bestScore = -INT_MAX;
	for (int64_t i = 0; i < stList_length(blockOrdering); i++) {
		directions[i] = -1;
		scores[i] = 0;
		stPinchEnd *end = stList_get(blockOrdering, i);
		stPinchBlock *block = stPinchEnd_getBlock(end);
		if (stSet_search(ignoredBlocks, block)) continue;

		//only true if negative pinches aren't allowed,
		//might have to be removed later
		assert(stPinchEnd_getOrientation(end) == 1);

		stHash_insert(blockIndex, block, (void*) i + 1);


		stSortedSet *threadsInBlock = getThreads(stPinchBlock_getFirst(block));
		stSortedSet *sharedThreads = stSortedSet_getIntersection(pathThreads, threadsInBlock);

		int64_t blockWeight = (2*stSortedSet_size(sharedThreads) - stSortedSet_size(pathThreads)) * stPinchBlock_getLength(block);
		stSortedSet_destruct(threadsInBlock);
		stSortedSet_destruct(sharedThreads);

		int64_t bestLeftScore = -INT_MAX;
		bool hasPredecessors = false;

		stSet *adjacentEnds = stPinchEnd_getConnectedPinchEnds(end);
		stSetIterator *adjacentEndsIt = stSet_getIterator(adjacentEnds);
		stPinchEnd *adjEnd;
		while((adjEnd = stSet_getNext(adjacentEndsIt)) != NULL) {
			stPinchBlock *leftBlock = stPinchEnd_getBlock(adjEnd);

			if (!stHash_search(blockIndex, leftBlock)) continue;
			int64_t leftBlockPosition = (int64_t) stHash_search(blockIndex, leftBlock) - 1;
			//this edge violates the ordering so it is ignored
			if (leftBlockPosition >= i) continue;
			hasPredecessors = true;

			assert(stPinchEnd_getOrientation(adjEnd) == 0);
			stPinchEnd *endInOrdering = stList_get(blockOrdering, leftBlockPosition);

			//Make sure this block is being traversed in the direction consistent with the ordering
			if (stPinchEnd_getOrientation(endInOrdering) == stPinchEnd_getOrientation(adjEnd)) continue;
			
			int64_t bestAdjLength = getMinAdjacencyLength(adjEnd, end);
			int64_t scoreFromLeftBlock = scores[leftBlockPosition] + blockWeight - gapPenalty * bestAdjLength;

			if (scoreFromLeftBlock > bestLeftScore) {
				bestLeftScore = scoreFromLeftBlock;
				directions[i] = leftBlockPosition;
			}
		}
		stSet_destructIterator(adjacentEndsIt);
		stSet_destruct(adjacentEnds);

		if (hasPredecessors) {
			scores[i] = bestLeftScore;
		}
		else {
			//base case for blocks with no predecessors
			scores[i] = blockWeight;
		}

		if (scores[i] > bestScore) {
			bestScore = scores[i];
			bestPathStart = i;
		}
	}

	assert(directions[0] == -1);
	stList *path = stList_construct();

	int64_t pos = bestPathStart;
	while ((pos >= 0) && (scores[pos] > 0)) {
		stPinchEnd *end = stList_get(blockOrdering, pos);
		stList_append(path, end);
		pos = directions[pos];
	}
	stList_reverse(path);
	*pathScore = bestScore;

	free(scores);
	free(directions);
	stHash_destruct(blockIndex);
	return path;
}

stPinchBlock *getHighestWeightBlock(stPinchThreadSet *graph) {
	stPinchBlock *bestBlock = NULL;
	int64_t highestWeight = 0;
	stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(graph);
	stPinchBlock *block;
	while((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
		int64_t blockWeight = stPinchBlock_getDegree(block) * stPinchBlock_getLength(block);
		if (blockWeight > highestWeight) {
			highestWeight = blockWeight;
			bestBlock = block;
		}
	}
	return bestBlock;
}

/*
stList *extendDensePath(stPinchThreadSet *graph) {
	stPinchBlock *startBlock = getHighestWeightBlock(graph);
	if (!startBlock) return NULL;

	stSet *seenBlocks = stSet_construct();
	stSet *seenThreads = stSet_construct();
	stHash *directions = stHash_construct();

	int64_t bestScore = 0;
	int64_t curPathScore = stPinchBlock_getDegree(startBlock) * stPinchBlock_getLength(startBlock);

	int64_t maxBadIterations = 5;
	int64_t badIterations = 0;

	stPinchBlock *lastBlock = startBlock;

	stSortedSet *threads = getThreads(startBlock);
	stSortedSet *pathWeight = stSortedSet_size(threads);

	//extend to the right
	stPinchEnd *end = stPinchEnd_construct(startBlock, _5PRIME);
	while (end && (badIterations < maxBadIterations)) {
		stPinchBlock *block = stPinchEnd_getBlock(end);
		stSet_insert(seenBlocks, block);

		int64_t bestNextBlockScore = 0;
		stPinchEnd *bestNextEnd = NULL;

		stPinchEnd *oppositeEnd = stPinchEnd_construct(block, !stPinchEnd_getOrientation(end));
		stSet *connectedEnds = stPinchEnd_getConnectedEnds(oppositeEnd);
		stSetIterator *connectedEndsIt = stSet_getIterator(connectedEnds);
		stPinchEnd *adjEnd = NULL;
		while((adjEnd = stSet_getNext(connectedEndsIt)) != NULL) {
			stPinchBlock *nextBlock = stPinchEnd_getBlock(adjEnd);
			if (stSet_search(seenBlocks, nextBlock)) continue;

			int64_t adjacencyLength = getMinAdjacencyLength(oppositeEnd, adjEnd);
			
			stSortedSet *nextBlockThreads = getThreads(nextBlock);
			stSortedSet *sharedThreads = stSortedSet_getIntersection(threads, nextBlockThreads);

			int64_t blockWeight = stSortedSet_size(sharedThreads) * stPinchBlock_getLength(nextBlock);

			int64_t nextBlockScore = blockWeight - adjacencyLength * pathWeight;
			if (nextBlockScore > bestNextBlockScore) {
				bestNextBlockScore = nextBlockScore;
				bestNextEnd = adjEnd;
			}
			
		}
		stSet_destructIterator(connectedEndsIt);
		stSet_destruct(connectedEnds);
		stPinchEnd_destruct(oppositeEnd);

		stPinchBlock *nextBlock = stPinchEnd_getBlock(nextEnd);

		curPathScore += bestNextBlockScore;
		stHash_insert(directions, block,)
		if (curPathScore > bestScore) {
			bestScore = curPathScore;
			lastBlock = nextBlock;
			badIterations = 0;
		}
		else {
			badIterations++;
		}

		end = bestNextEnd;
		
	}
}
*/