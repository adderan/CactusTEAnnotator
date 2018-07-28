import sys
import networkx
import re
import argparse
import multiset

def getPartitions(graphFile):
    partitions = multiset.Multiset()
    with open(graphFile, 'r') as graphRead:
        for line in graphRead:
            partition = set()
            threadSets = line.split(",")
            for threadSet in threadSets:
                threadSet = threadSet.split()
                threadSet = [int(thread) for thread in threadSet]
                partition.insert(frozenset(threadSet))
            partition = frozenset(partition)
            partitions.insert(partition)

    partitions = partitions.items()
    partitions.sort(key = lambda x: x[1], reverse=True)
    return partitions

def buildTree(graphFile):
    partitions = getPartitions(graphFile)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graphFile")
    #parser.add_argument("--output_categories", action="store_", default=True)

    args = parser.parse_args()
    graph = POGraph(args.graphFile)
    partitions = graph.getPartitions()
    tree = buildTree(graph.threads, partitions)
    assert networkx.is_tree(tree)
    leafPartitioning = getLeafPartitioning(tree)
    for partition in leafPartitioning:
        sys.stdout.write(" ".join(partition) + "\n")

if __name__ == "__main__":
    main()
