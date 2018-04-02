import sys
import networkx
import re
import argparse
import multiset

def buildTree(threads, partitions):
    #initialize star tree
    tree = networkx.DiGraph()
    tree.add_node("root")
    for seqName in threads:
        tree.add_node(seqName)
        tree.add_edge("root", seqName)
    assert networkx.is_tree(tree)

    for partition in partitions:
        applyPartition(tree, partition[0])
    return tree

def getParent(tree, node):
    assert len(tree.in_edges(node)) == 1
    parent = list(tree.in_edges(node))[0][0]
    return parent

def applyPartition(tree, partition):
    #Make sure all threads in thread set have common parent
    for threadSet in partition:
        parents = set()
        for thread in threadSet:
            parents.add(getParent(tree, thread))
        if len(parents) > 1:
            #Partition isn't compatible with the tree
            #print "Skipping partition %s because nodes have different parents" % partition
            return

    #Build dictionary of parent nodes to the thread sets that
    #descend from them, so we can split thread sets that share a parent
    parentToThreadSets = {}
    for threadSet in partition:
        parent = getParent(tree, threadSet[0])
        if parent not in parentToThreadSets:
            parentToThreadSets[parent] = []
        parentToThreadSets[parent].append(threadSet)

    #Create parent nodes for each thread set
    for parent in parentToThreadSets: 
        if len(parentToThreadSets[parent]) < 2:
            #nothing to do
            continue
        #print parentToThreadSets[parent]
        for threadSet in parentToThreadSets[parent]:
            if len(threadSet) < 2:
                continue
            print("Joining %s" % "_".join(threadSet))
            new_parent = "_".join(threadSet) + "_anc"
            tree.add_node(new_parent)
            tree.add_edge(parent, new_parent)
            for thread in threadSet:
                assert tree.has_edge(parent, thread)
                tree.remove_edge(parent, thread)
                tree.add_edge(new_parent, thread)
                if not networkx.is_tree(tree):
                    import pdb; pdb.set_trace()
        assert networkx.is_tree(tree)

class Node:
    def __init__(self, base, incidentNodes, incidentThreads, ringNodes):
        self.base = base
        self.incidentNodes = incidentNodes
        self.incidentThreads = incidentThreads
        self.ringNodes = ringNodes

class POGraph:
    def __init__(self, graphFile):
        self.threads = []
        self.nodes = []

        graphRead = open(graphFile, 'r')
        for line in graphRead:
            if line.startswith("SOURCENAME"):
                lineinfo = line.split("=")
                self.threads.append(lineinfo[1].rstrip())

            elif line.startswith("A:") or line.startswith("C:") or line.startswith("G:") or line.startswith("T:"):
                base, nodeinfo = line.split(":")
                labels = re.findall("[ASL]\d+", nodeinfo)

                incidentNodes = []
                incidentThreads = []
                ringNodes = []

                for label in labels:
                    attribute = label[0]
                    value = int(label[1:])
                    if attribute == 'A':
                        ringNodes.append(value)
                    elif attribute == 'L':
                        incidentNodes.append(value)
                    elif attribute == 'S':
                        incidentThreads.append(value)
                self.nodes.append(Node(base=base, incidentNodes=incidentNodes, incidentThreads=incidentThreads, ringNodes=ringNodes))
        graphRead.close()

    def getPartitions(self):
        right_partitions = [set() for i in range(len(self.nodes))]
        left_partitions = [set() for i in range(len(self.nodes))]
        for j in range(len(self.nodes)):
            #Nodes i with edges to j
            for i in self.nodes[j].incidentNodes:
                threads_i = set([self.threads[thread_id] for thread_id in self.nodes[i].incidentThreads])
                threads_j = set([self.threads[thread_id] for thread_id in self.nodes[j].incidentThreads])
                #threads connecting i and j
                threads_ij = threads_i.intersection(threads_j)
                left_partitions[j].add(tuple(threads_ij))
                right_partitions[i].add(tuple(threads_ij))

        partitions = left_partitions + right_partitions
        partition_set = multiset.Multiset()
        for partition in partitions:
            partition_set.add(tuple(partition))
        partitions = partition_set.items()
        #Sort by multiplicity of the partition
        partitions.sort(key=lambda x: x[1], reverse=True)
        return partitions

def getCategories(tree):
    categories = {}
    for node in tree.nodes():
        if tree.out_degree(node) > 0:
            #not a leaf
            continue
        parent = getParent(tree, node)
        if not parent in categories:
            categories[parent] = []
        categories[parent].append(node)
    return categories.values()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graphFile")
    #parser.add_argument("--output_categories", action="store_", default=True)

    args = parser.parse_args()
    graph = POGraph(args.graphFile)
    partitions = graph.getPartitions()
    tree = buildTree(graph.threads, partitions)
    categories = getCategories(tree)
    print categories

if __name__ == "__main__":
    main()
