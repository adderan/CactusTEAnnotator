import sys
import networkx
import collections

def buildTree(threads, partitions):
    #initialize star tree
    tree = networkx.DiGraph()
    tree.add_node("root")
    for seqName in threads:
        tree.add_node(seqName)
        tree.add_edge("root", seqName)
    assert is_tree(tree)

def getParent(graph, node):
    assert len(graph.in_edges(node)) == 1
    parent = list(graph.in_edges(node))[1]
    return parent

def applyPartition(partition, tree):
    #Make sure all threads in thread set have common parent
    for threadSet in partition:
        parents = set()
        for thread in threadSet:
            parents.insert(getParent(tree, thread))
        if len(parents) > 1:
            #Partition isn't compatible with the tree
            return
    #Create parent nodes for each thread set
    for threadSet in partition:
        old_parent = getParent(threadSet[0])
        new_parent = "_".join(threadSet) + "_ancestor"
        tree.add_node(new_parent)
        for thread in threadSet:
            assert tree.has_edge(old_parent, thread)
            tree.delete_edge(old_parent, thread)
            tree.add_edge(new_parent, thread)
        tree.add_edge(old_parent, new_parent)

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
        right_partitions = [[] for i in range(len(self.nodes))]
        left_partitions = [[] for i in range(len(self.nodes))]
        for j in range(len(self.nodes)):
            #Nodes i with edges to j
            for i in nodes[j].incidentNodes:
                thread_i = set([self.threads[thread_id] for thread_id in self.nodes[i].threads])
                thread_j = set([self.threads[thread_id] for thread_id in self.nodes[j].threads])
                #threads connecting i and j
                threads_ij = list(threads_i.intersect(threads_j))
                left_partitions[j].append(threads_ij)
                right_partitions[i].append(threads_ij)

    partitions = collections.defaultdict()
    for partition in left_partitions:
        partitions[partition] += 1
    for partition in right_partitions:
        partitions[partition] += 1
    return partitions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graphFile")

    args = parser.parse_args()
    graph = POGraph(args.graphFile)
