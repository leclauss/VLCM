#include "graph_utils.h"

AdjacencyList getGraph(const DistanceGraph &graph, int window) {
    AdjacencyList adjacencyList(boost::num_vertices(graph));
    BGL_FORALL_EDGES(e, graph, DistanceGraph) {
            if (edgeExists(graph[e], window)) {
                addEdge(adjacencyList, boost::source(e, graph), boost::target(e, graph));
            }
        }
    return adjacencyList;
}

AdjacencyList getUnionGraph(const DistanceGraph &graph, int windowMin, int windowMax) {
    AdjacencyList adjacencyList(boost::num_vertices(graph));
    BGL_FORALL_EDGES(e, graph, DistanceGraph) {
            if (edgeExists(graph[e], windowMin, windowMax)) {
                addEdge(adjacencyList, boost::source(e, graph), boost::target(e, graph));
            }
        }
    return adjacencyList;
}

AdjacencyList getFullUnionGraph(const DistanceGraph &graph) {
    AdjacencyList adjacencyList(boost::num_vertices(graph));
    BGL_FORALL_EDGES(e, graph, DistanceGraph) {
            addEdge(adjacencyList, boost::source(e, graph), boost::target(e, graph));
        }
    return adjacencyList;
}

AdjacencyList reduceGraph(const AdjacencyList &adjacencyList, int minDegree, int exclusion) {
    int minEdges = minDegree - 1;

    std::vector<bool> nodeActive(adjacencyList.size(), true);
    bool reduced;
    do {
        reduced = false;
        for (int i = 0; i < adjacencyList.size(); i++) {
            if (nodeActive[i]) {
                int lastNode = -exclusion;
                int trueEdges = 0;
                for (int node: adjacencyList[i]) {
                    if (nodeActive[node]) {
                        if (node >= lastNode + exclusion) {
                            trueEdges++;
                            lastNode = node;
                        }
                    }
                }
                if (trueEdges < minEdges) {
                    nodeActive[i] = false;
                    reduced = true;
                }
            }
        }
    } while (reduced);

    AdjacencyList smallAdjacencyList(adjacencyList.size());
    for (int i = 0; i < adjacencyList.size(); i++) {
        if (nodeActive[i]) {
            for (int edge: adjacencyList[i]) {
                if (nodeActive[edge]) {
                    smallAdjacencyList[i].push_back(edge);
                }
            }
        }
    }

    return smallAdjacencyList;
}

std::set<int> expandSubset(const std::vector<int> &set, int neg, int pos, int nodes) {
    std::set<int> result;
    for (int i: set) {
        for (int offset = neg; offset <= pos; offset++) {
            int newNode = i + offset;
            if (newNode >= 0 && newNode < nodes) {
                result.insert(i + offset);
            }
        }
    }
    return result;
}