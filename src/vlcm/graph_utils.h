#ifndef VLCM_GRAPH_UTILS_H
#define VLCM_GRAPH_UTILS_H

#include <vector>
#include <boost/graph/compressed_sparse_row_graph.hpp>


typedef std::pair<int, int> Edge;
typedef std::pair<int, int> WindowRange;
typedef std::vector<WindowRange> WindowRanges;
typedef boost::compressed_sparse_row_graph<boost::directedS, boost::no_property, WindowRanges> DistanceGraph;
typedef std::vector<std::vector<int>> AdjacencyList;


inline void addEdge(AdjacencyList &adjacencyList, int source, int target) {
    adjacencyList[source].push_back(target);
    adjacencyList[target].push_back(source);
}

inline bool edgeExists(const WindowRanges &ranges, int window) {
    return std::any_of(ranges.begin(), ranges.end(), [&](const WindowRange &range) {
        return range.first <= window && window <= range.second;
    });
}

inline bool edgeExists(const WindowRanges &ranges, int windowMin, int windowMax) {
    return std::any_of(ranges.begin(), ranges.end(), [&](const WindowRange &range) {
        return !((range.first > windowMin && range.first > windowMax) ||
                (range.second < windowMin && range.second < windowMax));
    });
}

AdjacencyList getGraph(const DistanceGraph &graph, int window);

template<typename Container>
AdjacencyList getSubgraph(const DistanceGraph &graph, int window, const Container &nodes) {
    AdjacencyList adjacencyList(boost::num_vertices(graph));
    for (int node: nodes) {
        BGL_FORALL_OUTEDGES(boost::vertex(node, graph), e, graph, DistanceGraph) {
                if (edgeExists(graph[e], window)) {
                    int target = boost::target(e, graph);
                    if (std::find(nodes.begin(), nodes.end(), target) != nodes.end()) { // if target in nodes
                        addEdge(adjacencyList, node, target);
                    }
                }
            }
    }
    return adjacencyList;
}

AdjacencyList getUnionGraph(const DistanceGraph &graph, int windowMin, int windowMax);

AdjacencyList getFullUnionGraph(const DistanceGraph &graph);

AdjacencyList reduceGraph(const AdjacencyList &adjacencyList, int minDegree, int exclusion);

std::set<int> expandSubset(const std::vector<int> &set, int neg, int pos, int nodes);

#endif //VLCM_GRAPH_UTILS_H
