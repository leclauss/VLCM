#ifndef VLCM_FAST_CLIQUES_H
#define VLCM_FAST_CLIQUES_H

#include <vector>
#include "utils.h"
#include "graph_utils.h"
#include "lmc.h"

inline void printClique(int window, const Clique &clique) {
    std::cout << window << ";";
    for (size_t i = 0; i < clique.size(); i++) {
        std::cout << clique[i];
        if (i < clique.size() - 1) {
            std::cout << ",";
        }
    }
    std::cout << std::endl;
}

inline Clique fastMaxClique(const AdjacencyList &adjacencyList, const Clique &initClique, int exclusion) {
    auto smallAdjacencyList = reduceGraph(adjacencyList, initClique.size() + 1, exclusion);
    auto maxClique = getMaximumClique(smallAdjacencyList);
    if (maxClique.size() > initClique.size()) {
        return maxClique;
    } else {
        return initClique;
    }
}

inline Clique fastMaxClique(const AdjacencyList &adjacencyList, int exclusion) {
    auto initClique = getInitClique(adjacencyList);
    return fastMaxClique(adjacencyList, initClique, exclusion);
}
std::vector<Clique> getMaximumCliques(const DistanceGraph& graph, int tsLen, int windowMin, int windowMax, int k, bool output);

#endif //VLCM_FAST_CLIQUES_H
