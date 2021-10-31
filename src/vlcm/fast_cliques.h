#ifndef VLCM_FAST_CLIQUES_H
#define VLCM_FAST_CLIQUES_H

#include <vector>
#include "utils.h"
#include "graph_utils.h"

std::vector<Clique> getMaximumCliques(const DistanceGraph& graph, int tsLen, int windowMin, int windowMax, bool output);

#endif //VLCM_FAST_CLIQUES_H
