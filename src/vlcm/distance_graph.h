#ifndef VLCM_DISTANCE_GRAPH_H
#define VLCM_DISTANCE_GRAPH_H

#include <vector>
#include "graph_utils.h"
#include "precompute.h"

typedef std::vector<Edge> EdgeList;
typedef std::vector<WindowRanges> WindowRangesList;

std::pair<EdgeList, WindowRangesList> computeDistanceGraph(const TimeSeriesInfo &ts, const PrecomputedInfo &info);

std::pair<EdgeList, WindowRangesList> computeDistanceGraphNaive(const TimeSeriesInfo &ts, const PrecomputedInfo &info);

std::pair<EdgeList, WindowRangesList>
mergeEdgeLists(const std::vector<EdgeList> &edgeLists, const std::vector<WindowRangesList> &windowRangesLists);

double getBestMultiplier(double correlation);

#endif //VLCM_DISTANCE_GRAPH_H
