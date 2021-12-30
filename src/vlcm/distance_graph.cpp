#include "distance_graph.h"

#include <cmath>
#include <algorithm>
#include <queue>

#include "precompute.h"

#define EXCLUSION_DIVISOR 1 // allowed overlap: 1=no overlapping; 2=subsequences may overlap in half; ...

struct DotProduct {
    double value;
    int window;
    bool isDP;
};


std::pair<EdgeList, WindowRangesList> computeDistanceGraph(const TimeSeriesInfo &ts, const PrecomputedInfo &info) {
    int exclusion = info.w / EXCLUSION_DIVISOR;

    EdgeList edgeList;
    WindowRangesList windowRanges;

    size_t rows = info.n - exclusion;
    auto correlations = new double[rows];
    auto correlationsIndex = new int[rows];
    auto dotProducts = new DotProduct[rows];
    auto ranges = new WindowRanges[rows];

    for (int row = 0; row < rows; row++) {
        size_t maxDiags = info.n - row;
        size_t calculatedDiags = maxDiags - exclusion;


        for (int diag = exclusion; diag < maxDiags; diag++) {
            int relativeDiag = diag - exclusion;
            int col = diag + row;

            // exact distance calculation (SCAMP)
            double corr = info.cov[diag] * info.norms[col] * info.norms[row]; // TODO cov[<diag] is never used

            // save correlation for pruning
            double corrPositive = std::max(corr, 0.0);
            correlations[relativeDiag] = corrPositive * corrPositive;
            correlationsIndex[relativeDiag] = relativeDiag;
            dotProducts[relativeDiag] = {info.cov[diag], info.w, false};

            // update for next row (SCAMP)
            info.cov[diag] += info.df[col] * info.dg[row] + info.dg[col] * info.df[row];

            if (corr >= ts.threshold) {
                // found edge (row, col) in graph info.w with correlation corr
                ranges[relativeDiag].emplace_back(info.w, info.w);
            }
        }

        // sort correlation index array (descending by correlation)
        std::sort(correlationsIndex, correlationsIndex + calculatedDiags, [&](const int &a, const int &b) {
            return (correlations[a] > correlations[b]);
        });

        // calculate distances for longer windows
        double sigBaseLenSquaredInverse = info.sigSquaredInverse[row];
        size_t maxW = std::min(static_cast<size_t>(info.maxW),
                               EXCLUSION_DIVISOR * (ts.len - row) / (EXCLUSION_DIVISOR + 1));
        double pruneThresholdBase = (1 - ts.threshold) * 2 / info.w * sigBaseLenSquaredInverse;
        for (int extendedW = info.w + 1; extendedW <= maxW; extendedW++) {
            int minCol = row + extendedW / EXCLUSION_DIVISOR;
            size_t maxCol = ts.len - extendedW;

            // calculate pruneThreshold
            double meanExtendedLen = (ts.sums[row + extendedW] - ts.sums[row]) / extendedW;
            double sigExtendedLenSquared =
                    (ts.squaredSums[row + extendedW] - ts.squaredSums[row]) / extendedW -
                    meanExtendedLen * meanExtendedLen;
            double pruneThreshold = 1 - pruneThresholdBase * extendedW * sigExtendedLenSquared;

            // find first index where correlation is lower than pruneThreshold
            auto it = std::upper_bound(correlationsIndex, correlationsIndex + calculatedDiags, pruneThreshold,
                                       [&](const double &value, const int &index) {
                                           return (value > correlations[index]);
                                       });
            size_t pruneIndex = it - correlationsIndex;

            // for each non-pruned index: update previous dot product and calculate exact distance
            for (int i = 0; i < pruneIndex; i++) {
                int index = correlationsIndex[i];
                int col = index + exclusion + row;
                if (minCol <= col && col <= maxCol) { // check overlapping
                    if (!dotProducts[index].isDP) { // first time: calculate dot product from stored cov
                        dotProducts[index].value += info.means[row] * info.means[col];
                        dotProducts[index].isDP = true;
                    }

                    // update dot product
                    for (int j = dotProducts[index].window; j < extendedW; j++) {
                        dotProducts[index].value += ts.data[row + j] * ts.data[col + j];
                    }
                    dotProducts[index].window = extendedW;

                    // calculate corr
                    double sumsCol = ts.sums[col + extendedW] - ts.sums[col];
                    double sigColSquaredTimesWSquared =
                            (ts.squaredSums[col + extendedW] - ts.squaredSums[col]) * extendedW -
                            sumsCol * sumsCol; //TODO check if near 0
                    double corr = (dotProducts[index].value - meanExtendedLen * sumsCol) /
                                  std::sqrt(sigExtendedLenSquared * sigColSquaredTimesWSquared);
                    if (corr >= ts.threshold) {
                        // found edge (row, col) in graph extendedW with correlation corr
                        if (ranges[index].empty() || ranges[index].back().second != extendedW - 1) {
                            ranges[index].emplace_back(extendedW, extendedW);
                        } else {
                            ranges[index].back().second = extendedW;
                        }
                    }
                }
            }
        }

        // fill edge list
        for (int index = 0; index < calculatedDiags; index++) {
            if (!ranges[index].empty()) {
                int col = index + exclusion + row;
                edgeList.emplace_back(row, col);
                windowRanges.push_back(ranges[index]);
                ranges[index].clear();
            }
        }

    }

    delete[] correlations;
    delete[] correlationsIndex;
    delete[] dotProducts;
    delete[] ranges;

    return std::make_pair(edgeList, windowRanges);
}


// only calculates edges for length info.w!
std::pair<EdgeList, WindowRangesList> computeDistanceGraphNaive(const TimeSeriesInfo &ts, const PrecomputedInfo &info) {
    int exclusion = info.w / EXCLUSION_DIVISOR;

    EdgeList edgeList;
    WindowRangesList windowRanges;
    for (int row = 0; row < info.n - exclusion; row++) {
        size_t maxDiags = info.n - row;

        for (int diag = exclusion; diag < maxDiags; diag++) {
            int col = diag + row;
            double corr = info.cov[diag] * info.norms[col] * info.norms[row];
            info.cov[diag] += info.df[col] * info.dg[row] + info.dg[col] * info.df[row];
            if (corr >= ts.threshold) {
                // found edge (row, col) in graph info.w with correlation corr
                edgeList.emplace_back(row, col);
                WindowRanges ranges;
                ranges.emplace_back(info.w, info.w);
                windowRanges.push_back(ranges);
            }
        }
    }

    return std::make_pair(edgeList, windowRanges);
}


std::pair<EdgeList, WindowRangesList>
mergeEdgeLists(const std::vector<EdgeList> &edgeLists, const std::vector<WindowRangesList> &windowRangesLists) {
    typedef std::tuple<int, int, int> EdgeEntry;
    EdgeList resultEdgeList;
    WindowRangesList resultWindowRanges;

    size_t numLists = edgeLists.size();
    std::vector<size_t> currentListIndex(numLists, 0);

    std::priority_queue<EdgeEntry, std::vector<EdgeEntry>, std::greater<>> minHeap;
    for (int list = 0; list < numLists; list++) {
        if (!edgeLists[list].empty()) {
            const Edge &edge = edgeLists[list][currentListIndex[list]];
            minHeap.emplace(edge.first, edge.second, list);
        }
    }

    while (!minHeap.empty()) {
        const EdgeEntry &entry = minHeap.top();
        Edge edge = std::make_pair(std::get<0>(entry), std::get<1>(entry));
        int list = std::get<2>(entry);
        minHeap.pop();
        if (!resultEdgeList.empty() && resultEdgeList.back() == edge) {
            const WindowRanges &newRanges = windowRangesLists[list][currentListIndex[list]];
            if (resultWindowRanges.back().back().second + 1 == newRanges.front().first) {
                resultWindowRanges.back().back().second = newRanges.front().second;
                for (int i = 1; i < newRanges.size(); i++) {
                    resultWindowRanges.back().push_back(newRanges[i]);
                }
            } else {
                for (const auto &newRange: newRanges) {
                    resultWindowRanges.back().push_back(newRange);
                }
            }
        } else {
            resultEdgeList.push_back(edge);
            resultWindowRanges.push_back(windowRangesLists[list][currentListIndex[list]]);
        }
        if (++currentListIndex[list] < edgeLists[list].size()) {
            const Edge &nextEdge = edgeLists[list][currentListIndex[list]];
            minHeap.emplace(nextEdge.first, nextEdge.second, list);
        }

    }

    return std::make_pair(resultEdgeList, resultWindowRanges);
}

double getBestMultiplier(double correlation) {
    if (correlation <= 0) return 1;
    if (correlation <= 0.55) return 1 + correlation / (0.55 - 0) * (1.1 - 1.0);
    if (correlation <= 0.6) return 1.1 + (correlation - 0.55) / (0.6 - 0.55) * (1.125 - 1.1);
    if (correlation <= 0.65) return 1.125 + (correlation - 0.6) / (0.65 - 0.6) * (1.15 - 1.125);
    if (correlation <= 0.7) return 1.15 + (correlation - 0.65) / (0.7 - 0.65) * (1.175 - 1.15);
    if (correlation <= 0.75) return 1.175 + (correlation - 0.7) / (0.75 - 0.7) * (1.225 - 1.175);
    if (correlation <= 0.8) return 1.225 + (correlation - 0.75) / (0.8 - 0.75) * (1.25 - 1.225);
    if (correlation <= 0.85) return 1.25 + (correlation - 0.8) / (0.85 - 0.8) * (1.3 - 1.25);
    if (correlation <= 0.9) return 1.3 + (correlation - 0.85) / (0.9 - 0.85) * (1.4 - 1.3);
    if (correlation <= 0.95) return 1.4 + (correlation - 0.9) / (0.95 - 0.9) * (1.6 - 1.4);
    if (correlation <= 0.99) return 1.6 + (correlation - 0.95) / (0.99 - 0.95) * (2.25 - 1.6);
    if (correlation >= 0.9999999) return 100000;
    return 0.0225 / (1 - correlation);

}