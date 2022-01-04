#include <iostream>
#include <vector>
#include <chrono>
#include <boost/graph/compressed_sparse_row_graph.hpp>

#include "utils.h"
#include "exception.h"
#include "precompute.h"
#include "distance_graph.h"
#include "graph_utils.h"
#include "fast_cliques.h"

bool verbose;
#define PRINT(msg) if(verbose) std::cout << msg << std::endl;
#define TIME_START if(verbose) {START}
#define TIME_END(msg) if(verbose) {END(msg)}

void runSingleVLCM(const TimeSeriesInfo &ts, int windowMin, int windowMax) {
    TIME_INIT

    // compute statistics for entire problem
    TIME_START
    PrecomputedInfo info(ts, windowMin, windowMax);
    TIME_END("Precomputation done")

    // compute distance graph
    TIME_START
    auto[edgeList, edgePropertyList] = computeDistanceGraph(ts, info);
    TIME_END("Calculated distance graph")

    // create csr graph
    TIME_START
    DistanceGraph graph(boost::edges_are_sorted, edgeList.begin(), edgeList.end(), edgePropertyList.begin(), info.n,
                        edgeList.size());
    TIME_END("Converted into CSR format")

    // calculate maximum cliques
    TIME_START
    getMaximumCliques(graph, ts.len, windowMin, windowMax, 16, true);
    TIME_END("Calculated maximum cliques")
}

void runVLCM(const TimeSeries &ts, int windowMin, int windowMax, double threshold) {
    TIME_INIT

    TIME_START
    TimeSeriesInfo tsInfo(ts, threshold);
    TIME_END("Preprocessed time series")

    double multiplier = getBestMultiplier(threshold);
    int currentW = windowMin;
    int currentMaxW = std::ceil(windowMin * multiplier);
    while (currentW <= windowMax) {
        int maxW = std::min(currentMaxW, windowMax);

        PRINT("Starting windows " << currentW << " .. " << maxW)
        TIME_START
        runSingleVLCM(tsInfo, currentW, maxW);
        TIME_END("Finished windows " << currentW << " .. " << maxW)

        currentW = currentMaxW + 1;
        currentMaxW = std::ceil(currentW * multiplier);
    }
}

void runVLCM_Naive(const TimeSeries &ts, int windowMin, int windowMax, double threshold) {
    TIME_INIT

    TIME_START
    TimeSeriesInfo tsInfo(ts, threshold);
    TIME_END("Preprocessed time series")

    for (int w = windowMin; w <= windowMax; w++) {
        PrecomputedInfo info(tsInfo, w, w);
        auto[edgeList, windowRangesList] = computeDistanceGraphNaive(tsInfo, info);
        DistanceGraph graph(boost::edges_are_sorted, edgeList.begin(), edgeList.end(), windowRangesList.begin(),
                            info.n, edgeList.size());
        auto clique = getMaximumClique(getFullUnionGraph(graph));
        printClique(w, clique);
    }
}

void verify(const TimeSeries &ts, int windowMin, int windowMax, double threshold) {
    if (windowMin < 3) throw VLCMException("Error: window size must be at least 3");
    if (windowMin > windowMax) throw VLCMException("Error: window size smaller (or equal) than maximum window");
    if (ts.size() / 2 < windowMax) throw VLCMException("Error: window size must be smaller than half the time series length");
    if (threshold < -1 || 1 < threshold) throw VLCMException("Error: correlation threshold must be in [-1,1]");
}

int main(int argc, char **argv) {
    // params
    if (argc < 5) {
        std::cout << "Usage: ./VLCM <time-series-path> <window-min> <window-max> <correlation> [vn]" << std::endl;
        return 1;
    }
    std::string tsPath = argv[1];
    int windowMin = std::stoi(argv[2]);
    int windowMax = std::stoi(argv[3]);
    double threshold = std::stod(argv[4]);
    verbose = false;
    bool naive = false;
    if (argc > 5) {
        std::string additionalArgs(argv[5]);
        verbose = additionalArgs.find('v') != std::string::npos;
        naive = additionalArgs.find('n') != std::string::npos;
    }

    TIME_INIT
    try {
        // read time series
        TIME_START
        TimeSeries ts;
        readFile(tsPath, ts);
        TIME_END("Time series with length " << ts.size() << " read from " << tsPath)

        // default: set window max to half ts length
        if (windowMax <= 0) {
            windowMax = static_cast<int>(ts.size() / 2);
        }

        // verify args
        verify(ts, windowMin, windowMax, threshold);
        PRINT("Args: windows=" << windowMin << ".." << windowMax << "; correlation=" << threshold)

        if (naive) {
            // run naive version (CliqueMotif for each length)
            TIME_START
            runVLCM_Naive(ts, windowMin, windowMax, threshold);
            TIME_END("VLCM_Naive finished")
        } else {
            // run Variable-Length CliqueMotif
            TIME_START
            runVLCM(ts, windowMin, windowMax, threshold);
            TIME_END("VLCM finished")
        }
    } catch (const VLCMException &e) {
        PRINT(e.what())
        return 1;
    }

    return 0;
}


