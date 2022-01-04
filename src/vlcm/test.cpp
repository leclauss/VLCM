#include <iostream>
#include <vector>
#include <chrono>

#include "utils.h"
#include "exception.h"
#include "precompute.h"
#include "distance_graph.h"
#include "lmc.h"
#include "fast_cliques.h"

using namespace std;

const string realDatasetsPathPrefix = "../../../data/";

DistanceGraph getCompleteGraph(const TimeSeries &ts, int windowMin, int windowMax, double threshold) {
    TimeSeriesInfo tsInfo(ts, threshold);

    vector<EdgeList> edgeLists;
    vector<WindowRangesList> windowRangesLists;
    {
        double multiplier = getBestMultiplier(threshold);
        int currentW = windowMin;
        int currentMaxW = std::ceil(windowMin * multiplier);
        while (currentW <= windowMax) {
            int maxW = std::min(currentMaxW, windowMax);
            PrecomputedInfo info(tsInfo, currentW, maxW);
            auto[edgeList, edgePropertyList] = computeDistanceGraph(tsInfo, info);
            edgeLists.push_back(edgeList);
            windowRangesLists.push_back(edgePropertyList);
            currentW = currentMaxW + 1;
            currentMaxW = std::ceil(currentW * multiplier);
        }
    }
    auto[edgeList, edgePropertyList] = mergeEdgeLists(edgeLists, windowRangesLists);
    DistanceGraph graph(boost::edges_are_sorted, edgeList.begin(), edgeList.end(), edgePropertyList.begin(),
                        tsInfo.len - windowMin + 1, edgeList.size());

    return graph;
}

int size(AdjacencyList &graph) {
    int size = 0;
    for (auto &edges: graph) {
        size += edges.size();
    }
    return size / 2;
}

void testGraphSimilarity() {
    cout << "[TEST - Graph Similarity]" << endl;
    // params
    string tsPath = realDatasetsPathPrefix + "nyc_taxi.csv";
    TimeSeries ts;
    readFile(tsPath, ts);

    double threshold = 0.95;
    int windowMin = 10;
    int windowMax = ts.size() / 2;

    auto graph = getCompleteGraph(ts, windowMin, windowMax, threshold);

    // output similarity

    cout << "window;size;next size;union size" << endl;

    AdjacencyList currentGraph = getGraph(graph, windowMin);
    for (int w = windowMin; w < windowMax; w++) {
        AdjacencyList unionGraph = getUnionGraph(graph, w, w + 1);
        int unionSize = size(unionGraph);
        int currentSize = size(currentGraph);
        currentGraph = getGraph(graph, w + 1);
        int nextSize = size(currentGraph);
        cout << w << ";" << currentSize << ";" << nextSize << ";" << unionSize << endl;
    }
}

void testUpperLowerBound() {
    cout << "[TEST - Upper and Lower Bounds]" << endl;
    // params
    string tsPath = realDatasetsPathPrefix + "nyc_taxi.csv";
    TimeSeries ts;
    readFile(tsPath, ts);

    double threshold = 0.95;
    int windowMin = 10;
    int windowMax = ts.size() / 2;

    auto graph = getCompleteGraph(ts, windowMin, windowMax, threshold);

    // output bounds
    cout << "window;lower bound;exact size;upper bound" << endl;

    Clique previousClique = fastMaxClique(getGraph(graph, windowMin), windowMin);

    for (int w = windowMin + 1; w <= windowMax; w++) {
        auto unionGraph = getUnionGraph(graph, w, windowMax);
        auto unionClique = fastMaxClique(unionGraph, windowMin);

        auto lbClique = getMaximumClique(
                getSubgraph(graph, w, expandSubset(previousClique, -1, 0, ts.size() - w + 1)));
        auto exactClique = fastMaxClique(getGraph(graph, w), lbClique, w);
        cout << w << ";" << lbClique.size() << ";" << exactClique.size() << ";" << unionClique.size() << endl;
        previousClique = exactClique;
    }
}

void testGraphReduction() {
    cout << "[TEST - Graph Reduction]" << endl;
    // params
    string tsPath = realDatasetsPathPrefix + "nyc_taxi.csv";
    TimeSeries ts;
    readFile(tsPath, ts);

    double threshold = 0.95;
    int windowMin = 10;
    int windowMax = ts.size() / 2;

    auto graph = getCompleteGraph(ts, windowMin, windowMax, threshold);

    // output graph size and reduced graph size
    cout << "window;graph size;reduced graph size" << endl;

    Clique previousClique = fastMaxClique(getGraph(graph, windowMin), windowMin);
    for (int w = windowMin + 1; w <= windowMax; w++) {
        auto wGraph = getGraph(graph, w);
        int graphSize = size(wGraph);

        auto lbClique = getMaximumClique(
                getSubgraph(graph, w, expandSubset(previousClique, -1, 0, ts.size() - w + 1)));
        auto reducedGraph = reduceGraph(wGraph, lbClique.size() + 1, w);
        int reducedSize = size(reducedGraph);

        cout << w << ";" << graphSize << ";" << reducedSize << endl;
    }
}

void testMultiplierOptimization() {
    cout << "[TEST - Multiplier Optimization]" << endl;
    cout << "dataset;correlation;multiplier;time" << endl;
    // params
    vector<double> corrs = {0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5};
    vector<string> files = {"exchange", "nyc_taxi", "Twitter_AMZN", "machine_temp", "insect_b"}; //"LSF5_10"
    TIME_INIT
    for (double threshold: corrs) {
        for (const string &fileName: files) {
            string tsPath = realDatasetsPathPrefix + fileName + ".csv";
            TimeSeries ts;
            readFile(tsPath, ts);

            int windowMin = 10;
            int windowMax = ts.size() / 2;

            // create graph
            TimeSeriesInfo tsInfo(ts, threshold);
            for (double multiplier = 1.1; multiplier < 2.51; multiplier += 0.05) {
                START
                int currentW = windowMin;
                int currentMaxW = std::ceil(windowMin * multiplier);
                while (currentW <= windowMax) {
                    int maxW = std::min(currentMaxW, windowMax);
                    PrecomputedInfo info(tsInfo, currentW, maxW);
                    auto[edgeList, edgePropertyList] = computeDistanceGraph(tsInfo, info);
                    DistanceGraph graph(boost::edges_are_sorted, edgeList.begin(), edgeList.end(),
                                        edgePropertyList.begin(),
                                        info.n, edgeList.size());
                    currentW = currentMaxW + 1;
                    currentMaxW = std::ceil(currentW * multiplier);
                }
                STOP
                cout << fileName << ";" << threshold << ";" << multiplier << ";" << TIME_MS << endl;
            }
        }
    }
}

void testUnionThresholdOptimization() {
    cout << "[TEST - Union Threshold Optimization]" << endl;
    cout << "dataset;correlation;threshold;time" << endl;
    // params
    vector<double> corrs = {0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5};
    vector<string> files = {"exchange", "nyc_taxi", "Twitter_AMZN", "machine_temp", "insect_b"}; //"LSF5_10"
    TIME_INIT
    for (double threshold: corrs) {
        for (const string &fileName: files) {
            string tsPath = realDatasetsPathPrefix + fileName + ".csv";
            TimeSeries ts;
            readFile(tsPath, ts);

            int windowMin = 10;
            int windowMax = ts.size() / 2;

            auto graph = getCompleteGraph(ts, windowMin, windowMax, threshold);
            // find maximum cliques
            for (int k = 1; k < 5000; k *= 2) {
                START
                getMaximumCliques(graph, ts.size(), windowMin, windowMax, k, false);
                STOP
                cout << fileName << ";" << threshold << ";" << k << ";" << TIME_MS << endl;
            }
        }
    }
}

int main() {
    // run all tests
    testGraphSimilarity();
    testUpperLowerBound();
    testGraphReduction();

    // the following tests might take multiple days/weeks. consider changing the parameters in the respective functions
    testMultiplierOptimization();
    testUnionThresholdOptimization();
    return 0;
}


