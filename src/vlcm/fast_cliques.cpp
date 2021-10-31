#include <iostream>
#include "fast_cliques.h"
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

std::vector<Clique>
getMaximumCliques(const DistanceGraph &graph, int tsLen, int windowMin, int windowMax, bool output) {
    std::vector<Clique> cliques(windowMax - windowMin + 1);

    std::vector<int> cliqueUpperBound(windowMax + 1, tsLen / windowMin);
    std::vector<int> cliqueLowerBound(windowMax + 1, 1);

    auto canSkip = [&cliqueLowerBound, &cliqueUpperBound](int w) {
        return cliqueLowerBound[w] == cliqueUpperBound[w];
    };

    int multiLength = 1;

    Clique previousClique = fastMaxClique(getGraph(graph, windowMin), windowMin);
    cliques[0] = previousClique;
    if (output) {
        printClique(windowMin, previousClique);
    }
    for (int w = windowMin + 1, i = 1; w <= windowMax; w++, i++) {

        // get large initial clique from previous clique
        auto currentClique = getMaximumClique(getSubgraph(graph, w,
                                                          expandSubset(previousClique, -1,
                                                                       0, tsLen - w + 1)));
        cliqueLowerBound[w] = currentClique.size();
        //cout << w << ": " << cliqueLowerBound[w] << " - " << cliqueUpperBound[w] << endl;

        if (!canSkip(w)) {

            if (multiLength >= 6) {
                auto unionClique = fastMaxClique(getUnionGraph(graph, w, w + multiLength - 1),
                                                 currentClique, w);
                for (int wi = w; wi < w + multiLength; wi++) {
                    cliqueUpperBound[wi] = std::min(cliqueUpperBound[wi], static_cast<int>(unionClique.size()));
                }
                //cout << "UB " << w << " - " << (w+multiLength-1) << " : " << unionClique.size() << endl;
                if (!canSkip(w)) {
                    currentClique = fastMaxClique(getGraph(graph, w), currentClique, w);
                    multiLength /= 2;
                    //cout << "  FAIL -> " << currentClique.size() << endl;
                }
            } else {
                currentClique = fastMaxClique(getGraph(graph, w), currentClique, w);
                //cout << "  C -> " << currentClique.size() << endl;
            }

        }


        if (currentClique.size() == previousClique.size()) {
            multiLength++;
        } else {
            //multiLength = max(multiLength / 2, 1);
        }
        multiLength = std::min(multiLength, windowMax - w);
        //cout << "  " << currentClique.size() << "  multiLength " << multiLength << endl;
        previousClique = currentClique;
        cliques[i] = currentClique;

        // output clique
        if (output) {
            printClique(w, currentClique);
        }
    }
    return cliques;
}