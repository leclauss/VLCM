#include <iostream>
#include "fast_cliques.h"

std::vector<Clique>
getMaximumCliques(const DistanceGraph &graph, int tsLen, int windowMin, int windowMax, int k, bool output) {
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

        if (!canSkip(w)) {

            if (multiLength >= k) {
                auto unionClique = fastMaxClique(getUnionGraph(graph, w, w + multiLength - 1),
                                                 currentClique, w);
                for (int wi = w; wi < w + multiLength; wi++) {
                    cliqueUpperBound[wi] = std::min(cliqueUpperBound[wi], static_cast<int>(unionClique.size()));
                }
                if (!canSkip(w)) {
                    currentClique = fastMaxClique(getGraph(graph, w), currentClique, w);
                    multiLength /= 2;
                }
            } else {
                currentClique = fastMaxClique(getGraph(graph, w), currentClique, w);
            }

        }


        if (currentClique.size() == previousClique.size()) {
            multiLength++;
        }
        multiLength = std::min(multiLength, windowMax - w);
        previousClique = currentClique;
        cliques[i] = currentClique;

        // output clique
        if (output) {
            printClique(w, currentClique);
        }
    }
    return cliques;
}