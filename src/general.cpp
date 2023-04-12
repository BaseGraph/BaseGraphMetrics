#include "BaseGraph/extensions/metrics/general.h"

#include <list>
#include <queue>
#include <unordered_map>
#include <vector>

#include "BaseGraph/algorithms/paths.hpp"
#include "BaseGraph/directed_graph.hpp"
#include "BaseGraph/types.h"
#include "BaseGraph/undirected_graph.hpp"

using namespace std;

namespace BaseGraph {
namespace metrics {

template <typename T>
static double getClosenessCentralityOfVertex(const T &graph,
                                             VertexIndex vertex) {
    vector<size_t> shortestPathLengths =
        getShortestPathLengthsFromVertex(graph, vertex);
    size_t componentSize = 0;
    unsigned long long int sum = 0;

    for (VertexIndex &vertex : graph) {
        if (shortestPathLengths[vertex] != algorithms::BASEGRAPH_VERTEX_MAX) {
            componentSize += 1;
            sum += shortestPathLengths[vertex];
        }
    }
    return sum > 0 ? ((double)componentSize - 1) / sum : 0;
}

template <typename T> vector<double> getClosenessCentralities(const T &graph) {
    vector<double> closenessCentralities(graph.getSize(), 0);

    for (VertexIndex vertex : graph)
        closenessCentralities[vertex] =
            getClosenessCentralityOfVertex(graph, vertex);

    return closenessCentralities;
}

template <typename T>
static double getHarmonicCentralityOfVertex(const T &graph,
                                            VertexIndex vertex) {
    vector<size_t> shortestPathLengths =
        getShortestPathLengthsFromVertex(graph, vertex);

    double harmonicSum = 0;
    for (VertexIndex &vertex : graph)
        if (shortestPathLengths[vertex] != 0 &&
            shortestPathLengths[vertex] != algorithms::BASEGRAPH_VERTEX_MAX)
            harmonicSum += 1.0 / shortestPathLengths[vertex];

    return harmonicSum;
}

template <typename T> vector<double> getHarmonicCentralities(const T &graph) {
    vector<double> harmonicCentralities(graph.getSize(), 0);

    for (VertexIndex vertex : graph)
        harmonicCentralities[vertex] =
            getHarmonicCentralityOfVertex(graph, vertex);

    return harmonicCentralities;
}

template <>
vector<double> getBetweennessCentralities(const DirectedGraph &graph,
                                          bool normalizeWithGeodesicNumber) {
    size_t verticesNumber = graph.getSize();
    vector<double> betweennesses;
    betweennesses.resize(verticesNumber, 0);

    algorithms::MultiplePredecessors distancesPredecessors;
    list<list<VertexIndex>> currentGeodesics;
    for (const VertexIndex &i : graph) {
        distancesPredecessors =
            algorithms::findAllPredecessorsOfVertex(graph, i);
        for (const VertexIndex &j : graph) {
            currentGeodesics =
                algorithms::findMultiplePathsToVertexFromPredecessors(
                    graph, i, j, distancesPredecessors);
            if (currentGeodesics.empty())
                continue; // vertices i and j are not in the same component

            for (auto &geodesic : currentGeodesics) {
                for (auto &vertexOnGeodesic : geodesic) {
                    if (vertexOnGeodesic != i && vertexOnGeodesic != j) {
                        if (normalizeWithGeodesicNumber)
                            betweennesses[vertexOnGeodesic] +=
                                1. / currentGeodesics.size();
                        else
                            betweennesses[vertexOnGeodesic] += 1;
                    }
                }
            }
        }
    }
    return betweennesses;
}

template <>
vector<double> getBetweennessCentralities(const UndirectedGraph &graph,
                                          bool normalizeWithGeodesicNumber) {
    size_t verticesNumber = graph.getSize();
    vector<double> betweennesses;
    betweennesses.resize(verticesNumber, 0);

    algorithms::MultiplePredecessors distancesPredecessors;
    list<list<VertexIndex>> currentGeodesics;
    for (const VertexIndex &i : graph) {
        distancesPredecessors =
            algorithms::findAllPredecessorsOfVertex(graph, i);
        for (const VertexIndex &j : graph) {
            if (i >= j)
                continue;

            currentGeodesics =
                algorithms::findMultiplePathsToVertexFromPredecessors(
                    graph, i, j, distancesPredecessors);
            if (currentGeodesics.empty())
                continue; // vertices i and j are not in the same component

            for (auto &geodesic : currentGeodesics) {
                for (auto &vertexOnGeodesic : geodesic) {
                    if (vertexOnGeodesic != i && vertexOnGeodesic != j) {
                        if (normalizeWithGeodesicNumber)
                            betweennesses[vertexOnGeodesic] +=
                                1. / currentGeodesics.size();
                        else
                            betweennesses[vertexOnGeodesic] += 1;
                    }
                }
            }
        }
    }
    return betweennesses;
}

template <typename T>
std::vector<size_t> getShortestPathLengthsFromVertex(const T &graph,
                                                     VertexIndex source) {
    return algorithms::findPredecessorsOfVertex(graph, source).first;
}

template <typename T> vector<size_t> getDiameters(const T &graph) {
    size_t verticesNumber = graph.getSize();
    vector<size_t> diameters(verticesNumber);

    vector<size_t> shortestPathLengths;
    size_t largestDistance;
    for (const VertexIndex &i : graph) {
        shortestPathLengths = getShortestPathLengthsFromVertex(graph, i);
        largestDistance = shortestPathLengths[0];

        for (const VertexIndex &j : graph)
            if (shortestPathLengths[j] > largestDistance &&
                shortestPathLengths[j] != algorithms::BASEGRAPH_VERTEX_MAX)
                largestDistance = shortestPathLengths[j];

        if (largestDistance == algorithms::BASEGRAPH_VERTEX_MAX)
            diameters[i] = 0;
        else
            diameters[i] = largestDistance;
    }
    return diameters;
}

template <typename T>
static double getShortestPathAverageOfVertex(const T &graph,
                                             VertexIndex vertex) {
    vector<size_t> shortestPathLengths =
        getShortestPathLengthsFromVertex(graph, vertex);

    size_t sum = 0;
    size_t componentSize = 1;

    for (VertexIndex &vertex : graph)
        if (shortestPathLengths[vertex] != 0 &&
            shortestPathLengths[vertex] != algorithms::BASEGRAPH_VERTEX_MAX) {
            sum += shortestPathLengths[vertex];
            componentSize++;
        }

    return componentSize > 1 ? (double)sum / (componentSize - 1) : 0;
}

template <typename T> vector<double> getShortestPathAverages(const T &graph) {
    vector<double> averageShortestPaths(graph.getSize(), 0);

    for (VertexIndex vertex : graph)
        averageShortestPaths[vertex] =
            getShortestPathAverageOfVertex(graph, vertex);

    return averageShortestPaths;
}

template <typename T>
vector<unordered_map<size_t, double>>
getShortestPathsDistribution(const T &graph) {
    auto connectedComponents = findConnectedComponents(graph);

    vector<size_t> shortestPathLengths;
    vector<unordered_map<size_t, double>> shortestPathDistribution(
        connectedComponents.size());
    size_t componentIndex = 0;

    for (auto component : connectedComponents) {
        auto &currentDistribution = shortestPathDistribution[componentIndex];

        if (component.size() > 1) {
            for (const VertexIndex &vertex : component) {
                shortestPathLengths =
                    getShortestPathLengthsFromVertex(graph, vertex);

                for (const size_t &pathLength : shortestPathLengths) {
                    if (pathLength != 0 &&
                        pathLength != algorithms::BASEGRAPH_VERTEX_MAX) {
                        if (currentDistribution.find(pathLength) ==
                            currentDistribution.end())
                            currentDistribution[pathLength] = 1;
                        else
                            currentDistribution[pathLength]++;
                    }
                }
            }
            for (auto &element : currentDistribution)
                element.second /= component.size();
        }

        componentIndex++;
    }
    return shortestPathDistribution;
}

template <typename T>
static double getShortestPathHarmonicAverageOfVertex(const T &graph,
                                                     VertexIndex vertex) {
    vector<size_t> shortestPathLengths =
        getShortestPathLengthsFromVertex(graph, vertex);
    size_t componentSize = 0;

    double sumOfInverse = 0;
    for (VertexIndex &vertex : graph) {
        if (shortestPathLengths[vertex] != 0 &&
            shortestPathLengths[vertex] != algorithms::BASEGRAPH_VERTEX_MAX) {
            componentSize += 1;
            sumOfInverse += 1.0 / shortestPathLengths[vertex];
        }
    }
    return componentSize > 0 ? sumOfInverse / componentSize : 0;
}

template <typename T>
vector<double> getShortestPathHarmonicAverages(const T &graph) {
    vector<double> harmonicAverages(graph.getSize(), 0);

    for (VertexIndex vertex : graph)
        harmonicAverages[vertex] =
            getShortestPathHarmonicAverageOfVertex(graph, vertex);

    return harmonicAverages;
}

template <typename T> list<Component> findConnectedComponents(const T &graph) {
    size_t verticesNumber = graph.getSize();
    if (verticesNumber == 0)
        throw logic_error("There are no vertices.");

    list<Component> connectedComponents;
    Component currentComponent;
    VertexIndex currentVertex, startVertex;

    queue<VertexIndex> verticesToProcess;
    vector<bool> processedVertices;
    bool allVerticesProcessed = false;
    processedVertices.resize(verticesNumber, false);

    list<VertexIndex>::const_iterator vertexNeighbour;

    while (!allVerticesProcessed) {
        allVerticesProcessed = true;
        for (VertexIndex i = 0; i < verticesNumber && allVerticesProcessed;
             ++i) {
            if (!processedVertices[i]) {
                allVerticesProcessed = false;
                startVertex = i;
            }
        }

        if (!allVerticesProcessed) {
            currentComponent.clear();
            verticesToProcess.push(startVertex);
            processedVertices[startVertex] = true;

            while (!verticesToProcess.empty()) {
                currentVertex = verticesToProcess.front();

                for (const VertexIndex &vertexNeighbour :
                     graph.getOutEdgesOf(currentVertex)) {
                    if (!processedVertices[vertexNeighbour]) {
                        verticesToProcess.push(vertexNeighbour);
                        processedVertices[vertexNeighbour] = true;
                    }
                }
                currentComponent.push_back(currentVertex);
                verticesToProcess.pop();
            }
            connectedComponents.push_back(currentComponent);
        }
    }
    return connectedComponents;
}

// Allowed classes for metrics

template vector<double> getClosenessCentralities(const DirectedGraph &);
template vector<double> getClosenessCentralities(const UndirectedGraph &);
template vector<double> getHarmonicCentralities(const DirectedGraph &);
template vector<double> getHarmonicCentralities(const UndirectedGraph &);

template vector<size_t> getShortestPathLengthsFromVertex(const DirectedGraph &,
                                                         VertexIndex);
template vector<size_t>
getShortestPathLengthsFromVertex(const UndirectedGraph &, VertexIndex);
template vector<size_t> getDiameters(const DirectedGraph &);
template vector<size_t> getDiameters(const UndirectedGraph &);
template vector<double> getShortestPathAverages(const DirectedGraph &);
template vector<double> getShortestPathAverages(const UndirectedGraph &);
template vector<double> getShortestPathHarmonicAverages(const DirectedGraph &);
template vector<double>
getShortestPathHarmonicAverages(const UndirectedGraph &);
template vector<unordered_map<size_t, double>>
getShortestPathsDistribution(const DirectedGraph &);
template vector<unordered_map<size_t, double>>
getShortestPathsDistribution(const UndirectedGraph &);

template list<Component> findConnectedComponents(const DirectedGraph &);
template list<Component> findConnectedComponents(const UndirectedGraph &);

} // namespace metrics
} // namespace BaseGraph
