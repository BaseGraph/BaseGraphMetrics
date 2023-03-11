#include "BaseGraph/extensions/metrics/undirected.h"
#include "BaseGraph/directedgraph.h"
#include "BaseGraph/types.h"
#include "BaseGraph/undirectedgraph.h"

#include <cstdio>
#include <vector>
#include <list>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <algorithm>


using namespace std;


namespace BaseGraph { namespace metrics {

template<typename T>
static std::list<T> intersectionOf(const std::list<T>& a, const std::list<T>& b);
template <typename T>
static double getAverage(const T& iterable);


size_t countTrianglesAroundVertex(const UndirectedGraph& graph, VertexIndex vertex1){
    size_t triangleNumber = 0;
    auto vertexNeighbourhood = graph.getNeighboursOf(vertex1);

    for(VertexIndex& vertex2: vertexNeighbourhood)
        triangleNumber += intersectionOf(vertexNeighbourhood, graph.getNeighboursOf(vertex2)).size();

    return triangleNumber/2;  // Triangles are all counted twice
}

list<array<VertexIndex, 3>> findAllTriangles(const UndirectedGraph& graph){
    list<array<VertexIndex, 3>> triangles;

    for(VertexIndex& vertex1: graph) {
        const list<VertexIndex>& vertex1Neighbours = graph.getNeighboursOf(vertex1);

        for (const VertexIndex& vertex2: vertex1Neighbours) {

            if (vertex1 < vertex2)
                for (const VertexIndex& vertex3: intersectionOf(vertex1Neighbours, graph.getNeighboursOf(vertex2)))
                    if (vertex2 < vertex3)
                        triangles.push_back({vertex1, vertex2, vertex3});
        }
    }
    return triangles;
}

size_t countTriangles(const UndirectedGraph& graph){
    size_t triangleTotal = 0;
    for(VertexIndex& vertex: graph)
        triangleTotal += countTrianglesAroundVertex(graph, vertex);

    return triangleTotal/3;
}

vector<double> getDegreeDistribution(const UndirectedGraph &graph) {
    vector<double> degreeDistribution(graph.getSize());

    size_t n=graph.getSize();
    for (VertexIndex& vertex: graph)
        degreeDistribution[vertex] = (double) graph.getDegreeOf(vertex)/n;
    return degreeDistribution;
}

double getGlobalClusteringCoefficient(const UndirectedGraph& graph) {
    vector<size_t> vertexTriangleNumbers;
    vertexTriangleNumbers.resize(graph.getSize());

    for (VertexIndex& vertex: graph)
        vertexTriangleNumbers[vertex] = countTrianglesAroundVertex(graph, vertex);
    return getGlobalClusteringCoefficient(graph, vertexTriangleNumbers);
}

double getGlobalClusteringCoefficient(const UndirectedGraph& graph, const vector<size_t>& vertexTriangleNumbers) {
    size_t vertexDegree;
    double globalTriangleNumber = 0;
    double globalWedgeNumber = 0;

    for (VertexIndex& vertex: graph) {
        vertexDegree = graph.getDegreeOf(vertex);

        if(vertexDegree > 1)
            globalWedgeNumber += (double) vertexDegree*(vertexDegree - 1) / 2;
        globalTriangleNumber += vertexTriangleNumbers[vertex];
    }
    return globalTriangleNumber/globalWedgeNumber;
}

vector<double> getRedundancy(const UndirectedGraph& graph) {
    vector<double> localClusteringCoefficients = getLocalClusteringCoefficients(graph);

    for (VertexIndex vertex: graph)
        localClusteringCoefficients[vertex] *= (double) (graph.getDegreeOf(vertex)-1);
    return localClusteringCoefficients;
}

vector<double> getLocalClusteringCoefficients(const UndirectedGraph& graph) {
    vector<double> localClusteringCoefficients;
    localClusteringCoefficients.resize(graph.getSize(), 0);

    int vertexDegree;
    double triangleNumber;

    for(VertexIndex& vertex: graph) {
        vertexDegree = graph.getDegreeOf(vertex);
        triangleNumber = countTrianglesAroundVertex(graph, vertex);

        if(vertexDegree > 1)
          localClusteringCoefficients[vertex] = 2.0*triangleNumber / vertexDegree / (vertexDegree - 1);
        else
          localClusteringCoefficients[vertex] = 0;
    }
    return localClusteringCoefficients;
}

unordered_map<size_t, double> getClusteringSpectrum(const UndirectedGraph& graph) {
    unordered_map<size_t, double> clusteringSpectrum;
    unordered_map<size_t, size_t> numberOfSummedValues;

    vector<double> localClusteringCoefficients = getLocalClusteringCoefficients(graph);

    size_t degree;

    for (VertexIndex& vertex: graph) {
        degree = graph.getDegreeOf(vertex);
        if (degree < 2) continue;

        if (clusteringSpectrum.find(degree) == clusteringSpectrum.end()) {
            clusteringSpectrum[degree] = localClusteringCoefficients[vertex];
            numberOfSummedValues[degree] = 1;
        }
        else {
            clusteringSpectrum[degree] += localClusteringCoefficients[vertex];
            numberOfSummedValues[degree]++;
        }
    }

    for (auto& degree_clustering: clusteringSpectrum)
        degree_clustering.second /= numberOfSummedValues[degree_clustering.first];

    return clusteringSpectrum;
}

vector<size_t> getKShells(const UndirectedGraph& graph){
    return getKShellsAndOnionLayers(graph).first;
}

list<VertexIndex> getKCore(const UndirectedGraph& graph, size_t k){
    return getKCore(k, getKShells(graph));
}

list<VertexIndex> getKCore(size_t k, const vector<size_t>& kshells){
    list<VertexIndex> kcore;

    for (size_t i=0; i<kshells.size(); i++)
        if (kshells[i] <= k)
            kcore.push_back(i);
    return kcore;
}

vector<size_t> getOnionLayers(const UndirectedGraph& graph){
    return getKShellsAndOnionLayers(graph).second;
}

pair<vector<size_t>, vector<size_t>> getKShellsAndOnionLayers(const UndirectedGraph& graph) {
    // Algorithm of Batagelj and Zaversnik modified by Hébert-Dufresne, Grochow and Allard.
    vector<size_t> verticesOnionLayer(graph.getSize());
    vector<size_t> verticesKShell(graph.getSize());

    vector<size_t> effectiveDegrees = graph.getDegrees();

    // set<effective degree, vertex index>: Sort vertices by degree
    set<pair<size_t, VertexIndex>> higherLayers;
    list<pair<size_t, VertexIndex>> currentLayer;

    for (VertexIndex& vertex: graph)
        higherLayers.insert({effectiveDegrees[vertex], vertex});


    size_t onionLayerDegree;
    size_t onionLayer = 0;

    while(!higherLayers.empty()) {
        onionLayer += 1;
        onionLayerDegree = higherLayers.begin()->first;

        for (auto it=higherLayers.begin(); it!=higherLayers.end() && it->first==onionLayerDegree; ) {
            const auto& vertex = it->second;
            verticesKShell[vertex] = onionLayerDegree;
            verticesOnionLayer[vertex] = onionLayer;

            currentLayer.push_back(*it);
            higherLayers.erase(it++);
        }

        // Ajust layers neighbours' effective degree
        while(!currentLayer.empty()) {
            const auto& vertex = currentLayer.begin()->second;

            for (const VertexIndex& neighbour: graph.getNeighboursOf(vertex)) {
                const auto& neighbourDegree = effectiveDegrees[neighbour];

                auto it = higherLayers.find({neighbourDegree, neighbour});
                if(it != higherLayers.end() && neighbourDegree > onionLayerDegree) {
                    higherLayers.erase(it);
                    higherLayers.insert( {neighbourDegree-1, neighbour} );
                    effectiveDegrees[neighbour]--;
                }
            }
            currentLayer.erase(currentLayer.begin());
        }
    }
    return {verticesKShell, verticesOnionLayer};
}

list<size_t> getNeighbourhoodDegreesOfVertex(const UndirectedGraph& graph, VertexIndex vertex) {
    list<size_t> neighbourDegrees;

    for (const VertexIndex& neighbour: graph.getNeighboursOf(vertex))
        neighbourDegrees.push_back( graph.getDegreeOf(neighbour) );

    return neighbourDegrees;
}

vector<double> getNeighbourDegreeSpectrum(const UndirectedGraph &graph, bool normalized) {
    vector<double> degreeSpectrum(graph.getSize());

    for (VertexIndex& vertex: graph)
        degreeSpectrum[vertex] = getAverage(getNeighbourhoodDegreesOfVertex(graph, vertex));

    if (normalized) {
        double firstMoment = 0;
        double secondMoment = 0;
        size_t degree;
        for (VertexIndex& vertex: graph) {
            degree = graph.getDegreeOf(vertex);
            firstMoment += degree;
            secondMoment += degree*degree;
        }
        for (VertexIndex& vertex: graph)
            degreeSpectrum[vertex] *= firstMoment/secondMoment;
    }
    return degreeSpectrum;
}

unordered_map<size_t, list<double>> getOnionSpectrum(const UndirectedGraph& graph) {
    auto kshells_onionLayers = getKShellsAndOnionLayers(graph);
    return getOnionSpectrum(graph, kshells_onionLayers.first, kshells_onionLayers.second);
}

unordered_map<size_t, list<double>> getOnionSpectrum(const UndirectedGraph& graph, const vector<size_t>& kshells,
        const vector<size_t>& onionLayers) {
    if (graph.getSize()!=kshells.size() || graph.getSize()!=onionLayers.size())
        throw logic_error("The onion layers vector and the k-shells vector must be the size of the graph");

    size_t onionLayersNumber = *max_element(onionLayers.begin(), onionLayers.end());
    unordered_map<size_t, list<double>> onionSpectrum;


    vector<size_t> onionLayerToKShell(onionLayersNumber + 1);
    vector<size_t> onionLayerSizes(onionLayersNumber + 1, 0);

    for(VertexIndex& vertex: graph) {
        const size_t& layer = onionLayers[vertex];
        onionLayerToKShell[layer] = kshells[vertex];
        onionLayerSizes[layer] += 1;
    }

    for(size_t layer=1; layer<=onionLayersNumber; layer++)
        onionSpectrum[onionLayerToKShell[layer]].push_back((double) onionLayerSizes[layer] / graph.getSize());

    return onionSpectrum;
}

double getDegreeCorrelation(const UndirectedGraph& graph) {
    return getDegreeCorrelation(graph, getAverage(graph.getDegrees()));
}

double getDegreeCorrelation(const UndirectedGraph& graph, double averageDegree) {
    size_t degree;
    size_t n = graph.getSize();

    vector<double> excessDegreeDistribution(1);

    double excessDegree;
    double firstMoment = 0;
    double secondMoment = 0;

    for (VertexIndex& vertex: graph) {
        degree = graph.getDegreeOf(vertex);
        if (degree == 0) continue;

        if (degree>excessDegreeDistribution.size())
            excessDegreeDistribution.resize(degree, 0);

        excessDegree = degree/(averageDegree*n);
        excessDegreeDistribution[degree-1] += excessDegree;

        firstMoment += (degree-1)*excessDegree;
        secondMoment += (degree-1)*(degree-1)*excessDegree;
    }
    double excessDegreeVariance = secondMoment - firstMoment*firstMoment;

    double degreeCorrelationCoefficient = 0;
    size_t neighbourDegree;
    size_t edgeNumber = graph.getEdgeNumber();
    for (VertexIndex& vertex: graph) {
        degree = graph.getDegreeOf(vertex);
        if (degree < 2) continue;

        for (const VertexIndex& neighbour: graph.getNeighboursOf(vertex)) {
            if (vertex > neighbour) {
                neighbourDegree = graph.getDegreeOf(neighbour);
                degreeCorrelationCoefficient += (double) (degree-1)*(neighbourDegree-1)/edgeNumber;
            }
        }
    }

    size_t maxDegree = excessDegreeDistribution.size();
    auto& edDistribution = excessDegreeDistribution;
    for (size_t degree=1; degree<maxDegree; degree++){
        degreeCorrelationCoefficient -= degree*degree * edDistribution[degree]*edDistribution[degree];

        for (size_t degree2=degree+1; degree2<maxDegree; degree2++)
            degreeCorrelationCoefficient -= 2*degree*degree2 * edDistribution[degree]*edDistribution[degree2];
    }

    return degreeCorrelationCoefficient/excessDegreeVariance;
}

double getModularity(const UndirectedGraph& graph, const vector<size_t>& vertexCommunities) {
    if (graph.getSize() == 0) throw std::logic_error("Graph is empty");
    if (vertexCommunities.size() != graph.getSize()) throw std::logic_error("Vertex communities vector must be the size of the graph");

    size_t communityNumber = *std::max_element(vertexCommunities.begin(), vertexCommunities.end());
    size_t intraCommunityStubs = 0;
    vector<size_t> communityDegreeSum(communityNumber+1, 0);

    for (VertexIndex& vertex: graph) {
        communityDegreeSum[vertexCommunities[vertex]] += graph.getDegreeOf(vertex);

        for (const VertexIndex& neighbour: graph.getNeighboursOf(vertex))
            if (vertexCommunities[vertex] == vertexCommunities[neighbour])
                intraCommunityStubs++;
    }
    double modularity = 0;
    size_t edgeNumber = graph.getEdgeNumber();

    modularity += (double) intraCommunityStubs/(2*edgeNumber);

    for (size_t& degreeSum: communityDegreeSum)
        modularity -= (degreeSum/(2.*edgeNumber)) * (degreeSum/(2.*edgeNumber));

    return modularity;
}


static size_t count4Cliques(const UndirectedGraph& graph, std::vector<int>& X, const std::unordered_set<VertexIndex>& triangles) {
    size_t cliqueCount = 0;

    for (auto w: triangles) {
        for (auto r: graph.getNeighboursOf(w))
            if (X[r]==2)
                cliqueCount++;
        X[w] = 0;
    }
    return cliqueCount;
}

static size_t countCycles(const UndirectedGraph& graph, std::vector<int>& X, const std::unordered_set<VertexIndex>& twoStars_u) {
    size_t cycleCount = 0;

    for (auto w: twoStars_u) {
        for (auto r: graph.getNeighboursOf(w))
            if (X[r]==3)
                cycleCount++;
        X[w] = 0;
    }
    return cycleCount;
}

static inline size_t choose2(size_t value) {
    return (double)value*(value-1)/2;
}

static inline size_t choose3(size_t value) {
    return (double)value*(value-1)*(value-2)/6;
}

static inline size_t choose4(size_t value) {
    return (double) value*(value-1)*(value-2)*(value-3)/24;
}

std::unordered_map<std::string, size_t> countMotifs(const UndirectedGraph& graph, size_t maxMotifSize) {
    if (maxMotifSize > 4)
        fprintf(stderr, "Warning: only motifs of size 3 and 4 are supported.\n");

    std::unordered_map<std::string, size_t> motifFrequencies = {
        {"3:clique", 0}, {"3:2-star", 0}, // connected
        {"3:edge", 0}, {"3:empty", 0}, // disconnected
        {"4:clique", 0}, {"4:chordalcycle", 0}, {"4:tailedtriangle", 0}, {"4:cycle", 0}, {"4:3-star", 0}, {"4:path", 0}, // connected
        {"4:triangle", 0}, {"4:2-star", 0}, {"4:2 edges", 0}, {"4:edge", 0}, {"4:empty", 0} // disconnected
    };

    std::unordered_set<VertexIndex> twoStars_u, twoStars_v, triangles, neighbourhood;
    std::vector<int> X(graph.getSize(), 0);
    auto edgeNumber = graph.getEdgeNumber();

    size_t N_TT(0), N_SuSv(0), N_TSuVSv(0), N_SS(0), N_TI(0), N_SuVSvI(0), N_II(0), N_II1(0);
    size_t N_edge_3;

    for (const auto& edge: graph.edges) {
        twoStars_u.clear(); twoStars_v.clear(); triangles.clear();

        auto& u = edge.first;
        auto& v = edge.second;
        neighbourhood = {u, v};

        for (auto w: graph.getNeighboursOf(u)) {
            if (w==v)
                continue;
            X[w] = 1;
            twoStars_u.insert(w);
            neighbourhood.insert(w);
        }

        for (auto w: graph.getNeighboursOf(v)) {
            if (w==u)
                continue;

            neighbourhood.insert(w);
            if (X[w]==1) {
                X[w] = 2;
                triangles.insert(w);
                twoStars_u.erase(w);
            }
            else {
                twoStars_v.insert(w);
                X[w] = 3;
            }
        }
        N_edge_3 = graph.getSize() - neighbourhood.size();

        motifFrequencies["3:clique"] += triangles.size();
        motifFrequencies["3:2-star"] += twoStars_u.size()+twoStars_v.size();
        motifFrequencies["3:edge"] += N_edge_3;

        if (maxMotifSize == 3)
            continue;

        motifFrequencies["4:clique"] += count4Cliques(graph, X, triangles);
        motifFrequencies["4:cycle"] += countCycles(graph, X, twoStars_u);

        N_TT += choose2(triangles.size());
        N_SuSv += twoStars_u.size()*twoStars_v.size();
        N_TSuVSv += triangles.size()*(twoStars_u.size() + twoStars_v.size());
        N_SS += choose2(twoStars_u.size()) + choose2(twoStars_v.size());

        N_TI += triangles.size()*N_edge_3;
        N_SuVSvI += N_edge_3*(twoStars_u.size() + twoStars_v.size());
        N_II += choose2(N_edge_3);
        N_II1 += edgeNumber - graph.getDegreeOf(u, false) - graph.getDegreeOf(v, false) + 1;

        for (auto w: graph.getNeighboursOf(v))
            X[w] = 0;
    }
    motifFrequencies["3:clique"] /= 3;
    motifFrequencies["3:2-star"] /= 2;
    motifFrequencies["3:empty"] = choose3(graph.getSize())
        -motifFrequencies["3:clique"]-motifFrequencies["3:2-star"]-motifFrequencies["3:edge"];

    if (maxMotifSize == 3)
        return motifFrequencies;

    motifFrequencies["4:clique"] /= 6;
    motifFrequencies["4:chordalcycle"] = N_TT - 6*motifFrequencies["4:clique"];
    motifFrequencies["4:tailedtriangle"] = N_TSuVSv-4*motifFrequencies["4:chordalcycle"];
    motifFrequencies["4:cycle"] /= 4;
    motifFrequencies["4:3-star"] = N_SS - motifFrequencies["4:tailedtriangle"];
    motifFrequencies["4:path"] = N_SuSv - 4*motifFrequencies["4:cycle"];
    motifFrequencies["4:triangle"] = N_TI - motifFrequencies["4:tailedtriangle"];
    motifFrequencies["4:2-star"] = N_SuVSvI - 2*motifFrequencies["4:path"];
    motifFrequencies["4:2 edges"] = N_II1;
    motifFrequencies["4:edge"] = N_II - 2*motifFrequencies["4:2 edges"];

    size_t nonEmptyMotifs4 = 0;
    for (auto motif_count: motifFrequencies)
        if (motif_count.first[0]=='4' && motif_count.first != "4:empty")
            nonEmptyMotifs4 += motif_count.second;

    motifFrequencies["4:empty"] = choose4(graph.getSize()) - nonEmptyMotifs4;

    return motifFrequencies;
}


// From https://stackoverflow.com/questions/38993415/how-to-apply-the-intersection-between-two-lists-in-c
template<typename T>
static std::list<T> intersectionOf(const std::list<T>& a, const std::list<T>& b){
    std::list<T> rtn;
    std::unordered_multiset<T> st;
    std::for_each(a.begin(), a.end(), [&st](const T& k){ st.insert(k); });
    std::for_each(b.begin(), b.end(),
        [&st, &rtn](const T& k){
            auto iter = st.find(k);
            if(iter != st.end()){
                rtn.push_back(k);
                st.erase(iter);
            }
        }
    );
    return rtn;
}

template <typename T>
static double getAverage(const T& iterable) {
    if (iterable.size() == 0)
        return 0;

    size_t sum=0;
    for (const auto& element: iterable)
        sum += element;
    return (double) sum/iterable.size();
}


}} // namespace BaseGraph::metrics
