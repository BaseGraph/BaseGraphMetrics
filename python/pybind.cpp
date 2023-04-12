#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "BaseGraph/algorithms/paths.hpp"
#include "BaseGraph/directed_graph.hpp"
#include "BaseGraph/extensions/metrics/directed.h"
#include "BaseGraph/extensions/metrics/general.h"
#include "BaseGraph/extensions/metrics/undirected.h"
#include "BaseGraph/undirected_graph.hpp"

namespace py = pybind11;
using namespace BaseGraph;

// From https://github.com/pybind/pybind11/issues/1042
template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence &&seq) {
    // Move entire object to heap (Ensure is moveable!). Memory handled via
    // Python capsule
    Sequence *seq_ptr = new Sequence(std::move(seq));
    auto capsule = py::capsule(
        seq_ptr, [](void *p) { delete reinterpret_cast<Sequence *>(p); });
    return py::array(seq_ptr->size(), // shape of array
                     seq_ptr->data(), // c-style contiguous strides for Sequence
                     capsule          // numpy array references this parent
    );
}

// Metrics that aren't validated with Networkx are tagged with /**/
PYBIND11_MODULE(_metrics, m) {
    // Required Python import the module to work
    py::module_::import("basegraph");

    // General metrics
    m.def("get_closeness_centralities", [&](const DirectedGraph &graph) {
        return as_pyarray(metrics::getClosenessCentralities(graph));
    });
    m.def("get_closeness_centralities", [&](const UndirectedGraph &graph) {
        return as_pyarray(metrics::getClosenessCentralities(graph));
    });
    m.def("get_harmonic_centralities", [&](const DirectedGraph &graph) {
        return as_pyarray(metrics::getHarmonicCentralities(graph));
    });
    m.def("get_harmonic_centralities", [&](const UndirectedGraph &graph) {
        return as_pyarray(metrics::getHarmonicCentralities(graph));
    });
    m.def("get_betweenness_centralities",
          [&](const DirectedGraph &graph, bool normalize) {
              return as_pyarray(
                  metrics::getBetweennessCentralities(graph, normalize));
          });
    m.def("get_betweenness_centralities",
          [&](const UndirectedGraph &graph, bool normalize) {
              return as_pyarray(
                  metrics::getBetweennessCentralities(graph, normalize));
          });
    m.def("get_shortest_path_lengths_from_vertex",
          [&](const DirectedGraph &graph, VertexIndex vertex) {
              return as_pyarray(
                  metrics::getShortestPathLengthsFromVertex(graph, vertex));
          });
    m.def("get_shortest_path_lengths_from_vertex",
          [&](const UndirectedGraph &graph, VertexIndex vertex) {
              return as_pyarray(
                  metrics::getShortestPathLengthsFromVertex(graph, vertex));
          });
/**/m.def("get_diameters", [](const DirectedGraph& graph) {
              return as_pyarray(metrics::getDiameters(graph));
          });
/**/m.def("get_diameters", [](const UndirectedGraph& graph) {
              return as_pyarray(metrics::getDiameters(graph));
          });
    m.def("get_shortest_path_averages", [](const DirectedGraph& graph) {
              return as_pyarray(metrics::getShortestPathAverages(graph));
          });
    m.def("get_shortest_path_averages", [](const UndirectedGraph& graph) {
              return as_pyarray(metrics::getShortestPathAverages(graph));
          });
    m.def("get_shortest_path_harmonic_averages",
          py::overload_cast<const DirectedGraph &>(
              &metrics::getShortestPathHarmonicAverages<DirectedGraph>));
    m.def("get_shortest_path_harmonic_averages",
          py::overload_cast<const UndirectedGraph &>(
              &metrics::getShortestPathHarmonicAverages<UndirectedGraph>));

    /**/ m.def("get_shortest_paths_distribution",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::getShortestPathsDistribution<DirectedGraph>));
    /**/ m.def("get_shortest_paths_distribution",
               py::overload_cast<const UndirectedGraph &>(
                   &metrics::getShortestPathsDistribution<UndirectedGraph>));
    /**/ m.def("find_connected_components",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::findConnectedComponents<DirectedGraph>));
    m.def("find_connected_components",
          py::overload_cast<const UndirectedGraph &>(
              &metrics::findConnectedComponents<UndirectedGraph>));

    // Undirected metrics
    m.def("get_degree_correlation", py::overload_cast<const UndirectedGraph &>(
                                        &metrics::getDegreeCorrelation));
    m.def("find_all_triangles", &metrics::findAllTriangles);
    m.def("count_triangles_around_vertex",
          &metrics::countTrianglesAroundVertex);
    m.def("count_triangles", &metrics::countTriangles);

    m.def("get_local_clustering_coefficients",
          py::overload_cast<const UndirectedGraph &>(
              &metrics::getLocalClusteringCoefficients));
    m.def("get_global_clustering_coefficient",
          py::overload_cast<const UndirectedGraph &>(
              &metrics::getGlobalClusteringCoefficient));
    m.def("get_clustering_spectrum", &metrics::getClusteringSpectrum);
    m.def("get_redundancy", &metrics::getRedundancy);

    m.def("get_kshells_and_onion_layers", &metrics::getKShellsAndOnionLayers);
    m.def("get_kshells", &metrics::getKShells);
    m.def("get_onion_layers", &metrics::getOnionLayers);
    m.def("get_onion_spectrum", py::overload_cast<const UndirectedGraph &>(
                                    &metrics::getOnionSpectrum));
    /**/ m.def("get_kcore", py::overload_cast<const UndirectedGraph &, size_t>(
                                &metrics::getKCore));
    m.def("get_neighbourhood_degrees_of_vertex",
          &metrics::getNeighbourhoodDegreesOfVertex);
    m.def("get_neighbourhood_degree_spectrum",
          &metrics::getNeighbourDegreeSpectrum);

    m.def("get_modularity", &metrics::getModularity);

    // Directed metrics
    m.def("get_density", &metrics::getDensity);
    /**/ m.def("find_all_directed_triangles",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::findAllDirectedTriangles));
    /**/ m.def("get_triangle_spectrum", &metrics::getTriangleSpectrum);
    m.def("get_undirected_local_clustering_coefficients",
          py::overload_cast<const DirectedGraph &>(
              &metrics::getUndirectedLocalClusteringCoefficients));
    m.def("get_undirected_global_clustering_coefficient",
          py::overload_cast<const DirectedGraph &>(
              &metrics::getUndirectedGlobalClusteringCoefficient));

    m.def("get_reciprocity", &metrics::getReciprocity);
    /**/ m.def("get_reciprocal_degrees", &metrics::getReciprocalDegrees);
    /**/ m.def("get_jaccard_reciprocities",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::getJaccardReciprocities));
    /**/ m.def("get_reciprocity_ratios",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::getReciprocityRatios));
    /**/ m.def("get_out_degree_histogram", &metrics::getOutDegreeHistogram);
    /**/ m.def("get_in_degree_histogram",
               py::overload_cast<const DirectedGraph &>(
                   &metrics::getInDegreeHistogram));
}
