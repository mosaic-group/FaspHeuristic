//
// Created by gonciarz on 2019-03-07.
//

#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include "graph.h"
#include "graphExt.h"
#include "tools/easylogging++.h"
#include "tools/dynamicBitset.h"

#include <cassert>
#include <type_traits>
#include <queue>
#include <random>
#include <cstdint>
#include <functional>


namespace Graph::Tools {
    /**
     * Finds all vertices accessible from start vertex using DFS method
     * @param aGraph input graph
     * @param aSrc start vertex
     * @return sorted set of vertices
     */
    template<typename VERTEX_TYPE>
    static typename Graph<VERTEX_TYPE>::VerticesSet
    depthFirstSearch(const Graph<VERTEX_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                     bool aReversedSearch = false) {

        std::size_t size = aGraph.getNumOfVertices();

        typename Graph<VERTEX_TYPE>::VerticesSet visited;
        typename Graph<VERTEX_TYPE>::Vertices stack;
        stack.reserve(size);
        stack.emplace_back(aSrc);

        while (!stack.empty()) {
            const auto currentVertex = stack.back();
            stack.pop_back();
            visited.emplace(currentVertex);

            // find all not visted outgoing vertices and add them to stack
            const auto &vertices = aReversedSearch ? aGraph.getInVertices(currentVertex) : aGraph.getOutVertices(
                    currentVertex);
            for (const auto &v : vertices) {
                if (visited.find(v) == visited.end()) stack.emplace_back(v);
            }
        }

        return visited;
    }

    template<typename VERTEX_TYPE>
    static bool pathExistsDFS(const Graph<VERTEX_TYPE> &aGraph,
                              const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                              const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                              bool aReversedSearch = false) {
        if (aSrc == aDst) return true;

        std::size_t size = aGraph.getNumOfVertices();

        std::unordered_set<VERTEX_TYPE> visited;
        visited.reserve(size);

        typename Graph<VERTEX_TYPE>::Vertices stack;
        stack.reserve(size);
        stack.emplace_back(aSrc);

        while (!stack.empty()) {
            const auto currentVertex = stack.back();
            stack.pop_back();
            visited.emplace(currentVertex);

            // find all not visted vertices and add them to stack
            const auto &vertices = aReversedSearch ? aGraph.getInVertices(currentVertex) : aGraph.getOutVertices(
                    currentVertex);
            for (const auto &v : vertices) {
                if (v == aDst) return true;
                else if (visited.find(v) == visited.end()) stack.emplace_back(v);
            }
        }

        return false;
    }


    /**
     * Finds all edges with cycles, that is if we have edge a->b there is a path from b to a
     * @param aGraph input graph
     * @return container with edges
     */
    template<typename VERTEX_TYPE>
    static typename Graph<VERTEX_TYPE>::Edges
    findEdgesWithCycles(const Graph<VERTEX_TYPE> &aGraph) {
        typename Graph<VERTEX_TYPE>::Edges edges;
        for (const auto &e : aGraph.getEdges()) {
            if (pathExistsDFS(aGraph, e.dst, e.src)) {
                edges.emplace_back(e);
            }
        }
        return edges;
    }

    /**
     * Finds all edges with cycles, that is if we have edge a->b there is a path from b to a
     * @param aGraph input graph
     * @return true if there are still cycles
    */
    template<typename VERTEX_TYPE>
    static bool isAcyclic(const Graph<VERTEX_TYPE> &aGraph) {
        for (const auto &e : aGraph.getEdges()) {
            if (pathExistsDFS(aGraph, e.dst, e.src)) return false;
        }

        return true;
    }

    template<typename VERTEX_TYPE, typename EDGE_PROP_TYPE>
    static auto findPathWithPositiveCapacity(const Graph<VERTEX_TYPE> &aGraph,
                                             const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                             const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                             const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {

        typename Graph<VERTEX_TYPE>::Vertices path;
        if (aGraph.hasVertex(aSrc) && aGraph.hasVertex(aDst)) {
            if (aSrc == aDst) {
                path.emplace_back(aSrc);
                return std::pair{true, path};
            }

            auto parents = Ext::getVertexProperties<VERTEX_TYPE>(aGraph);
            typename Graph<VERTEX_TYPE>::VerticesSet visited;

            std::queue<VERTEX_TYPE> stack;
            stack.push(aSrc);

            while (!stack.empty()) {
                auto currentVertex = stack.front();
                stack.pop();
                visited.emplace(currentVertex);

                if (currentVertex == aDst) {
                    // Traverse path back to source and build the path
                    path.emplace_back(currentVertex);
                    while (true) {
                        currentVertex = parents.at(currentVertex);
                        path.emplace_back(currentVertex);
                        if (currentVertex == aSrc) break;
                    }
                    // Finally reverse it to have path from src to dst
                    std::reverse(path.begin(), path.end());
                    return std::pair{true, path};
                };

                const auto &vertices = aGraph.getOutVertices(currentVertex);
                for (const auto &v : vertices) {
                    if (visited.find(v) == visited.end() && aWeights.at({currentVertex, v}) > 0) {
                        parents.at(v) = currentVertex;
                        stack.emplace(v);
                    }
                }
            }
        }
        return std::pair{false, path};
    }

    /**
     * Calculate minimum s-t cut via Ford-Fulkerson max flow - min cut metaGraphWorker
     * @param aGraph input graph
     * @param aSrc source vertex
     * @param aDst destination vertex
     * @param aWeights weights of edges
     * @return max flow value, graph and capacities after algorithm ends
     */
    template<typename VERTEX_TYPE, typename EDGE_PROP_TYPE>
    static auto minStCutFordFulkersonBase(const Graph<VERTEX_TYPE> &aGraph,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                          const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        // Copy graph and weights - they will be modified
        auto g{aGraph};
        auto c{aWeights};
        // Add edges in backwerd direction if not yet exists
        for (const auto &e : g.getEdges()) {
            auto backwardEdge = typename Graph<VERTEX_TYPE>::Edge{e.dst, e.src};
            if (!g.hasEdge(backwardEdge)) {
                g.addEdge(backwardEdge);
                c.emplace(std::move(backwardEdge), 0);
            }
        }

        EDGE_PROP_TYPE maxFlow = 0;

        while (true) {
            // Find new path from source to sink, if not found than full capacity of network is reached
            auto[pathExists, path] = findPathWithPositiveCapacity(g, aSrc, aDst, c);
            if (!pathExists) break;
            // Find edge wiht minimum capacity left
            EDGE_PROP_TYPE pathFlow = std::numeric_limits<EDGE_PROP_TYPE>::max();
            for (std::size_t i = 0; i < path.size() - 1; ++i) {
                pathFlow = std::min(pathFlow, c[{path[i], path[i + 1]}]);
            }

            // Update max flow and capacities of edges on the found path
            maxFlow += pathFlow;
            for (std::size_t i = 0; i < path.size() - 1; ++i) {
                c.at({path[i], path[i + 1]}) -= pathFlow;
                c.at({path[i + 1], path[i]}) += pathFlow;
            }
        }

        return std::tuple{maxFlow, g, c};
    }

    /**
     * Calculate minimum s-t cut via Ford-Fulkerson max flow - min cut metaGraphWorker
     * @param aGraph input graph
     * @param aSrc source vertex
     * @param aDst destination vertex
     * @param aWeights weights of edges
     * @return max flow value
     */
    template<typename VERTEX_TYPE, typename EDGE_PROP_TYPE>
    static auto minStCutFordFulkerson(const Graph<VERTEX_TYPE> &aGraph,
                                      const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                      const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                      const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        auto maxFlow = std::get<0>(minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights));
        return maxFlow;
    }

    /**
     * Calculate minimum s-t cut via Ford-Fulkerson max flow - min cut metaGraphWorker
     * @param aGraph input graph
     * @param aSrc source vertex
     * @param aDst destination vertex
     * @param aWeights weights of edges
     * @return max flow value, set of s-t cut edges
     */
    template<typename VERTEX_TYPE, typename EDGE_PROP_TYPE>
    static auto minStCutFordFulkersonEdges(const Graph<VERTEX_TYPE> &aGraph,
                                           const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                           const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                           const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        auto[maxFlow, g, c] = minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights);
        typename Graph<VERTEX_TYPE>::Edges edges;

        // Find all edges that are fully used in terms of capacity,
        // those witch are reachable from source to edge beginning
        // and not rachable from source to edge end are those that belong to min-cut set
        for (const auto &e : aGraph.getEdges()) {
            if (aWeights.at(e) > 0 &&
                c.at(e) == 0 &&
                findPathWithPositiveCapacity(g, aSrc, e.dst, c).first == false &&
                findPathWithPositiveCapacity(g, aSrc, e.src, c).first == true) {
                edges.emplace_back(e);
            }
        }
        return std::pair{maxFlow, edges};
    }

    /**
     * Removes up to requested number of edges with cycles (it may remove less).
     * TODO: Edges with cycles need to be found only first time, then the only cycles left might be a part of left edges only.
     *       So we could speedup procedure to search only within edgesWithCycles
     */
    template<typename VERTEX_TYPE>
    static auto getRandomSubgraph(const Graph<VERTEX_TYPE> &aGraph, int aNumEdgesToRemove) {
        auto g{aGraph};
        int edgesRemovedCnt = 0;
        typename Graph<VERTEX_TYPE>::Edges edges;
        while (edgesRemovedCnt < aNumEdgesToRemove) {
            auto edgesWithCycles = Tools::findEdgesWithCycles(g);
            auto n = edgesWithCycles.size();
            if (n == 0) {
                LOG(TRACE) << "Removed less edges than requested (" << edgesRemovedCnt << ", " << aNumEdgesToRemove
                           << ")";
                break;
            } else {
                auto rndEdge = edgesWithCycles[rand() % n];
                g.removeEdge(rndEdge);
                ++edgesRemovedCnt;
                edges.emplace_back(std::move(rndEdge));
            }
        }

        return std::pair{g, edges};
    }

    /**
     * Generates a graph using Erods-Renyi ranodm model
     * @tparam EDGE_PROP_TYPE
     * @tparam VERTEX_TYPE
     * @tparam GRAPH_TYPE
     * @param aNumOfVertices number of vertices in a graph
     * @param aProbabilityOfEdge probability in G(n,p) Erdos-Renyi model
     * @return generated graph
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
    static auto generateErdosRenyiGraph(int aNumOfVertices, float aProbabilityOfEdge) {
        Graph<VERTEX_TYPE> g;
        Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> c;

        // add vertices in range: 0..aNumOfVertices-1
        for (int i = 0; i < aNumOfVertices; ++i) g.addVertex(i);


        int numOfArcs = 0;

        std::random_device rd;
        std::mt19937 mt(rd());

        // Go through all possible arcs in graph and check if should be generated
        for (int i = 0; i < aNumOfVertices - 1; ++i) {
            for (int j = i + 1; j < aNumOfVertices; ++j) {
                typename Graph<VERTEX_TYPE>::Edge e1 = {i, j};
                if (std::uniform_real_distribution<>(0.0, 1.0)(mt) < aProbabilityOfEdge) {
                    c[e1] = 1;
                    g.addEdge(std::move(e1));
                    ++numOfArcs;
                }
                typename Graph<VERTEX_TYPE>::Edge e2 = {j, i};
                if (std::uniform_real_distribution<>(0.0, 1.0)(mt) < aProbabilityOfEdge) {
                    c[e2] = 1;
                    g.addEdge(std::move(e2));
                    ++numOfArcs;
                }
            }
        }
        LOG(INFO) << "Erdos-Renyi graph: #v=" << aNumOfVertices << " #e=" << numOfArcs << " #theoreticalEdges="
                  << ((aNumOfVertices) * (aNumOfVertices - 1) * aProbabilityOfEdge);
        return std::pair{g, c};
    }

    // https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    template<typename VERTEX_TYPE>
    static auto stronglyConnectedComponents(const Graph<VERTEX_TYPE> &aGraph) {
        int index_counter = 0;
        const int numOfV = aGraph.getNumOfVertices();
        Ext::VertexProperties<VERTEX_TYPE, VERTEX_TYPE> lowLinks; lowLinks.reserve(numOfV);
        Ext::VertexProperties<VERTEX_TYPE, VERTEX_TYPE> index;    index.reserve(numOfV);
        std::unordered_set<VERTEX_TYPE> stackHelper;    stackHelper.reserve(numOfV);
        typename Graph<VERTEX_TYPE>::Vertices stack;    stack.reserve(numOfV);
        std::vector<std::unordered_set<VERTEX_TYPE>> result; result.reserve(numOfV);

        // worker function that is called recursively
        std::function<void(const VERTEX_TYPE &)> strongconnect = [&](const VERTEX_TYPE &node) {
            index[node] = index_counter;
            lowLinks[node] = index_counter;
            ++index_counter;
            stack.emplace_back(node);
            stackHelper.emplace(node);

            for (const auto &successor : aGraph.getOutVertices(node)) {
                if (lowLinks.find(successor) == lowLinks.end()) {
                    // successor has not yet been visited; recurse on it
                    strongconnect(successor);
                    lowLinks[node] = std::min(lowLinks.at(node), lowLinks.at(successor));
                }
                else if (stackHelper.find(successor) != stackHelper.end()) {
                        // the successor is in the stack and hence in the current strongly connected component (SCC)
                        lowLinks[node] = std::min(lowLinks.at(node), index.at(successor));
                }
            }

            // If `node` is a root node, pop the stack and generate an SCC
            if (lowLinks.at(node) == index.at(node)) {
                std::unordered_set<VERTEX_TYPE> connectedComponent;
                connectedComponent.reserve(stack.size());

                while (true) {
                    auto successor = stack.back();
                    stack.pop_back();
                    stackHelper.erase(successor);

                    connectedComponent.emplace(successor);
                    if (successor == node) break;
                }
                result.emplace_back(std::move(connectedComponent));
            }
        };

        for (auto &node : aGraph.getVertices()) {
            if (lowLinks.find(node) == lowLinks.end()) {
                strongconnect(node);
            }
        }
//        std::cout << "SCC: " << result << std::endl;
        return result;
    }

    // https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    template<typename VERTEX_TYPE>
    static auto stronglyConnectedComponents2(const Graph<VERTEX_TYPE> &aGraph) {
        auto vertices = aGraph.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        auto maxV = (maxId == vertices.end() ? 0 : *maxId) + 1;
        DynamicBitset <uint64_t, VERTEX_TYPE> sh(maxV);

        int index_counter = 0;
        const int numOfV = aGraph.getNumOfVertices();
        std::vector<int16_t> lowLinks(maxV, -1);
        std::vector<int16_t> index(maxV, -1);
        typename Graph<VERTEX_TYPE>::Vertices stack;    stack.reserve(numOfV);
        std::vector<std::unordered_set<VERTEX_TYPE>> result; result.reserve(numOfV);

        std::function<void(const VERTEX_TYPE &)> strongconnect = [&](const VERTEX_TYPE &node) {
            typename Graph<VERTEX_TYPE>::Vertices sv; sv.reserve(numOfV);
            std::vector<VERTEX_TYPE> si; si.reserve(numOfV);
            std::vector<bool> init; init.reserve(numOfV);

            sv.push_back(node);
            si.push_back(0);
            init.push_back(true);
            while (sv.size() > 0) {

                VERTEX_TYPE currentNode = sv.back(); sv.pop_back();
                int ci = si.back(); si.pop_back();
                bool initRun = init.back(); init.pop_back();

            processSuccessor:

                if (initRun) {
                    index[currentNode] = index_counter;
                    lowLinks[currentNode] = index_counter;
                    ++index_counter;
                    stack.emplace_back(currentNode);
                    sh.set(currentNode);
                }
                auto ov = aGraph.getOutVertices(currentNode);
                for(int i = ci; i < ov.size(); ++i) {
                    auto successor = ov[i];
                    if (lowLinks[successor] == -1) {
                            // save the state (it would be recurrent call in default version of Trajan's algorithm)
                            sv.push_back(currentNode);
                            si.push_back(i + 1);
                            init.push_back(false);
                            // set values for successor and repeat from beginning ('goto' is bad... I know).
                            currentNode = successor;
                            ci = 0;
                            initRun = true;
                            goto processSuccessor;
                    }
                     else if (sh.test(successor)) {
                        // the successor is in the stack and hence in the current strongly connected component (SCC)
                        lowLinks[currentNode] = std::min(lowLinks.at(currentNode), index[successor]);
                    }
                }

                // If `node` is a root node, pop the stack and generate an SCC
                if (lowLinks.at(currentNode) == index[currentNode]) {
                    std::unordered_set<VERTEX_TYPE> connectedComponent; connectedComponent.reserve(stack.size());

                    while (true) {
                        auto successor = stack.back();
                        stack.pop_back();
                        sh.clear(successor);
                        connectedComponent.emplace(successor);
                        if (successor == currentNode) break;
                    }
                    result.emplace_back(std::move(connectedComponent));
                }
                if (sv.size() > 0) {
                    auto predecessor = sv.back();
                    lowLinks[predecessor] = std::min(lowLinks.at(predecessor), lowLinks.at(currentNode));
                }
            }
        };

        for (auto &node : aGraph.getVertices()) {
            if (lowLinks[node] == -1) {
                strongconnect(node);
            }
        }

        return result;
    }

    /**
     * Generates property container with integer random weights in range [min, max]
     * @tparam EDGE_PROP_TYPE
     * @tparam VERTEX_TYPE
     * @tparam GRAPH_TYPE
     * @param graph input graph
     * @param min minimum weight value
     * @param max maximum weight value
     * @return container with random weights for every edge
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
    auto getRandomWeights(const Graph<VERTEX_TYPE> &graph, const EDGE_PROP_TYPE min, const EDGE_PROP_TYPE max) {
        // helper for easier generation of random values in range [min, max]
        auto randInt = [](EDGE_PROP_TYPE min, EDGE_PROP_TYPE max) -> EDGE_PROP_TYPE {
            static std::mt19937 mt(std::random_device{}());
            return std::uniform_int_distribution<EDGE_PROP_TYPE>(min, max)(mt);
        };

        auto c = Ext::getEdgeProperties(graph, 0);
        for (auto& [k, v] : c) {
            v = randInt(min, max);
        }

        return c;
    }
}


#endif
