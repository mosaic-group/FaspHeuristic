#ifndef GRAPHFASPFAST_H
#define GRAPHFASPFAST_H


#include "graphToIgraph.h"
#include "graph.h"
#include "graphExt.h"
#include "graphTools.h"
#include "tools/tools.h"
#include "tools/easylogging++.h"
#include "tools/prettyprint.h"
#include "tools/timer.h"
#include "tools/dynamicBitset.h"
#include "tools/stack.h"
#include <algorithm>
#include <future>

namespace Graph::FaspFast {

    /**
     * Exception class thrown in case never-ending mincut loop
     */
    class MinCutFailed : public std::exception {};

    /**
     * PathHero... class is a approach to speed up things. For all algorithms it has implemented it uses
     * pre-allocated memory - that is why it need to be initialized with maximum number of vertices (or if not
     * in sequence 0, 1, 2 ... then with a maximum value of vertex id + 1.
     * @tparam VERTEX_TYPE
     */
    template <typename VERTEX_TYPE>
    class alignas(128) PathHero {

        // Allocation of structures/memory/containers shared by all algorithms from PathHero class
        DynamicBitset<uint32_t, uint16_t> iVisited;
        Stack<VERTEX_TYPE> stack;
        std::vector<VERTEX_TYPE> parents;
        typename Graph<VERTEX_TYPE>::Edges inEdges;
        typename Graph<VERTEX_TYPE>::Edges outEdges;
        typename Graph<VERTEX_TYPE>::Vertices path;



    public:
        inline static int cnt = 0;
        explicit PathHero(size_t aN) : iVisited(aN), stack(aN), parents(aN, 0) {
            // both can have at most aN in/outgoing edges
            inEdges.reserve(aN);
            outEdges.reserve(aN);
            path.reserve(aN);
        }

        /**
         * @return true if path from src to dst exists
         */
        template<template<typename> class GRAPH_TYPE>
        bool pathExistsDFS(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                           const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                           const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                           bool aReversedSearch = false) {
            cnt++;
            if (aSrc == aDst) return true;

            iVisited.clearAll();

            stack.clear();
            stack.push(aSrc);
            iVisited.set(aSrc);

            while (!stack.empty()) {
                const auto currentVertex = stack.pop();

                // find all not visted vertices and add them to stack
                const auto &vertices = aReversedSearch ? aGraph.getInVertices(currentVertex) : aGraph.getOutVertices(currentVertex);
                for (const auto &v : vertices) {
                    if (v == aDst) return true;
                    else if (iVisited.test(v) == 0) {
                        stack.push(v);
                        iVisited.set(v);
                    }
                }
            }

            return false;
        }

        /**
         * @return container with all verticess accessible from src vertex (if reverse search is set to true then it goes against edge direction)
         */
        template<template<typename> class GRAPH_TYPE>
        auto depthFirstSearch(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                              const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                              bool aReversedSearch = false) {
            iVisited.clearAll();

            stack.clear();
            stack.push(aSrc);
            iVisited.set(aSrc);

            while (!stack.empty()) {
                const auto currentVertex = stack.pop();

                // find all not visted outgoing vertices and add them to stack
                const auto &vertices = aReversedSearch ? aGraph.getInVertices(currentVertex) : aGraph.getOutVertices(currentVertex);
                for (const auto &v : vertices) {
                    if (iVisited.test(v) == 0) {
                        stack.push(v);
                        iVisited.set(v);
                    }
                }
            }

            return DynamicBitset{iVisited};
        }

        /**
         * Finds path with positive (>0) capacity from src to dst
         */
        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        auto findPathWithPositiveCapacity(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                          const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            path.clear();
            if (aGraph.hasVertex(aSrc) && aGraph.hasVertex(aDst)) {
                if (aSrc == aDst) {
                    path.emplace_back(aSrc);
                    return true;
                }

                iVisited.clearAll();

                stack.clear();
                stack.push(aSrc);
                iVisited.set(aSrc);

                while (!stack.empty()) {
                    auto currentVertex = stack.pop();
                    if (currentVertex == aDst) {
                        // Traverse path back to source and build the path
                        path.emplace_back(currentVertex);
                        while (true) {
                            currentVertex = parents[currentVertex];
                            path.emplace_back(currentVertex);
                            if (currentVertex == aSrc) break;
                        }
                        // Finally reverse it to have path from src to dst
                        std::reverse(path.begin(), path.end());
                        return true;
                    };

                    const auto &vertices = aGraph.getOutVertices(currentVertex);
                    for (const auto &v : vertices) {
                        if (iVisited.test(v) == 0 && aWeights.at({currentVertex, v}) > 0) {
                            parents[v] = currentVertex;
                            stack.push(v);
                            iVisited.set(v);
                        }
                    }
                }
            }
            return false;
        }

        /**
         * Calculate minimum s-t cut via Ford-Fulkerson max flow - min cut metaGraphWorker
         * @param aGraph input graph
         * @param aSrc source vertex
         * @param aDst destination vertex
         * @param aWeights weights of edges
         * @return max flow value, graph and capacities after algorithm ends
         */
        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        auto minStCutFordFulkersonBase(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
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
            int cnt = 0;
            while(true) {
                // Find new path from source to sink, if not found than full capacity of network is reached
                auto pathExists = findPathWithPositiveCapacity(g, aSrc, aDst, c);
                if (!pathExists) break;
                // Find edge wiht minimum capacity left
                EDGE_PROP_TYPE pathFlow = std::numeric_limits<EDGE_PROP_TYPE>::max();
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    pathFlow = std::min(pathFlow, c[{path[i], path[i + 1]}]);
                }

                // Update max flow and capacities of edges on the found path
                maxFlow += pathFlow;
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    c.at({path[i], path[i + 1]}) -= pathFlow;
                    c.at({path[i + 1], path[i]}) += pathFlow;
                }

                if (++cnt > 10000) throw MinCutFailed();
            }
            return std::tuple{maxFlow, g, c};
        }

        /**
         * Implementation of maximum flow / min cut algorithm:
         * https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
         * @return maxFlow value
         */
        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        auto minStCutFordFulkerson(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                   const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
            auto maxFlow = std::get<0>(minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights));
            return maxFlow;
        }

        /**
         * G - cleaning graph phase
         * @param aGraph - input graph (it will be modified!)
         * @param aEdge  - edge that will be used as a starting place for cleaning
         */
        template<template <typename> class GRAPH_TYPE>
        void G(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge aEdge) {
            // Remove outgoing edges from destination and ingoing to source.
            for (const auto &vo : aGraph.getOutVertices(aEdge.dst)) {
                aGraph.removeEdge({aEdge.dst, vo});
            }
            for (const auto &vi : aGraph.getInVertices(aEdge.src)) {
                aGraph.removeEdge({vi, aEdge.src});
            }

            // Remove all vertices which are not accessible from source (regular search) and not accessible from destination
            // (reverse search - against edge direction). This will leave as with a graph with all paths from srouce to dest.
            // We need at least 2 vertices: source, destination (source might be equal to destination in case looped edge)
            // and at least one more vertex to be verified, if less then nothing to be done here
            auto numOfVertices = aGraph.getNumOfVertices();
            if (numOfVertices >= 2) {
                auto A = depthFirstSearch(aGraph, aEdge.src);
                auto B = depthFirstSearch(aGraph, aEdge.dst, true /* reverseSearch */);

                const auto &vertices = aGraph.getVertices();

                for (auto v : vertices) {
                    if (!(A.test(v) && B.test(v))) {
                        // Removes all vertices which are in set V\(A∩B)
                        aGraph.removeVertex(v);
                    }
                }
            }
        }

        /**
         * G* - cleaning graph phase
         * @param aGraph - input graph (it will be modified!)
         * @param aEdge  - edge that will be used as a starting place for cleaning
         */
        template<template <typename> class GRAPH_TYPE>
        void GStar(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge) {
            auto vFrom = aEdge.dst;
            auto vTo = aEdge.src;

            while (true) {
                // Run cleanup phase 'G' which will remove not reachable edges/vertices.
                // NOTE: G() will remove aEdge from aGraph
                G(aGraph, {vFrom, vTo});
                bool wasGraphModified = false;

                // For every edge in a graph remove:
                // - all edges 'e' which after removing its input edges does not have a path
                //   from e.dst to vTo (so this are loops through which we would need to go back)
                // - all edges 'e' which after removing its output edges does not have a path
                //   from vFrom to e.src
                for (const auto &e : aGraph.getEdges()) {
                    bool wasCurrentEdgeRemoved = false;
                    const auto &inVertices = aGraph.getInVertices(e.src);
                    inEdges.clear();
                    for (const auto &v : inVertices) {
                        typename Graph<VERTEX_TYPE>::Edge inEdge{v, e.src};
                        if (e != inEdge) inEdges.emplace_back(std::move(inEdge));
                    }
                    aGraph.removeEdges(inEdges);
                    if (e.src == e.dst || !pathExistsDFS(aGraph, e.dst, vTo)) {
                        aGraph.removeEdge(e);
                        wasGraphModified = true;
                        wasCurrentEdgeRemoved = true;
                    }
                    aGraph.addEdges(inEdges);

                    if (!wasCurrentEdgeRemoved) {
                        const auto &outVertices = aGraph.getOutVertices(e.dst);
                        outEdges.clear();
                        for (const auto &v : outVertices) {
                            typename Graph<VERTEX_TYPE>::Edge outEdge{e.dst, v};
                            if (e != outEdge) outEdges.emplace_back(std::move(outEdge));
                        }
                        aGraph.removeEdges(outEdges);
                        if (e.src == e.dst || !pathExistsDFS(aGraph, vFrom, e.src)) {
                            aGraph.removeEdge(e);
                            wasGraphModified = true;
                        }
                        aGraph.addEdges(outEdges);
                    }
                }

                if (!wasGraphModified) break;
            }
        }

        /**
         * step2b - cleaning graph phase (Yes! we need better name for it).
         * Finds all edges which are part of other cycle (not only cycle going through aEdge).
         */
        template<template <typename> class GRAPH_TYPE>
        auto step2b(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                    const Graph<VERTEX_TYPE, GRAPH_TYPE> &aCleanedGraph,
                    const typename Graph<VERTEX_TYPE>::Edge aEdge) {
            // Remove edge to not find any path going through it
            aGraph.removeEdge(aEdge);

            typename Graph<VERTEX_TYPE>::Edges S;
            for (const auto &h : aCleanedGraph.getEdges()) {
                if (pathExistsDFS(aGraph, h.dst, h.src)) {
                    S.emplace_back(h);
                }
            }

            // revert
            aGraph.addEdge(aEdge);
            return S;
        }
    };

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithm(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        auto outGraph{aGraph}; // eventually acyclic graph
        while(true) {
            bool wasGraphModified = false;

            for (const auto &e : outGraph.getEdges()) {

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                path.GStar(workGraph, e);

                auto S = path.step2b(outGraph, workGraph, e);
                workGraph.removeEdges(S);
                path.GStar(workGraph, e);

                auto mc = path.minStCutFordFulkerson(workGraph, e.dst, e.src, aWeights);

                if (mc >= aWeights.at(e)) {
                    wasGraphModified = true;
                    outGraph.removeEdge(e);
                    removedEdges.emplace_back(e);
                }
            }

            if (!wasGraphModified) break;
        }

        return removedEdges;
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithm(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        auto vertices = aGraph.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

        return superAlgorithm(aGraph, aWeights, path);
    }


    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto deltaFASP(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        auto g{aGraph};
        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        auto vertices = aGraph.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

        auto edgesToRemove = superAlgorithm(g, aWeights, path);
        g.removeEdges(edgesToRemove);
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

        while (true) {
            if (Tools::isAcyclic(g)) break;

            EDGE_PROP_TYPE capacityToRemove = 0;
            typename Graph<VERTEX_TYPE>::Edge edgeToRemove;
            bool firstLoop = true;
            const auto &edges = Tools::findEdgesWithCycles(g);
            for (const auto &e : edges) {
                auto gs{g};
                path.GStar(gs, e);

                auto mc = path.minStCutFordFulkerson(gs, e.dst, e.src, aWeights);
                const auto delta = (mc - aWeights.at(e));
                if (firstLoop || delta > capacityToRemove) {
                    firstLoop = false;
                    capacityToRemove = delta;
                    edgeToRemove = e;
                }
            }

            g.removeEdge(edgeToRemove);
            removedEdges.emplace_back(std::move(edgeToRemove));

            auto etr = superAlgorithm(g, aWeights, path);
            if (etr.size() > 0) {
                g.removeEdges(etr);
                removedEdges.reserve(removedEdges.size() + etr.size());
                removedEdges.insert(removedEdges.end(), etr.begin(), etr.end());
            }
        }

        EDGE_PROP_TYPE capacity = 0;
        for (const auto &e : removedEdges) {
            capacity += aWeights.at(e);
        }
        LOG(DEBUG) << "PATH called: " << path.cnt ;
        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(DELTA) capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto randomFASP(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        constexpr int numOfReps = 20;
        auto g{aGraph};

        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        // find max vertex id in graph and setup PathHero
        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path{static_cast<size_t>(maxId == vertices.end() ? 1 : *maxId + 1)};

        // initial run of superAlgorithm (SA)
        auto edgesToRemove = superAlgorithm(g, aWeights, path);
        g.removeEdges(edgesToRemove);
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

//        srand(time(NULL));

        int numEdgesToRemove = 1;
        while (true) {
            if (Tools::isAcyclic(g)) break;

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;

            std::future<typename Graph<VERTEX_TYPE>::Edges> tasks[numOfReps];
            for (int i = 0; i < numOfReps; ++i) {
                tasks[i] = std::async(std::launch::async,
                                      [&]() {
                                          auto[randGraph, _] = Tools::getRandomSubgraph(g, numEdgesToRemove);
                                          return superAlgorithm(randGraph, aWeights);
                                      });
            }
            for (int i = 0; i < numOfReps; ++i) {
                auto edgesToRemove = tasks[i].get();
                for (auto &e : edgesToRemove) {
                    edgesCnt.try_emplace(e, 0).first->second++;
                }
            }
            if (edgesCnt.size() == 0) {
                numEdgesToRemove++;
                continue;
            }
            numEdgesToRemove = 1;

            if (edgesCnt.size() > 0) {
                auto maxCnt = std::max_element(edgesCnt.begin(), edgesCnt.end(), [](const auto &a, const auto &b) -> bool { return a.second < b.second; });

                // Remove best candidate
                g.removeEdge(maxCnt->first);
                removedEdges.emplace_back(std::move(maxCnt->first));

                // Try to solve the rest with superAlgoritm - maybe it will now succeed!
                auto edgesToRemove = superAlgorithm(g, aWeights, path);
                g.removeEdges(edgesToRemove);
                removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
                removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());
            }
            else {
                LOG(FATAL) << "SHOULD NOT HAPPEN!";
                break;
            }

        }

        EDGE_PROP_TYPE capacity = 0;
        for (const auto &e : removedEdges) {
            capacity += aWeights.at(e);
        }

        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }
}


#endif
