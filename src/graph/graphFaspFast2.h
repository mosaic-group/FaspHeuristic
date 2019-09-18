#ifndef GRAPHFASPFAST2_H
#define GRAPHFASPFAST2_H


#include "graphToIgraph.h"
#include "graph.h"
#include "graphExt.h"
#include "graphTools.h"
#include "graphFasp.h"
#include "tools/tools.h"
#include "tools/easylogging++.h"
#include "tools/prettyprint.h"
#include "tools/timer.h"
#include "tools/dynamicBitset.h"
#include "tools/stack.h"
#include <algorithm>
#include <future>

namespace Graph::FaspFast2 {

    /**
     * Exception class thrown in case never-ending mincut loop
     */
    class MinCutFailed : public std::exception {
    };

    template <typename VERTEX_TYPE>
    using EdgesSet = std::unordered_set<typename Graph<VERTEX_TYPE>::Edge, Ext::EdgeHasher<VERTEX_TYPE>>;

    /**
     * PathHero... class is a approach to speed up things. For all algorithms it has implemented it uses
     * pre-allocated memory - that is why it need to be initialized with maximum number of vertices (or if not
     * in sequence 0, 1, 2 ... then with a maximum value of vertex id + 1.
     * @tparam VERTEX_TYPE
     */
    template<typename VERTEX_TYPE>
    class PathHero {

        // Allocation of structures/memory/containers shared by all algorithms from PathHero class
        DynamicBitset<uint32_t, uint16_t> iVisited;
        Stack<VERTEX_TYPE> stack;
        std::vector<VERTEX_TYPE> parents;
        typename Graph<VERTEX_TYPE>::Edges inEdges;
        typename Graph<VERTEX_TYPE>::Edges outEdges;
        typename Graph<VERTEX_TYPE>::Vertices path;


    public:
        volatile inline static std::atomic_int saCnt = 0;

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
        bool pathExistsDFS(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                           const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                           const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                           bool aReversedSearch = false) {
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
         * Finds all edges with cycles, that is if we have edge a->b there is a path from b to a
         * @param aGraph input graph
         * @return true if there are still cycles
        */
        template<template<typename> class GRAPH_TYPE>
        bool isAcyclic(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph) {
            for (const auto &e : aGraph.getEdges()) {
                if (pathExistsDFS(aGraph, e.dst, e.src)) return false;
            }

            return true;
        }

        /**
         * Finds all edges with cycles, that is if we have edge a->b there is a path from b to a
         * @param aGraph input graph
         * @return container with edges
         */
        template<template<typename> class GRAPH_TYPE>
        auto findEdgesWithCycles(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph) {
            typename Graph<VERTEX_TYPE>::Edges edges;
            for (const auto &e : aGraph.getEdges()) {
                if (pathExistsDFS(aGraph, e.dst, e.src)) {
                    edges.emplace_back(e);
                }
            }
            return edges;
        }

        /**
         * Removes up to requested number of edges with cycles (it may remove less).
         * TODO: 1. Edges with cycles need to be found only first time, then the only cycles left might be a part of left edges only.
         *       So we could speedup procedure to search only within edgesWithCycles
         *       2. Should we use here SCC istead of DFS searches?
         */
        template<template<typename> class GRAPH_TYPE>
        auto getRandomSubgraph(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, int aNumEdgesToRemove) {
            int edgesRemovedCnt = 0;
            while (edgesRemovedCnt < aNumEdgesToRemove) {
                auto edgesWithCycles = findEdgesWithCycles(aGraph);
                auto n = edgesWithCycles.size();
                if (n == 0) {
//                    LOG(TRACE) << "Removed less edges than requested (" << edgesRemovedCnt << ", " << aNumEdgesToRemove << ")";
//                    LOG(TRACE) << edgesWithCycles;
                    return edgesRemovedCnt > 0 ? true : false;
                } else {
                    auto rndEdge = edgesWithCycles[rand() % n];
                    aGraph.removeEdge(rndEdge);
                    ++edgesRemovedCnt;
                }
            }

            return aNumEdgesToRemove > 0 ? true : false;
        }

        template<template <typename> class GRAPH_TYPE>
        auto findNotBlueEdges(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph) {
            using Edge = typename Graph<VERTEX_TYPE>::Edge;
            std::unordered_map<Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;

            for (const auto &e : aGraph.getEdges()) {
                // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
                //    They consist from edges which are cycles not belonging only to aEdge so remove them.
                aGraph.removeEdge(e);

                auto scc = Tools::stronglyConnectedComponents(aGraph);
                for (const auto &s : scc) {
                    if (s.size() == 1) continue;
                    // we get sets of vertices from SCC, find all edges connecting vertices in given SCC and remove
                    for (auto &v : s) {
                        const auto outV = aGraph.getOutVertices(v);
                        for (const auto &vo : outV) {
                            if (s.find(vo) != s.end()) {
                                edgesCnt.try_emplace(e, 0).first->second++;
                            }
                        }
                    }
                }

                aGraph.addEdge(e);
            }

            auto maxCnt = std::max_element(edgesCnt.begin(), edgesCnt.end(), [](const auto &a, const auto &b) -> bool { return a.second < b.second; });

            std::vector<Edge> emax;
            for (auto & [k, v] : edgesCnt) {
                if (v == maxCnt->second) emax.emplace_back(k);
            }

            return emax;
        }

        template<template<typename> class GRAPH_TYPE>
        auto getRandomSubgraphNotBlue(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, int aNumEdgesToRemove, const EdgesSet<VERTEX_TYPE> &blueEdges) {
            int edgesRemovedCnt = 0;
            while (edgesRemovedCnt < aNumEdgesToRemove) {
                auto edgesWithCycles = findEdgesWithCycles(aGraph);
//                auto edgesWithCycles = findNotBlueEdges(aGraph);
//                auto weights = Ext::getEdgeProperties(aGraph, 1);
//                auto [_, edgesWithCycles] = Fasp::GR(aGraph, weights);
                auto n = edgesWithCycles.size();
//                std::cout << "n: " << n << std::endl;
                if (n == 0) {
//                    LOG(TRACE) << "Removed less edges than requested (" << edgesRemovedCnt << ", " << aNumEdgesToRemove << ")";
//                    LOG(TRACE) << edgesWithCycles;
                    return edgesRemovedCnt > 0 ? true : false;
                } else {
                    auto rndEdge = edgesWithCycles[rand() % n];
                    while (blueEdges.find(rndEdge) != blueEdges.end()) {rndEdge = edgesWithCycles[rand() % n];}
                    aGraph.removeEdge(rndEdge);
                    ++edgesRemovedCnt;
                }
            }

            return aNumEdgesToRemove > 0 ? true : false;
        }

        /**
         * Finds path from src to dst
         */
        template<template<typename> class GRAPH_TYPE>
        auto findPathDfs(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                         const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                         const typename Graph<VERTEX_TYPE>::VertexId &aDst) {

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
                        if (iVisited.test(v) == 0) {
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
         * @return container with all verticess accessible from src vertex (if reverse search is set to true then it goes against edge direction)
         */
        template<template<typename> class GRAPH_TYPE>
        auto depthFirstSearch(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
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
        auto findPathWithPositiveCapacity(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                          const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                          const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {

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
        auto minStCutFordFulkersonBase(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                       const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                       const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                       const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

            // Copy graph and weights - they will be modified
            auto g{aGraph};
            auto c{aWeights};
//            std::cout << "MC1:" << aGraph.getStrRepresentationOfGraph() << std::endl;
//            std::cout << "MC1:" << aWeights << std::endl;
//            std::cout << "MC2:" << g.getStrRepresentationOfGraph() << std::endl;
//            std::cout << "MC2:" << c << std::endl;

            // Add edges in backwerd direction if not yet exists
            for (const auto &e : g.getEdges()) {
                auto backwardEdge = typename Graph<VERTEX_TYPE>::Edge{e.dst, e.src};
                if (!g.hasEdge(backwardEdge)) {
                    g.addEdge(backwardEdge);
                    c.insert_or_assign(std::move(backwardEdge), 0);
                }
            }

//            std::cout << "MC3:" << g.getStrRepresentationOfGraph() << std::endl;
//            std::cout << "MC3:" << c << std::endl;

            EDGE_PROP_TYPE maxFlow = 0;
            int cnt = 0;
            while (true) {
                // Find new path from source to sink, if not found than full capacity of network is reached
                auto pathExists = findPathWithPositiveCapacity(g, aSrc, aDst, c);
                if (!pathExists) break;
                // Find edge wiht minimum capacity left
                EDGE_PROP_TYPE pathFlow = std::numeric_limits<EDGE_PROP_TYPE>::max();
                for (size_t i = 0; i < path.size() - 1; ++i) {
//                    std::cout << path[i] << "-" << path[i + 1] << "    ";
                    pathFlow = std::min(pathFlow, c[{path[i], path[i + 1]}]);
                }
//                std::cout << "\n";

                // Update max flow and capacities of edges on the found path
                maxFlow += pathFlow;
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    c.at({path[i], path[i + 1]}) -= pathFlow;
                    c.at({path[i + 1], path[i]}) += pathFlow;
                }

                if (++cnt > 10000) throw MinCutFailed();
            }

//            std::cout << c << std::endl;
//            std::cout << g.getStrRepresentationOfGraph() << std::endl;
            return std::tuple{maxFlow, g, c};
        }

        /**
         * Implementation of maximum flow / min cut algorithm:
         * https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
         * @return maxFlow value
         */
        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        auto minStCutFordFulkerson(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                   const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
            auto maxFlow = std::get<0>(minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights));
            return maxFlow;
        }

        /**
         * G - cleaning graph phase
         * @param aGraph - input graph (it will be modified!)
         * @param aEdge  - edge that will be used as a starting place for cleaning
         */
        template<template<typename> class GRAPH_TYPE>
        void G(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge aEdge) {
            // Remove outgoing edges from destination and ingoing to source.
            auto outv = aGraph.getOutVertices(aEdge.dst);
            for (const auto &vo : outv) {
                aGraph.removeEdge({aEdge.dst, vo});
            }
            auto inv = aGraph.getInVertices(aEdge.src);
            for (const auto &vi : inv) {
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

        template<template<typename> class GRAPH_TYPE>
        void G2(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge aEdge) {
            // Remove outgoing edges from destination and ingoing to source.
            auto outv = aGraph.getOutVertices(aEdge.dst);
            for (const auto &vo : outv) {
                aGraph.removeEdge({aEdge.dst, vo});
            }
            auto inv = aGraph.getInVertices(aEdge.src);
            for (const auto &vi : inv) {
                aGraph.removeEdge({vi, aEdge.src});
            }
        }

        /**
         * G* - cleaning graph phase
         * @param aGraph - input graph (it will be modified!)
         * @param aEdge  - edge that will be used as a starting place for cleaning
         */
        template<template<typename> class GRAPH_TYPE>
        void GStar2(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge) {
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
         * G* - cleaning graph phase
         * @param aGraph - input graph (it will be modified!)
         * @param aEdge  - edge that will be used as a starting place for cleaning
         */
        template<template<typename> class GRAPH_TYPE>
        bool GStar(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge) {

            // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
            //    They consist from edges which are cycles not belonging only to aEdge so remove them.
            aGraph.removeEdge(aEdge);
            auto scc = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) continue;
                // we get sets of vertices from SCC, find all edges connecting vertices in given SCC and remove
                for (auto &v : s) {
                    const auto outV = aGraph.getOutVertices(v);
                    for (const auto &vo : outV) {
                        if (s.find(vo) != s.end()) {
                            aGraph.removeEdge({v, vo});
                        }
                    }
                }
            }

            // 2. There can be only one (or none) connected component with size > 1, if it is found
            // then it consist aEdge and all its isolated cycles
            aGraph.addEdge(aEdge);
            auto scc2 = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc2) {
                if (s.size() > 1) return true;
            }
            return false;
        }

        //, typename EDGE_PROP_TYPE>
        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        bool GStarBlue(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge, EdgesSet<VERTEX_TYPE> &aBlueEdges, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph2, typename Graph<VERTEX_TYPE>::Edge &redEdge, int &maxMcRedEdge) {


            auto pathExists = findPathDfs(aGraph, aEdge.dst, aEdge.src);
            auto somePath = path;

            // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
            //    They consist from edges which are cycles not belonging only to aEdge so remove them.
            aGraph.removeEdge(aEdge);
            auto scc = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) continue;
                // we get sets of vertices from SCC, find all edges connecting vertices in given SCC and remove
                for (auto &v : s) {
                    const auto outV = aGraph.getOutVertices(v);
                    for (const auto &vo : outV) {
                        if (s.find(vo) != s.end()) {
                            aGraph.removeEdge({v, vo});
//                            std::cout << typename Graph<VERTEX_TYPE>::Edge{v, vo} << std::endl;
                        }
                    }
                }
//                std::cout << std::endl;
            }

            // 1.5 update blue edges
            for (auto &e : aGraph.getEdges()) {
                aBlueEdges.emplace(e);
            }

            // 2. There can be only one (or none) connected component with size > 1, if it is found
            // then it consist aEdge and all its isolated cycles
            aGraph.addEdge(aEdge);
            auto scc2 = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc2) {
                if (s.size() > 1) return true;
            }

            // 3. Find 'red edges' when there is no isolated cycles
//            {
//                int cnt = 0;
//                for (size_t i = 0; i < somePath.size() - 1; ++i) {
//                    typename Graph<VERTEX_TYPE>::Edge e = {somePath[i], somePath[i + 1]};
//
//                    if (!aGraph.hasEdge(e)) {
//                        cnt++;
//                        auto d = minStCutFordFulkerson(aGraph2, e.dst, e.src, aWeights);
//                        if (d > maxMcRedEdge) {
//                            maxMcRedEdge = d;
//                            redEdge = e;
//                        }
//                    }
//                }
////                std::cout << aGraph.getStrRepresentationOfGraph() << std::endl;
////                std::cout << cnt << "/" << (somePath.size() - 1) << " " << aEdge << " Edges: " << p << std::endl;
////                std::cout << "-------\n";
//            }


            return false;
        }

        /**
         * step2b - cleaning graph phase (Yes! we need better name for it).
         * Finds all edges which are part of other cycle (not only cycle going through aEdge).
         */
        template<template<typename> class GRAPH_TYPE>
        auto step2b(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                    const Graph <VERTEX_TYPE, GRAPH_TYPE> &aCleanedGraph,
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

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto superAlgorithm(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path, bool aUseWeights = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdges;
        path.saCnt++;
        auto outGraph{aGraph}; // eventually acyclic graph

        while (true) {
            bool wasGraphModified = false;

            for (const auto &e : outGraph.getEdges()) {
                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                if (!path.GStar(workGraph, e)) continue;

                // If we have weights assigned to edges then we need to do min-cut, if not it is always safe
                // to remove current edge
                bool shouldRemoveCurrentEdge = true;
                if (aUseWeights) {
                    auto mc = path.minStCutFordFulkerson(workGraph, e.dst, e.src, aWeights);
                    if (mc < aWeights.at(e)) {
                        // Do not remove that edge
                        shouldRemoveCurrentEdge = false;
                    }
                }

                if (shouldRemoveCurrentEdge) {
                    wasGraphModified = true;
                    outGraph.removeEdge(e);
                    removedEdges.emplace_back(e);
                }
            }
            if (!wasGraphModified) break;
        }

        return removedEdges;
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto superAlgorithm(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, bool aUseWeights = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        auto vertices = aGraph.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

        return superAlgorithm(aGraph, aWeights, path);
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithmBlue(Graph<VERTEX_TYPE, GRAPH_TYPE> &outGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path, bool aUseWeights = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdges;
        typename Graph<VERTEX_TYPE>::Edge edgeDummy;
        int mcDummy;
        path.saCnt++;
//        auto outGraph{aGraph}; // eventually acyclic graph

        EdgesSet<VERTEX_TYPE> setOfEdges;

        while(true) {
            bool wasGraphModified = false;

            setOfEdges.clear();

            for (const auto &e : outGraph.getEdges()) {

                if (setOfEdges.find(e) != setOfEdges.end()) continue; // e is in 'blue edges' set

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                if (!path.GStarBlue(workGraph, e, setOfEdges, aWeights, outGraph, edgeDummy, mcDummy)) continue;

                // If we have weights assigned to edges then we need to do min-cut, if not it is always safe
                // to remove current edge
                bool shouldRemoveCurrentEdge = true;
                if (aUseWeights) {
                    auto mc = path.minStCutFordFulkerson(workGraph, e.dst, e.src, aWeights);
                    if (mc < aWeights.at(e)) {
                        // Do not remove that edge
                        shouldRemoveCurrentEdge = false;
                    }
                }

                if (shouldRemoveCurrentEdge) {
                    wasGraphModified = true;
                    outGraph.removeEdge(e);
                    removedEdges.emplace_back(e);
                }
            }
            if (!wasGraphModified) break;
        }
//        std::cout <<" BLUE: " << setOfEdges.size() << " G: " << aGraph.getNumOfEdges() << std::endl;
        return std::pair{removedEdges, setOfEdges};
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithmBlue2(Graph<VERTEX_TYPE, GRAPH_TYPE> &outGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path, bool aUseWeights = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdgesSA;
        typename Graph<VERTEX_TYPE>::Edges removedEdgesGR;
        path.saCnt++;
//        auto outGraph{aGraph}; // eventually acyclic graph

        typename Graph<VERTEX_TYPE>::Edge redEdge;
        int maxMcRedEdge = 0;
        EdgesSet<VERTEX_TYPE> setOfEdges;

        int cnt = 0;
        while(true) {
            cnt++;
            bool wasGraphModified = false;

            setOfEdges.clear();

//            auto ed = Fasp::GR(outGraph, aWeights).second;

            for (const auto &e : outGraph.getEdges()) {

                if (setOfEdges.find(e) != setOfEdges.end()) continue; // e is in 'blue edges' set

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                if (!path.GStarBlue(workGraph, e, setOfEdges, aWeights, outGraph, redEdge, maxMcRedEdge)) continue;
//                {
//                    const auto &itBegin = ed.cbegin();
//                    const auto &itEnd = ed.cend();
//                    if (std::find(itBegin, itEnd, e) == itEnd) continue; // e not in GR set
//                }

                // If we have weights assigned to edges then we need to do min-cut, if not it is always safe
                // to remove current edge
                bool shouldRemoveCurrentEdge = true;
                if (aUseWeights) {
                    auto mc = path.minStCutFordFulkerson(workGraph, e.dst, e.src, aWeights);
                    if (mc < aWeights.at(e)) {
                        // Do not remove that edge
                        shouldRemoveCurrentEdge = false;
                    }
                }

                if (shouldRemoveCurrentEdge) {
                    wasGraphModified = true;
                    outGraph.removeEdge(e);
                    removedEdgesSA.emplace_back(e);
//                    ed = Fasp::GR(outGraph, aWeights).second;
                }
            }
            if (!wasGraphModified) break;
        }
        if (removedEdgesSA.size() == 0) {
                auto ed = Fasp::GR(outGraph, aWeights).second;
                const auto &itBegin = ed.cbegin();
                const auto &itEnd = ed.cend();
                auto n = ed.size();
                if (maxMcRedEdge == -1) { //turn off
//                    std::cout << "MC EDGE: " << redEdge << " " << maxMcRedEdge << std::endl;
                    removedEdgesGR.push_back(redEdge);
                }
                else if (n > 0) removedEdgesGR.push_back(ed[0]);
        }
//        std::cout <<" BLUE: " << setOfEdges.size() << " G: " << aGraph.getNumOfEdges() << std::endl;
        return std::tuple{removedEdgesSA, setOfEdges, removedEdgesGR};
    }

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto deltaFASP(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
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
                path.GStar2(gs, e);

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

        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(DELTA) capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }

    /**
     * Modification to original idea. Instead of regenerating graphs from beginning, it will keep of removeing 1 edge
     * till any solution is found in any of random subgraphs.
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto randomFASP_sequential(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
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

        srand(time(NULL));

        while (true) {

            if (path.isAcyclic(g)) break;

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, double, Ext::EdgeHasher < VERTEX_TYPE>>
            edgesCnt;

            std::future<std::pair<typename Graph<VERTEX_TYPE>::Edges, int>> tasks[numOfReps];
            int edgesFound = 0;
            int workersDone = 0;
            int workersDone2 = 0;
            int randSuccessCnt = 0;
            std::mutex mutex;
            std::mutex mutex2;
            std::condition_variable cv;

            for (int i = 0; i < numOfReps; ++i) {
                tasks[i] = std::async(std::launch::async,
                                      [&]() {
                                          auto randGraph{g};

                                          auto vertices = randGraph.getVertices();
                                          auto maxId = std::max_element(vertices.begin(), vertices.end());
                                          PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

                                          int level = 1;
                                          while (true) {
                                              auto randSuccess = path.getRandomSubgraph(randGraph, 1 /*numEdgesToRemove*/);

                                              auto saEdges = superAlgorithm(randGraph, aWeights, path);

                                              { // ------------ Wait for results & update data for all worker threads
                                                  std::unique_lock<std::mutex> lock(mutex);
                                                  if (saEdges.size() > 0) ++edgesFound;
                                                  if (randSuccess) ++randSuccessCnt;
                                                  ++workersDone;
                                                  if (workersDone == numOfReps) {
                                                      workersDone2 = 0;
                                                      cv.notify_all();
                                                  }
                                                  else { cv.wait(lock, [&] { return workersDone == numOfReps; }); }
                                              }

                                              // -------------- If found any solution then quit (same will happen with all threads running)
                                              if (edgesFound > 0 || randSuccessCnt == 0) return std::pair{saEdges, level};

                                              { // ------------ Reinitialize counters to make first condition_variable ready for next loop
                                                  std::unique_lock<std::mutex> lock(mutex2);
                                                  ++workersDone2;
                                                  if (workersDone2 == numOfReps) {
                                                      workersDone = 0;
                                                      randSuccessCnt = 0;
                                                      cv.notify_all();
                                                  }
                                                  else { cv.wait(lock, [&] { return workersDone2 == numOfReps; }); }
                                              }

                                              level++;
                                          };
                                      });
            }
            for (int i = 0; i < numOfReps; ++i) {
                auto[edgesToRemove, _] = tasks[i].get();
                for (auto &e : edgesToRemove) {
                    edgesCnt.try_emplace(e, 0).first->second++;
                }
            }

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
            } else {
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

    /**
     * Working original idea of how random FASP heuristic should work
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto randomFASP_orig(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
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

        srand(time(NULL));

        int numEdgesToRemove = 1;
        int numCnt = 0;
        int emptyCnt = 0;
        while (true) {
            numCnt++;

            if (path.isAcyclic(g)) break;

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher < VERTEX_TYPE>>
            edgesCnt;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
            std::future<typename Graph<VERTEX_TYPE>::Edges> tasks[numOfReps];
            for (int i = 0; i < numOfReps; ++i) {
                tasks[i] = std::async(std::launch::async,
                                      [&]() {
                                          auto workGraph{g};
                                          auto vertices = workGraph.getVertices();
                                          auto maxId = std::max_element(vertices.begin(), vertices.end());
                                          PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

                                          path.getRandomSubgraph(workGraph, numEdgesToRemove);
                                          return superAlgorithm(workGraph, aWeights, path);
                                      });
            }
            for (int i = 0; i < numOfReps; ++i) {
                auto edgesToRemove = tasks[i].get();
                for (auto &e : edgesToRemove) {
                    edgesCnt.try_emplace(e, 0).first->second++;
                }
            }
            if (edgesCnt.size() == 0) {
                emptyCnt++;
                numEdgesToRemove++;
                continue; // reapeat loop - we have not found any solutions
            }
            numEdgesToRemove = 1;

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

        EDGE_PROP_TYPE capacity = 0;
        for (const auto &e : removedEdges) {
            capacity += aWeights.at(e);
        }

        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }


    /**
     * new idea - remove random edges till SA succeedes and continue till graph acyclic
     * when done take original graph and remove only edges removed by SA if still cyclic then
     * repeate whole procedure and continu till all removed edges are obtain thanks to SA only
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
    static auto randomFASP_new(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
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

        srand(time(NULL));

        while (true) {
            if (path.isAcyclic(g)) break;

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher < VERTEX_TYPE>>
            edgesCnt;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
            std::future<typename Graph<VERTEX_TYPE>::Edges> tasks[numOfReps];
            for (int i = 0; i < numOfReps; ++i) {
                tasks[i] = std::async(std::launch::async,
                                      [&]() {
                                          typename Graph<VERTEX_TYPE>::Edges setOfEdges;
                                          int numEdges2Remove = 1;
                                          bool randEdgeRemoved = false;
                                          auto inputGraph{g};
                                          auto workGraph{g};

                                          auto vertices = workGraph.getVertices();
                                          auto maxId = std::max_element(vertices.begin(), vertices.end());
                                          PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

                                          while (true) {
                                              // Still need to help of random edges
                                              auto randSuccess = path.getRandomSubgraph(workGraph, numEdges2Remove);
                                              if (randSuccess) randEdgeRemoved = true;

                                              if (randSuccess == false) {
                                                  // we could not remove any edge - graph is acyclic

                                                  if (randEdgeRemoved) {
                                                      // we have some randomg edges removed already
                                                      workGraph = inputGraph;
                                                      workGraph.removeEdges(setOfEdges); // remove so far discovered SA edges

                                                      randEdgeRemoved = false;
                                                  } else {
                                                      // If all edges are remove only via SA and graph is acyclic then return solution
                                                      // we have acyclic graph using SA edges only!
                                                      return setOfEdges;
                                                  }
                                              }

                                              auto sa = superAlgorithm(workGraph, aWeights, path);
                                              setOfEdges.insert(std::end(setOfEdges), std::begin(sa), std::end(sa));
                                              workGraph.removeEdges(sa);
                                          }
                                      });
            }

            {
                std::vector<typename Graph<VERTEX_TYPE>::Edges> collect;
                for (int i = 0; i < numOfReps; ++i) {
                    auto sa = tasks[i].get();

//                    std::cout << "SOL:: " << sa.size() << " " << sa << std::endl;
                    collect.emplace_back(std::move(sa));
                }

                bool first = true;
                int min = 0;
                size_t minIdx = 0;
                for (int i = 0; i < numOfReps; ++i) {
                    if (first || collect[i].size() < min) {
                        first = false;
                        min = collect[i].size();
                        minIdx = i;
                    }
                }

                g.removeEdges(collect[minIdx]);
                removedEdges.insert(std::end(removedEdges), std::begin(collect[minIdx]), std::end(collect[minIdx]));
                std::cout << "AFTER1: " << collect[minIdx].size() << " " << collect[minIdx] << std::endl;
                std::cout << "AFTER2: " << path.findEdgesWithCycles(g) << std::endl;
            }

            // Try to solve the rest with superAlgoritm - maybe it will now succeed!
            auto edgesToRemove = superAlgorithm(g, aWeights, path);
            g.removeEdges(edgesToRemove);
            removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
            removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());
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

    /**
     * Working original idea of how random FASP heuristic should work
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto randomFASP_blueEdges(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        constexpr int numOfReps = 20;
        int useGrThreshold = 8;
        auto g{aGraph};


        // ========= calc threshold ==========
        auto [_, grEdges] = Fasp::GR(g, aWeights);
        int faspSize = grEdges.size();
        double density = g.getNumOfEdges() / g.getNumOfVertices();
        if (faspSize > 0) useGrThreshold = 2 * g.getNumOfEdges() / faspSize;
        useGrThreshold = std::min(useGrThreshold, 2 * 5);
        std::cout << ":::::::::THRESHOLD = " << useGrThreshold << std::endl;
        // ===================================


        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        // find max vertex id in graph and setup PathHero
        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path{static_cast<size_t>(maxId == vertices.end() ? 1 : *maxId + 1)};
//        {
//            Timer<true, false> t(true, "SCC1");
//            t.start_timer("");
//            auto scc = Tools::stronglyConnectedComponents(g);
//            std::cout << scc.size() << "\n";
//            auto numOfScc = std::count_if(scc.begin(), scc.end(), [](auto &in) {
//                if (in.size() > 1) std::cout << in.size() << " ";
//                return in.size() > 1;
//            });
//            std::cout << "NUM OF SCC: " << numOfScc << std::endl;
//            t.stop_timer();
//
//        }
        // initial run of superAlgorithm (SA)
        auto [edgesToRemove, blueEdgesx] = superAlgorithmBlue(g, aWeights, path);
        auto blueEdges = blueEdgesx;
        g.removeEdges(edgesToRemove);
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

        srand(time(NULL));

        int numEdgesToRemoveInitVal = 1; //faspSize == 0 ? 1 : g.getNumOfEdges() / faspSize / 2;
        std::cout << "EDGES TO REMOVE INIT VAL = " << numEdgesToRemoveInitVal << std::endl;
        int numEdgesToRemove = numEdgesToRemoveInitVal;

        int numCnt = 0;
        int emptyCnt = 0;
        while (true) {
            numCnt++;


            if (path.isAcyclic(g)) break;

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCntGR;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
            std::future<std::pair<typename Graph<VERTEX_TYPE>::Edges, typename Graph<VERTEX_TYPE>::Edges>> tasks[numOfReps];
            for (int i = 0; i < numOfReps; ++i) {
                tasks[i] = std::async(std::launch::async,
                                      [&] () {
                                          auto workGraph{g};
                                          auto vertices = workGraph.getVertices();
                                          auto maxId = std::max_element(vertices.begin(), vertices.end());
                                          PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included

                                          path.getRandomSubgraphNotBlue(workGraph, numEdgesToRemove, blueEdges);
                                          auto [edgesSA, _, edgesGR] = superAlgorithmBlue2(workGraph, aWeights, path);
                                          return std::pair{edgesSA, edgesGR};
                                      });
            }
            for (int i = 0; i < numOfReps; ++i) {
                auto [edgesToRemove, edgesToRemoveGR] = tasks[i].get();
                for (auto &e : edgesToRemove) {
                    edgesCnt.try_emplace(e, 0).first->second++;
                }
                for (auto &e : edgesToRemoveGR) {
                    edgesCntGR.try_emplace(e, 0).first->second++;
                }
            }
            if (edgesCnt.size() == 0) {
                std::cout << "------ USING GR ------ " << edgesCntGR.size() <<"\n";
                edgesCnt.swap(edgesCntGR); // In case when SA didn't find any edges use these from GR heuristic
            }
            std::cout << "....\n";
//            std::cout << "ITER: " << edgesCnt.size() << " #2R: " << numEdgesToRemove << std::endl;
            if (edgesCnt.size() == 0) {
                emptyCnt++;
                numEdgesToRemove++;

                if (numEdgesToRemove >= useGrThreshold/ 2) {
                    useGrThreshold++;
                    std::cout << ">>>>>>>> " << useGrThreshold << std::endl;
                    auto [c, faspEdges] = Fasp::GR(g, aWeights);
                    int n = faspEdges.size();
                    auto e = faspEdges[0]; // rand() % n
                    g.removeEdge(e);
                    removedEdges.emplace_back(std::move(e));
                    numEdgesToRemove = numEdgesToRemoveInitVal;
                }
                else
                continue; // reapeat loop - we have not found any solutions
            }
            else {
                if (useGrThreshold > 6) useGrThreshold--;
            }
            numEdgesToRemove = numEdgesToRemoveInitVal;


            if (edgesCnt.size() > 0) {
                auto maxCnt = std::max_element(edgesCnt.begin(), edgesCnt.end(), [](const auto &a, const auto &b) -> bool { return a.second < b.second; });
//                std::cout << "CNT: " << *maxCnt << std::endl;
                // Remove best candidate
                g.removeEdge(maxCnt->first);
                removedEdges.emplace_back(std::move(maxCnt->first));
            }
            // Try to solve the rest with superAlgoritm - maybe it will now succeed!
            auto [edgesToRemove, blueEdges2] = superAlgorithmBlue(g, aWeights, path);
            blueEdges = blueEdges2;
            g.removeEdges(edgesToRemove);
            removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
            removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

//            Timer<true, false> t(true, "SCC");
//            t.start_timer("");
//            auto scc = Tools::stronglyConnectedComponents(g);
//            auto numOfScc = std::count_if(scc.begin(), scc.end(), [](auto &in) {
//                if (in.size() > 1) std::cout << in.size() << " ";
//                return in.size() > 1;
//            });
//            std::cout << "NUM OF SCC: " << numOfScc << std::endl;
//            t.stop_timer();
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
