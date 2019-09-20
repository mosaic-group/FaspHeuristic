#ifndef FASPHEURISTIC_GRAPHFASPFASTFINAL_H
#define FASPHEURISTIC_GRAPHFASPFASTFINAL_H

#include "graph/graph.h"
#include "graph/graphExt.h"
#include "graph/graphFasp.h"
#include "tools/dynamicBitset.h"
#include "tools/stack.h"
#include "tools/tools.h"
#include <future>
#include <cstddef>

namespace Graph::FaspFastFinal {

    /**
     * Exception class thrown in case never-ending mincut loop
     */
    class MinCutFailed : public std::exception {};

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
        DynamicBitset <uint32_t, uint16_t> iVisited;
        Stack <VERTEX_TYPE> stack;
        std::vector<VERTEX_TYPE> parents;
        typename Graph<VERTEX_TYPE>::Edges inEdges;
        typename Graph<VERTEX_TYPE>::Edges outEdges;
        typename Graph<VERTEX_TYPE>::Vertices path;


    public:
        explicit PathHero(std::size_t aN) : iVisited(aN), stack(aN), parents(aN, 0) {
            // both can have at most aN in/outgoing edges
            inEdges.reserve(aN);
            outEdges.reserve(aN);
            path.reserve(aN);
        }

        /**
         * @return true if path from src to dst exists
         */
        template<bool FORWARD_SEARCH=true, template<typename> class GRAPH_TYPE>
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
                const auto &vertices = FORWARD_SEARCH ? aGraph.getOutVertices(currentVertex) : aGraph.getInVertices(currentVertex);
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

        template<template<typename> class GRAPH_TYPE, typename EDGE_PROP_TYPE>
        auto getRedEdge(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            using Edge = typename Graph<VERTEX_TYPE>::Edge;
            Edge redEdge{};
            int maxMcRedEdge = 0;

            for (const auto &e : aGraph.getEdges()) {
                if (!pathExistsDFS(aGraph, e.dst, e.src)) continue; // optimization

                auto pathExists = findPathDfs(aGraph, e.dst, e.src);
                if (!pathExists) continue;
                auto somePath = path;

                aGraph.removeEdge(e);
                auto scc = Tools::stronglyConnectedComponents(aGraph);

                std::vector<Edge> redEdges;
                bool prevCandidate = false;
                for (std::size_t i = 0; i < somePath.size(); ++i) {
                    typename Graph<VERTEX_TYPE>::VertexId v = somePath[i];

                    bool currCandidate = false;
                    for (auto &s : scc) if (s.size() > 1 && s.find(v) != s.end()) {currCandidate = true; break;}
                    if (currCandidate && prevCandidate) redEdges.emplace_back(Edge{somePath[i - 1], somePath[i]});
                    prevCandidate = currCandidate;
                }

                aGraph.addEdge(e);
                for (const auto &e : redEdges) {
                    auto d = minStCutFordFulkerson(aGraph, e.dst, e.src, aWeights);
                    if (d > maxMcRedEdge) {
                        maxMcRedEdge = d;
                        redEdge = e;
                    }
                }

                if (maxMcRedEdge > 1) break;
            }

            return std::pair{maxMcRedEdge, redEdge};
        }

        template<template<typename> class GRAPH_TYPE>
        bool GStarBlue(Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge, EdgesSet<VERTEX_TYPE> &aBlueEdges) {
            // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
            //    They consist from edges which are cycles not belonging only to aEdge so remove them.
            aGraph.removeEdge(aEdge);
            const auto scc = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) continue;
                // we get sets of vertices from SCC, find all edges connecting vertices in given SCC and remove
                for (const auto &v : s) {
                    const auto outV = aGraph.getOutVertices(v);
                    for (const auto &vo : outV) {
                        if (s.find(vo) != s.end()) {
                            aGraph.removeEdge({v, vo});
                        }
                    }
                }
            }

            // 1.5 update blue edges
            for (auto &e : aGraph.getEdges()) {
                aBlueEdges.emplace(e);
            }

            // 2. There can be only one (or none) connected component with size > 1, if it is found
            // then it consist aEdge and all its isolated cycles
            if (pathExistsDFS(aGraph, aEdge.dst, aEdge.src)) return true;
//            aGraph.addEdge(aEdge);


//            auto scc2 = Tools::stronglyConnectedComponents(aGraph);
//            for (const auto &s : scc2) {
//                if (s.size() > 1) return true;
//            }

            return false;
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
                    }

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

            // Add edges in backwerd direction if not yet exists
            for (const auto &e : g.getEdges()) {
                auto backwardEdge = typename Graph<VERTEX_TYPE>::Edge{e.dst, e.src};
                if (!g.hasEdge(backwardEdge)) {
                    g.addEdge(backwardEdge);
                    c.insert_or_assign(std::move(backwardEdge), 0);
                }
            }

            EDGE_PROP_TYPE maxFlow = 0;
            int cnt = 0;
            while (true) {
                // Find new path from source to sink, if not found than full capacity of network is reached
                auto pathExists = findPathWithPositiveCapacity(g, aSrc, aDst, c);
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
        auto minStCutFordFulkerson(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                   const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
            auto maxFlow = std::get<0>(minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights));
            return maxFlow;
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
                    return (edgesRemovedCnt > 0);
                } else {
                    auto rndEdge = edgesWithCycles[::Tools::randInt(0, n - 1)];
                    while (blueEdges.find(rndEdge) != blueEdges.end()) {rndEdge = edgesWithCycles[::Tools::randInt(0, n - 1)];}
                    aGraph.removeEdge(rndEdge);
                    ++edgesRemovedCnt;
                }
            }

            return (aNumEdgesToRemove > 0);
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
                    }

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
    };

    // ------------------------------------------------------------------------

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithmBlue(Graph<VERTEX_TYPE, GRAPH_TYPE> &outGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path, bool aUseWeights = false, bool relaxSA = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdgesSA;
        typename Graph<VERTEX_TYPE>::Edges removedEdgesGR;

        EdgesSet<VERTEX_TYPE> setOfEdges{};

        int cnt = 0;
        Timer<false, false> t(true);
        while(true) {
            cnt++;
            bool wasGraphModified = false;

            setOfEdges.clear();

            for (const auto &e : outGraph.getEdges()) {
                if (setOfEdges.find(e) != setOfEdges.end()) continue; // e is in 'blue edges' set

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                if (!path.GStarBlue(workGraph, e, setOfEdges)) continue;

                // If we have weights assigned to edges then we need to do min-cut, if not it is always safe to remove current edge
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
                }
            }
            if (!wasGraphModified) break;
        }
        if (removedEdgesSA.size() == 0) {
            t.start_timer("RED");
            auto [maxMcRedEdge, redEdge] = path.getRedEdge(outGraph, aWeights);
            std::cout << t.stop_timer() << " ";
            if (relaxSA && maxMcRedEdge > 0) {
                std::cout << "xRE";
                removedEdgesGR.push_back(redEdge);
            }
            else {
                auto ed = Fasp::GR(outGraph, aWeights).second;
                auto n = ed.size();

                if (n > 0) {
                    std::cout << "xGR";
                    removedEdgesGR.push_back(ed[0]);
                }
            }

        }
        return std::tuple{removedEdgesSA, setOfEdges, removedEdgesGR};
    }

    /**
     * Working original idea of how random FASP heuristic should work
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto randomFASP(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {

        auto cleanGraphWithScc = [](Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph) {
            int cnt1 = 0, cntBig = 0;
            auto scc = Tools::stronglyConnectedComponents(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) {cnt1++; aGraph.removeVertex(*s.begin());}
                else cntBig++;
            }
            std::cout << "SCC  #1=" << cnt1 << " #BIG=" << cntBig << "\n";
        };

        Timer<true, false, true> t(true);
        constexpr int numOfReps = 20;

        Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> props[numOfReps];
        for (auto &p : props) p = aWeights;

        auto g{aGraph};
        t.start_timer("init scc");
        cleanGraphWithScc(g);
        t.stop_timer();

        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        // find max vertex id in graph and setup PathHero
        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path{static_cast<std::size_t>(maxId == vertices.end() ? 1 : *maxId + 1)};

        // initial run of superAlgorithm (SA)
        auto [edgesToRemove, blueEdgesx, dummy2] = superAlgorithmBlue(g, aWeights, path);
        auto blueEdges = blueEdgesx;
        g.removeEdges(edgesToRemove);
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

        int numEdgesToRemoveInitVal = 1;
        //auto [dummy1, grEdges] = Fasp::GR(g, aWeights); int faspSize = grEdges.size(); numEdgesToRemoveInitVal = faspSize == 0 ? 1 : g.getNumOfEdges() / faspSize / 2;
        std::cout << "EDGES TO REMOVE INIT VAL = " << numEdgesToRemoveInitVal << std::endl;
        int numEdgesToRemove = numEdgesToRemoveInitVal;
        int cnt = 1;
        t.start_timer("rand loop");
        while (true) {
            std::cout << "----- loop=" << cnt++ << "\n";
            if (path.isAcyclic(g)) break;

            cleanGraphWithScc(g);

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCntGR;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
            t.start_timer("random graphs");
            alignas(1024) std::future<std::pair<typename Graph<VERTEX_TYPE>::Edges, typename Graph<VERTEX_TYPE>::Edges>> tasks[numOfReps];
            int i = 0;
            for (auto &task : tasks) {

                task= std::async(std::launch::async,
                     [&, i] () {
                         Timer<true, false> tt(true);
                         if (i == 0) tt.start_timer("1 - copy graph");
                         auto workGraph{g};
                         if (i == 0) tt.stop_timer();

                         if (i == 0) tt.start_timer("2 - prepare path hero");
                         auto vertices = workGraph.getVertices();
                         auto maxId = std::max_element(vertices.begin(), vertices.end());
                         PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included
                         if (i == 0) tt.stop_timer();

                         if (i == 0) tt.start_timer("3 - random subgraph");
                         path.getRandomSubgraphNotBlue(workGraph, numEdgesToRemove, blueEdges);
                         if (i == 0) tt.stop_timer();

                         if (i == 0) tt.start_timer("4 - SA blue");
                         auto [edgesSA, _, edgesGR] = superAlgorithmBlue(workGraph, props[i], path, false, true);
                         if (i == 0) tt.stop_timer();
                         return std::pair{edgesSA, edgesGR};
                     });
                i++;
            }

            for (auto &task : tasks) {
                auto [edgesToRemove, edgesToRemoveGR] = task.get();
                for (auto &e : edgesToRemove) edgesCnt.try_emplace(e, 0).first->second++;
                for (auto &e : edgesToRemoveGR) edgesCntGR.try_emplace(e, 0).first->second++;
            }
            t.stop_timer();

            if (edgesCnt.size() == 0) {
                std::cout << "------ USING alternative (GR/RED_EDGE) ------ " << edgesCntGR.size() <<"\n";
                edgesCnt.swap(edgesCntGR); // In case when SA didn't find any edges use these from GR heuristic
                std::cout << "ER2: " << edgesCnt << "\n";
            }
            if (edgesCnt.size() == 0) {
                numEdgesToRemove++;
                continue; // reapeat loop - we have not found any solutions
            }
            numEdgesToRemove = numEdgesToRemoveInitVal;

            if (edgesCnt.size() > 0) {
                auto maxCnt = std::max_element(edgesCnt.begin(), edgesCnt.end(), [](const auto &a, const auto &b) -> bool { return a.second < b.second; });
                // Remove best candidate
                g.removeEdge(maxCnt->first);
                removedEdges.emplace_back(std::move(maxCnt->first));
            }

            t.start_timer("main SA");
            // Try to solve the rest with superAlgoritm - maybe it will now succeed!
            auto [edgesToRemove, blueEdges2, dummy3] = superAlgorithmBlue(g, aWeights, path);
            blueEdges = blueEdges2;
            g.removeEdges(edgesToRemove);
            removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
            removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());
            t.stop_timer();
        }
        t.stop_timer();

        EDGE_PROP_TYPE capacity = 0;
        for (const auto &e : removedEdges) { capacity += aWeights.at(e); }

        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }
}


#endif //FASPHEURISTIC_GRAPHFASPFASTFINAL_H
