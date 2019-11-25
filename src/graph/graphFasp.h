#ifndef FASPHEURISTIC_GRAPHFASP_H
#define FASPHEURISTIC_GRAPHFASP_H

#include "graph/graph.h"
#include "graph/graphExt.h"
#include "graph/graphFaspTools.h"
#include "tools/dynamicBitset.h"
#include "tools/stack.h"
#include "tools/tools.h"
#include <future>
#include <cstddef>
#include <vector>
#include "graph/hipr/hi_pr.h"

namespace Graph::FaspFastFinal {

    template <typename VERTEX_TYPE>
    using EdgesSet = std::unordered_set<typename Graph<VERTEX_TYPE>::Edge, Ext::EdgeHasher<VERTEX_TYPE>>;

    /**
     * PathHero... class is a approach to speed up things. For all algorithms it has implemented it uses
     * pre-allocated memory - that is why it need to be initialized with maximum number of vertices (or if not
     * in sequence 0, 1, 2 ... then with a maximum value of vertex id + 1.
     * @tparam VERTEX_TYPE
     */
    template<typename VERTEX_TYPE>
    class alignas (64) PathHero {

        // Allocation of structures/memory/containers shared by all algorithms from PathHero class
        DynamicBitset <uint32_t, uint16_t> iVisited;
        Stack <VERTEX_TYPE> stack;
        std::vector<VERTEX_TYPE> parents;
        typename Graph<VERTEX_TYPE>::Vertices path;
        std::vector<int16_t> lowLinks;
        std::vector<int16_t> index;

    public:
        explicit PathHero(std::size_t aN) : iVisited(aN), stack(aN), parents(aN, 0), lowLinks(aN, -1), index(aN, -1) {
            // both can have at most aN in/outgoing edges
            path.reserve(aN);
        }

        /**
         * @return true if path from src to dst exists
         */
        template<bool FORWARD_SEARCH=true>
        bool pathExistsDFS(const Graph <VERTEX_TYPE> &aGraph,
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
        bool isAcyclic(const Graph <VERTEX_TYPE> &aGraph) {
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
        auto findEdgesWithCycles(const Graph <VERTEX_TYPE> &aGraph) {
            typename Graph<VERTEX_TYPE>::Edges edges;
            for (const auto &e : aGraph.getEdges()) {
                if (pathExistsDFS(aGraph, e.dst, e.src)) {
                    edges.emplace_back(e);
                }
            }
            return edges;
        }

        template<typename EDGE_PROP_TYPE>
        auto getRedEdge(Graph <VERTEX_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            using Edge = typename Graph<VERTEX_TYPE>::Edge;
            Edge redEdge{};
            int maxMcRedEdge = 0;

            for (const auto &e : aGraph.getEdges()) {
                EDGE_PROP_TYPE eCapacity = aWeights.at(e);
                if (!pathExistsDFS(aGraph, e.dst, e.src)) continue; // optimization

                auto pathExists = findPathDfs(aGraph, e.dst, e.src);
                if (!pathExists) continue;
                auto somePath = path;

                aGraph.removeEdge(e);
                auto scc = stronglyConnectedComponents2(aGraph);

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
                for (auto &ee : redEdges) {
                    auto d = minStCutFordFulkerson(aGraph, ee.dst, ee.src, aWeights);
                    if (d > maxMcRedEdge) {
                        maxMcRedEdge = d;
                        redEdge = ee;
                    }
                }

                if (maxMcRedEdge > eCapacity) break;
            }

            return std::pair{maxMcRedEdge, redEdge};
        }

        bool GStarBlue(Graph <VERTEX_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge, EdgesSet<VERTEX_TYPE> &aBlueEdges, bool weighted) {
            // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
            //    They consist from edges which are cycles not belonging only to aEdge so remove them.
            aGraph.removeEdge(aEdge);
            const auto scc = stronglyConnectedComponents2(aGraph);//Tools::stronglyConnectedComponents2(aGraph);
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
            if (!weighted)
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
         * Implementation of maximum flow / min cut algorithm:
         * https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
         * @return maxFlow value
         */
        template<typename EDGE_PROP_TYPE>
        auto minStCutFordFulkerson(const Graph <VERTEX_TYPE> &aGraph,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                                   const typename Graph<VERTEX_TYPE>::VertexId &aDst,
                                   const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
//            auto maxFlow = std::get<0>(minStCutFordFulkersonBase(aGraph, aSrc, aDst, aWeights));

//            Graph<int, GraphMap> g;
//            std::unordered_map<int, int> em;
//            int nv = 0;
//            for (auto &v : aGraph.getVertices()) {
//                if (em.find(v) == em.end()) {
//                    em[v] = nv++;
//                }
//                g.addVertexSafe(em[v]);
//                for (auto &vo : aGraph.getOutVertices(v)) {
//                    if (em.find(vo) == em.end()) {
//                        em[vo] = nv++;
//                    }
//                    g.addVertexSafe(em[vo]);
//                    g.addEdge({em[v], em[vo]});
//                }
//            }
//            std::cout << "DINIC:" << g << " " << nv << std::endl;
//            auto c = Dinic().runDinic(g, em[aSrc], em[aDst], aWeights, em);
//            auto c = Dinic().runDinic(aGraph, aSrc, aDst, aWeights, stack.capacity());
            auto c = HIPR().runHipr(aGraph, aSrc, aDst, aWeights, stack.capacity(), parents);
//            std::cout << "VS: " << dinicc << " " << hiprc << std::endl;
            auto maxFlow = c;
        return maxFlow;
        }


        auto getRandomSubgraphNotBlue(Graph <VERTEX_TYPE> &aGraph, int aNumEdgesToRemove, const EdgesSet<VERTEX_TYPE> &blueEdges) {
            int edgesRemovedCnt = 0;
            typename Graph<VERTEX_TYPE>::Edge lastRndEdge{};
            while (edgesRemovedCnt < aNumEdgesToRemove) {
                auto edgesWithCycles = findEdgesWithCycles(aGraph);
//                auto edgesWithCycles = findNotBlueEdges(aGraph);
//                auto weights = Ext::getEdgeProperties(aGraph, 1);
//                auto [_, edgesWithCycles] = Fasp::GR(aGraph, weights);
                auto n = edgesWithCycles.size();
//                std::cout << "n: " << n << std::endl;
                if (n == 0) {
                    if (edgesRemovedCnt > 0) aGraph.addEdge(lastRndEdge);
                    else LOG(TRACE) << "Removed less edges than requested (" << edgesRemovedCnt << ", " << aNumEdgesToRemove << ")";
//                    LOG(TRACE) << edgesWithCycles;
                    return (edgesRemovedCnt > 0);
                } else {
                    lastRndEdge = edgesWithCycles[::Tools::randInt(0, n - 1)];
                    while (blueEdges.find(lastRndEdge) != blueEdges.end()) {lastRndEdge = edgesWithCycles[::Tools::randInt(0, n - 1)];}
                    aGraph.removeEdge(lastRndEdge);
                    ++edgesRemovedCnt;
                }
            }

            return (aNumEdgesToRemove > 0);
        }

        /**
         * Finds path from src to dst
         */
        auto findPathDfs(const Graph <VERTEX_TYPE> &aGraph,
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

        auto stronglyConnectedComponents2(const Graph<VERTEX_TYPE> &aGraph) {
            iVisited.clearAll();
            auto &sh = iVisited;
            stack.clear();

            int index_counter = 0;
            const int numOfV = aGraph.getNumOfVertices();
            for (std::size_t i = 0; i < lowLinks.size(); ++i) lowLinks[i] = -1;

            std::vector<std::unordered_set<VERTEX_TYPE>> result; result.reserve(numOfV);

            std::function<void(const VERTEX_TYPE &)> strongconnect = [&](const VERTEX_TYPE &node) {

                struct R {
                    const VERTEX_TYPE v;
                    const VERTEX_TYPE idx;
                    R(VERTEX_TYPE aV, VERTEX_TYPE aIdx) : v(aV), idx(aIdx) {}
                };

                std::vector<R> r; r.reserve(numOfV);
                r.push_back(R{node, 0});

                bool initRun = true;

                while (r.size() > 0) {

                    auto &b = r.back();
                    VERTEX_TYPE currentNode =  b.v;
                    int ci = b.idx;
                    r.pop_back();

                    processSuccessor:

                    if (initRun) {
                        index[currentNode] = index_counter;
                        lowLinks[currentNode] = index_counter;
                        ++index_counter;
                        stack.push(currentNode);
                        sh.set(currentNode);
                        initRun = false;
                    }
                    auto ov = aGraph.getOutVertices(currentNode);
                    for(std::size_t i = ci; i < ov.size(); ++i) {
                        auto successor = ov[i];
                        if (lowLinks[successor] == -1) {
                            // save the state (it would be recurrent call in default version of Trajan's algorithm)
                            r.push_back(R(currentNode, i + 1));
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
                            auto successor = stack.pop();
                            sh.clear(successor);
                            connectedComponent.emplace(successor);
                            if (successor == currentNode) break;
                        }
                        result.emplace_back(std::move(connectedComponent));
                    }
                    if (r.size() > 0) {
                        auto predecessor = r.back().v;
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
    };

    // ------------------------------------------------------------------------

    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
    static auto superAlgorithmBlue(Graph<VERTEX_TYPE> &outGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, PathHero<VERTEX_TYPE> &path, bool aUseWeights = false, bool relaxSA = false) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");
        typename Graph<VERTEX_TYPE>::Edges removedEdgesSA;
        typename Graph<VERTEX_TYPE>::Edges removedEdgesGR;

        EdgesSet<VERTEX_TYPE> setOfEdges{};

        int cnt = 0;
        while(true) {
            cnt++;
            bool wasGraphModified = false;

            setOfEdges.clear();

            for (const auto &e : outGraph.getEdges()) {
                if (setOfEdges.find(e) != setOfEdges.end()) continue; // e is in 'blue edges' set

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!path.pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                if (!path.GStarBlue(workGraph, e, setOfEdges, aUseWeights)) continue;

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
        if (removedEdgesSA.size() == 0 && relaxSA) {
            auto [maxMcRedEdge, redEdge] = path.getRedEdge(outGraph, aWeights);
            if (maxMcRedEdge > 0) {
                removedEdgesGR.push_back(redEdge);
            }
//            else {
//                auto ed = Fasp::GR(outGraph, aWeights).second;
//                auto n = ed.size();
//                if (n > 0) {
//                    removedEdgesGR.push_back(ed[0]);
//                }
//            }
        }
        return std::tuple{removedEdgesSA, setOfEdges, removedEdgesGR};
    }

    /**
     * Working original idea of how random FASP heuristic should work
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, bool WEIGHTED = false>
    static auto randomFASP(const Graph<VERTEX_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        if (WEIGHTED) std::cout << "Graph with WEIGHTS!\n";
        auto cleanGraphWithScc = [](Graph<VERTEX_TYPE> &aGraph, PathHero<VERTEX_TYPE> &path) {
            int cnt1 = 0, cntBig = 0;
            auto scc = path.stronglyConnectedComponents2(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) {cnt1++; aGraph.removeVertex(*s.begin());}
                else cntBig++;
            }
//            std::cout << "SCC  #1=" << cnt1 << " #BIG=" << cntBig << "\n";
        };

        constexpr int numOfReps = 20;



        auto g{aGraph};
        // find max vertex id in graph and setup PathHero
        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        PathHero<VERTEX_TYPE> path{static_cast<std::size_t>(maxId == vertices.end() ? 1 : *maxId + 1)};
//        std::cout << g << std::endl;
//        std::cout << "Max V = " << (maxId == vertices.end() ? 0 : *maxId) << std::endl;

        cleanGraphWithScc(g, path);

        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        int saEdgesCnt = 0;
        int saRndEdgesCnt = 0;
        int redRndEdgesCnt = 0;
        // initial run of superAlgorithm (SA)
        auto [edgesToRemove, blueEdgesx, dummy2] = superAlgorithmBlue(g, aWeights, path, WEIGHTED);
        auto blueEdges = blueEdgesx;
        g.removeEdges(edgesToRemove);
        saEdgesCnt += edgesToRemove.size();
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

        int numEdgesToRemoveInitVal = 3;
        //auto [dummy1, grEdges] = Fasp::GR(g, aWeights); int faspSize = grEdges.size(); numEdgesToRemoveInitVal = faspSize == 0 ? 1 : g.getNumOfEdges() / faspSize / 2;
//        std::cout << "EDGES TO REMOVE INIT VAL = " << numEdgesToRemoveInitVal << std::endl;
        int numEdgesToRemove = numEdgesToRemoveInitVal;
        [[maybe_unused]] int cnt = 1;
        while (true) {

//            std::cout << "----- loop=" << cnt++ << "\n";
            if (path.isAcyclic(g)) break;

            cleanGraphWithScc(g, path);

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCntGR;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
//            t.start_timer("random graphs");
//            alignas(64) Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> props[numOfReps];
//            for (auto &p : props) p = aWeights;
//            alignas(64) std::future<std::pair<typename Graph<VERTEX_TYPE>::Edges, typename Graph<VERTEX_TYPE>::Edges>> tasks[numOfReps];
//            int i = 0;
//            for (auto &task : tasks) {
//
//                task= std::async(std::launch::async,
//                     [&, i, numEdgesToRemove] () {
//                         Timer<false, false> tt(true);
//                         if (i == 0) tt.start_timer("1 - copy graph");
//                         auto workGraph{g};
//                         if (i == 0) {tt.stop_timer();}
//
//                         if (i == 0) tt.start_timer("2 - prepare path hero");
//                         auto vertices = workGraph.getVertices();
//                         auto maxId = std::max_element(vertices.begin(), vertices.end());
//                         PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included
//                         if (i == 0) tt.stop_timer();
//
//                         if (i == 0) tt.start_timer("3 - random subgraph");
//                         path.getRandomSubgraphNotBlue(workGraph, numEdgesToRemove, blueEdges);
//                         if (i == 0) tt.stop_timer();
//
//                         if (i == 0) tt.start_timer("4 - SA blue");
//                         auto [edgesSA, _, edgesGR] = superAlgorithmBlue(workGraph, props[i], path, false, true);
//                         if (i == 0) tt.stop_timer();
//                         return std::pair{edgesSA, edgesGR};
//                     });
//                i++;
//            }
//
//            for (auto &task : tasks) {
//                auto [edgesToRemove, edgesToRemoveGR] = task.get();
//                for (auto &e : edgesToRemove) edgesCnt.try_emplace(e, 0).first->second++;
//                for (auto &e : edgesToRemoveGR) edgesCntGR.try_emplace(e, 0).first->second++;
//            }
//            t.stop_timer();
            for (int i = 0; i < numOfReps; ++i) {
                 Timer<false, false, false> tt("RAND_GRAPHS");
                 if (i == 0) tt.start_timer("1 - copy graph");
                 auto workGraph{g};
                 if (i == 0) {tt.stop_timer();}

                 if (i == 0) tt.start_timer("2 - prepare path hero");
                 auto vertices = workGraph.getVertices();
                 auto maxId = std::max_element(vertices.begin(), vertices.end());
                 PathHero<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included
                 if (i == 0) tt.stop_timer();

                 if (i == 0) tt.start_timer("3 - random subgraph");
                 path.getRandomSubgraphNotBlue(workGraph, numEdgesToRemove, blueEdges);
                 if (i == 0) tt.stop_timer();

                 if (i == 0) tt.start_timer("4 - SA blue");
                 auto [edgesSA, _, edgesGR] = superAlgorithmBlue(workGraph, aWeights, path, WEIGHTED, true);
                 if (i == 0) tt.stop_timer();

                 if (i == 0) tt.start_timer("5 - Getting edges");
                 auto [edgesToRemove, edgesToRemoveGR] = std::pair{edgesSA, edgesGR};
                 for (auto &e : edgesToRemove) edgesCnt.try_emplace(e, 0).first->second++;
                 for (auto &e : edgesToRemoveGR) edgesCntGR.try_emplace(e, 0).first->second++;
                 if (i == 0) tt.stop_timer();
            }

            if (edgesCnt.size() > 0) saRndEdgesCnt++;
            else if (edgesCntGR.size() > 0) redRndEdgesCnt++;
            if (edgesCnt.size() == 0) {
//                std::cout << "------ USING alternative (GR/RED_EDGE) ------ " << edgesCntGR.size() <<"\n";
                edgesCnt.swap(edgesCntGR); // In case when SA didn't find any edges use these from GR heuristic
//                std::cout << "ER2: " << edgesCnt << "\n";
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

            // Try to solve the rest with superAlgoritm - maybe it will now succeed!
            auto [edgesToRemove, blueEdges2, dummy3] = superAlgorithmBlue(g, aWeights, path, WEIGHTED);
            blueEdges = blueEdges2;
            g.removeEdges(edgesToRemove);
            saEdgesCnt += edgesToRemove.size();
            removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
            removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());
        }

        EDGE_PROP_TYPE capacity = 0;
        for (const auto &e : removedEdges) { capacity += aWeights.at(e); }

        return std::tuple{capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt};
    }
}

#endif //FASPHEURISTIC_GRAPHFASP_H
