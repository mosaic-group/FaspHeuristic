#ifndef GRAPHFASP_H
#define GRAPHFASP_H


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

namespace Graph::Fasp {

    template <typename VERTEX_TYPE>
    using EdgesSet = std::unordered_set<typename Graph<VERTEX_TYPE>::Edge, Ext::EdgeHasher<VERTEX_TYPE>>;

    /**
     * GraphSpeedUtils... class is a approach to speed up things. For all algorithms it has implemented it uses
     * pre-allocated memory - that is why it need to be initialized with maximum number of vertices (or if not
     * in sequence 0, 1, 2 ... then with a maximum value of vertex id + 1.
     * @tparam VERTEX_TYPE
     */
    template<typename VERTEX_TYPE>
    class alignas (64) GraphSpeedUtils {

        // Allocation of structures/memory/containers shared by all algorithms from GraphSpeedUtils class
        DynamicBitset <uint32_t, VERTEX_TYPE> iVisited;
        Stack <VERTEX_TYPE> stack;
        std::vector<VERTEX_TYPE> parents;
        std::vector<VERTEX_TYPE> lowLinks;
        std::vector<VERTEX_TYPE> index;
        std::vector<VERTEX_TYPE> path;

    public:
        explicit GraphSpeedUtils(std::size_t aN) : iVisited(aN), stack(aN), parents(aN, 0), lowLinks(aN, -1), index(aN, -1) {
            path.reserve(aN);
        }

        /**
         * Checks if there is a path from src to dst in given graph
         * @tparam FORWARD_SEARCH if false search is done by going from head to tail of each edge instead
         *                         of normal direction
         * @param aGraph input graph
         * @param aSrc start vertex
         * @param aDst destination vertex
         * @return true if path from src to dst exists
         */
        template<bool FORWARD_SEARCH=true>
        bool pathExistsDFS(const Graph <VERTEX_TYPE> &aGraph,
                           const typename Graph<VERTEX_TYPE>::Vertex &aSrc,
                           const typename Graph<VERTEX_TYPE>::Vertex &aDst) {
            // If we are already in destination job is done
            if (aSrc == aDst) return true;

            // Initialize needed structures
            iVisited.clearAll();
            iVisited.set(aSrc);
            stack.clear();
            stack.push(aSrc);

            while (!stack.empty()) {
                const auto currentVertex = stack.pop();

                // find all not visted vertices and add them to stack
                const auto &vertices = FORWARD_SEARCH ? aGraph.getOutVertices(currentVertex) : aGraph.getInVertices(currentVertex);
                for (const auto &v : vertices) {
                    if (v == aDst) return true;
                    if (iVisited.test(v) == 0) {
                        stack.push(v);
                        iVisited.set(v);
                    }
                }
            }

            return false;
        }

        /**
         * Checks if there is a path from src to dst in given graph
         * @param aGraph input graph
         * @param aSrc start vertex
         * @param aDst destination vertex
         * @return pair [pathExist, path] if (pathExist) then path contain vertices of path from src to dst, otherwise path may contain garbage
         */
        auto findPathDFS(const Graph <VERTEX_TYPE> &aGraph,
                         const typename Graph<VERTEX_TYPE>::Vertex &aSrc,
                         const typename Graph<VERTEX_TYPE>::Vertex &aDst) {

            path.clear();
            if (aGraph.hasVertex(aSrc) && aGraph.hasVertex(aDst)) {
                if (aSrc == aDst) {
                    path.emplace_back(aSrc);
                    return std::pair{true, path};
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
                        return std::pair{true, path};
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
            return std::pair{false, path};
        }

        /**
         * Checks if graph is acyclic
         * @param aGraph input graph
         * @return true if graph is acyclic
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

        /**
         * Find maximum flow / min cut using HIPR implementation from:
         * http://www.avglab.com/andrew/soft.html
         * @return maxFlow value
         */
        template<typename EDGE_PROP_TYPE>
        auto minStCut(const Graph <VERTEX_TYPE> &aGraph,
                      const typename Graph<VERTEX_TYPE>::Vertex &aSrc,
                      const typename Graph<VERTEX_TYPE>::Vertex &aDst,
                      const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
            auto maxFlow = HIPR().runHipr(aGraph, aSrc, aDst, aWeights, parents);
            return maxFlow;
        }

        /**
         * Finds strongly connected components in the provided graph. It is an implementation
         * of Trajan's algorithm converted from recursive to iterative version (it is
         * about 2x faster and much more cache friendly):
         * <a href=https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm>Trajan's algorithm on wiki</a>
         *
         * NOTE: It is using max value of VERTEX_TYPE to indicate unprocessed vertex
         *
         * @param aGraph - input graph
         * @return vector of sets, each element of vector is one SCC
         */
        auto stronglyConnectedComponents(const Graph<VERTEX_TYPE> &aGraph, bool aSkipSingleVertexSCCs = false) {
            stack.clear();
            // we will use iVisited to keep track of all elements in the stack - this is O(1) check instead of going
            // through all the stack elements
            iVisited.clearAll();

            int index_counter = 0;
            const int numOfV = aGraph.getNumOfVertices();

            const auto unprocessedVertex = std::numeric_limits<VERTEX_TYPE>::max();

            for (auto &lowLink : lowLinks) lowLink = unprocessedVertex;

            std::vector<std::unordered_set<VERTEX_TYPE>> result; result.reserve(numOfV);

            std::function<void(const VERTEX_TYPE &)> strongconnect = [&](const VERTEX_TYPE &node) {

                // state vector and 'S' structure are for storing current index and vertex to roll back to it after processing
                // children of this vertex. In regular Trajan algorithm this is not needed since we just call strongconnect recurently so
                // when we are back everything is there. In iterative approach we need to somehow rollback to previous values.
                struct S {
                    VERTEX_TYPE vertex;
                    VERTEX_TYPE childIndex;
                    S(VERTEX_TYPE aV, VERTEX_TYPE aIdx) : vertex(aV), childIndex(aIdx) {}
                };
                std::vector<S> state; state.reserve(numOfV);
                state.push_back(S{node, 0});

                bool initRun = true;

                while (state.size() > 0) {
                    // we have process all the children so roll back to previous state
                    auto [currentNode, currChildIndex] = state.back();
                    state.pop_back();

                    // Yes, this is label for 'goto'. You cannot blame me if you do not try to write
                    // iterative approach by yourself.
                    processSuccessor:

                    if (initRun) {
                        index[currentNode] = index_counter;
                        lowLinks[currentNode] = index_counter;
                        ++index_counter;
                        stack.push(currentNode);
                        iVisited.set(currentNode);
                        initRun = false;
                    }
                    auto ov = aGraph.getOutVertices(currentNode);
                    for(std::size_t i = currChildIndex; i < ov.size(); ++i) {
                        auto successor = ov[i];
                        if (lowLinks[successor] == unprocessedVertex) {
                            // save the state (it would be recurrent call in default version of Trajan's algorithm)
                            state.push_back(S(currentNode, i + 1));
                            // set values for successor and repeat from beginning ('goto' is bad... I know).
                            currentNode = successor;
                            currChildIndex = 0;
                            initRun = true;
                            goto processSuccessor;
                        }
                        else if (iVisited.test(successor)) {
                            // the successor is in the stack and hence in the current strongly connected component (SCC)
                            lowLinks[currentNode] = std::min(lowLinks[currentNode], index[successor]);
                        }
                    }

                    // If `node` is a root node, pop the stack and generate an SCC
                    if (lowLinks[currentNode] == index[currentNode]) {
                        std::unordered_set<VERTEX_TYPE> connectedComponent; connectedComponent.reserve(stack.size());

                        while (true) {
                            auto successor = stack.pop();
                            iVisited.clear(successor);
                            connectedComponent.emplace(successor);
                            if (successor == currentNode) break;
                        }
                        if (!aSkipSingleVertexSCCs || connectedComponent.size() > 1) result.emplace_back(std::move(connectedComponent));
                    }
                    if (state.size() > 0) {
                        auto predecessor = state.back().vertex;
                        lowLinks[predecessor] = std::min(lowLinks[predecessor], lowLinks[currentNode]);
                    }
                }
            };

            for (auto &node : aGraph.getVertices()) {
                if (lowLinks[node] == unprocessedVertex) {
                    strongconnect(node);
                }
            }

            return result;
        }

        /**
         * Generates subgraph of input graph (modifies it) by removing wanted number of edges with cycles.
         * There will be one cycle left in a input graph (if there was at least one cycle in input graph).
         * @param[in, out] aGraph input/output graph
         * @param[in] aNumEdgesToRemove number of edges to be removed (migh be less if not sufficient number of edges with cycles).
         * @param[in]aBlueEdges edges that should not be removed
         */
        void getRandomSubgraph(Graph <VERTEX_TYPE> &aGraph, int aNumEdgesToRemove, const EdgesSet<VERTEX_TYPE> &aBlueEdges) {
            int edgesRemovedCnt = 0;
            typename Graph<VERTEX_TYPE>::Edge lastRemovedRndEdge{};

            while (edgesRemovedCnt < aNumEdgesToRemove) {
                // We need to find edges with cycles each time (previously removed edge could kill a cycle with a lot of edges).
                // From found edges remove 'blue edges' which should not be processed
                auto edgesWithCycles = findEdgesWithCycles(aGraph);
                edgesWithCycles.erase(
                        std::remove_if(edgesWithCycles.begin(), edgesWithCycles.end(), [&aBlueEdges] (const typename EdgesSet<VERTEX_TYPE>::value_type &edge) { return aBlueEdges.find(edge) != aBlueEdges.end(); }),
                        edgesWithCycles.end());

                auto n = edgesWithCycles.size();

                if (n == 0) {
                    // If there is no more cycles revert last removed edge - we still want ISO-CUT algorithm to have something to cut.
                    if (edgesRemovedCnt > 0) aGraph.addEdge(lastRemovedRndEdge);
                    return;
                } else {
                    lastRemovedRndEdge = edgesWithCycles[::Tools::randInt(0, n - 1)];
                    aGraph.removeEdge(lastRemovedRndEdge);
                    ++edgesRemovedCnt;
                }
            }
        }

        /**
         * Finds 'red edge' in a graph (edge which is a part of cycle of some edge 'e' and at the same time is
         * a part of another cycle of G\e). After finding first 'red edge' which is good enough (mincut of such edge is
         * bigger then capacity of edge 'e' we consider that good 'red edge' was found.
         * @tparam EDGE_PROP_TYPE
         * @param aGraph - input graph
         * @param aWeights - weights of edges
         * @param aEdges - set of edges to search
         * @return pair of mincut and red edge
         */
        template<typename EDGE_PROP_TYPE>
        auto getRedEdge(Graph<VERTEX_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, const typename Graph<VERTEX_TYPE>::Edges &aEdges) {
            using Edge = typename Graph<VERTEX_TYPE>::Edge;

            Edge maxRedEdge{};
            int maxMincutOfRedEdge = 0;

            for (const auto &e : aEdges) {
                // find any path being a part of cycle for this edge
                const auto [pathExists, path] = findPathDFS(aGraph, e.dst, e.src);
                if (!pathExists) continue;

                // find SCCs not including edge 'e'
                aGraph.removeEdge(e);
                auto scc = stronglyConnectedComponents(aGraph, true);
                aGraph.addEdge(e);

                // find 'red edges' in a found path
                std::vector<Edge> redEdges;
                for (std::size_t i = 1; i < path.size(); ++i) {
                    for (auto &s : scc) {
                        // if edge is a part of any connected component it is a 'red edge'!
                        if (s.find(path[i - 1]) != s.end() && s.find(path[i]) != s.end()) {
                            redEdges.emplace_back(Edge{path[i - 1], path[i]});
                            break;
                        }
                    }
                }

                // find the 'red edge' with highest mincut
                for (auto &redEdge : redEdges) {
                    auto mc = minStCut(aGraph, redEdge.dst, redEdge.src, aWeights);
                    if (mc > maxMincutOfRedEdge) {
                        maxMincutOfRedEdge = mc;
                        maxRedEdge = redEdge;
                    }
                }

                // if found 'red edge' has mincut higher then weight of edge 'e' we found good candidate!
                auto eCapacity = aWeights.at(e);
                if (maxMincutOfRedEdge > eCapacity) break;
            }

            return std::pair{maxMincutOfRedEdge, maxRedEdge};
        }

        /**
         * Takes input graph aGraph and removes all edges which are not isolated w.r.t. aEdge
         * Also updates 'blue edges' container with all edges that are left after removing SCCs. These edges
         * don't need to be processed since they would have same isolated cycles as aEdge.
         * @param aGraph - input graph (it will be modified).
         * @param aEdge - edge to work with
         * @param aBlueEdges - container to be updated with newly found 'blue edges'
         * @param weighted - set to true if graph is weighted (skips 'blue edges' optimization)
         * @return true if isolated cycles found, false otherwise
         */
        bool findIsolatedCycles(Graph <VERTEX_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge, EdgesSet<VERTEX_TYPE> &aBlueEdges, bool weighted) {
            // 1. Remove an edge of interest 'aEdge' and find all connected components bigger than 1
            //    They consist from edges which are cycles not belonging only to aEdge so remove them.
            aGraph.removeEdge(aEdge);
            const auto scc = stronglyConnectedComponents(aGraph, true);
            for (const auto &s : scc) {
                // we get sets of vertices from SCC, find all edges connecting vertices in given SCC and remove
                for (const auto &v : s) {
                    const auto outVertices = aGraph.getOutVertices(v);
                    for (const auto &vo : outVertices) {
                        if (s.find(vo) != s.end()) {
                            aGraph.removeEdge({v, vo});
                        }
                    }
                }
            }

            // 2 update blue edges
            if (!weighted)
                for (auto &e : aGraph.getEdges()) {
                    aBlueEdges.emplace(e);
                }

            // 3. Check if there is stil a path from head to tail of aEdge, if yes we have
            //    isolated cycles.
            return pathExistsDFS(aGraph, aEdge.dst, aEdge.src);
        }
    };

    // ------------------------------------------------------------------------

    template<bool aUseWeights = true, typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
    static auto isoCut(Graph<VERTEX_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights, GraphSpeedUtils<VERTEX_TYPE> &aUtils, bool aUseRelaxApproach = false) {
        typename Graph<VERTEX_TYPE>::Edges removedEdges;
        EdgesSet<VERTEX_TYPE> blueEdges;

        while(true) {
            bool wasGraphModified = false;

            blueEdges.clear();

            for (const auto &e : aGraph.getEdges()) {
                if (blueEdges.find(e) != blueEdges.end()) continue; // e is in 'blue edges' set

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!aUtils.pathExistsDFS(aGraph, e.dst, e.src)) continue;

                // Find isolated cycles of edge 'e'. If there are no iso-cycles then continue with next edge.
                auto isoCyclesOfCurrentEdge{aGraph};
                if (!aUtils.findIsolatedCycles(isoCyclesOfCurrentEdge, e, blueEdges, aUseWeights)) continue;

                // If we have weights assigned to edges then we need to do min-cut, if not it is always safe to remove current edge
                bool shouldRemoveCurrentEdge = true; // default for non weighted graphs
                if (aUseWeights) {
                    // Check if edge should be removed
                    auto mc = aUtils.minStCut(isoCyclesOfCurrentEdge, e.dst, e.src, aWeights);
                    if (mc < aWeights.at(e)) shouldRemoveCurrentEdge = false;
                }

                // If there is edge to remove we should also set flag to repeat whole edges loop.
                // Removing one edge can cause that other edges have isolated cycles.
                if (shouldRemoveCurrentEdge) {
                    wasGraphModified = true;
                    aGraph.removeEdge(e);
                    removedEdges.emplace_back(e);
                }
            }
            if (!wasGraphModified) break;
        }

        // If we have not found any edge we will use guess to find best candidate in relax mode.
        typename Graph<VERTEX_TYPE>::Edges removedEdgesRelaxed;
        if (removedEdges.size() == 0 && aUseRelaxApproach) {
            auto [maxMcRedEdge, redEdge] = aUtils.getRedEdge(aGraph, aWeights, aGraph.getEdges());
            if (maxMcRedEdge > 0) {
                removedEdgesRelaxed.push_back(redEdge);
            }
        }

        return std::tuple{removedEdges, blueEdges, removedEdgesRelaxed};
    }

    template<typename VERTEX_TYPE>
    static auto isoCut(Graph<VERTEX_TYPE> &aGraph, GraphSpeedUtils<VERTEX_TYPE> &aUtils, bool aUseRelaxApproach = false) {
        return isoCut<false>(aGraph, Ext::getEdgeProperties(aGraph, 1), aUtils, aUseRelaxApproach);
    }

    // ------------------------------------------------------------------------

    /**
     * Working original idea of how random FASP heuristic should work
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, bool WEIGHTED = false, bool PARALLELIZED=false>
    static auto randomFASP(const Graph<VERTEX_TYPE> &aGraph, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        if (WEIGHTED) std::cout << "Graph with WEIGHTS!\n";
        auto cleanGraphWithScc = [](Graph<VERTEX_TYPE> &aGraph, GraphSpeedUtils<VERTEX_TYPE> &path) {
            int cnt1 = 0, cntBig = 0;
            auto scc = path.stronglyConnectedComponents(aGraph);
            for (const auto &s : scc) {
                if (s.size() == 1) {cnt1++; aGraph.removeVertex(*s.begin());}
                else cntBig++;
            }
        };

        constexpr int numOfReps = 20;



        auto g{aGraph};
        // find max vertex id in graph and setup GraphSpeedUtils
        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        GraphSpeedUtils<VERTEX_TYPE> path{static_cast<std::size_t>(maxId == vertices.end() ? 1 : *maxId + 1)};

        cleanGraphWithScc(g, path);

        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        int saEdgesCnt = 0;
        int saRndEdgesCnt = 0;
        int redRndEdgesCnt = 0;
        // initial run of superAlgorithm (SA)
        auto [edgesToRemove, blueEdgesx, dummy2] = WEIGHTED ? isoCut(g, aWeights, path) : isoCut(g, path);
        auto blueEdges = blueEdgesx;
        g.removeEdges(edgesToRemove);
        saEdgesCnt += edgesToRemove.size();
        removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
        removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());

        int numEdgesToRemoveInitVal = 3;
        //auto [dummy1, grEdges] = Fasp::GR(g, aWeights); int faspSize = grEdges.size(); numEdgesToRemoveInitVal = faspSize == 0 ? 1 : g.getNumOfEdges() / faspSize / 2;
//        std::cout << "EDGES TO REMOVE INIT VAL = " << numEdgesToRemoveInitVal << std::endl;
        int numEdgesToRemove = numEdgesToRemoveInitVal;
        while (true) {

//            std::cout << "----- loop=" << cnt++ << "\n";
            if (path.isAcyclic(g)) break;

            cleanGraphWithScc(g, path);

            // if we are here there are still cycles not handled by SA
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCnt;
            std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, int, Ext::EdgeHasher<VERTEX_TYPE>> edgesCntGR;

            // Run each randomly generated graph in seperate thread and later collect all solutions found
            if (PARALLELIZED) {
                alignas(64) Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> props[numOfReps];
                for (auto &p : props) p = aWeights;
                alignas(64) std::future<std::pair<typename Graph<VERTEX_TYPE>::Edges, typename Graph<VERTEX_TYPE>::Edges>> tasks[numOfReps];
                int i = 0;
                for (auto &task : tasks) {

                    task= std::async(std::launch::async,
                         [&, i, numEdgesToRemove] () {
                             Timer<false, false> tt{};
                             if (i == 0) tt.start_timer("1 - copy graph");
                             auto workGraph{g};
                             if (i == 0) {tt.stop_timer();}

                             if (i == 0) tt.start_timer("2 - prepare path hero");
                             auto vertices = workGraph.getVertices();
                             auto maxId = std::max_element(vertices.begin(), vertices.end());
                             GraphSpeedUtils<VERTEX_TYPE> path(maxId == vertices.end() ? 1 : *maxId + 1); // maxId included
                             if (i == 0) tt.stop_timer();

                             if (i == 0) tt.start_timer("3 - random subgraph");
                             path.getRandomSubgraph(workGraph, numEdgesToRemove, blueEdges);
                             if (i == 0) tt.stop_timer();

                             if (i == 0) tt.start_timer("4 - SA blue");
                             auto [edgesSA, _, edgesGR] = WEIGHTED ? isoCut(workGraph, props[i], path, true) : isoCut(workGraph, path, true);
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
            }
            else {
                for (int i = 0; i < numOfReps; ++i) {
                    Timer<false, false, false> tt("RAND_GRAPHS");
                    if (i == 0) tt.start_timer("1 - copy graph");
                    auto workGraph{g};
                    if (i == 0) { tt.stop_timer(); }

                    if (i == 0) tt.start_timer("2 - prepare path hero");
                    auto vertices2 = workGraph.getVertices();
                    auto maxId2 = std::max_element(vertices2.begin(), vertices2.end());
                    GraphSpeedUtils<VERTEX_TYPE> path(maxId2 == vertices2.end() ? 1 : *maxId2 + 1); // maxId included
                    if (i == 0) tt.stop_timer();

                    if (i == 0) tt.start_timer("3 - random subgraph");
                    path.getRandomSubgraph(workGraph, numEdgesToRemove, blueEdges);
                    if (i == 0) tt.stop_timer();

                    if (i == 0) tt.start_timer("4 - SA blue");
                    auto[edgesSA, _, edgesGR] = WEIGHTED ? isoCut(workGraph, aWeights, path, true) : isoCut(workGraph, path, true);
                    if (i == 0) tt.stop_timer();

                    if (i == 0) tt.start_timer("5 - Getting edges");
                    auto[edgesToRemove, edgesToRemoveGR] = std::pair{edgesSA, edgesGR};
                    for (auto &e : edgesToRemove) edgesCnt.try_emplace(e, 0).first->second++;
                    for (auto &e : edgesToRemoveGR) edgesCntGR.try_emplace(e, 0).first->second++;
                    if (i == 0) tt.stop_timer();
                }
            }

            if (edgesCnt.size() > 0) saRndEdgesCnt++;
            else if (edgesCntGR.size() > 0) redRndEdgesCnt++;
            if (edgesCnt.size() == 0) {
                edgesCnt.swap(edgesCntGR); // In case when SA didn't find any edges use these from GR heuristic
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
            auto [edgesToRemove, blueEdges2, dummy3] = WEIGHTED ? isoCut(g, aWeights, path) : isoCut(g, path);
            blueEdges = blueEdges2;
            g.removeEdges(edgesToRemove);
            saEdgesCnt += edgesToRemove.size();
            removedEdges.reserve(removedEdges.size() + edgesToRemove.size());
            removedEdges.insert(removedEdges.begin(), edgesToRemove.begin(), edgesToRemove.end());
        }

        EDGE_PROP_TYPE capacity = 0;
        std::cout << aWeights << std::endl;
        for (const auto &e : removedEdges) { std::cout << " " << e ; capacity += aWeights.at(e); }

        return std::tuple{capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt};
    }
}

#endif
