//
// Created by gonciarz on 2019-03-18.
//

#ifndef GRAPHFASP_H
#define GRAPHFASP_H


#include "graphToIgraph.h"
#include "graph.h"
#include "graphExt.h"
#include "graphTools.h"
#include "tools/tools.h"
#include "tools/easylogging++.h"
#include "tools/prettyprint.h"
#include "tools/timer.h"
#include <algorithm>
#include <random>
#include <chrono>


namespace Graph::Fasp {

    /**
     * FASP heuristic - implementation of "A fast and effective heuristic for the feedback arc set problem" Eades, Lin, Smyth 1993
     * @tparam EDGE_PROP_TYPE type of edge property - must be a signed number with weight of each edge in a graph
     * @tparam VERTEX_TYPE
     * @tparam GRAPH_TYPE
     * @param aGraph input graph
     * @param aWeights input weights
     * @return capacity of cut edges
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static EDGE_PROP_TYPE GR(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        typename Graph<VERTEX_TYPE>::Vertices s1, s2;
        s1.reserve(aGraph.getNumOfVertices());
        s2.reserve(aGraph.getNumOfVertices());
        Graph<VERTEX_TYPE, GRAPH_TYPE> g(aGraph);

        while (true) {
            // Handle sinks
            bool changed = false;
            do {
                changed = false;
                for (const auto &v : g.getVertices()) {
                    if (g.getOutVertices(v).size() == 0) {
                        changed = true;
                        g.removeVertex(v);
                        s2.insert(s2.begin(), v);
                    }
                }
            } while (changed);

            // Handle sources
            do {
                changed = false;
                for (const auto &v : g.getVertices()) {
                    if (g.getInVertices(v).size() == 0) {
                        changed = true;
                        g.removeVertex(v);
                        s1.emplace_back(v);
                    }
                }
            } while (changed);

            // Handle deltas
            const auto &vertices = g.getVertices();
            if (vertices.empty()) break;
            EDGE_PROP_TYPE maxDelta = std::numeric_limits<EDGE_PROP_TYPE>().lowest();
            typename Graph<VERTEX_TYPE>::VertexId maxDeltaVertex;
            for (const auto &v : vertices) {
                EDGE_PROP_TYPE temp = 0;
                for (const auto &vo : g.getOutVertices(v)){ temp += aWeights.at({v, vo});}
                for (const auto &vi : g.getInVertices(v)) { temp -= aWeights.at({vi, v});}
                if (temp > maxDelta) {
                    maxDelta = temp;
                    maxDeltaVertex = v;
                }
            }
            s1.emplace_back(maxDeltaVertex);
            g.removeVertex(maxDeltaVertex);
        }

        // Calculate capacity, and find all edges to cut
        s1.insert(s1.end(), s2.begin(), s2.end());
        size_t idx = 1;
        EDGE_PROP_TYPE capacity = 0;
        typename Graph<VERTEX_TYPE>::Edges edges;
        for (const auto &v : s1) {
            for (const auto &vi : aGraph.getInVertices(v)) {
                if (std::find(s1.begin(), s1.begin() + idx, vi) == s1.begin() + idx) {
                    capacity += aWeights.at({vi, v});
                    edges.emplace_back(typename Graph<VERTEX_TYPE>::Edge{vi ,v});
                }
            }
            ++idx;
        }

        // Debug printout and check of solution.
        LOG(DEBUG) << "FASP(GR)    capacity = " << capacity << " edgeCnt = " << edges.size() << " edgeList = " << edges;
        LOG(TRACE) << "Edges with cycles: " << Tools::findEdgesWithCycles(g);

        return capacity;
    }

    template<typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
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
            auto A = Tools::depthFirstSearch(aGraph, aEdge.src);
            auto B = Tools::depthFirstSearch(aGraph, aEdge.dst, true /* reverseSearch */);
            typename Graph<VERTEX_TYPE>::VerticesSet AandB;
            std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), std::inserter(AandB, AandB.begin()));
            const auto &vertices = aGraph.getVertices();
            if (vertices.size() == AandB.size()) return; // Optimization: V == A∩B
            for (auto v : vertices) {
                if (AandB.find(v) == AandB.end()) {
                    // Removes all vertices which are in set V\(A∩B)
                    aGraph.removeVertex(v);
                }
            }
        }
    }

    template<typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    void GStar(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Edge &aEdge) {
//        LOG(DEBUG) << aGraph << " " << aEdge;
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
                typename Graph<VERTEX_TYPE>::Edges inEdges;
                for (const auto &v : aGraph.getInVertices(e.src)) {
                    typename Graph<VERTEX_TYPE>::Edge inEdge{v, e.src};
                    if (e != inEdge) inEdges.emplace_back(std::move(inEdge));
                }
                aGraph.removeEdges(inEdges);
                if (e.src == e.dst || !Tools::pathExistsDFS(aGraph, e.dst, vTo)) {
                    aGraph.removeEdge(e);
                    wasGraphModified = true;
                    wasCurrentEdgeRemoved = true;
                }
                aGraph.addEdges(inEdges);

                if (!wasCurrentEdgeRemoved) {
                    typename Graph<VERTEX_TYPE>::Edges outEdges;
                    for (const auto &v : aGraph.getOutVertices(e.dst)) {
                        typename Graph<VERTEX_TYPE>::Edge outEdge{e.dst, v};
                        if (e != outEdge) outEdges.emplace_back(std::move(outEdge));
                    }
                    aGraph.removeEdges(outEdges);
                    if (e.src == e.dst || !Tools::pathExistsDFS(aGraph, vFrom, e.src)) {
                        aGraph.removeEdge(e);
                        wasGraphModified = true;
                    }
                    aGraph.addEdges(outEdges);
                }
            }

            if (!wasGraphModified) break;
        }
    }

    template<typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    auto step2b(Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Graph<VERTEX_TYPE, GRAPH_TYPE> &aCleanedGraph, const typename Graph<VERTEX_TYPE>::Edge aEdge) {
        // Remove edge to not find any path going through it
        aGraph.removeEdge(aEdge);

        typename Graph<VERTEX_TYPE>::Edges S;
        for (const auto &h : aCleanedGraph.getEdges()) {
            if (Tools::pathExistsDFS(aGraph, h.dst, h.src)) {
                S.emplace_back(h);
            }
        }

        // revert
        aGraph.addEdge(aEdge);
        return S;
    }

    /**
     * FASP heuristic - implementation of "A fast and effective heuristic for the feedback arc set problem" Eades, Lin, Smyth 1993
     * @tparam EDGE_PROP_TYPE type of edge property - must be a signed number with weight of each edge in a graph
     * @tparam VERTEX_TYPE
     * @tparam GRAPH_TYPE
     * @param aGraph input graph
     * @param aWeights input weights
     * @return capacity of cut edges
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto superAlgorithm(const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph, const Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> &aWeights) {
        assert(std::is_signed<EDGE_PROP_TYPE>::value && "Weights are expected to be signed type");

        typename Graph<VERTEX_TYPE>::Edges removedEdges;

        std::cout << "superAlgorithm START" << std::endl;
        auto outGraph{aGraph}; // eventually acyclic graph
        while(true) {
            bool wasGraphModified = false;

            for (const auto &e : outGraph.getEdges()) {

                // Optimization, if there is no path back from dst to src, then edge has no cycles.
                if (!Tools::pathExistsDFS(outGraph, e.dst, e.src)) continue;

                auto workGraph{outGraph};
                GStar(workGraph, e);

                auto S = step2b(outGraph, workGraph, e);

                workGraph.removeEdges(S);
                GStar(workGraph, e);

                auto mc = Tools::minStCutFordFulkerson(workGraph, e.dst, e.src, aWeights);

                if (mc >= aWeights.at(e)) {
                    std::cout << "    REMOVING EDGE: " << e << std::endl;
                    wasGraphModified = true;
                    outGraph.removeEdge(e);
                    removedEdges.emplace_back(e);
                }
            }

            if (!wasGraphModified) break;
        }
        std::cout << "superAlgorithm DONE" << std::endl;

        return removedEdges;
    }

    /**
     * Generates a graph with known FASP solution (minimal capacity size of edges to cut).
     *
     * @tparam EDGE_PROP_TYPE
     * @tparam VERTEX_TYPE
     * @tparam GRAPH_TYPE
     * @param aNumOfVertices number of vertices in a graph
     * @param aFaspCapacity wanted size of fasp solution (capacity of edges to cut)
     * @param aLowerBondOfNumOfEdges minimum number of edges (more edges might be added if needed to achieve wanted FASP)
     * @param aOnlyNewEdges if true it adds new edges only (filler edges), if false then if edge exists increases it weight.
     *
     * @return generated graph
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto generateGraphWithKnownFasp(int aNumOfVertices, int aFaspCapacity, int aLowerBondOfNumOfEdges, int aMaxRandomFaspValue = 10, bool aAddRandomFaspWeights = false, bool aOnlyNewEdges = true) {
        Graph<VERTEX_TYPE, GRAPH_TYPE> g;
        Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> c;

        // helper for easier generation of random ints - gives rand int in range [min, max]
        auto randInt = [](int min, int max) -> int {
            static std::mt19937 mt(std::random_device{}());
            return std::uniform_int_distribution<>(min, max)(mt);
        };

        auto randWeights = [&randInt](unsigned int numOfEdges, int capacity) {
            std::vector<EDGE_PROP_TYPE> v(numOfEdges);
            unsigned int edgeIdx = 0;
            while (capacity > 0) {
                if (edgeIdx == numOfEdges - 1) { v[edgeIdx] = capacity; break; }

                int c = randInt(0, capacity);
                v[edgeIdx++] = c;
                capacity -= c;
            }
            return v;
        };

        // add vertices in range: 0..aNumOfVertices-1
        for (int i = 0; i < aNumOfVertices; ++i) g.addVertex(i);

        // generate mapping to random order of generated vertices
        std::vector<VERTEX_TYPE> verticesShuffle(aNumOfVertices);
        for (size_t i = 0; i < verticesShuffle.size(); ++i) verticesShuffle[i] = i;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle (verticesShuffle.begin(), verticesShuffle.end(), std::default_random_engine(seed));


        int numOfArcs = 0;
        EDGE_PROP_TYPE finalCapacity = aFaspCapacity;

        // add requested number of feedback arcs
        for (int f = 0; f < aFaspCapacity; ++f) {
            // add random leftward arc (arcs to be cut) from i to j
            int i, j;
            int randomValueDelta = 0;

            while (true) {
                i = randInt(1, aNumOfVertices - 1);
                j = randInt(0, i - 1);

                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge newEdge = {verticesShuffle[i], verticesShuffle[j]};
                if (g.hasEdge(newEdge)) continue; // make sure we add new arc
                c[newEdge] = 1;
                g.addEdge(std::move(newEdge));
                numOfArcs++;
                if (aAddRandomFaspWeights) {
                    // increase randomly weight of edge
                    randomValueDelta = randInt(0, aMaxRandomFaspValue - 1);
                    c[newEdge] += randomValueDelta;
                    finalCapacity += randomValueDelta;
                }
                break;
            }


            // finish the cycle by adding rightward arc(s) from j to i
            std::vector<typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge> backEdges;
            while (j != i) {
                int k = randInt(j + 1, i);
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e = {verticesShuffle[j], verticesShuffle[k]};
                if (aAddRandomFaspWeights) backEdges.push_back(e);
                if (g.hasEdge(e)) {
                    c.at(e) += 1;
                }
                else {
                    c[e] = 1;
                    g.addEdge(std::move(e));
                    numOfArcs++;
                }
                j = k;
            }

            if (aAddRandomFaspWeights) {
                // we need to increase weights of righward arcs to be sure that cutting them is not
                // better option than FASP arcs which is our goal.
                for (const auto &e : backEdges) c[e] += randomValueDelta;
            }
        }

        // add additional rightward arcs (till lower bond is reached)
        while (numOfArcs < aLowerBondOfNumOfEdges) {
            int i = randInt(0, aNumOfVertices - 2);
            int j = randInt(i + 1, aNumOfVertices - 1);

            typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e = {verticesShuffle[i], verticesShuffle[j]};
            if (g.hasEdge(e)) {
                if (aOnlyNewEdges) continue;
                c.at(e) += 1;
            }
            else {
                c[e] = 1;
                g.addEdge(std::move(e));
                numOfArcs++;
            }
        }

        return std::tuple{g, c, finalCapacity};
    }


    /**
      * Generates a graph with known FASP solution - all edges will have same (=1) capacity.
      *
      * @tparam EDGE_PROP_TYPE
      * @tparam VERTEX_TYPE
      * @tparam GRAPH_TYPE
      * @param aNumOfVertices number of vertices in a graph
      * @param aFaspCapacity wanted size of fasp solution (capacity of edges to cut)
      * @param aNumOfEdges number of edges
      * @param aOnlyNewEdges if true it adds new edges only (filler edges), if false then if edge exists increases it weight.
      *
      * @return generated graph or empty graph in case when input parameters (#V/#E/#FASP) are wrong
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto generateGraphWithKnownFaspAndSameWeights(int aNumOfVertices, int aFaspCapacity, int aNumOfEdges) {
        // Create graph and property container for capacity
        Graph<VERTEX_TYPE, GRAPH_TYPE> g;
        Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> c;

        // If requested parameters are wrong, return empty graph.
        const int maxNumOfEdges = aNumOfVertices * (aNumOfVertices - 1);
        bool wrongMaxNumberOfLeftwardEdges = aFaspCapacity > maxNumOfEdges/2; // max number of leftward edges
        bool wrongNumOfRequestedEdgesTooLow = aNumOfEdges < 2 * aFaspCapacity; // requested num of edges lower than needed to generate fasp
        bool wrongNumOfRequestedEdgesTooHigh = aNumOfEdges > maxNumOfEdges/2 + aFaspCapacity; // requested num of edges higher than sum of leftward edges and max number of rightward edges
        if (wrongMaxNumberOfLeftwardEdges || wrongNumOfRequestedEdgesTooLow || wrongNumOfRequestedEdgesTooHigh) {
                LOG(ERROR) << "Graph not generated!!! " << wrongMaxNumberOfLeftwardEdges << "/" << wrongNumOfRequestedEdgesTooLow << "/" << wrongNumOfRequestedEdgesTooHigh;
                return std::pair{g, c};
        }

        // add vertices in range: 0..aNumOfVertices-1
        for (int i = 0; i < aNumOfVertices; ++i) g.addVertex(i);

        // generate mapping to random order of generated vertices
        std::vector<VERTEX_TYPE> verticesShuffle(aNumOfVertices);
        for (size_t i = 0; i < verticesShuffle.size(); ++i) verticesShuffle[i] = i;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle (verticesShuffle.begin(), verticesShuffle.end(), std::default_random_engine(seed));

        // helper for easier generation of random ints - gives rand int in range [min, max]
        auto randInt = [](int min, int max) -> int {
            static std::mt19937 mt(std::random_device{}());
            return std::uniform_int_distribution<>(min, max)(mt);
        };

        // Let's keep a log of all fasp arc created it will be useful later
        std::vector<std::pair<int, int>> faspSet;
        int numOfArcs = 0;

        // Generate random fasp set arcs with backward arcs forming nice cycles (i ---> j (fasp arc leftward) and j ---> i (cycle arc rightward))
        for (int f = 0; f < aFaspCapacity; ++f) {
            while (true) {
                int i = randInt(1, aNumOfVertices - 1);
                int j = randInt(0, i - 1);

                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge faspEdge = {verticesShuffle[i], verticesShuffle[j]};
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge cycleEdge = {verticesShuffle[j], verticesShuffle[i]};

                if (g.hasEdge(faspEdge)) continue; // make sure we add new arc

                // Add edges to graph and update capacities
                c[faspEdge] = 1;
                c[cycleEdge] = 1;
                g.addEdge(std::move(faspEdge));
                g.addEdge(std::move(cycleEdge));
                numOfArcs += 2;
                faspSet.emplace_back(std::pair{i, j});
                break;
            }
        }

        // Try to randomize cycle arcs a litte bit. For cyccle edge from j ---> i try to find random number
        // k and divide this edge into two  j ---> k and k ---> i. Then repeat process. If not lucky in finding new
        // edges just leave it as it is.
        for (int f = 0; f < aFaspCapacity; ++f) {
            int randIdx = randInt(0, faspSet.size() - 1);
            //  j ---> i
            auto [i, j] = faspSet[randIdx];
            faspSet.erase(faspSet.begin() + randIdx);

            while (j != i) {
                // Control number of edges
                if (numOfArcs >= aNumOfEdges) break;

                if (j + 1 == i) break; // no more intermediate points

                int k = randInt(j + 1, i - 1);

                // j =====> k =====> i
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e1 = {verticesShuffle[j], verticesShuffle[k]};
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e2 = {verticesShuffle[k], verticesShuffle[i]};

                if (g.hasEdge(e1) || g.hasEdge(e2)) break; // no luck - give  up

                // Add two new edges j->k and k->i and remove old one j->i
                c[e1] = 1;
                c[e2] = 1;
                g.addEdge(std::move(e1));
                g.addEdge(std::move(e2));
                g.removeEdge({verticesShuffle[j], verticesShuffle[i]});
                c.erase({verticesShuffle[j], verticesShuffle[i]});
                ++numOfArcs; // one edge removed and two added so +1
                j = k;
            }
        }

        // add additional rightward arcs (till lower bond is reached)
        while (numOfArcs < aNumOfEdges) {
            int i = randInt(0, aNumOfVertices - 2);
            int j = randInt(i + 1, aNumOfVertices - 1);

            typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e = {verticesShuffle[i], verticesShuffle[j]};

            if (g.hasEdge(e)) continue;

            c[e] = 1;
            g.addEdge(std::move(e));
            ++numOfArcs;
        }

        return std::pair{g, c};
    }

    /**
  * Generates a graph with known FASP solution - all edges will have same (=1) capacity.
  *
  * @tparam EDGE_PROP_TYPE
  * @tparam VERTEX_TYPE
  * @tparam GRAPH_TYPE
  * @param aNumOfVertices number of vertices in a graph
  * @param aFaspCapacity wanted size of fasp solution (capacity of edges to cut)
  * @param aNumOfEdges number of edges
  * @param aOnlyNewEdges if true it adds new edges only (filler edges), if false then if edge exists increases it weight.
  *
  * @return generated graph or empty graph in case when input parameters (#V/#E/#FASP) are wrong
 */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template <typename> class GRAPH_TYPE>
    static auto generateGraphWithKnownFaspAndRandWeights(int aNumOfVertices, int aFaspCapacity, int aNumOfEdges) {
        // Create graph and property container for capacity
        Graph<VERTEX_TYPE, GRAPH_TYPE> g;
        Ext::EdgeProperties<EDGE_PROP_TYPE, VERTEX_TYPE> c;

        // If requested parameters are wrong, return empty graph.
        const int maxNumOfEdges = aNumOfVertices * (aNumOfVertices - 1);
        bool wrongMaxNumberOfLeftwardEdges = aFaspCapacity > maxNumOfEdges/2; // max number of leftward edges
        bool wrongNumOfRequestedEdgesTooLow = aNumOfEdges < 2 * aFaspCapacity; // requested num of edges lower than needed to generate fasp
        bool wrongNumOfRequestedEdgesTooHigh = aNumOfEdges > maxNumOfEdges/2 + 1; // requested num of edges higher than sum of leftward edges and max number of rightward edges
        if (wrongMaxNumberOfLeftwardEdges || wrongNumOfRequestedEdgesTooLow || wrongNumOfRequestedEdgesTooHigh) {
            LOG(ERROR) << "Graph not generated!!! " << wrongMaxNumberOfLeftwardEdges << "/" << wrongNumOfRequestedEdgesTooLow << "/" << wrongNumOfRequestedEdgesTooHigh;
            return std::pair{g, c};
        }

        // add vertices in range: 0..aNumOfVertices-1
        for (int i = 0; i < aNumOfVertices; ++i) g.addVertex(i);

        // generate mapping to random order of generated vertices
        std::vector<VERTEX_TYPE> verticesShuffle(aNumOfVertices);
        for (size_t i = 0; i < verticesShuffle.size(); ++i) verticesShuffle[i] = i;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle (verticesShuffle.begin(), verticesShuffle.end(), std::default_random_engine(seed));

        // helper for easier generation of random ints - gives rand int in range [min, max]
        auto randInt = [](int min, int max) -> int {
            static std::mt19937 mt(std::random_device{}());
            return std::uniform_int_distribution<>(min, max)(mt);
        };

        // Let's keep a log of all fasp arc created it will be useful later
        std::vector<std::pair<int, int>> faspSet;
        int numOfArcs = 0;

        // Generate random fasp set arcs with backward arcs forming nice cycles (i ---> j (fasp arc leftward) and j ---> i (cycle arc rightward))
        for (int f = 0; f < aFaspCapacity; ++f) {
            while (true) {
                int i = randInt(1, aNumOfVertices - 1);
                int j = randInt(0, i - 1);

                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge faspEdge = {verticesShuffle[i], verticesShuffle[j]};
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge cycleEdge = {verticesShuffle[j], verticesShuffle[i]};

                if (g.hasEdge(faspEdge)) continue; // make sure we add new arc

                // Add edges to graph and update capacities
                c[faspEdge] = 1;
                c[cycleEdge] = 1;
                g.addEdge(std::move(faspEdge));
                g.addEdge(std::move(cycleEdge));
                numOfArcs += 2;
                faspSet.emplace_back(std::pair{i, j});
                break;
            }
        }

        // Try to randomize cycle arcs a litte bit. For cyccle edge from j ---> i try to find random number
        // k and divide this edge into two  j ---> k and k ---> i. Then repeat process. If not lucky in finding new
        // edges just leave it as it is.
        for (int f = 0; f < aFaspCapacity; ++f) {
            int randIdx = randInt(0, faspSet.size() - 1);
            //  j ---> i
            auto [i, j] = faspSet[randIdx];
            faspSet.erase(faspSet.begin() + randIdx);

            while (j != i) {
                // Control number of edges
                if (numOfArcs >= aNumOfEdges) break;

                if (j + 1 == i) break; // no more intermediate points

                int k = randInt(j + 1, i - 1);

                // j =====> k =====> i
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e1 = {verticesShuffle[j], verticesShuffle[k]};
                typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e2 = {verticesShuffle[k], verticesShuffle[i]};

                if (g.hasEdge(e1) || g.hasEdge(e2)) break; // no luck - give  up

                // Add two new edges j->k and k->i and remove old one j->i
                c[e1] = 1;
                c[e2] = 1;
                g.addEdge(std::move(e1));
                g.addEdge(std::move(e2));
                g.removeEdge({verticesShuffle[j], verticesShuffle[i]});
                c.erase({verticesShuffle[j], verticesShuffle[i]});
                ++numOfArcs; // one edge removed and two added so +1
                j = k;
            }
        }

        // add additional rightward arcs (till lower bond is reached)
        while (numOfArcs < aNumOfEdges) {
            int i = randInt(0, aNumOfVertices - 2);
            int j = randInt(i + 1, aNumOfVertices - 1);

            typename Graph<VERTEX_TYPE, GRAPH_TYPE>::Edge e = {verticesShuffle[i], verticesShuffle[j]};

            if (g.hasEdge(e)) continue;

            c[e] = 1;
            g.addEdge(std::move(e));
            ++numOfArcs;
        }

        return std::pair{g, c};
    }
}


#endif
