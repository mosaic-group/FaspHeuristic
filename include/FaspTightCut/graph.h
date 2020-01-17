#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <functional>


namespace Graph {

    template <typename VERTEX_TYPE>
    class GraphBase {
        static_assert(std::is_integral<VERTEX_TYPE>::value && "VERTEX_TYPE is expected to be a integral number");

    protected:
        // -------- Base Type defnitions ---------------------------------------
        using Vertex = VERTEX_TYPE;
        class Edge {
        public:
            Vertex src;
            Vertex dst;
            Edge() : src(-1), dst(-1) {}
            Edge(VERTEX_TYPE aSrc, VERTEX_TYPE aDst) : src(aSrc), dst(aDst) {}
            friend std::ostream &operator<<(std::ostream &os, const Edge &obj) {
                os << "{" << obj.src << ", " << obj.dst << "}";
                return os;
            }
        };

        friend bool operator==(const Edge &l, const Edge &r) {
            return (l.src == r.src) && (l.dst == r.dst);
        }

        friend bool operator!=(const Edge &l, const Edge &r) {
            return (l.src != r.src) || (l.dst != r.dst);
        }

        friend bool operator<(const Edge &l, const Edge &r) {
            return (l.src < r.src) || (l.src == r.src && l.dst < r.dst);
        }

        using Edges = std::vector<Edge>;
        using Vertices = std::vector<Vertex>;
        using VerticesSet = std::set<Vertex>;
        struct NeighborVertices {
            Vertices to;
            Vertices from;
        };
    };

    template <typename VERTEX_TYPE>
    class Graph : public GraphBase<VERTEX_TYPE> {
    public:

        using MyType = Graph<VERTEX_TYPE>;

        // Base type aliases
        using typename GraphBase<VERTEX_TYPE>::Vertex;
        using typename GraphBase<VERTEX_TYPE>::Edge;
        using typename GraphBase<VERTEX_TYPE>::Edges;
        using typename GraphBase<VERTEX_TYPE>::Vertices;
        using typename GraphBase<VERTEX_TYPE>::VerticesSet;
        using typename GraphBase<VERTEX_TYPE>::NeighborVertices;

        // -------- Type defnitions -----------------------------------------------
        using GraphDef = std::unordered_map<Vertex, NeighborVertices>;

    private:
        GraphDef graph;

    public:
        Graph() {}

        Graph(const Graph &obj) : graph(obj.graph) {}

        Graph(const Graph &&obj) noexcept {
            graph = std::move(obj.graph);
        }

        Graph &operator=(const Graph &obj) {
            graph = obj.graph;
            return *this;
        }

        Graph &operator=(Graph &&obj) noexcept {
            graph = std::move(obj.graph);
            return *this;
        }

        ~Graph() {}

        /**
         * Adds vertex to graph
         * @param aId id of a vertex
         */
        void addVertex(const Vertex &aId) {
            assert(graph.count(aId) == 0 && "Vertex already exists!");
            graph.emplace(aId, NeighborVertices());
        }

        /**
          * Adds vertex to graph or ignores if already there
          * @param aId id of a vertex
          */
        void addVertexSafe(const Vertex &aId) {
            if (graph.find(aId) == graph.end()) {
                graph.emplace(aId, NeighborVertices());
            }
        }

        /**
         * Removes vertex from graph with all edges that involve that vertex
         * @param aId id of a vertex
         */
        void removeVertex(const Vertex &aId) {
            assert(graph.count(aId) == 1 && "Vertex does not exists!");

            for (const auto &src : graph[aId].from) {
                auto &vTo = graph[src].to;
                vTo.erase(std::remove(vTo.begin(), vTo.end(), aId), vTo.end());
            }
            for (const auto &dst : graph[aId].to) {
                auto &vFrom = graph[dst].from;
                vFrom.erase(std::remove(vFrom.begin(), vFrom.end(), aId), vFrom.end());
            }
            graph.erase(aId);
        }

        /**
         * @param aV vertex
         * @return true if vertex exists
         */
        bool hasVertex(const Vertex &aV) const {
            if (graph.find(aV) == graph.end()) return false;
            return true;
        }

        /**
         * Adds edge to the graph
         * @param aE edge to be added
         */
        void addEdge(const Edge &aE) {
            // Check if vertices exist
            assert(graph.count(aE.src) == 1 && "Source vertex does not exists!");
            assert(graph.count(aE.dst) == 1 && "Destination vertex does not exists!");

            Vertices &to = graph[aE.src].to;
            Vertices &from = graph[aE.dst].from;

            // Check if there is no such a edge already (no duplication allowed)
            assert(std::find(to.begin(), to.end(), aE.dst) == to.end() && "Destination vertex exists!");
            assert(std::find(from.begin(), from.end(), aE.src) == from.end() && "Source vertex exists!");

            to.emplace_back(aE.dst);
            from.emplace_back(aE.src);
        }

        /**
         * Add set of edges to graph
         * @param aEdges
         */
        void addEdges(const Edges &aEdges) {
            for (const auto &e : aEdges) addEdge(e);
        }

        /**
         * Adds edge to the graph
         * @param aSrc source vertex
         * @param aDst destination vertex
         */
        void addEdge(const Vertex &aSrc, const Vertex &aDst) {addEdge({aSrc, aDst});}

        /**
         * Check if there is an edge in graph
         * @param aE edge to check
         * @return true if edge exists
         */
        bool hasEdge(const Edge &aE) {
            if (graph.find(aE.src) == graph.end()) return false;
            Vertices &to = graph[aE.src].to;
            if (std::find(to.begin(), to.end(), aE.dst) == to.end()) return false;
            
            return true;
        }
        /**
         * Removes edge from a graph
         * @param aE edge to be removed
         */
        void removeEdge(const Edge &aE) {
            assert(graph.count(aE.src) == 1 && "Source vertex does not exists!");
            assert(graph.count(aE.dst) == 1 && "Destination vertex does not exists!");

            auto &vTo = graph[aE.src].to;
            vTo.erase(std::remove(vTo.begin(), vTo.end(), aE.dst), vTo.end());
            auto &vFrom = graph[aE.dst].from;
            vFrom.erase(std::remove(vFrom.begin(), vFrom.end(), aE.src), vFrom.end());
        }

        /**
         * Remove set of edges from graph
         * @param aEdges
         */
        void removeEdges(const Edges &aEdges) {
            for (const auto &e : aEdges) removeEdge(e);
        }

        /**
         * @param aId vertex of interest
         * @return vertices outgoing from vertex aId
         */
        const Vertices& getOutVertices(const Vertex &aId) const {
            assert(graph.count(aId) == 1 && "Vertex does not exists!");

            return graph.at(aId).to;
        }

        /**
         * @param aId vertex of interest
         * @return vertices incoming to vertex aId
         */
        const Vertices& getInVertices(const Vertex &aId) const {
            assert(graph.count(aId) == 1 && "Vertex does not exists!");

            return graph.at(aId).from;
        }

        /**
         * Produces container with vertices
         * @return container with vertices
         */
        Vertices getVertices() const {
            Vertices v;
            for (const auto &e : graph) v.emplace_back(e.first);
            return v;
        }

        /**
         * Produces container with edges
         * @return container with edges
         */
        Edges getEdges() const {
            Edges edges;
            for (const auto &entry : graph) {
                for (const auto &dst : entry.second.to) {
                    edges.emplace_back(Edge{entry.first, dst});
                }
            }
            return edges;
        }

        /**
         * @return String representation of graph
         */
        std::string toString() const {
            std::ostringstream os;
            for (const auto &entry : graph) {
                os << entry.first << "->{";
                bool first = true;
                for (const auto &dstIdx : entry.second.to) {
                    if (!first) os << ",";
                    os << dstIdx;
                    first = false;
                }
                os << "}  ";
            }

            return os.str();
        }

        /**
         * @return number of vertices in a graph
         */
        int getNumOfVertices() const {
            return graph.size();
        }

        /**
         * @return number of edges in a graph
         */
        int getNumOfEdges() const {
            int cnt = 0;
            for (const auto &entry : graph) {
                cnt += entry.second.to.size();
            }
            return cnt;
        }

        /**
         * Outputs graph info to std::cout
         */
        friend std::ostream &operator<<(std::ostream &os, const Graph &obj) {
            os << "{Graph V/E: " << obj.getNumOfVertices() << "/" << obj.getNumOfEdges() << "}";
            return os;
        }
    };
}

#endif
