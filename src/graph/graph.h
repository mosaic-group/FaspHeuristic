//
// Created by gonciarz on 2019-03-04.
//

#ifndef GRAPH_H
#define GRAPH_H

#include "tools/prettyprint.h"
#include "tools/tools.h"
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
    protected:
        // -------- Base Type defnitions ---------------------------------------
        using VertexId = VERTEX_TYPE;
        class Edge {
        public:
            VertexId src;
            VertexId dst;

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
        using Vertices = std::vector<VertexId>;
        using VerticesSet = std::set<VertexId>;
        struct NeighborVertices {
            Vertices to;
            Vertices from;
        };
    };

    template <typename VERTEX_TYPE>
    class GraphMap : public GraphBase<VERTEX_TYPE> {
    public:
        // Base type aliases
        using typename GraphBase<VERTEX_TYPE>::VertexId;
        using typename GraphBase<VERTEX_TYPE>::Edge;
        using typename GraphBase<VERTEX_TYPE>::Edges;
        using typename GraphBase<VERTEX_TYPE>::Vertices;
        using typename GraphBase<VERTEX_TYPE>::VerticesSet;
        using typename GraphBase<VERTEX_TYPE>::NeighborVertices;

        // -------- Type defnitions -----------------------------------------------
        using GraphDef = std::unordered_map<VERTEX_TYPE, NeighborVertices>;

    private:
        // -------- Internal data -------------------------------------------------
        inline static int iGlobalId = 0;
        inline static int iGlobalCreations = 0;
        int iMyId;
        int iBornId;
        GraphDef graph;
        static constexpr bool showConstructorMsg = false;

    public:
        GraphMap() : iMyId(++iGlobalId), iBornId(iMyId) {
            ++iGlobalCreations;
            if (showConstructorMsg) std::cout << "Graph() " << *this << std::endl;
        }

        GraphMap(const GraphMap &obj) : iMyId(++iGlobalId), iBornId(iMyId), graph(obj.graph) {
            ++iGlobalCreations;
            if (showConstructorMsg) std::cout << "Graph(Graph&) " << *this << " from " << obj << std::endl;
        }

        GraphMap(const GraphMap &&obj) noexcept : iMyId(++iGlobalId), iBornId(iMyId) {
            iMyId = obj.iMyId;
            graph = std::move(obj.graph);
            if (showConstructorMsg) std::cout << "Graph(Graph&&) " << *this << " from " << obj << std::endl;
        }

        GraphMap &operator=(const GraphMap &obj) {
            iMyId = obj.iMyId;
            graph = obj.graph;
            if (showConstructorMsg) std::cout << "Graph=Graph& " << *this << " from " << obj << std::endl;
            return *this;
        }

        GraphMap &operator=(GraphMap &&obj) {
            iMyId = obj.iMyId;
            graph = std::move(obj.graph);
            if (showConstructorMsg) std::cout << "Graph=Graph&& " << *this << std::endl;
            return *this;
        }

        ~GraphMap() {
            if (showConstructorMsg) std::cout << "~Graph() " << *this << std::endl;
        }

        /**
         * Adds vertex to graph
         * @param aId id of a vertex
         */
        void addVertex(const VERTEX_TYPE &aId) {
            assert(graph.count(aId) == 0 && "Vertex already exists!");
            graph.emplace(aId, NeighborVertices());
        }

        /**
          * Adds vertex to graph or ignores if already there
          * @param aId id of a vertex
          */
        void addVertexSafe(const VERTEX_TYPE &aId) {
            if (graph.find(aId) == graph.end()) {
                graph.emplace(aId, NeighborVertices());
            }
        }

        /**
         * Removes vertex from graph with all edges that involve that vertex
         * @param aId id of a vertex
         */
        void removeVertex(const VERTEX_TYPE &aId) {
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
        bool hasVertex(const VERTEX_TYPE &aV) const {
            if (graph.find(aV) == graph.end()) return false;
            return true;
        }

        /**
         * Adds edge to the graph
         * @param aE edge to be added
         */
        void addEdge(const Edge &aE) {
            // Check if vertices exist
            assert(graph.count(aE.src) == 1 && "Vertex does not exists!");
            assert(graph.count(aE.dst) == 1 && "Vertex does not exists!");

            Vertices &to = graph[aE.src].to;
            Vertices &from = graph[aE.dst].from;

            // Check if there is no such a edge already (no dupplication allowed)
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
        void addEdge(const VERTEX_TYPE &aSrc, const VERTEX_TYPE &aDst) {addEdge({aSrc, aDst});}

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
            assert(graph.count(aE.src) == 1 && "Vertex does not exists!");
            assert(graph.count(aE.dst) == 1 && "Vertex does not exists!");

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
        const Vertices& getOutVertices(const VERTEX_TYPE &aId) const {
            assert(graph.count(aId) == 1 && "Vertex does not exists!");

            return graph.at(aId).to;
        }

        /**
         * @param aId vertex of interest
         * @return vertices incoming to vertex aId
         */
        const Vertices& getInVertices(const VERTEX_TYPE &aId) const {
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
                for(const auto &dst : entry.second.to) {
                    edges.emplace_back(Edge{entry.first, dst});
                }
            }
            return edges;
        }

        /**
         * @return String representation of graph
         */
        std::string getStrRepresentationOfGraph() const {
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
        friend std::ostream &operator<<(std::ostream &os, const GraphMap &obj) {
            os << "{Graph (" << obj.iMyId << "/" << obj.iBornId << "/" << GraphMap::iGlobalId << "/"
               << GraphMap::iGlobalCreations << ")" << obj.getNumOfVertices() << "/" << obj.getNumOfEdges() << "}";
            if (false) os << obj.getStrRepresentationOfGraph();
            return os;
        }
    };

    template <typename VERTEX_TYPE>
    class GraphVector : public GraphBase<VERTEX_TYPE>  {
    public:
        // Base type aliases
        using typename GraphBase<VERTEX_TYPE>::VertexId;
        using typename GraphBase<VERTEX_TYPE>::Edge;
        using typename GraphBase<VERTEX_TYPE>::Edges;
        using typename GraphBase<VERTEX_TYPE>::Vertices;
        using typename GraphBase<VERTEX_TYPE>::VerticesSet;
        using typename GraphBase<VERTEX_TYPE>::NeighborVertices;
        
        // -------- Type defnitions -----------------------------------------------
        struct GraphDef {
            std::vector<VERTEX_TYPE> vertices;
            std::vector<NeighborVertices> edges;
            size_t getIdx(const VERTEX_TYPE &id) const {
                const auto &itBegin = vertices.cbegin();
                const auto &itEnd = vertices.cend();
                return std::find(itBegin, itEnd, id) - itBegin;
            }
        };

    private:
        // -------- Internal data -------------------------------------------------
        inline static int iGlobalId = 0;
        inline static int iGlobalCreations = 0;
        int iMyId;
        int iBornId;
        GraphDef graph;

        static constexpr bool showConstructorMsg = false;

    public:
        GraphVector() : iMyId(++iGlobalId), iBornId(iMyId) {
            ++iGlobalCreations;
            if (showConstructorMsg) std::cout << "Graph() " << *this << std::endl;
        }

        GraphVector(const GraphVector &obj) : iMyId(++iGlobalId), iBornId(iMyId), graph(obj.graph) {
            ++iGlobalCreations;
            if (showConstructorMsg) std::cout << "Graph(Graph&) " << *this << " from " << obj << std::endl;
        }

        GraphVector(const GraphVector &&obj) noexcept : iMyId(++iGlobalId), iBornId(iMyId) {
            iMyId = obj.iMyId;
            graph = std::move(obj.graph);
            if (showConstructorMsg) std::cout << "Graph(Graph&&) " << *this << " from " << obj << std::endl;
        }

        GraphVector &operator=(const GraphVector &obj) {
            iMyId = obj.iMyId;
            graph = obj.graph;
            if (showConstructorMsg) std::cout << "Graph=Graph& " << *this << " from " << obj << std::endl;
            return *this;
        }

        GraphVector &operator=(GraphVector &&obj) {
            iMyId = obj.iMyId;
            graph = std::move(obj.graph);
            if (showConstructorMsg) std::cout << "Graph=Graph&& " << *this << std::endl;
            return *this;
        }

        ~GraphVector() {
            if (showConstructorMsg) std::cout << "~Graph() " << *this << std::endl;
        }

        /**
         * Adds vertex to graph
         * @param aId id of a vertex
         */
        void addVertex(const VERTEX_TYPE &aId) {
            assert(graph.getIdx(aId) == graph.vertices.size() && "Vertex already exists!");

            graph.vertices.emplace_back(aId);
            graph.edges.emplace_back(NeighborVertices());
        }

        /**
          * Adds vertex to graph or ignores if already there
          * @param aId id of a vertex
          */
        void addVertexSafe(const VERTEX_TYPE &aId) {
            if (std::find(graph.vertices.begin(), graph.vertices.end(), aId) == graph.vertices.end()) {
                addVertex(aId);
            }
        }

        /**
         * Removes vertex from graph with all edges that involve that vertex
         * @param aId id of a vertex
         */
        void removeVertex(const VERTEX_TYPE &aId) {
            assert(graph.getIdx(aId) < graph.vertices.size() && "Vertex does not exists!");

            auto idx = graph.getIdx(aId);

            // Remove connection from other vertices that are pointing to or being pointed to aId
            for (const auto &src : graph.edges[idx].from) {
                auto idxSrc = graph.getIdx(src);
                auto &vTo = graph.edges[idxSrc].to;
                vTo.erase(std::remove(vTo.begin(), vTo.end(), aId), vTo.end());
            }
            for (const auto &dst : graph.edges[idx].to) {
                auto idxDst = graph.getIdx(dst);
                auto &vFrom = graph.edges[idxDst].from;
                vFrom.erase(std::remove(vFrom.begin(), vFrom.end(), aId), vFrom.end());
            }

            // Remove aId from structures
            graph.edges.erase(graph.edges.begin() + idx);
            graph.vertices.erase(graph.vertices.begin() + idx);
        }

        /**
         * @param aV vertex
         * @return true if vertex exists
         */
        bool hasVertex(const VERTEX_TYPE &aV) const {
            if (graph.getIdx(aV) < graph.vertices.size()) return true;
            return false;
        }

        /**
         * Adds edge to the graph
         * @param aE edge to be added
         */
        void addEdge(const Edge &aE) {
            // Check if vertices exist
            assert(graph.getIdx(aE.src) < graph.vertices.size() && "Vertex does not exists!");
            assert(graph.getIdx(aE.dst) < graph.vertices.size() && "Vertex does not exists!");

            size_t srcIdx = graph.getIdx(aE.src);
            size_t dstIdx = graph.getIdx(aE.dst);
            Vertices &to = graph.edges[srcIdx].to;
            Vertices &from = graph.edges[dstIdx].from;

            // Check if there is no such a edge already (no dupplication allowed)
            assert(std::find(to.begin(), to.end(), aE.dst) == to.end() && "Destination vertex exists!");
            assert(std::find(from.begin(), from.end(), aE.src) == from.end() && "Source vertex exists!");

            to.emplace_back(aE.dst);
            from.emplace_back(aE.src);
        }

        /**
         * Add edge
         * @param aSrc source vertex
         * @param aDst destination vertex
         */
        void addEdge(const VERTEX_TYPE &aSrc, const VERTEX_TYPE &aDst) {addEdge({aSrc, aDst});}

        /**
         * Add set of edges to graph
         * @param aEdges
         */
        void addEdges(const Edges &aEdges) {
            for (const auto &e : aEdges) addEdge(e);
        }

        /**
         * Checks if there is an edge in graph
         * @param aE edge to find
         * @return true if edge exists, false otherwise
         */
        bool hasEdge(const Edge &aE) {
            size_t srcIdx = graph.getIdx(aE.src);
            if (srcIdx == graph.vertices.size()) return false;
            Vertices &to = graph.edges[srcIdx].to;
            if (std::find(to.begin(), to.end(), aE.dst) == to.end()) return false;

            return true;
        }

        /**
         * Removes edge from a graph
         * @param aE edge to be removed
         */
        void removeEdge(const Edge &aE) {
            assert(graph.getIdx(aE.src) < graph.vertices.size() && "Vertex does not exists!");
            assert(graph.getIdx(aE.dst) < graph.vertices.size() && "Vertex does not exists!");

            size_t srcIdx = graph.getIdx(aE.src);
            size_t dstIdx = graph.getIdx(aE.dst);

            auto &vTo = graph.edges[srcIdx].to;
            vTo.erase(std::remove(vTo.begin(), vTo.end(), aE.dst), vTo.end());
            auto &vFrom = graph.edges[dstIdx].from;
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
        const Vertices& getOutVertices(const VERTEX_TYPE &aId) const {
            assert(graph.getIdx(aId) < graph.vertices.size() && "Vertex does not exists!");

            size_t idx = graph.getIdx(aId);
            return graph.edges[idx].to;
        }

        /**
         * @param aId vertex of interest
         * @return vertices incoming to vertex aId
         */
        const Vertices& getInVertices(const VERTEX_TYPE &aId) const {
            assert(graph.getIdx(aId) < graph.vertices.size() && "Vertex does not exists!");

            size_t idx = graph.getIdx(aId);
            return graph.edges[idx].from;
        }

        /**
         * Produces container with vertices
         * @return container with vertices
         */
        Vertices getVertices() const {
            return graph.vertices;
        }

        /**
         * Produces container with edges
         * @return contianer with edges
         */
        Edges getEdges() const {
            Edges edges;
            for (const auto &entry : graph.vertices) {
                size_t idx = graph.getIdx(entry);
                for(const auto &dst : graph.edges[idx].to) {
                    edges.emplace_back(Edge{entry, dst});
                }
            }
            return edges;
        }

        /**
         * @return String representation of graph
         */
        std::string getStrRepresentationOfGraph() const {
            std::ostringstream os;

            for (const auto &entry : graph.vertices) {
                os << entry << "->{";
                bool first = true;
                size_t idx = graph.getIdx(entry);
                for (const auto &dstIdx : graph.edges[idx].to) {
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
            return graph.vertices.size();
        }

        /**
         * @return number of edges in a graph
         */
        int getNumOfEdges() const {
            int cnt = 0;
            for (const auto &entry : graph.edges) {
                cnt += entry.to.size();
            }
            return cnt;
        }

        /**
         * Outputs graph info to std::cout
         */
        friend std::ostream &operator<<(std::ostream &os, const GraphVector &obj) {
            os << "{Graph (" << obj.iMyId << "/" << obj.iBornId << "/" << GraphVector::iGlobalId << "/"
               << GraphVector::iGlobalCreations << ")" << obj.getNumOfVertices() << "/" << obj.getNumOfEdges() << "}";
            if (false) os << obj.getStrRepresentationOfGraph();
            return os;
        }
    };

    template <typename VERTEX_TYPE = uint16_t, template <typename> class GRAPH_TYPE = GraphVector>
    class Graph : public GRAPH_TYPE<VERTEX_TYPE> {};
}

#endif
