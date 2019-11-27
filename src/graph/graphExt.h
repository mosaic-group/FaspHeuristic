//
// Created by gonciarz on 2019-03-18.
//

#ifndef GRAPHEXT_H
#define GRAPHEXT_H

#include "graph.h"
#include <cstddef>
#include <map>
#include <unordered_map>

namespace Graph::Ext {
    template <typename VERTEX_TYPE>
    class EdgeHasher {
    public:
        std::size_t operator() (const typename Graph<VERTEX_TYPE>::Edge &key) const {
            typedef std::size_t result_type;
            const result_type s = std::hash<typename Graph<VERTEX_TYPE>::Vertex >{}(key.src);
            const result_type d = std::hash<typename Graph<VERTEX_TYPE>::Vertex >{}(key.dst);
            return (17 * 31 + s) * 31 + d;
        }
    };

    template<typename VERTEX_TYPE, typename EDGE_PROP_TYPE>
    using EdgeProperties = std::unordered_map<typename Graph<VERTEX_TYPE>::Edge, EDGE_PROP_TYPE, EdgeHasher<VERTEX_TYPE>>;

    template<typename VERTEX_TYPE, typename VERTEX_PROP_TYPE>
    using VertexProperties = std::unordered_map<typename Graph<VERTEX_TYPE>::Vertex, VERTEX_PROP_TYPE>;

    /**
     * Generates map to store properties with initial value for each edge.
     * @tparam EDGE_PROP_TYPE type of property
     * @param aInitValue init value for property
     * @return map wiht K: Edge and V: property of type T
     */
    template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
    auto getEdgeProperties(const Graph<VERTEX_TYPE> &graph, const EDGE_PROP_TYPE &aInitValue) {
        EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> p;
        for (const auto &srcVertex : graph.getVertices()) {
            for (const auto &dstVertex : graph.getOutVertices(srcVertex)) {
                p.emplace(typename Graph<VERTEX_TYPE>::Edge{srcVertex, dstVertex}, aInitValue);
            }
        }
        return p;
    }

    /**
     * Generates map to store properties with initial value for each vertex.
     * @tparam VERTEX_TYPE type of property
     * @param aInitValue init value for property
     * @return map wiht K: VertexId and V: property of type T
     */
    template<typename VERTEX_PROP_TYPE, typename VERTEX_TYPE>
    auto getVertexProperties(const Graph<VERTEX_TYPE> &graph, const VERTEX_PROP_TYPE &aInitValue = VERTEX_PROP_TYPE()) {
        VertexProperties <VERTEX_TYPE, VERTEX_PROP_TYPE> p;
        for (auto &v : graph.getVertices()) p.emplace(v, aInitValue);
        return p;
    }
}

#endif
