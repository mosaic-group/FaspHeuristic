//
// Created by gonciarz on 2019-03-29.
//

#ifndef TEST_GRAPHTOIGRAPH_H
#define TEST_GRAPHTOIGRAPH_H


#include "graph.h"
#include <iostream>


template<typename VERTEX_TYPE>
void printIGraphCommands(const Graph::Graph<VERTEX_TYPE> &aGraph) {
    std::cout << "g = igraph.Graph(directed=True)" << std::endl;
    std::string vertexLabel = "g.vs[\"label\"] = [";
    std::string vertexColor = "g.vs[\"color\"] = [";
    for (const auto &v : aGraph.getVertices()) {
        std::cout << "g.add_vertex(str(" << v << "))\n";
        vertexLabel += "\"";
        vertexLabel += std::to_string(v);
        vertexLabel += "\", ";
        vertexColor += "\"#FFFACD\", ";

    }
    vertexLabel += "]";
    vertexColor += "]";
    std::cout << vertexLabel << std::endl;
    std::cout << vertexColor << std::endl;
    for (const auto &e : aGraph.getEdges()) {
        std::cout << "g.add_edge(str(" << e.src << "), str(" << e.dst << "))\n";
    }
    std::cout << "igraph.plot(g, layout = g.layout(\"circle\"), vertex_size = 30)" << std::endl;
}

#endif //TEST_GRAPHTOIGRAPH_H
