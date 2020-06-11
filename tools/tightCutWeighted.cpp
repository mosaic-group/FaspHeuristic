/*
 * Simple tool that reads from stdin unweighted simple graph data and calculates FASP size (number of edges to cut to make a graph acyclic)
 *
 * Graph format (srcVertex dstVertex weightOfEdge):
 *
 * 1 2 5
 * 2 3 2
 * 3 4 3
 * 2 4 1
 * 4 1 10
 *
 * As a output (for above example) it will print out:
 *
 * Input graph: #vertices=4 #edges=5
 * FASP size (#edges to cut): 2, capacity of edges: 3
 *
 * NOTE: This is simple implementation and there are no syntax/errors checks done.
 */


#include <iostream>
#include <string>
#include <FaspTightCut/graphFasp.h>


int main() {
    using VERTEX_TYPE = int;
    using EDGE_PROP_TYPE = double;

    // Read graph from input
    Graph::Graph<VERTEX_TYPE> graph;
    Graph::Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> weights;
    int lineCnt = 1;
    for (std::string line; std::getline(std::cin, line); lineCnt++) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

        // We expect 1 (just vertex) or 3 numbers in line (srcVertex dstVertex edgeWeight)
        if (tokens.size() == 1 || tokens.size() == 3) {
            try {
                // Read vertices and create edge
                typename Graph::Graph<VERTEX_TYPE>::Vertex src = std::stoi(tokens[0]);
                graph.addVertexSafe(src);

                if (tokens.size() == 3) {
                        typename Graph::Graph<VERTEX_TYPE>::Vertex dst = std::stoi(tokens[1]);
                        graph.addVertexSafe(dst);
                        graph.addEdge(src, dst);

                        // Read weight and update properites of graph
                        EDGE_PROP_TYPE w = static_cast<EDGE_PROP_TYPE>(std::stod(tokens[2]));
                        weights[{src, dst}] = w;
                }

            }
            catch (const std::invalid_argument &e) {
                std::cout << "Invalid line #" << lineCnt << std::endl;
                std::cout << "\"" << line << "\"" << std::endl;
                std::cout<< "Error: '" << e.what() << "'" << std::endl;
                exit(1);
            }
        }
        else if (tokens.size() > 0) {
            std::cout << "Invalid line #" << lineCnt << std::endl;
            std::cout << "\"" << line << "\"" << std::endl;
            std::cout<< "Wrong number of elements (" << tokens.size() << ") in line (expected 1 or 3)." << std::endl;
            exit(1);
        }
    }


    // Print info about read graph
    std::cout << "# Input graph: #vertices=" << graph.getNumOfVertices() << " #edges=" << graph.getNumOfEdges() << std::endl;

    // Run TIGHT-CUT heuristic
    auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut<true, true, EDGE_PROP_TYPE, VERTEX_TYPE>(graph, weights);

    // Print result
    std::cout << "# FASP size (#edges to cut): " << removedEdges.size() << ", capacity of edges: " << capacity << std::endl;
    std::cout << "# Feedback arcs:" << std::endl;
    for (const auto &edge : removedEdges)
        std::cout << edge.src << ' ' << edge.dst << '\n';

    return 0;
}
