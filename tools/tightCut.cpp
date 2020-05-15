/*
 * Simple tool that reads from stdin unweighted simple graph data and calculates FASP size (number of edges to cut to make a graph acyclic)
 *
 * Graph format (first number is src vertex, next numbers are dst vertices of src vertex):
 *
 *     1	3    5
 *     3	5
 *     5    1
 *
 *     As a output (for above example) it will print out:
 *
 *     Input graph: #vertices=3 #edges=4
 *     FASP size (#edges to cut): 1
 */

#include <iostream>
#include <string>
#include <FaspTightCut/graphFasp.h>

int main() {
    using VERTEX_TYPE = int;

    // Read graph from input
    Graph::Graph<VERTEX_TYPE> graph;
    int lineCnt = 1;
    for (std::string line; std::getline(std::cin, line); lineCnt++) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

        // We expect at least vertex without any outgoing connections
        if (tokens.size() >= 1) {
            try {
                typename Graph::Graph<VERTEX_TYPE>::Vertex src = std::stoi(tokens[0]);
                graph.addVertexSafe(src);

                // Add edges to destination vertices if defined
                for (size_t i = 1; i < tokens.size(); ++i) {
                    typename Graph::Graph<VERTEX_TYPE>::Vertex dst = std::stoi(tokens[i]);
                    graph.addVertexSafe(dst);
                    graph.addEdge(src, dst);
                }
            }
            catch (const std::invalid_argument &e) {
                std::cout << "Invalid line #" << lineCnt << std::endl;
                std::cout << "\"" << line << "\"" << std::endl;
                std::cout<< "Error: '" << e.what() << "'" << std::endl;
                exit(1);
            }
        }
    }

    // Print info about read graph
    std::cout << "Input graph: #vertices=" << graph.getNumOfVertices() << " #edges=" << graph.getNumOfEdges() << std::endl;

    // Run TIGHT-CUT heuristic
    auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut<true, VERTEX_TYPE>(graph);

    // Print result
    std::cout << "FASP size (#edges to cut): " << capacity << std::endl;

    return 0;
}
