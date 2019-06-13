
#ifndef FASPHEURISTIC_BENCHBRUIJN_H
#define FASPHEURISTIC_BENCHBRUIJN_H

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"

#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphTools.h"
#include "graph/graphFasp.h"
#include "graph/graphFaspFast.h"
#include "graph/graphFaspFast2.h"
#include "graph/graphIO.h"



void benchBruijnGraphs(const std::string &outputDir) {

//    LOG(TRACE) << "Running benchBruijnGraphs";
//
//    std::string dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/Fasp/src/benchmarks/benchErdosRenyiGraphs/data";
//    std::string file = "erdos_renyi_n_100_c_5_seed_1.edges";
//
//    std::cout << " ======================= Processing [" << file << "] " << std::endl;

    if (false) {
        auto [gg, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(18, 13, 112);
        Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", gg);
        Graph::Graph g = Graph::IO::graphFromFile<int, Graph::GraphMap>("/tmp/graph.txt");
        Graph::Graph g2 = Graph::IO::graphFromFile<int, Graph::GraphMap>("/tmp/graph.txt");

        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        Graph::FaspFast::PathHero<int> path(maxId == vertices.end() ? 1 : *maxId + 1);
        Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);

        auto c = Graph::Ext::getEdgeProperties<int>(g, 1);
        auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);

        std::cout << " ================================= OLD +====================\n";
        std::cout << " OLD: " << Graph::FaspFast::superAlgorithm(g, c, path) << std::endl;
        std::cout << " ================================= NEW +====================\n";
        std::cout << " NEW: " << Graph::FaspFast2::superAlgorithm(g2, c2, path2) << std::endl;

//        Graph::FaspFast::randomFASP(g, c);
//        Graph::FaspFast2::randomFASP(g2, c2);
    }


    if (true) {
        auto [gg, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(30, 13, 150);
        std::string outF = "/Users/gonciarz/Documents/MOSAIC/work/repo/Fasp/src/benchmarks/benchErdosRenyiGraphs/data/erdos_renyi_n_100_c_5_seed_1.edges";

//        Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", gg); outF = "/tmp/graph.txt";
        Graph::Graph g = Graph::IO::graphFromFile<int, Graph::GraphMap>(outF);
        Graph::Graph g2 = Graph::IO::graphFromFile<int, Graph::GraphMap>(outF);


        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        Graph::FaspFast::PathHero<int> path(maxId == vertices.end() ? 1 : *maxId + 1);
        Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);


        auto c = Graph::Ext::getEdgeProperties<int>(g, 1);
        auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);

        Timer<true, false> t("");
        t.start_timer("old");
        Graph::FaspFast::deltaFASP(g, c);
        t.stop_timer();
        t.start_timer("new");
        Graph::FaspFast2::deltaFASP(g2, c2);
        t.stop_timer();
//        Graph::FaspFast::randomFASP(g, c);
//        Graph::FaspFast2::randomFASP(g2, c2);
    }

    if (false) {
        auto g = Graph::Graph<int, Graph::GraphMap>();
        for (int i = 1; i <= 8; ++i) g.addVertex(i);
        g.addEdge({1, 2});
        g.addEdge({2, 3});
        g.addEdge({3, 1});
        g.addEdge({2, 4});
        g.addEdge({4, 5});
        g.addEdge({5, 6});
        g.addEdge({6, 1});
        g.addEdge({4, 7});
        g.addEdge({7, 4});
        g.addEdge({6, 8});
        g.addEdge({8, 5});
        g.addEdge({3, 5});

        auto g2 = Graph::Graph(g);
//        auto g2 = Graph::Graph<int, Graph::GraphMap>();
//        for (int i = 1; i <= 8; ++i) g2.addVertex(i);
//        g2.addEdge({1, 2});
//        g2.addEdge({2, 3});
//        g2.addEdge({3, 1});
//        g2.addEdge({2, 4});
//        g2.addEdge({4, 5});
//        g2.addEdge({5, 6});
//        g2.addEdge({6, 1});
//        g2.addEdge({4, 7});
//        g2.addEdge({7, 4});
//        g2.addEdge({6, 8});
//        g2.addEdge({8, 5});
//        g2.addEdge({3, 5});



        std::cout << g.getStrRepresentationOfGraph() << std::endl;
        std::cout << g2.getStrRepresentationOfGraph() << std::endl;

        auto vertices = g.getVertices();
        auto maxId = std::max_element(vertices.begin(), vertices.end());
        Graph::FaspFast::PathHero<int> path(maxId == vertices.end() ? 1 : *maxId + 1);
        Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);

        path.GStar(g, {1, 2});
        std::cout << g.getStrRepresentationOfGraph() << std::endl;

        path2.GStar(g2, {1, 2});
        std::cout << g2.getStrRepresentationOfGraph() << std::endl;


    }
}


#endif
