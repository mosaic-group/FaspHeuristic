
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
#include "graph/graphIO.h"
#include "graph/graphFaspFast2.h"



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

        int cntA = 0;
        int cntB = 0;
        int cntAB = 0;
        std::vector<double> timesGS;
        std::vector<double> timesGS2;
        for (int i = 50; i <= 50; i += 30) {
            Timer<true, false> t("");
            int rep = 25;
            std::vector<double> tsa;
            std::vector<int> sa;
            double ct =0;
            int cn = 0;
            for (int r = 0; r < rep; ++r) {
//                auto [ge, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(i, 15, 4 * i);
                auto[ge, cc] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(i, 3.5 * (double) i / (i * (i-1)));

                Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", ge);
                Graph::Graph gg = Graph::IO::graphFromFile<int, Graph::GraphMap>("/tmp/graph.txt");

                auto vertices = gg.getVertices();
                auto maxId = std::max_element(vertices.begin(), vertices.end());
                Graph::FaspFast::PathHero<int> path1(maxId == vertices.end() ? 1 : *maxId + 1);
                Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);
                auto g1{gg};
                auto g2{gg};
                auto c1 = Graph::Ext::getEdgeProperties<int>(g1, 1);
                auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);


//                t.start_timer("G");
//                Graph::FaspFast::randomFASP(g1, c1);
//                std::cout << "CNT1: " << path1.cnt << std::endl;
//                t.stop_timer();

                    int b = 0;
                t.start_timer("G2 orig");
                b = Graph::FaspFast2::randomFASP_orig(g2, c2);
                std::cout << "CNT SA:" << path2.saCnt << std::endl;
                tsa.push_back(t.stop_timer());
                ct += tsa.back();
                sa.push_back(path2.saCnt);
                cn += path2.saCnt;
                path2.saCnt = 0;

                auto [ca, ed] = Graph::Fasp::GR(g2, c2);
                std::cout << "GR CAPACITY = " << ca << std::endl;

//
//
//                t.start_timer("G2 parallel");
//                int a = Graph::FaspFast2::randomFASP(g2, c2);
//                std::cout << "CNT SA: " << path2.saCnt << std::endl;
//                tsa.push_back(t.stop_timer());
//                ct += tsa.back();
//                sa.push_back(path2.saCnt);
//                cn += path2.saCnt;
//                path2.saCnt = 0;

                t.start_timer("G2 new");
                int a = Graph::FaspFast2::randomFASP_blueEdges(g1, c1);
                std::cout << "CNT SA: " << path2.saCnt << std::endl;
                tsa.push_back(t.stop_timer());
                ct += tsa.back();
                sa.push_back(path2.saCnt);
                cn += path2.saCnt;
                path2.saCnt = 0;

                if (a > b) cntB++;
                else if (b > a) cntA++;
                else cntAB++;
//
                std::cout << "======================>      " << cntA << " " << cntB << " " << cntAB << std::endl;
            }
            std::cout << "t = " << tsa << ";\n";
            std::cout << "n = " << sa << ";\n";
            std::cout << "SA time = " << ct / cn << std::endl;


        }
//        std::cout << timesSCC << std::endl;

        if (true) return;
            std::vector<double> times;
            std::vector<int> paths;
        std::vector<double> times2;
        std::vector<int> paths2;

        for (int i = 0; i < 100; ++i ) {
            auto [gg, cc] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(75, 0.03);
            auto g1{gg};
            auto g2{gg};

            auto vertices = g1.getVertices();
            auto maxId = std::max_element(vertices.begin(), vertices.end());
            Graph::FaspFast::PathHero<int> path1(maxId == vertices.end() ? 1 : *maxId + 1);
            Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);

            auto c1 = Graph::Ext::getEdgeProperties<int>(g1, 1);
            auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);

            Timer<true, false> t("");

            t.start_timer("old rand");
            Graph::FaspFast::randomFASP(g1, c1);
            times.push_back(t.stop_timer());
            std::cout << "PPPPPPPP : " << path1.cnt << std::endl;
            int x = path1.cnt;
            paths.push_back(x);
            path1.cnt = 0;

            t.start_timer("new rand");
            Graph::FaspFast2::randomFASP_sequential(g2, c2);
            times2.push_back(t.stop_timer());
            paths2.push_back(x);

        }

        std::cout << "times = " << times << ";" << std::endl;
        std::cout << "paths = " << paths << ";" << std::endl;
        std::cout << "times2 = " << times2 << ";" << std::endl;
        std::cout << "paths2 = " << paths2 << ";" << std::endl;

    }

    if (false) {
        auto [gg, cc] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(50, 0.041);
//        auto [gg, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(50, 13, 100);
//        std::string outF = "/Users/gonciarz/Documents/MOSAIC/work/repo/Fasp/src/benchmarks/benchErdosRenyiGraphs/data/erdos_renyi_n_70_c_7_seed_18.edges";

//        Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", gg); outF = "/tmp/graph.txt";
//        Graph::Graph g = Graph::IO::graphFromFile<int, Graph::GraphMap>(outF);
//        Graph::Graph g2 = Graph::IO::graphFromFile<int, Graph::GraphMap>(outF);

        std::string outF = "/tmp/graph.txt";
        Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", gg);
        Graph::Graph gf = Graph::IO::graphFromFile<int, Graph::GraphMap>(outF);
        std::cout << gf << std::endl;

        auto g{gf};
        auto g2{gf};

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

        t.start_timer("old rand");
        Graph::FaspFast::randomFASP(g, c);
        t.stop_timer();
        t.start_timer("new rand");
        Graph::FaspFast2::randomFASP_sequential(g2, c2);
        t.stop_timer();

//        t.start_timer("old rand");
//        std::cout << Graph::FaspFast::superAlgorithm(g, c) << std::endl;
//        t.stop_timer();
//        t.start_timer("new rand");
//        std::cout << Graph::FaspFast2::superAlgorithm(g2, c2) << std::endl;
//        t.stop_timer();
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

//        path.GStar(g, {1, 2});
//        std::cout << g.getStrRepresentationOfGraph() << std::endl;
//
//        path2.GStar(g2, {1, 2});
//        std::cout << g2.getStrRepresentationOfGraph() << std::endl;


    }
}


#endif
