//
// Created by gonciarz on 2019-03-05.
//

#include "graph.h"
#include "graphTools.h"
#include "graphFasp.h"
#include "graphFaspFast.h"
#include "graphIO.h"

#include "tools/tools.h"
#include "tools/timer.h"
#include "tools/prettyprint.h"
#include "tools/easylogging++.h"

#include "hdf5/dataHdf5.h"

#include <iostream>
#include <vector>

void testSuperAlgorithm() {
    Timer<true, false> t(true);



    int numOfVertices = 30;
    int reps = 60;

    DataHdf5<double> f("/tmp/out.h5");

    for (int nv = 1; nv < 20 ; nv += 1) {
        std::cout << "FASP SIZE: " << nv << std::endl;
        for (int r = 0; r < reps; r++) {
            using EDGE_PROP_TYPE = int;
            auto[g, c, _] = Graph::Fasp::generateGraphWithKnownFasp<EDGE_PROP_TYPE, int, Graph::GraphMap>(numOfVertices, nv, 60, true);

            Tools::printType<decltype(g)>();
            std::cout << Tools::demangle(typeid(g).name()) << std::endl;
            std::cout << typeid(g).name() << std::endl;

            f.put("vertices", g.getNumOfVertices());
            f.put("edges", g.getNumOfEdges());

            auto removedEdges = Graph::FaspFast::superAlgorithm(g, c);
            EDGE_PROP_TYPE capacity = 0;
            for (const auto &e : removedEdges) {
                capacity += c.at(e);
            }

            f.put("superAlgorithmCapacity", capacity);
            f.put("exactCapacity", nv);
        }
    }

    f.save();
//    t.start_timer("1");
//    t.stop_timer();

}


std::string getFilenameOfBenchmarkFaspWithConstVE(int v, int e, int f) {
    return std::string("knownFaspSolution") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(f) +
           ".graph.txt";
}

void testIncreasingSizeOfFaspWithConstVE(int v, int e) {
    Timer<true, false> t(true);

    for (int i = 0 ; i < 1; ++i) {
        auto[g, c, fc] = Graph::Fasp::generateGraphWithKnownFasp<int, int, Graph::GraphMap>(139, 4, 133, 4, true, false);
        std::cout << g << std::endl;
        std::cout << g.getStrRepresentationOfGraph() << std::endl;
        std::cout << c << std::endl;
        std::cout << fc << std::endl;
        Graph::Fasp::GR(g, c);
        Graph::FaspFast::deltaFASP(g, c);
        Graph::FaspFast::randomFASP(g, c);

        Graph::IO::graphWithWeightsToFile("/tmp/" + getFilenameOfBenchmarkFaspWithConstVE(v, e, fc), g, c);
        auto [gn, cn] = Graph::IO::graphWithWeightsFromFile<int, Graph::GraphMap, int>("/tmp/" + getFilenameOfBenchmarkFaspWithConstVE(v, e, fc));
        std::cout << gn << std::endl;
        std::cout << gn.getStrRepresentationOfGraph() << std::endl;
        std::cout << cn << std::endl;
    }
    int numOfVertices = v;
    int numOfEdges = e;
    int reps = 1;
//
//    DataHdf5<double> f1("/tmp/out.h5");
//    DataHdf5<double> f2("/tmp/outr.h5");
//
//    for (int faspSize = 1; faspSize < 2 ; faspSize += 1) {
//        std::cout << "FASP SIZE: " << faspSize << std::endl;
//        for (int r = 0; r < reps; r++) {
//            auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(numOfVertices, faspSize, numOfEdges);
//            auto rc = Graph::Tools::getRandomWeights(g, 1, 10);
//
//            f1.put("vertices", g.getNumOfVertices());
//            f1.put("edges", g.getNumOfEdges());
//            f1.put("gr", Graph::Fasp::GR(g, c));
//            f1.put("delta", Graph::FaspFast::deltaFASP(g, c));
//            f1.put("random", Graph::FaspFast::randomFASP(g, c));
//            f1.put("exact", faspSize);
//
//            f2.put("vertices", g.getNumOfVertices());
//            f2.put("edges", g.getNumOfEdges());
//            f2.put("gr", Graph::Fasp::GR(g, rc));
//            f2.put("delta", Graph::FaspFast::deltaFASP(g, rc));
//            f2.put("random", Graph::FaspFast::randomFASP(g, rc));
//            f2.put("exact", faspSize);
//        }
//    }
//
//    f1.save();
//    f2.save();
}

void test() {
    using EDGE_PROP_TYPE = int;
    int v = 50;
    int f = 10;
    int e = 120;

    int maxEdges = v * (v - 1) / 2 + f;
    e = e > maxEdges ? maxEdges : e;
    auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<EDGE_PROP_TYPE, int, Graph::GraphMap>(v, f, e);
    std::cout << c << std::endl;
    std::cout << g << std::endl;
    std::cout << g.getStrRepresentationOfGraph() << std::endl;
    std::cout << "V/E/F: " << v << " " << e << " " << f << std::endl;

    Graph::Fasp::GR(g, c);
    Graph::FaspFast::deltaFASP(g, c);
    Graph::FaspFast::randomFASP(g, c);
}

void runGraph(int argc, char **argv) {


//    test();
//    testSuperAlgorithm();

    testIncreasingSizeOfFaspWithConstVE(50, 120);

    //Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/elo.txt");
//    Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/home/gonciarz/fasp/test/graphDelaun4/planar_73_213_336.graph.txt");
//    Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
    //Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
//    Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_340.graph.txt");
    //Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
//Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/krzysg/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
  //  Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
// Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/Users/krzysg/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_73_213_336.graph.txt");
//    Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphVector>("/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4/planar_12_30_31.graph.txt");
//Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>("/home/gonciarz/fasp/test/graphDelaun4/planar_73_213_336.graph.txt");
//    LOG(INFO) << "Edges with cycles: " << Graph::Tools::findEdgesWithCycles(gv) << std::endl;
//
//    std::cout << gv << std::endl;
//    std::cout << Graph::Tools::getRandomSubgraph(gv, 9) << std::endl;



    // Randomize a little weights
//    for (auto &w : weights) {
//        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5;
//        r /= 10000;
//        w.second += r;
//        std::cout << w << std::endl;
//    }

    Timer<true, false> t(true);


//    std::string dir = "/home/gonciarz/fasp/test/graphDelaun4/";
    std::string dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/GraphsStuff/test/graphDelaun4";

//
//    int limitCnt = 1;
//
//    int cnt = 0;
//    int cntGr = 0;
//    int cntDelta = 0;
//    int cntRand = 0;
//
//    DataHdf5<double> f("/tmp/out.h5");
//
//    t.start_timer("ALL_GRAPHS");
//    for (auto &file : Graph::IO::getFilesInDir(dir)) {
//        if (!Tools::endsWith(file, "graph.txt")) continue;
//
////        if (limitCnt-- == 0) break;
//
//        cnt++;
//        std::cout << cnt << " ======================= Processing [" << file << "] " << cnt << std::endl;
//
//        Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>(dir + "/" + file);

//        Graph::IO::graphToFile("/tmp/myGraph.txt", gv);
//        Tools::replace(file, "graph.txt", "fas.txt");
//        int solution = Graph::IO::solutionFromFile(dir + "/" + file);
//        f.put("exact", solution);
//        f.put("edges", gv.getNumOfEdges());
//        f.put("vertices", gv.getNumOfVertices());
//        std::cout << "SOLUTION=" << solution << std::endl;
//
//        auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
//
//        t.start_timer("GR");
//        auto c = Graph::Fasp::GR(gv, weights);
//        f.put("GR", c);
//        if (solution - c == 0) cntGr++;
//        t.stop_timer();
//        t.start_timer("DELTA");
//        try {
//            auto c = Graph::FaspFast::deltaFASP(gv, weights);
//            f.put("delta", c);
//            if (solution - c == 0) cntDelta++;
//        }
//        catch (const Graph::FaspFast::MinCutFailed &e) {
//            std::cout << "Failed to calculate FASP due to mincut non-terminating loop." << std::endl;
//        }
//        t.stop_timer();
//        t.start_timer("RAND");
//        try {
//            auto c = Graph::FaspFast::randomFASP(gv, weights);
//            f.put("random", c);
//            if (solution - c == 0) cntRand++;
//        }
//        catch (const Graph::FaspFast::MinCutFailed &e) {
//            std::cout << "Failed to calculate FASP due to mincut non-terminating loop." << std::endl;
//        }
//        t.stop_timer();
//
//        f.put("randomDelta", 0); //not existing but added to be compatible with other matlab files
//    }
//    t.stop_timer();
//    std::cout << cnt << std::endl;
//    std::cout << cntGr << std::endl;
//    std::cout << cntDelta << std::endl;
//    std::cout << cntRand << std::endl;
//
//    f.save();
//    std::cout << f << std::endl;


//    t.start_timer("GENERATE GRAPH");
//    {
//        auto g = Graph::Graph<int, Graph::GraphMap>();
//        for (int i = 1; i <= 4; ++i) g.addVertex(i);
//        g.addEdge({1, 2});
//        g.addEdge({2, 3});
//        g.addEdge({3, 1});
//        g.addEdge({3, 4});
//        Graph::Tools::stronglyConnectedComponents(g);
//    }
//
//
//    auto [g, c] = Graph::Fasp::generateGraphWithKnownFasp<int, int, Graph::GraphMap>(20, 2, 10, false);
//
//
//    int n = 30;
//    int cn = 5;
//    float p = (float)cn/n;
////    auto [g, c] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(n, p);
//    std::cout << g << std::endl;
//    Graph::Fasp::GR(g, c);
//    Graph::FaspFast::deltaFASP(g, c);
//    Graph::FaspFast::randomFASP(g, c);
//    t.stop_timer();
//
//    Graph::Tools::testLambda();


//    std::vector<Graph::Graph<int, Graph::GraphMap>> v;
//    int noOfReadFiles = 0;
//    for (auto &file : Graph::IO::getFilesInDir(dir)) {
//        if (!Tools::endsWith(file,"graph.txt")) continue;
//        v.emplace_back(Graph::IO::graphFromFile<int, Graph::GraphMap>(dir + "/" + file));
//        noOfReadFiles++;
//    }
//    std::cout << "Number of graphs: " << noOfReadFiles << std::endl;
//
//    t.start_timer("1");
//    for (int i = 0; i < 1; ++i)
//        for (const auto &g : v) {
//        auto weights = Graph::Ext::getEdgeProperties<int>(g, 1);
//            auto c = Graph::FaspFast::randomFASP(g, weights);
//
////            auto r = Graph::Tools::stronglyConnectedComponents(g);
////            int cnt = 0;
////            for (auto x : r) if (r.size() > 1) cnt++;
////            if (Graph::Tools::stronglyConnectedComponents(g).size() > 1) std::cout << cnt << "\n";
//        }
//    t.stop_timer();


}

