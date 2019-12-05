
#ifndef FASPHEURISTIC_BENCHBRUIJN_H
#define FASPHEURISTIC_BENCHBRUIJN_H

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"

#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphFaspTools.h"
#include "graph/graphFasp.h"
#include "graph/graphIO.h"


void benchBruijnGraphs(const std::string &outputDir) {


    LOG(TRACE) << "Running benchBruijnGraphs";

    std::string dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/Fasp/src/benchmarks/benchErdosRenyiGraphs/data";
    std::string file = "erdos_renyi_n_100_c_5_seed_1.edges";

    std::cout << " ======================= Processing [" << file << "] " << std::endl;

        std::string inF = "/Users/gonciarz/ws/FaspRemote/src/benchmarks/benchErdosRenyiGraphs/data/erdos_renyi_n_75_c_9_seed_43.edges";
        Graph::Graph gg = Graph::IO::graphFromFile<int>(inF);
        std::cout << "INF: " << gg << std::endl;

        Graph::IO::graphToFile<int>("/tmp/graph.txt", gg);
        std::string outF = "/tmp/graph.txt";
        Graph::Graph gf = Graph::IO::graphFromFile<int>(outF);
        std::cout << gf << std::endl;

        auto g{gf};
        auto g2{gf};

        auto c = Graph::Ext::getEdgeProperties<int>(g, 1);
        auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);

        Timer<true, false> t("");

        t.start_timer("final");
    Graph::Fasp::tightCut(g, c);
        t.stop_timer();
}


#endif
