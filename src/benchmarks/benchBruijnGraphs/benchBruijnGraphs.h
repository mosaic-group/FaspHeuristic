
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



void benchBruijnGraphs(const std::string &outputDir) {

    LOG(TRACE) << "Running benchBruijnGraphs";

    std::string dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/Fasp/src/benchmarks/benchErdosRenyiGraphs/data";
    std::string file = "erdos_renyi_n_35_c_7_seed_18.edges";

    std::cout << " ======================= Processing [" << file << "] " << std::endl;


    Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>(dir + "/" + file);

    auto c = Graph::Ext::getEdgeProperties<int>(gv, 1);

    Graph::FaspFast::deltaFASP(gv, c);
//    Graph::FaspFast::randomFASP(gv, c);
}


#endif
