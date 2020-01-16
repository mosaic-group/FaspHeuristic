#ifndef BENCHISOCUT_H
#define BENCHISOCUT_H


#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "../../../include/FaspTightCut/graph.h"
#include "graph/graphFaspTools.h"
#include "../../../include/FaspTightCut/graphFasp.h"
#include <string>


/**
 * This benchmark generates graphs with const V/E and range of known FASP sizes. Runs ISO-CUT on these graphs and saves all data.
 *
 * Example command line:
 * benchmark benchIsoCutConstWeightVarFaspConstVE -d . --v 50 --e 500 --fmin 1 --fmax 10 --steps 6 --reps 3 --logScale
 */


std::string getFilenameOfBenchmarkIsoCutFaspWithConstVE(int v, int e, int fmin, int fmax, int steps, int reps, bool logDistr) {
    return std::string("IsoCutConstWeightVarFaspConstVE") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log" : "_lin") +
           "_r_" + std::to_string(reps) +
           ".h5";
}

void benchIsoCutConstWeightVarFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchIsoCut. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkIsoCutFaspWithConstVE(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> file(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    // Generate steps for FASP size depending on chosen distribution (linear or logarithmic).
    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);

    Graph::Fasp::GraphSpeedUtils<int> tools(numOfVertices + 1); // maxId included

    for (auto &faspSize : faspValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::FaspTools::generateGraphWithKnownFaspAndSameWeights<int, int>(numOfVertices, faspSize, numOfEdges);

            // Run ISO-CUT on generated graph
            auto [cutEdges, dummy1, dummy2] = isoCut(g, tools, false);

            file.put("vertices", g.getNumOfVertices());
            file.put("edges", g.getNumOfEdges());
            file.put("sa", cutEdges.size());
            file.put("exact", faspSize);
        }
    }

    file.save();
}


#endif
