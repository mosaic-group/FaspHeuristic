//
// Created by gonciarz on 2019-05-28.
//


#ifndef FASPHEURISTIC_BENCHSUPERALGORITHM_H
#define FASPHEURISTIC_BENCHSUPERALGORITHM_H


#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphFaspTools.h"
#include "graph/graphFasp.h"
#include <string>

std::string getFilenameOfBenchmarkSAFaspWithConstVE(int v, int e, int fmin, int fmax, int steps, int reps, bool logDistr) {
    return std::string("SuperAlgorithmConstWeightVarFaspConstVE") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log_" : "_lin_") +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchSuperAlgorithmConstWeightVarFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchSuperAlgorithm. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkSAFaspWithConstVE(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);

    Graph::FaspFastFinal::PathHero<int> path(numOfVertices + 1); // maxId included

    for (auto &faspSize : faspValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::FaspTools::generateGraphWithKnownFaspAndSameWeights<int, int>(numOfVertices, faspSize, numOfEdges);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());
            auto [edgesSA, _, edgesGR] = superAlgorithmBlue(g, c, path, false, false);
            f1.put("sa", edgesSA.size());
            f1.put("exact", faspSize);
        }
    }

    f1.save();
}


#endif //FASPHEURISTIC_BENCHSUPERALGORITHM_H
