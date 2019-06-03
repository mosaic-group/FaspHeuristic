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
#include "graph/graphTools.h"
#include "graph/graphFasp.h"
#include "graph/graphFaspFast.h"


std::string getFilenameOfBenchmarkFaspWithConstVE(int v, int e, int fmin, int fmax, int steps, int reps) {
    return std::string("SuperAlgorithmConstWeightVarFaspConstVE") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchSuperAlgorithmConstWeightVarFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps) {

    LOG(TRACE) << "Running benchSuperAlgorithm. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkFaspWithConstVE(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto faspValues =  Tools::linspace(minFasp, maxFasp, numOfSteps);

    Graph::FaspFast::PathHero<int> path(numOfVertices + 1); // maxId included

    for (auto &faspSize : faspValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(numOfVertices, faspSize, numOfEdges);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());
            f1.put("sa", Graph::FaspFast::superAlgorithm(g, c, path).size());
//            f1.put("gr", Graph::Fasp::GR(g, c));
//            f1.put("delta", Graph::FaspFast::deltaFASP(g, c));
//            f1.put("random", Graph::FaspFast::randomFASP(g, c));
            f1.put("exact", faspSize);
        }
    }

    f1.save();
}


#endif //FASPHEURISTIC_BENCHSUPERALGORITHM_H