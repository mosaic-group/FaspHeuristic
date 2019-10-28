//
// Created by gonciarz on 2019-05-28.
//


#ifndef BENCH_RESULTS_DISTR
#define BENCH_RESULTS_DISTR

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphTools.h"
#include "graph/graphFasp.h"
#include "graph/graphFaspFastFinal.h"


std::string getFilenameOfBenchResultsDistr(int v, int e, int fmin, int fmax, int steps, int reps, bool logDistr) {
    return std::string("BenchResultsDistr") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log_" : "lin") +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchResultsDistr(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchResultsDistr. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkHeuristicsFaspWithConstVE(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);

    Graph::FaspFastFinal::PathHero<int> path(numOfVertices + 1); // maxId included

    for (auto &faspSize : faspValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int>(numOfVertices, faspSize, numOfEdges);

            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::FaspFastFinal::randomFASP(g, c);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());
            f1.put("gr", Graph::Fasp::GR(g, c).first);
            f1.put("random", capacity);
            f1.put("exact", faspSize);
            f1.put("saEdges", saEdgesCnt);
            f1.put("saRndEdges", saRndEdgesCnt);
            f1.put("redRndEdges", redRndEdgesCnt);
        }
    }

    f1.save();
}


#endif
