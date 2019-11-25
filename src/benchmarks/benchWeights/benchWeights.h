#ifndef BENCH_WEIGHTS_H
#define BENCH_WEIGHTS_H

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphFaspTools.h"
#include "graph/graphFasp.h"
#include "graph/graphIO.h"

static std::string getFilenameOfBenchmarkTimingVarWeightVarFaspConstVE(int v, int e, int fmin, int fmax, int steps, int reps, bool logDistr) {
    return std::string("TimingVarWeightVarFaspConstVE") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log_" : "_lin_") +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchTimingVarWeightVarFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchTimingConstWeightVarFaspConstVE. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkTimingVarWeightVarFaspConstVE(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);


    Timer<true, false> t("benchTiming");

    int fidx = 0;
    for (auto &faspSize : faspValues) {
        Timer<true, false> t("");
        for (int r = 0; r < numOfReps; ++r) {
            fidx++;
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c, graphcapacity] = Graph::FaspTools::generateGraphWithKnownFasp<int, int>(numOfVertices, faspSize, numOfEdges, 10, true, false);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());

            t.start_timer("gr");
            auto gr = Graph::FaspTools::GR(g, c);
            auto grTime = t.stop_timer();

            t.start_timer("random");
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::FaspFastFinal::randomFASP<int, int, true>(g, c);
            auto randomTime = t.stop_timer();

            Graph::IO::graphWithWeightsToFile(outputDir + "/" + "xl-" + Tools::convertToStrWithLeadingZeros(fidx) + "-" + std::to_string(numOfVertices) + "-" + std::to_string(numOfEdges)+"-"+std::to_string(faspSize)+".al", g, c);

            f1.put("gr", gr.first);
            f1.put("grEdges", gr.second.size());
            f1.put("grTime", grTime);
            f1.put("random", capacity);
            f1.put("randomEdges", removedEdges.size());
            f1.put("randomTime", randomTime);
            f1.put("exact", graphcapacity);
            f1.put("exactTime", 1); // fake it
            f1.put("saEdges", saEdgesCnt);
            f1.put("saRndEdges", saRndEdgesCnt);
            f1.put("redRndEdges", redRndEdgesCnt);
        }
    }
    f1.save();
}

#endif
