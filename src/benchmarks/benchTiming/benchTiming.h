#ifndef BENCH_TIMING_H
#define BENCH_TIMING_H

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphTools.h"
#include "graph/graphFasp.h"
#include "graph/graphFaspFastFinal.h"
#include "graph/graphIO.h"

static std::string getFilenameOfBenchmarkTiming(int v, int e, int fmin, int fmax, int steps, int reps, bool logDistr) {
    return std::string("NewTimingConstWeightVarFaspConstVE") +
           "_v_" + std::to_string(v) +
           "_e_" + std::to_string(e) +
           "_f_" + std::to_string(fmin) + "-" + std::to_string(fmax) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log_" : "_lin_") +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchTimingConstWeightVarFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchTimingConstWeightVarFaspConstVE. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkTiming(numOfVertices, numOfEdges, minFasp, maxFasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);


    Timer<true, false> t("benchTiming");

    int fidx = 0;
    for (auto &faspSize : faspValues) {
        Timer<true, false> t("");
        for (int r = 0; r < numOfReps; ++r) {
            fidx++;
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int>(numOfVertices, faspSize, numOfEdges);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());

            t.start_timer("gr");
            auto gr = Graph::Fasp::GR(g, c);
            auto grTime = t.stop_timer();

//            t.start_timer("delta");
//            auto delta = Graph::FaspFastFinal::  deltaFASP(g, c);
//            auto deltaTime = t.stop_timer();

            t.start_timer("random");
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::FaspFastFinal::randomFASP(g, c);
            auto randomTime = t.stop_timer();

            Graph::IO::graphToFile(outputDir + "/" + "xl-" + Tools::convertToStrWithLeadingZeros(fidx) + "-" + std::to_string(numOfVertices) + "-" + std::to_string(numOfEdges)+"-"+std::to_string(faspSize)+".al", g);

            f1.put("gr", gr.first);
            f1.put("grTime", grTime);
            f1.put("random", capacity);
            f1.put("randomTime", randomTime);
            f1.put("exact", faspSize);
            f1.put("exactTime", 1); // fake it
            f1.put("saEdges", saEdgesCnt);
            f1.put("saRndEdges", saRndEdgesCnt);
            f1.put("redRndEdges", redRndEdgesCnt);
        }
    }
    f1.save();
}


std::string getFilenameOfBenchmarkConstDensityAndFaspGrowingsize(int vmin, int vmax, double d, int fasp, int steps, int reps, bool logDistr) {
    std::ostringstream streamObj3;
    streamObj3 << std::fixed;
    streamObj3 << std::setprecision(2);
    streamObj3 << d;
    std::string dens = streamObj3.str();

    std::replace( dens.begin(), dens.end(), '.', '_');
    return std::string("TimingConstWeightVarFaspConstVE") +
           "_v_" + std::to_string(vmin) + "-" + std::to_string(vmax) +
           "_d_" + dens +
           "_f_" + std::to_string(fasp) +
           "_s_" + std::to_string(steps) + (logDistr ? "_log_" : "lin") +
           "_r_" + std::to_string(reps) +
           ".h5";
}


void benchTimingConstDensityAndFaspGrowingsize(const std::string &outputDir, int minNumOfVertices, int maxNumOfVertices, double density, int fasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchTimingConstDensityAndFaspGrowingsize. Params: vmin=" << minNumOfVertices << " vmax=" << maxNumOfVertices << " d=" << density << " f=" << fasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    auto outputFile = outputDir + "/" + getFilenameOfBenchmarkConstDensityAndFaspGrowingsize(minNumOfVertices, maxNumOfVertices, density, fasp, numOfSteps, numOfReps, logDistribution);
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto vValues =  logDistribution ? Tools::logspace(minNumOfVertices, maxNumOfVertices, numOfSteps) : Tools::linspace(minNumOfVertices, maxNumOfVertices, numOfSteps);

    Timer<true, false> t("benchTiming");

    for (auto &cv : vValues) {
        Timer<true, false> t("");
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- V size=" << cv  << "/" << maxNumOfVertices << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int>(cv, fasp, cv * density);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());

            t.start_timer("gr");
            auto gr = Graph::Fasp::GR(g, c);
            auto grTime = t.stop_timer();

            t.start_timer("random");
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::FaspFastFinal::randomFASP(g, c);
            auto randomTime = t.stop_timer();

            f1.put("gr", gr.first);
            f1.put("grTime", grTime);
            f1.put("random", capacity);
            f1.put("randomTime", randomTime);
            f1.put("exact", fasp);
            f1.put("exactTime", 1);
            f1.put("saEdges", saEdgesCnt);
            f1.put("saRndEdges", saRndEdgesCnt);
            f1.put("redRndEdges", redRndEdgesCnt);
        }
    }
    f1.save();}



#endif
