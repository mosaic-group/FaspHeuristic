#ifndef BENCHHEURISTICS_H
#define BENCHHEURISTICS_H

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphFaspTools.h"
#include "graph/graphFasp.h"
#include "graph/graphIO.h"


/**
 * Benchmark for testing heuristics using graphs with  const weights, const number of vertices and edges, and changing fasp size according provided data.
 * FASP range can increase lineary or using logscale.
 * @param outputDir - where to save output files (h5 + generated graphs)
 * @param numOfVertices
 * @param numOfEdges
 * @param minFasp
 * @param maxFasp
 * @param numOfSteps - number of data points that FASP range is divided into, minFasp and maxFasp are always included
 * @param numOfReps - number of repetitions for given FASP size
 * @param logDistribution - if true log distribution is used
 *
 * Example cmd line:
 * benchmark benchVarFaspConstWeightVE --v 20 --e 100 --fmin 1 --fmax 11 --steps 6 --reps 6 --outputDirectory result1
 */
void benchVarFaspConstWeightVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {
    LOG(TRACE) << "Running benchVarFaspConstWeightVE. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;

    auto outputFile = outputDir + "/" +
                      "BenchVarFaspConstWeightVE" +
                      "_v_" + std::to_string(numOfVertices) +
                      "_e_" + std::to_string(numOfEdges) +
                      "_f_" + std::to_string(minFasp) + "-" + std::to_string(maxFasp) +
                      "_s_" + std::to_string(numOfSteps) + (logDistribution ? "_log_" : "_lin_") +
                      "_r_" + std::to_string(numOfReps) +
                      ".h5";
    DataHdf5<double> resultFile(outputFile);

    auto faspValues =  logDistribution ? Tools::logspace(minFasp, maxFasp, numOfSteps) : Tools::linspace(minFasp, maxFasp, numOfSteps);

    Timer<true, false> t("benchVarFaspConstWeightVE");
    int graphIndex = 0;

    for (auto &faspSize : faspValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- Fasp size=" << faspSize  << "/" << maxFasp << " Reps=" << r + 1 << "/" << numOfReps << "";

            // --- generate benchmarked graph
            auto[g, c] = Graph::FaspTools::generateGraphWithKnownFaspAndSameWeights<int, int>(numOfVertices, faspSize, numOfEdges);

            // --- run GR heuristic
            t.start_timer("gr");
            auto gr = Graph::FaspTools::GR(g, c);
            auto grTime = t.stop_timer();

            // --- run TIGHT-CUT heuristic
            t.start_timer("tightCut");
            auto [tightCutCapacity, tightCutEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut(g);
            auto tightCutTime = t.stop_timer();

            // save generated graph and results
            graphIndex++;
            Graph::IO::graphToFile(outputDir + "/" + "graph-" + Tools::convertToStrWithLeadingZeros(graphIndex) + "-" + std::to_string(numOfVertices) + "-" + std::to_string(numOfEdges) + "-" + std::to_string(faspSize) + ".al", g);

            resultFile.put("vertices", g.getNumOfVertices());
            resultFile.put("edges", g.getNumOfEdges());
            resultFile.put("vertices", g.getNumOfVertices());
            resultFile.put("edges", g.getNumOfEdges());
            resultFile.put("gr", gr.first);
            resultFile.put("grTime", grTime);
            resultFile.put("random", tightCutCapacity);
            resultFile.put("randomTime", tightCutTime);
            resultFile.put("exact", faspSize);
            resultFile.put("exactTime", 1); // fake it for matlab
            resultFile.put("saEdges", saEdgesCnt);
            resultFile.put("saRndEdges", saRndEdgesCnt);
            resultFile.put("redRndEdges", redRndEdgesCnt);
        }
    }

    resultFile.save();
}


void benchTimingConstDensityAndFaspGrowingsize(const std::string &outputDir, int minNumOfVertices, int maxNumOfVertices, double density, int fasp, int numOfSteps, int numOfReps, bool logDistribution) {

    LOG(TRACE) << "Running benchTimingConstDensityAndFaspGrowingsize. Params: vmin=" << minNumOfVertices << " vmax=" << maxNumOfVertices << " d=" << density << " f=" << fasp << " steps=" << numOfSteps << " reps=" << numOfReps;


    std::string result;
    std::ostringstream streamObj3;
    streamObj3 << std::fixed;
    streamObj3 << std::setprecision(2);
    streamObj3 << density;
    std::string dens = streamObj3.str();

    std::replace(dens.begin(), dens.end(), '.', '_');

    auto outputFile = outputDir + "/" +
                      "TimingConstWeightVarFaspConstVE" +
                      "_v_" + std::to_string(minNumOfVertices) + "-" + std::to_string(maxNumOfVertices) +
                      "_d_" + dens +
                      "_f_" + std::to_string(fasp) +
                      "_s_" + std::to_string(numOfSteps) + (logDistribution ? "_log_" : "lin") +
                      "_r_" + std::to_string(numOfReps) +
                      ".h5";

    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    auto vValues =  logDistribution ? Tools::logspace(minNumOfVertices, maxNumOfVertices, numOfSteps) : Tools::linspace(minNumOfVertices, maxNumOfVertices, numOfSteps);

    Timer<true, false> t("benchTiming");

    for (auto &cv : vValues) {
        Timer<true, false> t("");
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- V size=" << cv  << "/" << maxNumOfVertices << " Reps=" << r + 1 << "/" << numOfReps << "";
            auto[g, c] = Graph::FaspTools::generateGraphWithKnownFaspAndSameWeights<int, int>(cv, fasp, cv * density);

            f1.put("vertices", g.getNumOfVertices());
            f1.put("edges", g.getNumOfEdges());

            t.start_timer("gr");
            auto gr = Graph::FaspTools::GR(g, c);
            auto grTime = t.stop_timer();

            t.start_timer("random");
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut(g, c);
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
    f1.save();
}

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
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut<true, false, int, int>(g, c);
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
