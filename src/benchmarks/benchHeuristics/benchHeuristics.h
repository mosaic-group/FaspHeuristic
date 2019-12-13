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
#include <sstream>

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
                      "_s_" + std::to_string(numOfSteps) + (logDistribution ? "_log" : "_lin") +
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
            std::string graphFile{"graph-" + Tools::convertToStrWithLeadingZeros(graphIndex) + "-" + std::to_string(numOfVertices) + "-" + std::to_string(numOfEdges) + "-" + std::to_string(faspSize) + ".al"};
            Graph::IO::graphToFile(outputDir + "/" + graphFile, g);

            resultFile.put("fileName", graphFile);
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

/**
 * Benchmark for testing heuristic using graphs with const density (defined as a |e|/|v|), const fasp size and changin in provided
 * range number of vertices.
 * @param outputDir - where to save output files
 * @param minNumOfVertices
 * @param maxNumOfVertices
 * @param density - density of graph
 * @param faspSize - fasp size
 * @param numOfSteps - how many data points for provided vertices range
 * @param numOfReps - how many repetions for each data point
 * @param logDistribution true to use log distribution of vertices
 *
 * Example cmd line:
 * benchmark benchVarVConstDensityFaspWeight --vmin 20 --vmax 30 --d 3 --f 4 --steps 4 --reps 5 --outputDirectory result2
 */
void benchVarVConstDensityFaspWeight(const std::string &outputDir, int minNumOfVertices, int maxNumOfVertices, double density, int faspSize, int numOfSteps, int numOfReps, bool logDistribution) {
    LOG(TRACE) << "Running benchVarVConstDensityFaspWeight. Params: vmin=" << minNumOfVertices << " vmax=" << maxNumOfVertices << " d=" << density << " f=" << faspSize << " steps=" << numOfSteps << " reps=" << numOfReps;

    // get density string with limited precision and underscore instead of dot (4.3751 becomes "4_37")
    std::ostringstream tmpStr; tmpStr << std::fixed << std::setprecision(2) << density;
    std::string densityStr = tmpStr.str();
    std::replace(densityStr.begin(), densityStr.end(), '.', '_');

    auto outputFile = outputDir + "/" +
                      "BenchVarVConstDensityFaspWeight" +
                      "_v_" + std::to_string(minNumOfVertices) + "-" + std::to_string(maxNumOfVertices) +
                      "_d_" + densityStr +
                      "_f_" + std::to_string(faspSize) +
                      "_s_" + std::to_string(numOfSteps) + (logDistribution ? "_log" : "_lin") +
                      "_r_" + std::to_string(numOfReps) +
                      ".h5";
    DataHdf5<double> resultFile(outputFile);

    auto vValues =  logDistribution ? Tools::logspace(minNumOfVertices, maxNumOfVertices, numOfSteps) : Tools::linspace(minNumOfVertices, maxNumOfVertices, numOfSteps);

    Timer<true, false> t("benchVarVConstDensityFaspWeight");
    int graphIndex = 0;

    for (auto &cv : vValues) {
        for (int r = 0; r < numOfReps; ++r) {
            LOG(TRACE) << "--- Progress --- V size=" << cv  << "/" << maxNumOfVertices << " Reps=" << r + 1 << "/" << numOfReps << "";

            // generate benchmarked graph
            auto[g, c] = Graph::FaspTools::generateGraphWithKnownFaspAndSameWeights<int, int>(cv, faspSize, cv * density);

            // --- run GR heuristic
            t.start_timer("gr");
            auto gr = Graph::FaspTools::GR(g, c);
            auto grTime = t.stop_timer();

            // --- run TIGHT-CUT heuristic
            t.start_timer("tightCut");
            auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut(g, c);
            auto randomTime = t.stop_timer();

            // save generated graph and results
            graphIndex++;
            std::string graphFile{"graph-" + Tools::convertToStrWithLeadingZeros(graphIndex) + "-" + std::to_string(g.getNumOfVertices()) + "-" + std::to_string(g.getNumOfEdges()) + "-" + std::to_string(faspSize) + ".al"};
            Graph::IO::graphToFile(outputDir + "/" + graphFile, g);

            resultFile.put("fileName", graphFile);
            resultFile.put("vertices", g.getNumOfVertices());
            resultFile.put("edges", g.getNumOfEdges());
            resultFile.put("gr", gr.first);
            resultFile.put("grTime", grTime);
            resultFile.put("random", capacity);
            resultFile.put("randomTime", randomTime);
            resultFile.put("exact", faspSize);
            resultFile.put("exactTime", 1); // fake it for matlab
            resultFile.put("saEdges", saEdgesCnt);
            resultFile.put("saRndEdges", saRndEdgesCnt);
            resultFile.put("redRndEdges", redRndEdgesCnt);
        }
    }
    resultFile.save();
}

/**
 * Benchmark for testing heuristics using graphs with  variable weights, const number of vertices and edges.
 * Changing fasp size according provided data and then generating random values for cycles in range to 1-10. So
 * final capacity may vary lot.
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
 * benchmark benchVarWeightFaspConstVE --v 20 --e 100 --fmin 1 --fmax 11 --steps 6 --reps 6 --outputDirectory result3
 */
void benchVarWeightFaspConstVE(const std::string &outputDir, int numOfVertices, int numOfEdges, int minFasp, int maxFasp, int numOfSteps, int numOfReps, bool logDistribution) {
    LOG(TRACE) << "Running benchVarFaspConstWeightVE. Params: v=" << numOfVertices << " e=" << numOfEdges << " fmin=" << minFasp << " fmax=" << maxFasp << " steps=" << numOfSteps << " reps=" << numOfReps;

    auto outputFile = outputDir + "/" +
                      "BenchVarWeightFaspConstVE" +
                      "_v_" + std::to_string(numOfVertices) +
                      "_e_" + std::to_string(numOfEdges) +
                      "_f_" + std::to_string(minFasp) + "-" + std::to_string(maxFasp) +
                      "_s_" + std::to_string(numOfSteps) + (logDistribution ? "_log" : "_lin") +
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
            auto[g, c, graphcapacity] = Graph::FaspTools::generateGraphWithKnownFasp<int, int>(numOfVertices, faspSize, numOfEdges, 10, true, false);

            // --- run GR heuristic
            t.start_timer("gr");
            auto gr = Graph::FaspTools::GR(g, c);
            auto grTime = t.stop_timer();

            // --- run TIGHT-CUT heuristic
            t.start_timer("tightCut");
            auto [tightCutCapacity, tightCutEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut(g, c);
            auto tightCutTime = t.stop_timer();

            // save generated graph and results
            graphIndex++;
            std::string graphFile{"graph-" + Tools::convertToStrWithLeadingZeros(graphIndex) + "-" + std::to_string(numOfVertices) + "-" + std::to_string(numOfEdges)+"-"+std::to_string(faspSize)+".al"};
            Graph::IO::graphWithWeightsToFile(outputDir + "/" + graphFile, g, c);

            resultFile.put("fileName", graphFile);
            resultFile.put("vertices", g.getNumOfVertices());
            resultFile.put("edges", g.getNumOfEdges());
            resultFile.put("gr", gr.first);
            resultFile.put("grTime", grTime);
            resultFile.put("random", tightCutCapacity);
            resultFile.put("randomTime", tightCutTime);
            resultFile.put("exact", graphcapacity);
            resultFile.put("exactTime", 1); // fake it for matlab
            resultFile.put("saEdges", saEdgesCnt);
            resultFile.put("saRndEdges", saRndEdgesCnt);
            resultFile.put("redRndEdges", redRndEdgesCnt);
        }
    }

    resultFile.save();

}

#endif
