#ifndef GRAPHS_FROM_PAPER_1
#define GRAPHS_FROM_PAPER_1

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphFaspTools.h"
#include "graph/graphFasp.h"
#include "graph/graphIO.h"
#include <regex>

/*
 * Benchmarks of graphs from this paper:
 * "An exact method for the minimum feedback arc set problem" Ali Baharev et al.
 */

std::string getFilenameOfBenchGraphsFromPaper1() {
    return std::string("BenchGraphsFromPaper1.h5");
}


void benchGraphsFromPaper1(const std::string &outputDir, const std::string &inputDir) {
    LOG(TRACE) << "Running benchGraphsFromPaper1";

    // create file for results
    auto outputFile = outputDir + "/" + getFilenameOfBenchGraphsFromPaper1();
    DataHdf5<double> resultsFile(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    // Simple statistics
    int numOfCorrect = 0;
    int numOfGraphs = 0;
    int largestError = 0;

    for (auto &file : Graph::IO::getFilesInDir(inputDir)) {
        if (!Tools::endsWith(file, "graph.txt")) continue;

        LOG(INFO) << " ======================= Processing [" << file << "] ";

        // Extract exact FASP size from filename
        std::regex rx(R"(.*opt=([0-9]+).*)");
        std::smatch match;
        std::regex_match(file, match, rx);
        int exactSolution = 0;
        if (match.size() != 2) {
            LOG(ERROR) << "Cannot extract exact solution from file name! [" << file << "]";
            continue;
        }
        else {
            std::ssub_match sub_match = match[1];
            std::string num = sub_match.str();
            exactSolution = std::stoi(num);
        }

        // Read graph from file and run ISO-CUT
        Graph::Graph gv = Graph::IO::graphFromFile<int>(inputDir + "/" + file);
        auto [capacity, removedEdges, tightCutEdgesCnt, tightCutRndEdgesCnt, redEdgesRndEdgesCnt] = Graph::Fasp::tightCut(gv);

        // Save data
        resultsFile.put("random", capacity);
        resultsFile.put("exact", exactSolution);
        resultsFile.put("saEdges", tightCutEdgesCnt);
        resultsFile.put("saRndEdges", tightCutRndEdgesCnt);
        resultsFile.put("redRndEdges", redEdgesRndEdgesCnt);

        // Some statistics...
        if (capacity == exactSolution) ++numOfCorrect;
        else {
            LOG(ERROR) << "Not exact soltion found!";
            if (capacity - exactSolution > largestError) largestError = capacity - exactSolution;
        }
        ++numOfGraphs;
    }

    resultsFile.save();

    LOG(INFO) << "Results - #graphs exactly solved / #graphs: " << numOfCorrect << "/" << numOfGraphs;
    LOG(INFO) << "Largest error: " << largestError;
}


#endif
