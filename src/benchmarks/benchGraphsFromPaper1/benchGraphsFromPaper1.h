//
// Created by gonciarz on 2019-05-28.
//


#ifndef GRAPHS_FROM_PAPER_1
#define GRAPHS_FROM_PAPER_1

#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"
#include "hdf5/dataHdf5.h"
#include "graph/graph.h"
#include "graph/graphTools.h"
#include "graph/graphFasp.h"
#include "graph/graphFaspFast.h"
#include <regex>

/*
 * Benchmarks of graphs from this paper:
 *
 * An exact method for the minimum feedback arc set problem
 */

std::string getFilenameOfBenchGraphsFromPaper1() {
    return std::string("BenchGraphsFromPaper1.h5");
}


void benchGraphsFromPaper1(const std::string &outputDir, const std::string &inputDir) {

    LOG(TRACE) << "Running benchGraphsFromPaper1";

    std::regex rx(R"(.*opt=([0-9]+).*)");
    std::smatch match;

    auto outputFile = outputDir + "/" + getFilenameOfBenchGraphsFromPaper1();
    DataHdf5<double> f1(outputFile, /* create output file (dummy run) */ (outputDir == "" ? true : false));

    int numOfCorrect = 0;
    int numOfGraphs = 0;
    int largestError = 0;
    for (auto &file : Graph::IO::getFilesInDir(inputDir)) {
        if (!Tools::endsWith(file, "graph.txt")) continue;

        LOG(INFO) << " ======================= Processing [" << file << "] ";

        std::regex_match(file, match, rx);
        int correct = 0;
        if (match.size() != 2) LOG(ERROR) << "Cannot extract correct solution from file name! [" << file << "]";
        else {
            std::ssub_match sub_match = match[1];
            std::string num = sub_match.str();
            correct = std::stoi(num);
        }

        Graph::Graph gv = Graph::IO::graphFromFile<int, Graph::GraphMap>(inputDir + "/" + file);
        auto c = Graph::Ext::getEdgeProperties<int>(gv, 1);

        int faspCapacity = Graph::FaspFast2::randomFASP_sequential(gv, c);
        f1.put("random", faspCapacity);
        f1.put("exact", correct);

        if (faspCapacity == correct) ++numOfCorrect;
        else {
            LOG(ERROR) << "Not exact soltion found!";
            if (faspCapacity - correct > largestError) largestError = faspCapacity - correct;
        }
        ++numOfGraphs;
    }

    f1.save();

    LOG(INFO) << "Results - #graphs exactly solved / #graphs: " << numOfCorrect << "/" << numOfGraphs;
    LOG(INFO) << "Largest error: " << largestError;
}


#endif
