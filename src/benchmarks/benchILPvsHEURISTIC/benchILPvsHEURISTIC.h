#ifndef BENCHILPVSHEURISTIC_H
#define BENCHILPVSHEURISTIC_H

/**
 * Compares TIGHT-CUT and GR results to precomputed exact ILP solver results.
 *
 * This benchmark as a input takes input directory which should containg graph description '.al' files, solution '.mfas' files
 * and time of computing by exact solver '.timing' files.
 * @param inputDir - directory with input files
 * @param outputDir - directory for output h5 file wtih results
 * @param filename - name of h5 file (with extenstion)
 *
 * Example command line to execute:
 * benchmark benchILPvsHEURISTIC --inputDirectory /Users/gonciarz/ws/repo/FASP-benchmarks/data/random/ --outputDirectory . --fileName xyz.h5
 */
void benchILPvsHEURISTIC(const std::string &inputDir, const std::string &outputDir, const std::string &filename) {

    DataHdf5<double> f(outputDir + "/" + filename);

    auto checkIfExist = [](const std::vector<std::string> &files, const std::string &file) -> bool {return std::find(files.begin(), files.end(), file) != files.end();};

    double timeExact=0.0;
    double timeRandom=0.0;


    Timer<true, false> t{};
    auto files = Graph::IO::getFilesInDir(inputDir);
    std::sort(files.begin(), files.end());
    int cnt = 0;
    for (auto &graphFile : files) {
        // we are search for a graph files extenstion ".al"
        if (!Tools::endsWith(graphFile, ".al")) continue;

        // generate names of solution, timing and timeout indication files
        auto solutionFile{graphFile}; Tools::replace(solutionFile, ".al", ".mfas");
        auto timeoutFile{graphFile}; Tools::replace(timeoutFile, ".al", ".timeout");
        auto timingFile{graphFile}; Tools::replace(timingFile, ".al", ".timing");

        // if files are not needed files or there was timeout indication skip that graph
        if (!checkIfExist(files, solutionFile) || checkIfExist(files, timeoutFile) || !checkIfExist(files, timingFile)) {
            continue;
        }

        // Read graph and generate 'no weighted' weights (all are set to 1).
        Graph::Graph g = Graph::IO::graphFromFile<int16_t>(inputDir + "/" + graphFile);
        auto weights = Graph::Ext::getEdgeProperties<int16_t>(g, 1);

        // Read FASP solution and time of exact solver
        int solution = Graph::IO::solutionFromFile(inputDir + "/" + solutionFile);
        double timeExactOfGraph = Graph::IO::timingFromFile(inputDir + "/" + timingFile);

        std::cout << ++cnt << " ======================= Processing [" << graphFile << "] " << g << std::endl;

        // ---- EXACT ILP
        std::cout << "Exact solution=" << solution << " Time=" << timeExactOfGraph << std::endl;
        timeExact += timeExactOfGraph;

        // ---- TIGHT-CUT
        t.start_timer("-------- tightCut");
        auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut<false, false>(g, weights);
        auto thisStep = t.stop_timer();
        timeRandom+=thisStep;
        LOG(DEBUG) << "SA / SA RND / RED RND edges: " << saEdgesCnt << " / " << saRndEdgesCnt << " / " << redRndEdgesCnt;
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;

        // ---- GR
        t.start_timer("gr");
        auto gr = Graph::FaspTools::GR(g, weights);
        auto grTime = t.stop_timer();

        // Save all data
        f.put("gr", gr.first);
        f.put("grTime", grTime);
        f.put("random", capacity);
        f.put("randomTime", thisStep);
        f.put("exact", solution);
        f.put("exactTime", timeExactOfGraph);
        f.put("vertices", g.getNumOfVertices());
        f.put("edges", g.getNumOfEdges());
        f.put("fileName", graphFile);
        f.put("saEdges", saEdgesCnt);
        f.put("saRndEdges", saRndEdgesCnt);
        f.put("redRndEdges", redRndEdgesCnt);

        std::cout << "TIME exact/random: " << timeExact << "/" << timeRandom << std::endl;
    }

    f.save();
}


#endif
