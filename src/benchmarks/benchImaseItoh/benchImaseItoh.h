//
// Created by gonciarz on 10/30/19.
//

#ifndef FASPHEURISTIC_BENCHIMASEITOH_H
#define FASPHEURISTIC_BENCHIMASEITOH_H


void benchImaseItoh(const std::string &outputDir, const std::string &inputDir) {
    DataHdf5<double> f(outputDir + "/imaseItoh.h5");
    std::string dir = inputDir;

    Timer<true, false> t{};
    int cnt = 0;


    auto checkIfExist = [](const std::vector<std::string> &files, const std::string &file) -> bool {return std::find(files.begin(), files.end(), file) != files.end();};

    double timeRandom=0.0;

    auto files = Graph::IO::getFilesInDir(std::string(dir));

    std::sort(files.begin(), files.end());
    std::cout << files.size() << std::endl;
    for (auto &graphFile : files) {
        if (!Tools::endsWith(graphFile, ".al")) continue;


        auto solutionFile{graphFile}; Tools::replace(solutionFile, ".al", ".mfsize");

        if (!checkIfExist(files, solutionFile)) {
            continue;
        }

        Graph::Graph g = Graph::IO::graphFromFile<int16_t>(dir + "/" + graphFile);
        auto c = Graph::Ext::getEdgeProperties<int16_t>(g, 1);
        Graph::IO::graphToFile("/tmp/myGraph.txt", g);
        double solution = Graph::IO::timingFromFile(dir + "/" + solutionFile);

        std::cout << ++cnt << " ======================= Processing [" << graphFile << "] " << g << std::endl;
        std::cout << "Exact solution=" << solution << std::endl;

        t.start_timer("--------RANDOM new");
        auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut(g, c);
        auto thisStep = t.stop_timer();
        LOG(DEBUG) << "SA / SA RND / RED RND edges: " << saEdgesCnt << " / " << saRndEdgesCnt << " / " << redRndEdgesCnt;
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        timeRandom+=thisStep;

        t.start_timer("gr");
        auto gr = Graph::FaspTools::GR(g, c);
        auto grTime = t.stop_timer();

        f.put("gr", gr.first);
        f.put("grTime", grTime);
        f.put("random", capacity);
        f.put("randomTime", thisStep);
        f.put("exact", solution);
        f.put("exactTime", 1); // fake it
        f.put("vertices", g.getNumOfVertices());
        f.put("edges", g.getNumOfEdges());
        f.put("fileName", graphFile);
        f.put("saEdges", saEdgesCnt);
        f.put("saRndEdges", saRndEdgesCnt);
        f.put("redRndEdges", redRndEdgesCnt);
    }

    f.save();
}


#endif //FASPHEURISTIC_BENCHIMASEITOH_H
