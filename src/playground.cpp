#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

#include "graph/graph.h"
#include "graph/graphIO.h"
#include "hdf5/dataHdf5.h"
#include "graph/graphFasp.h"

#include <string>
using std::size_t;

void configureLogger() {
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
    defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Info, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Error, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Fatal, el::ConfigurationType::Enabled, "true");
    defaultConf.set(el::Level::Global, el::ConfigurationType::SubsecondPrecision, "1");
    constexpr const char *format = "[%datetime, %fbase:%line] %msg";
    defaultConf.set(el::Level::Global, el::ConfigurationType::Format, format);
    el::Loggers::reconfigureLogger("default", defaultConf);
}

// Easy helper class for comparing results from different solvers (just shows number of wins/draws etc.). Really bad design - just temporary ;-)
class FaspSolutionResult {
    std::vector<std::string> names;
    std::vector<int> results;
    std::vector<int> statistics;
    int eqStat = 0;
    int runCnt = 0;

    void print() {
        std::cout << "#WINS(";
        for (std::size_t i = 0; i < names.size(); ++i) {
            std::cout << names[i] << (i == names.size() - 1 ? "" : "/");
        }
        std::cout << ")=";

        for (std::size_t i = 0; i < names.size(); ++i) {
            std::cout << statistics[i] << (i == names.size() - 1 ? "" : "/");
        }
        std::cout << "     #DRAWS=" << eqStat << "     #RUNS=" << runCnt << std::endl;
    }

public:
    auto& getCnt(const std::string &aCounterName) {
        names.push_back(aCounterName);
        results.push_back(0);
        return results.back();
    }

    void calculateAndPrint() {
        size_t n = results.size();
        if (statistics.size() < n) statistics.resize(n);

        bool eqOnce = true;
        for (size_t i = 0; i < n; ++i) {
            auto a = results[i];
            bool best = true;
            bool eq = true;
            // Stats are really small so O(n^2) is OK ;-)
            for (size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                int b = results[j];
                if (a > b) {best = false; eq = false;}
                else if (a < b) {eq = false;}
            }

            if (best & !eq) statistics[i]++;
            else if (eq && eqOnce) { eqStat++; eqOnce = false; } // inc. only one time
        }

        runCnt++;
        print();

        names.clear();
        results.clear();
    }
};

void test(const char *inputDir) {

//    Graph::Graph<int, Graph::GraphMap> g;
//    for (int i = 0; i < 3; ++i) g.addVertex(i);
//    g.addEdge({0, 1});
//    g.addEdge({0,2});
//    g.addEdge({2,1});
//    auto c = Graph::Ext::getEdgeProperties<int16_t>(g, 2);
//    std::vector<int> mapv(3, 0);
//    std::cout << Graph::HIPR().runHipr(g, 0, 1, c, 2, mapv) << std::endl;
//    if (true) return;

    DataHdf5<double> f("/tmp/out.h5");
    std::string dir;
    dir = "/Users/gonciarz/ws2/repo/FASP-benchmarks/data/random/";
//    dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/FASP-benchmarks/data/de-bruijn/";
//    dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/FASP-benchmarks/data/tournaments";

    if (strcmp(inputDir, "") != 0) dir = inputDir;

    FaspSolutionResult fsr;
    Timer<true, false> t{};
    [[maybe_unused]] int limitCnt = 23;
    int cnt = 0;


    auto checkIfExist = [](const std::vector<std::string> &files, const std::string &file) -> bool {return std::find(files.begin(), files.end(), file) != files.end();};

    double timeExact=0.0;
    double timeRandom=0.0;

    auto files = Graph::IO::getFilesInDir(std::string(dir));

    std::sort(files.begin(), files.end());
    std::cout << files.size() << std::endl;
    for (auto &graphFile : files) {
        if (!Tools::endsWith(graphFile, ".al")) continue;

//        graphFile = "random-1463-410-533.al"; // 0.1s
//        graphFile = "random-1833-500-700.al"; // 1s
//        graphFile = "random-1224-350-700.al"; // 15s

        auto solutionFile{graphFile}; Tools::replace(solutionFile, ".al", ".mfas");
        auto timeoutFile{graphFile}; Tools::replace(timeoutFile, ".al", ".timeout");
        auto timingFile{graphFile}; Tools::replace(timingFile, ".al", ".timing");

        if (!checkIfExist(files, solutionFile) || checkIfExist(files, timeoutFile) || !checkIfExist(files, timingFile)) {
//            std::cout << "Skipping [" << graphFile << "] - solution/timing not found or there was a timing\n";
            continue;
        }

//        if (limitCnt-- == 0) break;

        Graph::Graph g = Graph::IO::graphFromFile<int16_t>(dir + "/" + graphFile);
        auto c = Graph::Ext::getEdgeProperties<int16_t>(g, 1);
        Graph::IO::graphToFile("/tmp/myGraph.txt", g);
        int solution = Graph::IO::solutionFromFile(dir + "/" + solutionFile);
        double timeExactOfGraph = Graph::IO::timingFromFile(dir + "/" + timingFile);

//        if (timeExact >= 1) continue;
//        auto density = (double)g.getNumOfEdges() / g.getNumOfVertices();
//        if ( density > 2.5 || density < 2) continue;


        std::cout << ++cnt << " ======================= Processing [" << graphFile << "] " << g << std::endl;
        std::cout << "Exact solution=" << solution << " Time=" << timeExactOfGraph << std::endl;
        timeExact += timeExactOfGraph;

        t.start_timer("--------RANDOM new");
        auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::FaspFastFinal::randomFASP(g, c);
        auto thisStep = t.stop_timer();
        LOG(DEBUG) << "SA / SA RND / RED RND edges: " << saEdgesCnt << " / " << saRndEdgesCnt << " / " << redRndEdgesCnt;
        LOG(DEBUG) << "FASP(RAND)  capacity = " << capacity << " edgeCnt = " << removedEdges.size() << " edgeList = " << removedEdges;
        fsr.getCnt("NewRandom") = capacity;
        timeRandom+=thisStep;

        t.start_timer("gr");
        auto gr = Graph::FaspTools::GR(g, c);
        auto grTime = t.stop_timer();

//        t.start_timer("---------Orig");
//        fsr.getCnt("OrigRandom") = Graph::FaspFast2::randomFASP_orig(g, c);
//        t.stop_timer();
//
//        t.start_timer("---------GR");
//        auto [grc, dummy] = Graph::Fasp::GR(g, c);
//        std::cout << "GR CAPACITY = " << grc << std::endl;
//        fsr.getCnt("GR") = grc;
//        t.stop_timer();

        fsr.getCnt("Exact") = solution;

        std::cout << "===========> ";
        fsr.calculateAndPrint();


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
    }

    f.save();
    std::cout << "TIME exact/random: " << timeExact << "/" << timeRandom << std::endl;
}


int main(int argc, char *argv[]) {
    configureLogger();
    test((argc > 1 ? argv[1] : ""));

    return 0;
}
