#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"

#include "graph/graphFaspFastFinal.h"
#include "graph/graph.h"
#include "graph/graphFaspFast2.h"
#include "graph/graphIO.h"
#include "hdf5/dataHdf5.h"

#include <string>
#include <cstddef>
using std::size_t;

INITIALIZE_EASYLOGGINGPP

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

void test() {
    DataHdf5<double> f("/tmp/out.h5");
    std::string dir;
    dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/FASP-benchmarks/data/random/";
//    dir = "/Users/gonciarz/Documents/MOSAIC/work/repo/FASP-benchmarks/data/de-bruijn/";
    FaspSolutionResult fsr;
Timer<true, false> t(true);
    int limitCnt = 1;
    int cnt = 0;

    auto checkIfExist = [](const std::vector<std::string> &files, const std::string &file) -> bool {return std::find(files.begin(), files.end(), file) != files.end();};

    auto files = Graph::IO::getFilesInDir(std::string(dir));
    std::sort(files.begin(), files.end());
    for (auto &graphFile : files) {
        if (!Tools::endsWith(graphFile, ".al")) continue;

//        graphFile = "random-1463-410-533.al"; // 0.1s
//        graphFile = "random-1833-500-700.al"; // 1s
        graphFile = "random-1224-350-700.al"; // 15s

        auto solutionFile{graphFile}; Tools::replace(solutionFile, ".al", ".mfas");
        auto timeoutFile{graphFile}; Tools::replace(timeoutFile, ".al", ".timeout");
        auto timingFile{graphFile}; Tools::replace(timingFile, ".al", ".timing");

        if (!checkIfExist(files, solutionFile) || checkIfExist(files, timeoutFile) || !checkIfExist(files, timingFile)) {
//            std::cout << "Skipping [" << graphFile << "] - solution/timing not found or there was a timing\n";
            continue;
        }

        if (limitCnt-- == 0) break;

        Graph::Graph g = Graph::IO::graphFromFile<int, Graph::GraphMap>(dir + "/" + graphFile);
        auto c = Graph::Ext::getEdgeProperties<int>(g, 1);
        Graph::IO::graphToFile("/tmp/myGraph.txt", g);
        int solution = Graph::IO::solutionFromFile(dir + "/" + solutionFile);
        double timeExact = Graph::IO::timingFromFile(dir + "/" + timingFile);

//        if (timeExact > 0.05) continue;
//        auto density = (double)g.getNumOfEdges() / g.getNumOfVertices();
//        if ( density > 2.5 || density < 2) continue;




        std::cout << ++cnt << " ======================= Processing [" << graphFile << "] " << g << std::endl;

        std::cout << "Exact solution=" << solution << " Time=" << timeExact << std::endl;

        t.start_timer("G2 new");
//        Graph::FaspFast2::randomFASP_blueEdges(g, c);

        fsr.getCnt("NewRandom") = Graph::FaspFastFinal::randomFASP(g, c);
        t.stop_timer();
        t.start_timer("Orig");
        fsr.getCnt("OrigRandom") = Graph::FaspFast2::randomFASP_orig(g, c);
        t.stop_timer();

        fsr.getCnt("Exact") = solution;
                    std::cout << "===========> ";
            fsr.calculateAndPrint();
//        f.put("exact", solution);
//        f.put("edges", gv.getNumOfEdges());
//        f.put("vertices", gv.getNumOfVertices());
    }
}

//void testBlueEdges() {
//    FaspSolutionResult fsr;
//    Timer<true, false> t("");
//
//    for (int i = 77; i <= 77; i += 30) {
//        int rep = 10;
//
//        for (int r = 0; r < rep; ++r) {
//            std::cout << "\n\n\n";
//
////                auto [ge, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(i, 15, 4 * i);
//            auto[ge, cc] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(i, 1.9 * (double) i / (i * (i - 1)));
//            Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", ge);
//            Graph::Graph gg = Graph::IO::graphFromFile<int, Graph::GraphMap>("/tmp/graph.txt");
//
//            auto g2{gg};
//            auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);
//
//            t.start_timer("G2 orig");
//            fsr.getCnt("OrigRandom") = Graph::FaspFast2::randomFASP_orig(g2, c2);
//            t.stop_timer();
//
//            t.start_timer("GR");
//            auto[ca, ed] = Graph::Fasp::GR(g2, c2);
//            std::cout << "GR CAPACITY = " << ca << std::endl;
//            fsr.getCnt("GR") = ca;
//            t.stop_timer();
//
//            t.start_timer("G2 new");
//            fsr.getCnt("NewRandom") = Graph::FaspFast2::randomFASP_blueEdges(g2, c2);
//            t.stop_timer();
//
//            std::cout << "===========> ";
//            fsr.calculateAndPrint();
//        }
//    }
//}

int main() {
    configureLogger();
//    testBlueEdges();
    test();

    return 0;
}
