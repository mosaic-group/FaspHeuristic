#include "graph/runGraph.h"
#include "tools/easylogging++.h"
#include "tools/tclap/CmdLine.h"
#include "tools/prettyprint.h"


#include "benchmarks/benchSuperAlgorithm/benchSuperAlgorithm.h"

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

void playgound(int v) {
    std::cout << "Playground " << v << "\n";
}

auto PrintAppArgs = [](int argc, char **argv) {std::cout << argc; for (int i = 0; i < argc; ++i) std::cout << " [" << argv[i] << "]"; std::cout << "\n";};

int main(int argc, char **argv) {
    configureLogger();

    std::cout << Tools::linspace(1, 45, 8) << std::endl;
    std::cout << Tools::logspace(1, 100, 20) << std::endl;
    if (true) return 0;

    try {
        TCLAP::CmdLine cmd("FASP heuristic benchamarks", ' ', "1.0", true);

        // Mandatory field
        std::vector<std::string> allowedBenchmarks;
        allowedBenchmarks.push_back("playground");
        allowedBenchmarks.push_back("benchSuperAlgorithmConstWeightVarFaspConstVE");
        allowedBenchmarks.push_back("benchConstVEvarFASP");
        TCLAP::ValuesConstraint<std::string> allowedVals( allowedBenchmarks );
        TCLAP::UnlabeledValueArg<std::string>  benchmarkName("benchmarkName", "name of benchmark to run", true, "", &allowedVals);
        cmd.add(benchmarkName);

        // Helper fields - they will be checked later on per-benchmark basis (each benchmark may req. different set of those fields)
        TCLAP::ValueArg<std::string> dirArg("d", "outputDirectory", "directory where output files will be saved", false, "", "outputDirectory");

        TCLAP::ValueArg<int> vArg("y", "v", "number of vertices in graph", false, 0, "#vertices");
        TCLAP::ValueArg<int> vMinArg("b", "vmin", "min (begin) number of vertices in graph", false, 0, "#minNumberOfVertices");
        TCLAP::ValueArg<int> vMaxArg("e", "vmax", "max (end) number of vertices in graph", false, 0, "#maxNumberOfVertices");

        TCLAP::ValueArg<int> eArg("w", "e", "number of vertices in graph", false, 0, "#edges");
        TCLAP::ValueArg<int> eMinArg("s", "emin", "min (begin) number of edges in graph", false, 0, "#minNumberOfEdges");
        TCLAP::ValueArg<int> eMaxArg("t", "emax", "max (end) number of edges in graph", false, 0, "#maxNumberOfEdges");

        TCLAP::ValueArg<int> fArg("f", "f", "size of FASP solution in graph", false, 0, "#faspSize");
        TCLAP::ValueArg<int> fMinArg("m", "fmin", "min (begin) size of FASP solution in graph", false, 0, "#minFaspSize");
        TCLAP::ValueArg<int> fMaxArg("x", "fmax", "max (end) size of FASP solution in graph", false, 0, "#maxFaspSize");

        TCLAP::ValueArg<int> stepsArg("z", "steps", "number of steps to be taken in given range of v/e/f", false, 0, "#steps");
        TCLAP::ValueArg<int> repsArg("r", "reps", "number of repetitions for each step", false, 0, "#repetitions");

        cmd.add(dirArg);

        cmd.add(vArg);
        cmd.add(vMinArg);
        cmd.add(vMaxArg);

        cmd.add(eArg);
        cmd.add(eMinArg);
        cmd.add(eMaxArg);

        cmd.add(fArg);
        cmd.add(fMinArg);
        cmd.add(fMaxArg);

        cmd.add(stepsArg);
        cmd.add(repsArg);

        cmd.parse( argc, argv );

        auto reqArgHdl = [] (TCLAP::ValueArg<int> &arg) {
            if (!arg.isSet()) {
                LOG(ERROR) << "Argument: [" << arg.longID("") << " " << arg.getDescription() << "] is required!";
                exit(-1);
            }
            return arg.getValue();
        };
        auto dirArgHdl = [] (TCLAP::ValueArg<std::string> &arg) {
            return arg.getValue();
        };


        if (benchmarkName.getValue() == allowedBenchmarks[0]) { // playground
            playgound(reqArgHdl(vArg));
        }
        else if (benchmarkName.getValue() == allowedBenchmarks[1]) {
            benchSuperAlgorithmConstWeightVarFaspConstVE(dirArgHdl(dirArg), reqArgHdl(vArg), reqArgHdl(eArg), reqArgHdl(fMinArg), reqArgHdl(fMaxArg), reqArgHdl(stepsArg), reqArgHdl(repsArg));
        }
        else {
            LOG(ERROR) << "Unknown benchmark [" << benchmarkName.getValue() << "]";
        }
    }
    catch (TCLAP::ArgException &e) {
        LOG(ERROR) << "Cmd line error: " << e.error() << " for arg " << e.argId();
    }

    return 0;
}
