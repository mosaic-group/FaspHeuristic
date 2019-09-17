#include "tools.h"
#include "prettyprint.h"
#include "tools/easylogging++.h"

#include "graph/graph.h"
#include "graph/graphFaspFast2.h"
#include "graph/graphIO.h"

INITIALIZE_EASYLOGGINGPP

int main() {
    int cntA = 0;
    int cntB = 0;
    int cntAB = 0;
    std::vector<double> timesGS;
    std::vector<double> timesGS2;
    for (int i = 800; i <= 800; i += 30) {
        Timer<true, false> t("");
        int rep = 1;
        std::vector<double> tsa;
        std::vector<int> sa;
        double ct = 0;
        int cn = 0;
        for (int r = 0; r < rep; ++r) {
//                auto [ge, cc] = Graph::Fasp::generateGraphWithKnownFaspAndSameWeights<int, int, Graph::GraphMap>(i, 15, 4 * i);
            auto[ge, cc] = Graph::Tools::generateErdosRenyiGraph<int, int, Graph::GraphMap>(i, 1.4 * (double) i / (i * (i - 1)));

            Graph::IO::graphToFile<int, Graph::GraphMap>("/tmp/graph.txt", ge);
            Graph::Graph gg = Graph::IO::graphFromFile<int, Graph::GraphMap>("/tmp/graph.txt");


            auto vertices = gg.getVertices();
            auto maxId = std::max_element(vertices.begin(), vertices.end());
            Graph::FaspFast2::PathHero<int> path2(maxId == vertices.end() ? 1 : *maxId + 1);
            auto g1{gg};
            auto g2{gg};
            auto c1 = Graph::Ext::getEdgeProperties<int>(g1, 1);
            auto c2 = Graph::Ext::getEdgeProperties<int>(g2, 1);


//                t.start_timer("G");
//                Graph::FaspFast::randomFASP(g1, c1);
//                std::cout << "CNT1: " << path1.cnt << std::endl;
//                t.stop_timer();

            int b = 0;
//                t.start_timer("G2 orig");
//                b = Graph::FaspFast2::randomFASP_orig(g2, c2);
//                std::cout << "CNT SA:" << path2.saCnt << std::endl;
//                tsa.push_back(t.stop_timer());
//                ct += tsa.back();
//                sa.push_back(path2.saCnt);
//                cn += path2.saCnt;
//                path2.saCnt = 0;

            auto[ca, ed] = Graph::Fasp::GR(g2, c2);
            std::cout << "GR CAPACITY = " << ca << std::endl;

//
//
//                t.start_timer("G2 parallel");
//                int a = Graph::FaspFast2::randomFASP(g2, c2);
//                std::cout << "CNT SA: " << path2.saCnt << std::endl;
//                tsa.push_back(t.stop_timer());
//                ct += tsa.back();
//                sa.push_back(path2.saCnt);
//                cn += path2.saCnt;
//                path2.saCnt = 0;

            t.start_timer("G2 new");
            int a = Graph::FaspFast2::randomFASP_blueEdges(g1, c1);
            std::cout << "CNT SA: " << path2.saCnt << std::endl;
            tsa.push_back(t.stop_timer());
            ct += tsa.back();
            sa.push_back(path2.saCnt);
            cn += path2.saCnt;
            path2.saCnt = 0;

            if (a > b) cntB++;
            else if (b > a) cntA++;
            else cntAB++;
//
            std::cout << "======================>      " << cntA << " " << cntB << " " << cntAB << std::endl;
        }
        std::cout << "t = " << tsa << ";\n";
        std::cout << "n = " << sa << ";\n";
        std::cout << "SA time = " << ct / cn << std::endl;

    }
    return 0;
}
