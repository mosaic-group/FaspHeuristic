#include "../include/FaspTightCut/graph.h"
#include "graph/graphIO.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <iostream>


/**
 * This test is using c++17 features included in 'filesystem' header. On some older compiler this is not yet implemented.
 * So it skips these tests... 
 */
#if __has_include(<filesystem>)
#include <filesystem>

namespace {

    using ::testing::UnorderedElementsAreArray;

    TEST(GraphIO, testCreateReadUnweighted) {
        // Simple test for saving and loading graph without weights
        // It creates small graph, saves it, and then loads it
        // Finally compares loaded graph with created graph and expect they
        // are same.

        // Prepare graph to test
        auto g = Graph::Graph<uint16_t>{};
        g.addVertex(1);
        g.addVertex(666);
        g.addVertex(3);
        g.addVertex(4);

        g.addEdge({1, 3});
        g.addEdge({3, 4});
        g.addEdge({4, 3});

        // Generate temporary file name
        auto tmpPath = std::string{std::filesystem::temp_directory_path()};
        std::string testGraphName{tmpPath + "/testGraph.txt"};

        // TESTED METHOD
        Graph::IO::graphToFile(testGraphName, g);

        // TESTED METHOD
        auto g2 = Graph::IO::graphFromFile(testGraphName);
        
        ASSERT_EQ(g.getNumOfEdges(), g2.getNumOfEdges());
        ASSERT_EQ(g.getNumOfVertices(), g2.getNumOfVertices());

        // Return elements might be in different order than in original graph
        ASSERT_THAT(g.getEdges(), UnorderedElementsAreArray(g2.getEdges()));
        ASSERT_THAT(g.getVertices(), UnorderedElementsAreArray(g2.getVertices()));

        // Clean up
        std::filesystem::remove(testGraphName);
    }

    TEST(GraphIO, testCreateReadWeighted) {
        // Simple test for saving and loading graph with weights
        // It creates small graph, saves it, and then loads it
        // Finally compares loaded graph with created graph and expect they
        // are same.

        // Prepare graph to test
        using VERTEX_TYPE = uint16_t;
        using EDGE_PROPERTY_TYPE = float;
        using Edge = Graph::Graph<VERTEX_TYPE>::Edge;

        auto g = Graph::Graph<VERTEX_TYPE>{};
        g.addVertex(1);
        g.addVertex(666);
        g.addVertex(3);
        g.addVertex(4);

        auto ep = Graph::Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROPERTY_TYPE>{};
        g.addEdge({1, 3}); ep.emplace(Edge{1, 3}, 4.5);
        g.addEdge({3, 4}); ep.emplace(Edge{3, 4}, 1.0);
        g.addEdge({4, 3}); ep.emplace(Edge{4, 3}, 2.1);

        // Generate temporary file name
        auto tmpPath = std::string{std::filesystem::temp_directory_path()};
        std::string testGraphName{tmpPath + "/testGraph.txt"};

        // TESTED METHOD
        Graph::IO::graphWithWeightsToFile(testGraphName, g, ep);

        // TESTED METHOD
        auto [g2, ep2] = Graph::IO::graphWithWeightsFromFile<VERTEX_TYPE, EDGE_PROPERTY_TYPE>(testGraphName);

        ASSERT_EQ(g.getNumOfEdges(), g2.getNumOfEdges());
        ASSERT_EQ(g.getNumOfVertices(), g2.getNumOfVertices());

        // Return elements might be in different order than in original graph
        ASSERT_THAT(g.getEdges(), UnorderedElementsAreArray(g2.getEdges()));
        ASSERT_THAT(g.getVertices(), UnorderedElementsAreArray(g2.getVertices()));
        ASSERT_THAT(ep, UnorderedElementsAreArray(ep2));

        // Clean up
        std::filesystem::remove(testGraphName);
    }

    TEST(GraphIO, testSolutionFromFile) {
        // Create simple solution file (list of FASP edges)
        auto tmpPath = std::string{std::filesystem::temp_directory_path()};
        std::string testSolutionFileName{tmpPath + "/testSolution.txt"};

        std::ofstream solFile(testSolutionFileName, std::ios_base::out);
        solFile << "2 8\n";
        solFile << "1 10\n";
        solFile << "2 3\n";
        solFile << "333 1\n";
        solFile.close();

        // TESTED METHOD
        int faspSize = Graph::IO::solutionFromFile(testSolutionFileName);

        ASSERT_EQ(faspSize, 4);

        // Clean up
        std::filesystem::remove(testSolutionFileName);
    }

    TEST(GraphIO, testTimingFromFile) {
        // Create simple timing file (with one number representing seconds)
        auto tmpPath = std::string{std::filesystem::temp_directory_path()};
        std::string testTimingFileName{tmpPath + "/testTiming.txt"};

        std::ofstream timeFile(testTimingFileName, std::ios_base::out);
        double timeExpected  = 1.234;
        timeFile << timeExpected << "\n";
        timeFile.close();

        // TESTED METHOD
        double timeRead = Graph::IO::timingFromFile(testTimingFileName);

        ASSERT_EQ(timeRead, timeExpected);

        // Clean up
        std::filesystem::remove(testTimingFileName);
    }

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif

