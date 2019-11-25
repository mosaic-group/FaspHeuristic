#include "graph/graph.h"
#include "graph/graphFaspTools.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <random>

using ::testing::ElementsAre;
using ::testing::UnorderedElementsAre;

namespace {


    template <typename T>
    class GraphTest : public ::testing::Test {
    public:
        using EdgeType =  typename Graph::Graph<T>::Edge;
    };

    using GraphTypes = ::testing::Types<int, int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t>;
    TYPED_TEST_SUITE(GraphTest, GraphTypes,);

    template <typename VT>
    void testCreateEmpty(Graph::Graph<VT> &g) {
        ASSERT_EQ(g.getNumOfVertices(), 0);
        ASSERT_EQ(g.getNumOfEdges(), 0);
        ASSERT_EQ(g.getVertices().size(), 0);
        ASSERT_EQ(g.getEdges().size(), 0);
    }

    TYPED_TEST(GraphTest, createEmptyGraph) {
        Graph::Graph<TypeParam> gm;
        testCreateEmpty(gm);
    }

    template <typename VT>
    void testCreate(Graph::Graph<VT> &g) {
        // Tested methods
        g.addVertex(1);
        g.addVertex(2);
        g.addVertex(3);
        g.addEdge({1, 2});
        g.addEdge({2, 3});
        g.addEdge({3, 1});
        g.addEdge({1, 3});

        ASSERT_EQ(g.getNumOfVertices(), 3);
        ASSERT_EQ(g.getNumOfEdges(), 4);
        ASSERT_EQ(g.getVertices().size(), 3);
        ASSERT_EQ(g.getEdges().size(), 4);
        ASSERT_TRUE(g.hasEdge({1, 2}));
        ASSERT_TRUE(g.hasEdge({2, 3}));
        ASSERT_TRUE(g.hasEdge({3, 1}));
        ASSERT_TRUE(g.hasEdge({1, 3}));
        ASSERT_FALSE(g.hasEdge({2, 1}));
        ASSERT_FALSE(g.hasEdge({3, 2}));
        ASSERT_FALSE(g.hasEdge({1, 4}));
        ASSERT_FALSE(g.hasEdge({4, 1}));
    }

    TYPED_TEST(GraphTest, createGraph) {
        Graph::Graph<TypeParam> gm;
        testCreate(gm);
    }

    template <typename VT>
    void testRemoveVertex(Graph::Graph<VT> &g) {
        g.addVertex(1);
        g.addVertex(2);
        g.addVertex(3);
        g.addEdge({1, 2});
        g.addEdge({2, 3});
        g.addEdge({3, 1});
        g.addEdge({1, 3});
        ASSERT_EQ(g.getNumOfVertices(), 3);
        ASSERT_EQ(g.getNumOfEdges(), 4);

        ASSERT_TRUE(g.hasVertex(2));

        // Tested method
        g.removeVertex(2);

        ASSERT_FALSE(g.hasVertex(2));

        ASSERT_EQ(g.getVertices().size(), 2);
        ASSERT_EQ(g.getEdges().size(), 2);
        ASSERT_FALSE(g.hasEdge({1, 2}));
        ASSERT_FALSE(g.hasEdge({2, 3}));
        ASSERT_TRUE(g.hasEdge({3, 1}));
        ASSERT_TRUE(g.hasEdge({1, 3}));
        ASSERT_FALSE(g.hasEdge({2, 1}));
        ASSERT_FALSE(g.hasEdge({3, 2}));
    }

    TYPED_TEST(GraphTest, removeVertex) {
        Graph::Graph<TypeParam> gm;
        testRemoveVertex(gm);
    }

    template <typename VT>
    void testGraphInfo(Graph::Graph<VT> &g) {
        g.addVertex(1);
        g.addVertex(2);

        // add/remove vertex
        ASSERT_FALSE(g.hasVertex(3));
        ASSERT_EQ(g.getNumOfVertices(), 2);
        g.addVertex(3);
        ASSERT_TRUE(g.hasVertex(3));
        ASSERT_EQ(g.getNumOfVertices(), 3);
        g.removeVertex(3);
        ASSERT_FALSE(g.hasVertex(3));
        ASSERT_EQ(g.getNumOfVertices(), 2);

        // add/remove edge
        ASSERT_FALSE(g.hasEdge({1, 2}));
        ASSERT_EQ(g.getNumOfEdges(), 0);
        g.addEdge({1, 2});
        ASSERT_TRUE(g.hasEdge({1, 2}));
        ASSERT_EQ(g.getNumOfEdges(), 1);
        g.removeEdge({1, 2});
        ASSERT_FALSE(g.hasEdge({1, 2}));
        ASSERT_EQ(g.getNumOfEdges(), 0);
    }

    TYPED_TEST(GraphTest, graphInfo) {
        Graph::Graph<TypeParam> gm;
        testGraphInfo(gm);
    }

    template <typename VT>
    void testGraphGetters(Graph::Graph<VT> &g) {
        using Edge = typename Graph::Graph<VT>::Edge;

        g.addVertex(1);
        g.addVertex(2);
        g.addVertex(4);

        // we cannot assume the order of elements
        ASSERT_THAT(g.getVertices(), UnorderedElementsAre(1, 2, 4));

        g.addEdge({1, 2});
        g.addEdge({1, 4});
        g.addEdge({4, 2});
        g.addEdge({2, 1});

        // we cannot assume the order of elements
        ASSERT_THAT(g.getEdges(), UnorderedElementsAre(Edge{1, 2}, Edge{1, 4}, Edge{4, 2}, Edge{2, 1}));
    }

    TYPED_TEST(GraphTest, graphGetters) {
        Graph::Graph<TypeParam> gm;
        testGraphGetters(gm);
    }

// Testing asserts - these tests can be run only in debug mode
#ifndef NDEBUG
    template <typename VT>
    void testGraphAsserts(Graph::Graph<VT> &g) {
        //Add same vertex twice
        g.addVertex(2);
        ASSERT_DEATH(g.addVertex(2), "Vertex already exists!");

        // Remove same vertex twice
        g.removeVertex(2);
        ASSERT_DEATH(g.removeVertex(2), "Vertex does not exists!");

        // Add edge - vertices not existing
        ASSERT_DEATH(g.addEdge({1, 3}), "Source vertex does not exists!");
        g.addVertex(1);
        ASSERT_DEATH(g.addEdge({1, 3}), "Destination vertex does not exists!");

        // Remove edge - vertices not existing
        ASSERT_DEATH(g.removeEdge({2, 3}), "Source vertex does not exists!");
        g.addVertex(2);
        ASSERT_DEATH(g.removeEdge({2, 3}), "Destination vertex does not exists!");

        // Get in/out vetices
        ASSERT_DEATH(g.getInVertices(5), "Vertex does not exists!");
        ASSERT_DEATH(g.getOutVertices(6), "Vertex does not exists!");
    }

    TYPED_TEST(GraphTest, graphAsserts) {
        Graph::Graph<TypeParam> gm;
        testGraphAsserts(gm);
    }
#endif

}

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
