#include "graph/graph.h"
#include "graph/graphTools.h"
#include "graph/graphFasp.h"

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

        // Tested method
        g.removeVertex(2);

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
    void testFindPathWithPositiveCapacity(Graph::Graph<VT> &g) {
        Graph::Graph<int> gv;
        for (int i = 1; i <= 4; ++i) gv.addVertex(i);
        gv.addEdge({1, 2});
        gv.addEdge({1, 3});
        gv.addEdge({2, 3});
        gv.addEdge({3, 4});
        gv.addEdge({2, 4});

        //       ____> 3_____
        //      |      ^     v
        //      1      |     4
        //      |____> 2_____^

        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 3}) = 0;
            weights.at({2, 4}) = -1;
            auto[pathExists, path] = Graph::Tools::findPathWithPositiveCapacity(gv, 1, 4, weights);
            ASSERT_TRUE(pathExists);
            ASSERT_THAT(path, ElementsAre(1, 2, 3, 4));
        }
        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 2}) = 0;
            auto[pathExists, path] = Graph::Tools::findPathWithPositiveCapacity(gv, 1, 4, weights);
            ASSERT_TRUE(pathExists);
            ASSERT_THAT(path, ElementsAre(1, 3, 4));
        }
        { // no path
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 2}) = 0;
            weights.at({1, 3}) = 0;
            auto[pathExists, path] = Graph::Tools::findPathWithPositiveCapacity(gv, 1, 4, weights);
            ASSERT_FALSE(pathExists);
        }
    }
}

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
