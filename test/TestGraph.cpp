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

    template <typename VT, template <typename> class G>
    void testCreateEmpty(Graph::Graph<VT, G> &g) {
        ASSERT_EQ(g.getNumOfVertices(), 0);
        ASSERT_EQ(g.getNumOfEdges(), 0);
        ASSERT_EQ(g.getVertices().size(), 0);
        ASSERT_EQ(g.getEdges().size(), 0);
    }

    TYPED_TEST(GraphTest, createEmptyGraph) {
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testCreateEmpty(gm);
        testCreateEmpty(gv);
    }

    template <typename VT, template <typename> class G>
    void testCreate(Graph::Graph<VT, G> &g) {
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
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testCreate(gm);
        testCreate(gv);
    }

    template <typename VT, template <typename> class G>
    void testRemoveVertex(Graph::Graph<VT, G> &g) {
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
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testRemoveVertex(gm);
        testRemoveVertex(gv);
    }

    template <typename VT, template <typename> class G>
    void testFaspG(Graph::Graph<VT, G> &g) {
        for (VT v = 0; v < 5; ++v) g.addVertex(v);
        g.addEdge({0, 1});
        g.addEdge({0, 2});
        g.addEdge({0, 3});
        g.addEdge({3, 2});
        g.addEdge({2, 1});
        g.addEdge({3, 1});
        g.addEdge({1, 4});
        g.addEdge({3, 4});
        g.addEdge({4, 0});

        ASSERT_EQ(g.getNumOfVertices(), 5);
        ASSERT_EQ(g.getNumOfEdges(), 9);

        // Tested method
        Graph::Fasp::G(g, {0, 1});

        ASSERT_EQ(g.getNumOfVertices(), 4);
        ASSERT_EQ(g.getNumOfEdges(), 6);
        ASSERT_TRUE(g.hasEdge({0, 1}));
        ASSERT_TRUE(g.hasEdge({0, 2}));
        ASSERT_TRUE(g.hasEdge({0, 3}));
        ASSERT_TRUE(g.hasEdge({3, 2}));
        ASSERT_TRUE(g.hasEdge({2, 1}));
        ASSERT_TRUE(g.hasEdge({3, 1}));
    }

    TYPED_TEST(GraphTest, faspG) {
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testFaspG(gm);
        testFaspG(gv);
    }

    template <typename VT, template <typename> class G>
    void testFaspGStar(Graph::Graph<VT, G> &g) {
        for (VT v = 0; v < 7; ++v) g.addVertex(v);
        g.addEdge({0, 1});
        g.addEdge({1, 2});
        g.addEdge({2, 3});
        g.addEdge({3, 0});

        g.addEdge({2, 4});
        g.addEdge({4, 2});

        g.addEdge({3, 5});
        g.addEdge({5, 6});
        g.addEdge({6, 5});
        g.addEdge({6, 2});

        g.addEdge({3, 3});

        // Tested method
        Graph::Fasp::GStar(g, {0, 1});

        ASSERT_EQ(g.getNumOfVertices(), 4);
        ASSERT_EQ(g.getNumOfEdges(), 3);
        ASSERT_TRUE(g.hasEdge({1, 2}));
        ASSERT_TRUE(g.hasEdge({2, 3}));
        ASSERT_TRUE(g.hasEdge({3, 0}));
    }

    TYPED_TEST(GraphTest, faspGStar) {
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testFaspGStar(gm);
        testFaspGStar(gv);
    }

    template <typename VT, template <typename> class G>
    void testFindPathWithPositiveCapacity(Graph::Graph<VT, G> &g) {
        Graph::Graph<int, Graph::GraphVector> gv;
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

    TYPED_TEST(GraphTest, findPathWithPositiveCapacity) {
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testFindPathWithPositiveCapacity(gm);
        testFindPathWithPositiveCapacity(gv);
    }

    template <typename VT, template <typename> class G>
    void testMinStCutFordFulkerson(Graph::Graph<VT, G> &gv) {
        using Edge = typename Graph::Graph<VT>::Edge;

        for (int i = 1; i <= 4; ++i) gv.addVertex(i);
        gv.addEdge({1, 2});
        gv.addEdge({1, 3});
        gv.addEdge({2, 4});
        gv.addEdge({3, 4});
        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 2}) = 4;
            weights.at({1, 3}) = 3;
            ASSERT_EQ(2, Graph::Tools::minStCutFordFulkerson(gv, 1, 4, weights));
            auto [mc, edges] = Graph::Tools::minStCutFordFulkersonEdges(gv, 1, 4, weights);
            ASSERT_EQ(2, mc);
            ASSERT_THAT(edges, UnorderedElementsAre(Edge{2, 4}, Edge{3, 4}));
        }
        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 2}) = 4;
            weights.at({3, 4}) = 3;
            ASSERT_EQ(2, Graph::Tools::minStCutFordFulkerson(gv, 1, 4, weights));
            auto [mc, edges] = Graph::Tools::minStCutFordFulkersonEdges(gv, 1, 4, weights);
            ASSERT_EQ(2, mc);
            ASSERT_THAT(edges, UnorderedElementsAre(Edge{2, 4}, Edge{1, 3}));
        }
        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({1, 3}) = 4;
            weights.at({2, 4}) = 3;
            ASSERT_EQ(2, Graph::Tools::minStCutFordFulkerson(gv, 1, 4, weights));
            auto [mc, edges] = Graph::Tools::minStCutFordFulkersonEdges(gv, 1, 4, weights);
            ASSERT_EQ(2, mc);
            ASSERT_THAT(edges, UnorderedElementsAre(Edge{1, 2}, Edge{3, 4}));
        }
        {
            auto weights = Graph::Ext::getEdgeProperties<int>(gv, 1);
            weights.at({2, 4}) = 4;
            weights.at({3, 4}) = 3;
            ASSERT_EQ(2, Graph::Tools::minStCutFordFulkerson(gv, 1, 4, weights));
            auto [mc, edges] = Graph::Tools::minStCutFordFulkersonEdges(gv, 1, 4, weights);
            ASSERT_EQ(2, mc);
            ASSERT_THAT(edges, UnorderedElementsAre(Edge{1, 2}, Edge{1, 3}));
        }
    }

    TYPED_TEST(GraphTest, minStCutFordFulkerson) {
        Graph::Graph<TypeParam, Graph::GraphMap> gm;
        Graph::Graph<TypeParam, Graph::GraphVector> gv;
        testMinStCutFordFulkerson(gm);
        testMinStCutFordFulkerson(gv);
    }


    // ====================== playground =========================================

    TEST(GraphTest, removeEdge) {
        Graph::Graph<int, Graph::GraphVector> gv;
        gv.addVertex(1);gv.addVertex(2);gv.addVertex(3);
        gv.addEdge({1, 2});
        gv.addEdge({2, 3});
        gv.addEdge({3, 1});

        auto c{gv};
//        std::cout << c.getStrRepresentationOfGraph() << std::endl;
        Graph::Fasp::G(c, {2, 1});

        Graph::Fasp::GStar(gv, {2, 1});
    }
}

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
