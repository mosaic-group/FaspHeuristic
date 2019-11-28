#include "graph/graphFasp.h"
#include "graph/graph.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>

using ::testing::UnorderedElementsAre;
using ::testing::ElementsAre;


namespace {

    TEST(TestGraphFasp, testPathHero) {
	using VERTEX_TYPE = uint16_t;
        Graph::Fasp::GraphSpeedUtils<VERTEX_TYPE> u{7};
        // Test graph:
        //                  /---> 3
        //                 /---> 4
        //  0 ---> 1 ---> 3 ---> 5 ---> 6
        //                ^-------------

        using Edge = typename Graph::Graph<VERTEX_TYPE >::Edge;

        Graph::Graph<VERTEX_TYPE> g;
        for (int i = 0; i < 7; ++i) g.addVertex(i);
        g.addEdge({0, 1});
        g.addEdge({1, 2});
        g.addEdge({1, 3});
        g.addEdge({1, 4});
        g.addEdge({3, 5});
        g.addEdge({5, 6});
        g.addEdge({6, 3});

        {   // 'pathExistsDFS'
            ASSERT_TRUE(u.pathExistsDFS(g, 0, 6));
            ASSERT_TRUE(u.pathExistsDFS<false>(g, 6, 0));

            ASSERT_FALSE(u.pathExistsDFS(g, 6, 0));
            ASSERT_FALSE(u.pathExistsDFS<false>(g, 0, 6));

            ASSERT_FALSE(u.pathExistsDFS(g, 5, 4));
            ASSERT_FALSE(u.pathExistsDFS<false>(g, 4, 5));

            ASSERT_TRUE(u.pathExistsDFS(g, 3, 6));
            ASSERT_TRUE(u.pathExistsDFS(g, 6, 3));

            ASSERT_TRUE(u.pathExistsDFS(g, 3, 3));
        }

        {   // 'findPathDFS'
            auto [pathExist, path] = u.findPathDFS(g, 1, 6);
            ASSERT_TRUE(pathExist);
            ASSERT_THAT(path, ElementsAre(1, 3, 5, 6));

            auto [pathExist2, path2] = u.findPathDFS(g, 6, 1);
            ASSERT_FALSE(pathExist2);
        }

        {   // 'isAcyclic'
            auto gg = g; // copy 'g'
            ASSERT_FALSE(u.isAcyclic(gg));
            gg.removeEdge({5, 6});
            ASSERT_TRUE(u.isAcyclic(gg));
        }

        {   // 'findEdgesWithCycles'
            auto edges = u.findEdgesWithCycles(g);

            ASSERT_EQ(edges.size(), 3);
            ASSERT_THAT(edges, UnorderedElementsAre(Edge{3, 5}, Edge{5, 6}, Edge{6, 3}));

            auto gg = g; // copy 'g'
            gg.removeEdge({5, 6});
            ASSERT_EQ(u.findEdgesWithCycles(gg).size(), 0);
        }

    }
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
