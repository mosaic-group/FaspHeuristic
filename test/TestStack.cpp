#include "FaspTightCut/tools/stack.h"

#include "gtest/gtest.h"

#include <iostream>

namespace {

    TEST(TestStack, testStack) {
        Stack s{4};

        ASSERT_EQ(s.capacity(), 4);
        ASSERT_EQ(s.size(), 0);
        ASSERT_TRUE(s.empty());

        s.push(5);
        ASSERT_EQ(s.capacity(), 4);
        ASSERT_EQ(s.size(), 1);
        ASSERT_FALSE(s.empty());

        s.push(4);
        s.push(3);
        s.push(2);
        ASSERT_EQ(s.size(), 4);
        ASSERT_FALSE(s.empty());
        ASSERT_EQ(2, s.pop());
        ASSERT_EQ(s.size(), 3);
        ASSERT_EQ(3, s.pop());
        ASSERT_EQ(s.size(), 2);
        ASSERT_EQ(4, s.pop());
        ASSERT_EQ(s.size(), 1);
        ASSERT_EQ(5, s.pop());
        ASSERT_EQ(s.size(), 0);
        ASSERT_TRUE(s.empty());
    }

    // In release mode there are no assertions...
#ifndef NDEBUG
    TEST(TestStack, assertionsTest) {
        Stack s{4};

        ASSERT_DEATH(s.pop(), "All elements were pop from a stack!");
        s.push(1);
        s.push(2);
        s.push(3);
        s.push(4);
        ASSERT_DEATH(s.push(5), "Too many elements already in a stack!");

    }

#endif
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
