//
// Created by gonciarz on 2019-04-02.
//

#include "tools/dynamicBitset.h"

#include "gtest/gtest.h"

template <typename T>
class BitSetTest : public ::testing::Test{};

using BitsetElementTypes = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;
TYPED_TEST_SUITE(BitSetTest, BitsetElementTypes,);

TYPED_TEST(BitSetTest, easyTest) {
    for (int bitsetSize = 0; bitsetSize <= 129; ++bitsetSize) {
        // create object to test
        DynamicBitset<TypeParam> bitset(bitsetSize);

        // check size
        ASSERT_EQ(bitsetSize, bitset.getSize());

        // initially bitset should be zeroed
        for (int i = 0; i < bitsetSize; ++i) {
            ASSERT_FALSE(bitset.test(i));
        }

        // set only one bit and check if works then clear it and set next
        for (int i = 0; i < bitsetSize; ++i) {
            bitset.set(i);
            for (int b = 0; b < bitsetSize; ++b) {
                if (b != i) ASSERT_FALSE(bitset.test(b));
                else
                    ASSERT_TRUE(bitset.test(b));
            }
            bitset.clear(i);
        }

        // set all bits...
        for (int i = 0; i < bitsetSize; ++i) {
            bitset.set(i);
        }

        // ...and test if clearAll works properly
        bitset.clearAll();

        for (int i = 0; i < bitsetSize; ++i) {
            ASSERT_FALSE(bitset.test(i));
        }
    }
}

// In release mode there are no assertions...
#ifndef NDEBUG

TYPED_TEST(BitSetTest, assertionsTest) {
    int bitsetSize = 65;
    DynamicBitset<TypeParam> bitset(bitsetSize);

    // Beyond element capacity of bitset
    ASSERT_DEATH(bitset.test(200), ".*Wrong index.*");

    // beyond offset (such element exists but its bit are not fully used)
    ASSERT_DEATH(bitset.test(66), ".*Wrong offset.*");

    // negative indices are not allowed
    ASSERT_DEATH(bitset.test(-1), ".*Wrong index.*");
}

#endif

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
