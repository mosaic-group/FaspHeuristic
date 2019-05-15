//
// Created by gonciarz on 2019-04-02.
//

#include "tools/dynamicBitset.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

template <typename T>
class BitSetTest : public ::testing::Test{};

using BitsetElementTypes = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;
TYPED_TEST_SUITE(BitSetTest, BitsetElementTypes,);

TYPED_TEST(BitSetTest, easyTest) {
    for (int bitsetSize = 0; bitsetSize <= 129; ++bitsetSize) {
        DynamicBitset<TypeParam> bitset(bitsetSize);

        // check size
        ASSERT_EQ(bitsetSize, bitset.getSize());

        // check if all 0
        for (int i = 0; i < bitsetSize; ++i) {
            ASSERT_FALSE(bitset.test(i));
        }

        // set only one and check then clear and set next
        for (int i = 0; i < bitsetSize; ++i) {
            bitset.set(i);
            for (int b = 0; b < bitsetSize; ++b) {
                if (b != i) ASSERT_FALSE(bitset.test(b));
                else ASSERT_TRUE(bitset.test(b));
            }
            bitset.clear(i);
        }


        bitset.clearAll();
        // check if all 0
        for (int i = 0; i < bitsetSize; ++i) {
            ASSERT_FALSE(bitset.test(i));
        }

        // TODO: Check for trying to use bit larger than initial size
        //       check asserts etc.
    }
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}