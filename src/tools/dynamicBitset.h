//
// Created by gonciarz on 2019-04-02.
//

#ifndef TEST_DYNAMICBITSET_H
#define TEST_DYNAMICBITSET_H

#include <sstream>
#include <iostream>
#include <cassert>
#include <bitset>
#include <memory>
#include <cstring>


template <typename ELEMENT_TYPE = uint8_t, typename IDX=uint32_t>
class DynamicBitset {
    // We use negative index to indicate empty bitset
    static_assert(!std::is_signed<ELEMENT_TYPE>::value && "Element type of bitset must be unsinged");

    static constexpr int BitsPerElement = sizeof(ELEMENT_TYPE) * 8;
    static constexpr ELEMENT_TYPE Bit = 1;

    const IDX iSize; // size of bitset in bits
    size_t iNumOfElements; // size of bitset in data elements
    std::unique_ptr<ELEMENT_TYPE[]> iData;

    auto get(IDX aBitNum) {
        IDX idx = aBitNum / BitsPerElement;
        IDX off = aBitNum % BitsPerElement;
        assert(idx >= 0 && idx < iNumOfElements && "Wrong index (bit number too big)");
        assert(idx * BitsPerElement + off < iSize && "Wrong offset (bit number too big)");
        return std::pair{std::ref(iData[idx]), Bit << off};
    }

public:
    /**
     * Create bitset
     * @param aSize size of bitset in number of bits
     */
    DynamicBitset(IDX aSize) : iSize(aSize) {
        iNumOfElements = (iSize + BitsPerElement - 1) / BitsPerElement;
        iData.reset(new ELEMENT_TYPE[iNumOfElements]);
        clearAll();
    }

    IDX getSize() const {return iSize;}

    /**
     * Set bit to 1
     * @param aBitNum
     */
    void set(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        d |= m;
    }

    /**
     * clear bit to 0
     * @param aBitNum
     */
    void clear(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        d &= ~m;
    }

    /**
     * Check if bit is set
     * @param aBitNum
     * @return true if set, false otherwise
     */
    bool test(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        return d & m;
    }

    /**
     * Set all bits to 0
     */
    void clearAll() {
        memset(iData.get(), 0, sizeof(ELEMENT_TYPE) * iNumOfElements);
    }

    friend std::ostream &operator<<(std::ostream &os, const DynamicBitset &obj) {
        os << "DynamicBitset<" << BitsPerElement << ">{";
        for (size_t i = 0; i < obj.iNumOfElements; ++i) {
            std::bitset<BitsPerElement> bs(obj.iData[i]);
            os << bs;
            if (i < obj.iData.size() - 1) os << ", ";
        }
        os << "}";
        return os;
    }
};


#endif
