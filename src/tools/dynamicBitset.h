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
    static_assert(!std::is_signed<ELEMENT_TYPE>::value && "Element type of bitset must be unsinged");

    static constexpr int BitsPerElement = sizeof(ELEMENT_TYPE) * 8;
    static constexpr ELEMENT_TYPE Bit = 1;

    const IDX iSize; // size of bitset in bits
    size_t iNumOfElements; // size of bitset in data elements
    std::unique_ptr<ELEMENT_TYPE[]> iData;

public:
    DynamicBitset(IDX aSize) : iSize(aSize) {
        iNumOfElements = (iSize + BitsPerElement - 1) / BitsPerElement;
        iData.reset(new ELEMENT_TYPE[iNumOfElements]);
        clearAll();
    }
    DynamicBitset(DynamicBitset &obj) : iSize(obj.iSize) {
        iNumOfElements = obj.iNumOfElements;
        iData.reset(new ELEMENT_TYPE[iNumOfElements]);
        memcpy(iData.get(), obj.iData.get(), sizeof(ELEMENT_TYPE) * iNumOfElements);
    }

    DynamicBitset(DynamicBitset&&) = default;
    DynamicBitset& operator=(DynamicBitset&&) = default;

    IDX getSize() const {return iSize;}

    auto get(IDX aBitNum) {
        IDX idx = aBitNum / BitsPerElement;
        IDX off = aBitNum % BitsPerElement;
        assert(idx >= 0 && idx < iNumOfElements && "Wrong index (bit number too big)");
        assert(idx * BitsPerElement + off < iSize && "Wrong offset (bit number too big)");
        return std::pair{std::ref(iData[idx]), Bit << off};
    }

    void set(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        d |= m;
    }

    void clear(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        d &= ~m;
    }

    bool test(IDX aBitNum) {
        auto [d, m] = get(aBitNum);
        return d & m;
    }

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
