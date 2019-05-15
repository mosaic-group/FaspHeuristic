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

    const IDX iSize;
    size_t numOfElements;
    std::unique_ptr<ELEMENT_TYPE[]> iData;

public:
    DynamicBitset(IDX aSize) : iSize(aSize) {
        numOfElements = (iSize + BitsPerElement - 1) / BitsPerElement;
        iData.reset(new ELEMENT_TYPE[numOfElements]);
        clearAll();
    }
    DynamicBitset(DynamicBitset &obj) : iSize(obj.iSize) {
        numOfElements = obj.numOfElements;
        iData.reset(new ELEMENT_TYPE[numOfElements]);
        memcpy(iData.get(), obj.iData.get(), sizeof(ELEMENT_TYPE) * numOfElements);
    }

    DynamicBitset(DynamicBitset&&) = default;
    DynamicBitset& operator=(DynamicBitset&&) = default;

    IDX getSize() const {return iSize;}

    auto get(IDX aBitNum) {
        IDX idx = aBitNum / BitsPerElement;
        IDX off = aBitNum % BitsPerElement;
        assert(idx < iSize && "Wrong index (bit number too big)");
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
        memset(iData.get(), 0, sizeof(ELEMENT_TYPE) * numOfElements);
    }

    friend std::ostream &operator<<(std::ostream &os, const DynamicBitset &obj) {
        os << "DynamicBitset<" << BitsPerElement << ">{";
        for (size_t i = 0; i < obj.numOfElements; ++i) {
            std::bitset<BitsPerElement> bs(obj.iData[i]);
            os << bs;
            if (i < obj.iData.size() - 1) os << ", ";
        }
        os << "}";
        return os;
    }
};


#endif
