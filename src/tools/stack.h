//
// Created by gonciarz on 2019-04-02.
//

#ifndef STACK_H
#define STACK_H

#include <memory>
#include <cassert>

/**
 * Simple implementation of fixed size stack.
 * @tparam ELEMENT_TYPE - type of stored element in a stack
 * @tparam IDX - index type (for few elements stacks it can be int8_t (it must be signed type)
 */
template <typename ELEMENT_TYPE = uint8_t, typename IDX=int>
class Stack {
    static_assert(std::is_signed<IDX>::value);

    IDX iIdx = -1;
    const IDX iSize;
    std::unique_ptr<ELEMENT_TYPE[]> iData;

public:
    /**
     * Constructor
     * @param aSize max size of a stack
     */
    explicit Stack(IDX aSize) : iSize(aSize), iData(new ELEMENT_TYPE[aSize]) {}

    /**
     * @return true if stack is empty
     */
    bool empty() const {
        return iIdx < 0;
    }

    /**
     * @return number of elements in stack
     */
    IDX size() const {
        return iIdx + 1;
    }

    /**
     * clears stack
     */
    void clear() {
        iIdx = -1;
    }

    /**
     * Add element to the stack
     * @param e - reference to element
     */
    void push(const ELEMENT_TYPE &e) {
        assert(iIdx < (iSize - 1) && "Too many elements already in a stack!");
        iData[++iIdx] = e;
    }

    /**
     * Pops and returns last element from stack
     */
    ELEMENT_TYPE pop() {
        assert(iIdx >= 0 && "All elements were pop from a stack!");
        return iData[iIdx--];
    }

    /**
     *  Returns max capacity of stack
     */
    IDX capacity() const {
        return iSize;
    }

    /**
     * string stream output
     */
    friend std::ostream &operator<<(std::ostream &os, const Stack &obj) {
        os << "Stack<" << sizeof(ELEMENT_TYPE) << ">[" << obj.size() << "]{";
        for (IDX i = 0; i < obj.iSize; ++i) {
            os << obj.iData[i];
            if (i < obj.iSize - 1) os << ", ";
        }
        os << "}";
        return os;
    }
};


#endif
