#ifndef RANGE_H
#define RANGE_H

#include <iostream>

class Range
{
public:
    Range();
    Range(int rank, int nDecomposedElements, int nElementsTotal);

    int firstElement() {
        return m_firstElement;
    }
    int lastElement() {
        return m_lastElement;
    }
    int nElements() {
        return m_nElements;
    }

    friend std::ostream& operator<<(std::ostream& os, const Range& dt);

private:
    int m_firstElement;
    int m_lastElement;
    int m_nElements;
};

#endif // RANGE_H
