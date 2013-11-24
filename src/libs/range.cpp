#include "range.h"

/*!
  * Dummy range object with no data.
  */
Range::Range() :
    m_firstElement(0),
    m_lastElement(0),
    m_nElements(0)
{
}

/*!
  * Sets up a range object by decomposing the number of elements to the
  * given number of decomposed elements (number of processors).
  */
Range::Range(int rank, int nDecomposedElements, int nElementsTotal)
{
    // Divide the rows evenly
    m_firstElement = (int)((rank * nElementsTotal) / nDecomposedElements);
    m_lastElement = (int)(((rank + 1) * nElementsTotal) / nDecomposedElements - 1); // our last pixel is the one before the next process' first pixel

    m_nElements = m_lastElement - m_firstElement + 1;
}

std::ostream& operator <<(std::ostream &out, const Range &range)
{
    out << "From " << range.m_firstElement << " to " << range.m_lastElement << " (total: " << range.m_nElements << ")";
    return out;
}
