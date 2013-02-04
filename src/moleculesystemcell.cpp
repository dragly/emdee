#include "moleculesystemcell.h"

MoleculeSystemCell::MoleculeSystemCell() :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions))
{
    cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
}

void MoleculeSystemCell::setBoundaries(mat boundaries)
{
    m_boundaries = boundaries;
    irowvec counters = zeros<irowvec>(m_nDimensions);
    for(int i = 0; i < pow3nDimensions; i++) {
        irowvec direction = counters - ones<irowvec>(m_nDimensions);
        rowvec lengths = m_boundaries.row(1) - m_boundaries.row(0);
        rowvec directionVec = conv_to<rowvec>::from(direction);
        cellShiftVectors.row(i) = lengths % directionVec;
        counters(0) += 1;
        for(uint idim = 1; idim < counters.size(); idim++) {
            if(counters(idim - 1) > 2) {
                counters(idim - 1) = 0;
                counters(idim) += 1;
            }
        }
    }
    //    cout << cellShiftVectors << endl;
}

ostream& operator<<(ostream& os, MoleculeSystemCell* dt)
{
    os << dt->m_boundaries << endl;
    return os;
}

ostream& operator<<(ostream& os, const MoleculeSystemCell& dt)
{
    os << dt.m_boundaries << endl;
    return os;
}
