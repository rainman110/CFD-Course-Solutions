#include "poisson.hpp"

#include <cfdlib/blas.hpp>

namespace  {
    inline real square(real v)
    {
        return v*v;
    }
}

Poisson::Poisson(real lenx, real leny, size_t imax, size_t jmax)
    : m_lenx(lenx)
    , m_leny(leny)
    , m_imax(imax)
    , m_jmax (jmax)
{
    
}

DMatrix Poisson::createMatrix() const
{
    // The resulting matrix is a block matrix of (imax * jmax) x (imax * jmax) elements
    // There are two characteristic blocks, each of the size jmax x jmax
    
    // First is the one on the diagonal
    DMatrix center_mat = matrix::diag(m_jmax, 2./square(dx()) + 2./square(dy()));

    real off_diag_value = - 1. / square(dy());
    for (int i = 1; i < center_mat.m(); ++i) {
        center_mat(i-1, i) = off_diag_value;
        center_mat(i, i-1) = off_diag_value;
    }

    // Second are two next to the diaginal
    DMatrix off_center_mat = matrix::diag(m_jmax, -1./square(dx()));

    // Assemble the final matrix
    DMatrix result = DMatrix(m_imax * m_jmax, m_imax * m_jmax);
    for (size_t i = 0; i < result.m(); i+=m_jmax)
    {
        matrix::Block(result, i, i, m_jmax, m_jmax) = center_mat;
        if (i + m_jmax < result.m()) {
            matrix::Block(result, i + m_jmax, i, m_jmax, m_jmax) = off_center_mat;
            matrix::Block(result, i, i + m_jmax, m_jmax, m_jmax) = off_center_mat;
        }
    }

    return result;
}

real Poisson::dx() const
{
    return m_lenx / static_cast<real>(m_imax + 1);
}

real Poisson::dy() const
{
    return m_leny / static_cast<real>(m_jmax + 1);
}


