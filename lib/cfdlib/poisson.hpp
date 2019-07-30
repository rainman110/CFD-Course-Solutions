#ifndef POISSON_H
#define POISSON_H

#include <cfdlib/real.hpp>
#include <cfdlib/fields.hpp>

class Poisson
{
public:
    Poisson(real lenx, real leny, size_t imax, size_t jmax);
    
    /// Creates the system matrix of the poisson problem with Dirichlect
    /// boundary conditions p = 0.
    DMatrix createMatrix() const;
    
    template <typename Function2P>
    DMatrix sampleFunction(Function2P func) const
    {
        DMatrix result(m_imax, m_jmax);
        for (size_t i = 1; i <= m_imax; ++i) {
            real x = static_cast<real>(i)*dx();
            for (size_t j = 1; j <= m_jmax; ++j) {
                real y = static_cast<real>(j)*dy();
                result(i - 1, j - 1) = func(x, y);
            }
        }
        
        return result;
    }
    
    real dx() const;
    real dy() const;
    
    
    
private:
    real m_lenx, m_leny;
    size_t m_imax, m_jmax;
};

#endif // POISSON_H
