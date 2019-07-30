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
    
    real dx() const;
    real dy() const;
    
    
    
private:
    real m_lenx, m_leny;
    size_t m_imax, m_jmax;
};

#endif // POISSON_H
