#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <exception>
#include <cmath>
#include <iostream>

#include <cfdlib/real.hpp>
#include <cfdlib/fields.hpp>
#include <cfdlib/blas.hpp>

inline void notify_iteration_finished(unsigned iter, real error)
{
    std::cout << "Iter " << iter << ": " << error << std::endl;
}

inline void solve_sor(const DMatrix& A, const DVector& b, DVector& x, real omega, real epsilon, unsigned int itermax)
{
    if (A.m() != A.n()) {
        throw std::runtime_error("Matrix not square in solve_sor");
    }
    
    if (A.m() != b.n() || b.n() != x.n()) {
        throw std::runtime_error("Matrix and vector sizes don't match");
    }

    unsigned int iter = 0;
    real error = 0.;
    do {
        error = 0.;
        for (size_t k = 0; k < A.m(); ++k) {
            real lhs_k = blas::dot(FieldRow<real>(A, k), x);
            const real dx = omega/A(k, k) * (b(k) - lhs_k);
            
            x(k) += dx;
            error += dx*dx;

        }
        error = sqrt(error);
        notify_iteration_finished(iter, error);

        ++iter;
    } while (error > epsilon && iter < itermax);
}

#endif // SOLVER_HPP
