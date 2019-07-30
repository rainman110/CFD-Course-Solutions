#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <exception>
#include <cmath>
#include <iostream>

#include <cfdlib/real.hpp>
#include <cfdlib/fields.hpp>

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
    DVector xnext(x.n(), 0.);
    do {
        error = 0.;
        for (size_t k = 0; k < A.m(); ++k) {
            real t1 = 0.;
            for (size_t i = 0; i < k; ++i) {
                t1 += A(k, i)*xnext(i);
            }
            real t2 = 0.;
            for (size_t i = k+1; i < A.m(); ++i) {
                t2 += A(k, i)*x(i);
            }

            xnext(k) = (1. - omega)*x(k) + omega/A(k, k) * (b(k) - t1 - t2);
            error = std::max(error, fabs(xnext(k) - x(k)));
        }
        x = xnext;

        notify_iteration_finished(iter, error);

        ++iter;
    } while (error > epsilon && iter < itermax);
}

#endif // SOLVER_HPP
