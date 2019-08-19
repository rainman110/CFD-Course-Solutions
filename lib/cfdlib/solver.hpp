#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <exception>
#include <cmath>
#include <iostream>

#include <cfdlib/real.hpp>
#include <cfdlib/fields.hpp>
#include <cfdlib/blas.hpp>
#include <cfdlib/timer.hpp>

inline void notify_iteration_finished(unsigned iter, real error)
{
    //std::cout << "Iter " << iter << ": " << error << std::endl;
}

inline void notify_finished(unsigned iter, real error, real milliseconds)
{
    std::cout << "Finished after " << iter << " iterations with error: " << error << std::endl;
    std::cout << "Runtime [ms]: " << milliseconds << "\n" << std::endl;
}

/**
 * @brief Solves the linear A*x = b system using the SOR (Successive Over-Relaxation) method
 * 
 * @param A The matrix of the linear system
 * @param b The right hand side
 * @param x The solution vector. Requires a initial solution.
 * @param omega Relaxation parameter of the SOR method in the range (0,2)
 * @param epsilon Error tolerance ot the residuum
 * @param itermax Maximum number of iterations
 * 
 * @returns L2-Norm of the residuum of the solution
 */
inline real solve_sor(const DMatrix& A, const DVector& b, DVector& x, real omega, real epsilon, unsigned int itermax)
{
    if (A.m() != A.n()) {
        throw std::runtime_error("Matrix not square in solve_sor");
    }
    
    if (A.m() != b.n() || b.n() != x.n()) {
        throw std::runtime_error("Matrix and vector sizes don't match");
    }
    
    Timer t;
    t.start();
    unsigned int iter = 0;
    real error = 0.;
    do {
        error = 0.;
        for (size_t k = 0; k < A.m(); ++k) {
            real lhs_k = blas::dot(FieldRow<real>(A, k), x);
            real res_k = b(k) - lhs_k ;

            x(k) += omega/A(k, k) * res_k;
            // Compute l2 norm
            error += res_k*res_k;

        }
        // note: this is the residuum of the previous iteration
        error = sqrt(error);
        notify_iteration_finished(iter, error);

        ++iter;
    } while (error > epsilon && iter < itermax);

    t.stop();
    
    // Compute residuum
    DVector residuum(b.n(), 0.);
    blas::gemv(1., A, x, 0., residuum);
    blas::axpy(-1., b, residuum);
    real norm = blas::nrm2(residuum);

    notify_finished(iter, norm, t.milliseconds());

    return norm;
}

#endif // SOLVER_HPP
