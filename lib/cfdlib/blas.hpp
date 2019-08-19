#ifndef BLAS_H
#define BLAS_H

#include <cmath>

#include <cfdlib/fields.hpp>

#include <exception>

/**
 * Simple blas-like linear algebra routines
 */
namespace blas
{

/**
 * Vector addition
 * 
 * y <- a*x + y
 */
template <typename T>
void axpy(const T& alpha, const Field1D<T>& x, Field1D<T>& y)
{
    if (x.n() != y.n()) {
        throw std::length_error("Fields x and y have different size in blas::axpy.");
    }

    for (size_t i = 0; i < x.n(); ++i) {
        y(i) += alpha * x(i);
    }

}

/**
 * Matrix addition
 * 
 * Y <- a*X + Y
 */
template <typename T>
void axpy(const T& alpha, const Field2D<T>& x, Field2D<T>& y)
{
    if (x.m() != y.m() || x.n() != y.n()) {
        throw std::length_error("Fields x and y have different size in blas::axpy.");
    }

    for (size_t i = 0; i < x.m(); ++i) {
        for (size_t j = 0; j < x.n(); ++j) {
            y(i, j) += alpha * x(i, j);
        }
    }
}

/**
 * Scalar (Dot) product of x and y
 * 
 * x and y must have same size 
 */
template <typename VectorLike1, typename VectorLike2>
real dot(const VectorLike1& x, const VectorLike2& y)
{
    if (x.n() != y.n()) {
        throw std::length_error("Fields x and y have different size in blas::dot.");
    }

    real result = 0.;
    for (size_t i = 0; i < x.n(); ++i) {
        result += x(i) * y(i);
    }

    return result;
}


/**
 * Vector 2-norm of vector x
 */
template <typename T>
T nrm2(const Field1D<T>& x)
{
    return static_cast<T>(sqrt(dot(x, x)));
}

/**
 * Copies the value of x to y
 * 
 * Note: y must have the same size as x
 * 
 * y <- x
 */
template <typename T>
void copy(const Field1D<T>& x, Field1D<T>& y)
{
    if (x.n() != y.n()) {
        throw std::length_error("Fields x and y have different size in blas::copy.");
    }

    y = x;
}

/**
 * Scales the 1d field x by the factor alpha
 * 
 * x <- a * x
 */
template <typename T>
void scal(const T& alpha, Field1D<T>& x)
{
    for (size_t i = 0; i < x.n(); ++i) {
        x(i) *= alpha;
    }
}

/**
 * Scales the 2d field x by the factor alpha
 * 
 * x <- a * x
 */
template <typename T>
void scal(const T& alpha, Field2D<T>& x)
{
    for (size_t i = 0; i < x.m(); ++i) {
        for (size_t j = 0; j < x.n(); ++j) {
            x(i, j) *= alpha;
        }
    }
}

/**
 * General Matrix Vector Product
 * 
 * Computes y <- a A*x + b y
 */
template <typename T>
void gemv(const T& alpha, const Field2D<T>& A, const Field1D<T>& x, const T& beta, Field1D<T>& y)
{
    if (A.m() != y.n()) {
        throw std::length_error("A.n and y have different size in blas::gemv.");
    }
    
    if (A.n() != x.n()) {
        throw std::length_error("A.m() and x  have different size in blas::gemv.");
    }
    
    for (size_t i=0; i < y.n(); ++i) {
        y(i) *= beta;
        for (size_t j = 0; j < A.n(); ++j) {
            y(i) += alpha*A(i,j)*x(j);
        }
    }
}

} // namespace blas

#endif // BLAS_H
