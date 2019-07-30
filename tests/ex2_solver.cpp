#define _USE_MATH_DEFINES
#include <cmath>

#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/solver.hpp>

TEST(Solver, sor_simple)
{
    DMatrix A = matrix::diag(3, 2.);
    A(0,1) = 5.;
    A(1,2) = 4.;
    
    DVector b(3, 0);
    b(0) = 7.;
    b(1) = 6.;
    b(2) = 2.;
    DVector x(3, 0.);
    
    solve_sor(A, b, x, 1.5, 1e-15, 1000);
    EXPECT_NEAR(1., x(0), 1e-10);
    EXPECT_NEAR(1., x(1), 1e-10);
    EXPECT_NEAR(1., x(2), 1e-10);
}
