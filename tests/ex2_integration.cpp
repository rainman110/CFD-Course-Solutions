#define _USE_MATH_DEFINES
#include <cmath>

#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/blas.hpp>
#include <cfdlib/poisson.hpp>
#include <cfdlib/solver.hpp>

TEST(Ex2, 1_createPoissonMatrix)
{
    Poisson p(0.5, 1.5, 3, 2);
    DMatrix pm = p.createMatrix();
    
    real tol = 1e-10;
    // 1. row
    EXPECT_NEAR(136., pm(0,0), tol);
    EXPECT_NEAR( -4., pm(0,1), tol);
    EXPECT_NEAR(-64., pm(0,2), tol);
    EXPECT_NEAR(  0., pm(0,3), tol);
    EXPECT_NEAR(  0., pm(0,4), tol);
    EXPECT_NEAR(  0., pm(0,5), tol);
    
    // 2. row
    EXPECT_NEAR( -4., pm(1,0), tol);
    EXPECT_NEAR(136., pm(1,1), tol);
    EXPECT_NEAR(  0., pm(1,2), tol);
    EXPECT_NEAR(-64., pm(1,3), tol);
    EXPECT_NEAR(  0., pm(1,4), tol);
    EXPECT_NEAR(  0., pm(1,5), tol);
    
    // 3. row
    EXPECT_NEAR(-64., pm(2,0), tol);
    EXPECT_NEAR(  0., pm(2,1), tol);
    EXPECT_NEAR(136., pm(2,2), tol);
    EXPECT_NEAR( -4., pm(2,3), tol);
    EXPECT_NEAR(-64., pm(2,4), tol);
    EXPECT_NEAR(  0., pm(2,5), tol);
    
    // 4. row
    EXPECT_NEAR(  0., pm(3,0), tol);
    EXPECT_NEAR(-64., pm(3,1), tol);
    EXPECT_NEAR( -4., pm(3,2), tol);
    EXPECT_NEAR(136., pm(3,3), tol);
    EXPECT_NEAR(  0., pm(3,4), tol);
    EXPECT_NEAR(-64., pm(3,5), tol);
    
    // 5. row
    EXPECT_NEAR(  0., pm(4,0), tol);
    EXPECT_NEAR(  0., pm(4,1), tol);
    EXPECT_NEAR(-64., pm(4,2), tol);
    EXPECT_NEAR(  0., pm(4,3), tol);
    EXPECT_NEAR(136., pm(4,4), tol);
    EXPECT_NEAR( -4., pm(4,5), tol);
    
    // 6. row
    EXPECT_NEAR(  0., pm(5,0), tol);
    EXPECT_NEAR(  0., pm(5,1), tol);
    EXPECT_NEAR(  0., pm(5,2), tol);
    EXPECT_NEAR(-64., pm(5,3), tol);
    EXPECT_NEAR( -4., pm(5,4), tol);
    EXPECT_NEAR(136., pm(5,5), tol);
}

TEST(Ex2, 2_solveExampleProblem)
{
    Poisson p(1., 1., 20, 20);
    
    DMatrix m = p.createMatrix();
    DVector b(m.m(), 0.);
    DVector x(m.m(), 1.);
    
    solve_sor(m, b, x, 1.5, 1e-10, 1000);
    
    real tol = 1e-8;
    for (size_t i = 0; i < x.n(); ++i ) {
        EXPECT_NEAR(0.0, x(i), tol);
    }
}

TEST(Ex2, 3_solveExampleProblemWithRHS)
{
    std::vector<size_t> grid_sizes({10, 20, 30, 40});
    
    for (auto grid_size : grid_sizes) {

        Poisson p(1., 1., grid_size, grid_size);
    
        DMatrix m = p.createMatrix();
        DVector x(m.m(), 1.);
        
        // Lets define the right hand side function and sample it
        DMatrix rhs = p.sampleFunction([](real x, real y) {
            return 8 * M_PI * M_PI * sin(2. * M_PI * x)*sin(2. * M_PI * y);
        });
    
        // convert into vector
        DVector b = matrix::flatten(rhs);
    
        // solve the problem
        solve_sor(m, b, x, 1.5, 1e-10, 1000);
    
        // compute the analytical solution
        DMatrix anal_solution = p.sampleFunction([](real x, real y) {
            return sin(2. * M_PI * x)*sin(2. * M_PI * y);
        });
     
        // convert into vector
        DVector x_anal = matrix::flatten(anal_solution);
    
        // l2 norm
        blas::axpy(-1., x, x_anal);
        double norm_error = blas::nrm2(x_anal);
        
        std::cout << "nrm2(num_sol - ana_sol) for grid size " << grid_size << ": " << norm_error << std::endl;
        EXPECT_LE(norm_error, 0.2);
    }
}
