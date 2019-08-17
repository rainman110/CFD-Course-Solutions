#define _USE_MATH_DEFINES
#include <cmath>

#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/blas.hpp>
#include <cfdlib/poisson.hpp>
#include <cfdlib/solver.hpp>

TEST(Ex2, 1_createPoissonMatrix)
{
    RegularGrid2d g(0.5, 1.5, 3, 2);
    DMatrix pm = Poisson::createSystemMatrix(g);
    
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
    RegularGrid2d g(1., 1., 20, 20);
    
    DMatrix m = Poisson::createSystemMatrix(g);
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

        RegularGrid2d g(1., 1., grid_size, grid_size);
        
        DMatrix m = Poisson::createSystemMatrix(g);
        DVector x(m.m(), 1.);
        
        // Lets define the right hand side function and sample it
        DMatrix rhs = g.sampleFunction([](real x, real y) {
            return 8 * M_PI * M_PI * sin(2. * M_PI * x)*sin(2. * M_PI * y);
        });
    
        // convert into vector
        DVector b = matrix::flatten(rhs);
    
        // solve the problem
        solve_sor(m, b, x, 1.5, 1e-10, 1000);
    
        // compute the analytical solution
        DMatrix anal_solution = g.sampleFunction([](real x, real y) {
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

TEST(Ex2, compareSolvers)
{
    // Here we use a grid of 30x30 as this will demonstrate
    // the already large runtime difference of both solvers
    RegularGrid2d g(1., 1., 30, 30);
    DMatrix rhs = g.sampleFunction([](real x, real y) {
        return 8 * M_PI * M_PI * sin(2. * M_PI * x)*sin(2. * M_PI * y);
    });
    // compute the analytical solution
    DMatrix anal_solution = g.sampleFunction([](real x, real y) {
        return -sin(2. * M_PI * x)*sin(2. * M_PI * y);
    });
    
    DMatrix m = Poisson::createSystemMatrix(g);
    // we need to negate the matrix since the effienct poisson solver does not
    // use the minus sign (as demanded in the excercise)
    blas::scal(-1., m);
    DVector x(m.m(), -1.);

    // solve the problem
    std::cout << "Solving using generic (slow) sor solver" << std::endl;
    solve_sor(m, matrix::flatten(rhs), x, 1.5, 1e-10, 1000);
    DVector diff_solution_2 = matrix::flatten(anal_solution);
    blas::axpy(-1., x, diff_solution_2);
    std::cout << "nrm2(num_sol - ana_sol) for generic solver: " << blas::nrm2(diff_solution_2) << std::endl;

    DMatrix sol(anal_solution.m(), anal_solution.n(), -1.);

    std::cout << "Solving using efficient (fast) sor solver" << std::endl;
    Poisson::solve_sor(g, rhs, sol, HomogenousDirichletGridCorner(), 1.5, 1e-10, 1000);

    // l2 norm
    DMatrix diff_solution_1 = anal_solution;
    blas::axpy(-1., sol, diff_solution_1);

    std::cout << "nrm2(num_sol - ana_sol) for efficient solver: " << blas::nrm2(matrix::flatten(diff_solution_1)) << std::endl;

    blas::axpy(-1., matrix::flatten(sol), x);
    double norm_error =  blas::nrm2(x);
    std::cout << "nrm2(num_sol1 - num_sol2): " << norm_error << std::endl;

    // check this out... we have truely zero difference!
    EXPECT_EQ(norm_error, 0.);
}
