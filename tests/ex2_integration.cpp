#define _USE_MATH_DEFINES
#include <cmath>

#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/blas.hpp>
#include <cfdlib/poisson.hpp>

TEST(Poisson, createMatrix)
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
