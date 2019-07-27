#define _USE_MATH_DEFINES
#include <cmath>

#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/functions.hpp>
#include <cfdlib/blas.hpp>
#include <cfdlib/vtk-io.hpp>


TEST(Ex1, Integration_a)
{
    DMatrix a(10, 10);
    
    for (size_t i = 0; i < a.m(); ++i) {
        for (size_t j = 0; j < a.n(); ++j) {
            a(i, j) = static_cast<double>((i+1)*(j+1));
        }
    }
    
    DVector x(10, 1.0);
    DVector y(10, 0.0);
    
    std::cout << "A:" << std::endl << a << std::endl;
    
    // compute y = A*x
    blas::gemv(1.0, a, x, 0.0, y);

    std::cout << "nrm2(A*x): " << blas::nrm2(y) << std::endl;
    EXPECT_NEAR(1079.177928, blas::nrm2(y), 1e-6);
    
    // Compute (5*x + A*x)*(1.5*x) using our blas routines
    blas::axpy(5.0, x, y);
    blas::scal(1.5, x);
    double result = blas::dot(x,y);
    
    std::cout << "(5*x + A*x)*(1.5*x): " << result << std::endl;
    EXPECT_NEAR(4612.5, result, 1e-6);
}

TEST(Ex1, Integration_b)
{
    DMatrix f(60, 40, 0.);
    
    for (size_t i = 0; i < f.m(); ++i) {
        for (size_t j = 0; j < f.n(); ++j) {
            f(i, j) = M_PI / 60. * static_cast<double>(i);
        }
    }

    write_vtk_scalar("../data/gradientField.vtk", "gradient", f, 0.1, 0.1);

    apply_function_to_field(f, [](double v) {
        return sin(v);
    });

    write_vtk_scalar("../data/sinusField.vtk", "sine", f, 0.1, 0.1);
}

TEST(Ex1, Integration_c)
{
    DMatrix fieldU = read_field2d<double>("../data/fieldU.dat");
    DMatrix fieldV = read_field2d<double>("../data/fieldV.dat");

    EXPECT_EQ(60, fieldU.m());
    EXPECT_EQ(60, fieldU.n());
    EXPECT_EQ(60, fieldV.m());
    EXPECT_EQ(60, fieldV.n());

    write_vtk_vector("../data/vectorField.vtk", "vectorfield", fieldU, fieldV, 0.1, 0.1);
}
