#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/functions.hpp>

class Field1DTests : public testing::Test
{
public:
    Field1DTests()
        : x(3, 1.0), y(3, 2.0)
    {}
    
    Field1D<double> x, y;
};

TEST_F(Field1DTests, Constructor)
{
    EXPECT_NEAR(1.0, x(0), 1e-15);
    EXPECT_NEAR(1.0, x(1), 1e-15);
    EXPECT_NEAR(1.0, x(2), 1e-15);
    
    EXPECT_NEAR(2.0, y(0), 1e-15);
    EXPECT_NEAR(2.0, y(1), 1e-15);
    EXPECT_NEAR(2.0, y(2), 1e-15);
}

TEST_F(Field1DTests, Fill)
{
    x.fill(10.);
    
    EXPECT_NEAR(10.0, x(0), 1e-15);
    EXPECT_NEAR(10.0, x(1), 1e-15);
    EXPECT_NEAR(10.0, x(2), 1e-15);
}

TEST_F(Field1DTests, SetValue)
{
    x(1) = 10.;

    EXPECT_NEAR(1.0, x(0), 1e-15);
    EXPECT_NEAR(10.0, x(1), 1e-15);
    EXPECT_NEAR(1.0, x(2), 1e-15);
}

TEST_F(Field1DTests, Size)
{
    EXPECT_EQ(3, x.n());
}

TEST_F(Field1DTests, ApplyFunctionToField)
{
    apply_function_to_field(x, [](double v) {
        return v+10.; 
    });
    
    EXPECT_NEAR(11.0, x(0), 1e-15);
    EXPECT_NEAR(11.0, x(1), 1e-15);
    EXPECT_NEAR(11.0, x(2), 1e-15);
    
    apply_function_to_field(x, [](double v) {
        return v*v; 
    });
    
    EXPECT_NEAR(121.0, x(0), 1e-15);
    EXPECT_NEAR(121.0, x(1), 1e-15);
    EXPECT_NEAR(121.0, x(2), 1e-15);
}
