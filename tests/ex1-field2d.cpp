#include <gtest/gtest.h>

#include <cfdlib/fields.hpp>
#include <cfdlib/functions.hpp>

class Field2DTests : public testing::Test
{
public:
    Field2DTests()
        : x(2,3, 1.0), y(2,3, 2.0)
    {}
    
    Field2D<double> x, y;
};

TEST_F(Field2DTests, SetValue)
{
    x(1,2) = 10.;
    
    EXPECT_NEAR(1.0, x(0,0), 1e-15);
    EXPECT_NEAR(1.0, x(0,1), 1e-15);
    EXPECT_NEAR(1.0, x(0,2), 1e-15);
    
    EXPECT_NEAR(1.0, x(1,0), 1e-15);
    EXPECT_NEAR(1.0, x(1,1), 1e-15);
    EXPECT_NEAR(10., x(1,2), 1e-15);
}


TEST_F(Field2DTests, Size)
{
    EXPECT_EQ(2, x.m());
    EXPECT_EQ(3, x.n());
}

TEST_F(Field2DTests, ApplyFunctionToField)
{
    x(1,2) = 2.;
    apply_function_to_field(x, [](double v) {
        return v+10.; 
    });
    
    EXPECT_NEAR(11.0, x(0,0), 1e-15);
    EXPECT_NEAR(11.0, x(0,1), 1e-15);
    EXPECT_NEAR(11.0, x(0,2), 1e-15);
    
    EXPECT_NEAR(11.0, x(1,0), 1e-15);
    EXPECT_NEAR(11.0, x(1,1), 1e-15);
    EXPECT_NEAR(12.0, x(1,2), 1e-15);
}
