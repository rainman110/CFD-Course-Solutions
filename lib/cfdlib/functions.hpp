#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <cfdlib/fields.hpp>

/// Elementwise application of a 1-parm function to the field
template<typename T, typename Function1P>
void apply_function_to_field(Field1D<T>& field, Function1P func)
{
    for (size_t i = 0; i < field.n(); ++i) {
        field(i) = func(field(i));
    }
}

/// Elementwise application of a 1-parm function to the field
template<typename T, typename Function1P>
void apply_function_to_field(Field2D<T>& field, Function1P func)
{
    for (size_t i = 0; i < field.m(); ++i) {
        for (size_t j = 0; j < field.n(); ++j) {
            field(i,j) = func(field(i,j));
        }
    }
}

#endif // FUNCTIONS_HPP
