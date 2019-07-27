#ifndef VTKIO_HPP
#define VTKIO_HPP

#include <cfdlib/fields.hpp>
#include <iostream>

/**
 * @brief Writes a 2d-field/matrix into a scalar vtk file
 * @param filename    Filename to write the data to
 * @param description Name of the field
 * @param mat         The matrix to write
 * @param dx          X-Cell size
 * @param dy          Y-Cell size
 */
void write_vtk_scalar(const std::string& filename, const std::string& description,
                      const DMatrix& mat, double dx, double dy);

/**
 * @brief Writes a 2d-field/matrix into a vector vtk file
 * @param filename    Filename to write the data to
 * @param description Name of the field
 * @param matU        U Values of the vector field
 * @param matV        V Values of the vector field
 * @param dx          X-Cell size
 * @param dy          Y-Cell size
 */
void write_vtk_vector(const std::string& filename, const std::string& description,
                      const DMatrix& matU, const DMatrix& matV, double dx, double dy);

#endif // VTKIO_HPP
