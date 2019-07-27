
#include <cfdlib/vtk-io.hpp>

void write_vtk_scalar(const std::string& filename, const std::string& description,
                      const DMatrix& mat, double dx, double dy)
{
        std::ofstream fout(filename, std::ios::out);
        
        fout <<
            "# vtk DataFile Version 3.0\n"
            "Scalar Field\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n";
        fout <<
            "DIMENSIONS " << mat.m() << " " << mat.n() << " 1\n";
        fout <<
            "X_COORDINATES " << mat.m() << " double\n";
        for (int i = 0; i < mat.m(); ++i) {
            fout << dx*i << " ";
        }
        fout << "\n";
        fout <<
            "Y_COORDINATES " << mat.n() << " double\n";
        for (int j = 0; j < mat.n(); ++j) {
            fout << dy*j << " ";
        }
        fout << "\n";
        fout <<
            "Z_COORDINATES 1 double\n"
            "0.0\n";
        
        fout <<
            "POINT_DATA " << mat.m()*mat.n() << "\n"
            "SCALARS " << description.c_str() << " double 1\n"
            "LOOKUP_TABLE default\n";
        
        // Write data
        for (size_t j = 0; j < mat.n(); ++j) {
            for (size_t i = 0; i < mat.m(); ++i) {
                fout << mat(i,j) << std::endl;
            }
        }
        
        fout.close();

}


void write_vtk_vector(const std::string &filename, const std::string &description, const DMatrix &matU, const DMatrix &matV, double dx, double dy)
{
    if (matU.m() != matV.m() || matU.n() != matV.n()) {
        throw std::runtime_error("Sizes of matU and matV don't match.");
    }
    
    std::ofstream fout(filename, std::ios::out);
    
    fout <<
        "# vtk DataFile Version 3.0\n"
        "Vector Field\n"
        "ASCII\n"
        "DATASET RECTILINEAR_GRID\n";
    fout <<
        "DIMENSIONS " << matU.m() << " " << matU.n() << " 1\n";
    fout <<
        "X_COORDINATES " << matU.m() << " double\n";
    for (int i = 0; i < matU.m(); ++i) {
        fout << dx*i << " ";
    }
    fout << "\n";
    fout <<
        "Y_COORDINATES " << matU.n() << " double\n";
    for (int j = 0; j < matU.n(); ++j) {
        fout << dy*j << " ";
    }
    fout << "\n";
    fout <<
        "Z_COORDINATES 1 double\n"
        "0.0\n";
    
    fout <<
        "POINT_DATA " << matU.m()*matU.n() << "\n"
        "VECTORS " << description.c_str() << " double\n";
    
    // Write data
    for (size_t j = 0; j < matU.n(); ++j) {
        for (size_t i = 0; i < matU.m(); ++i) {
            fout << matU(i,j) << " " << matV(i,j) << " 0.0" << std::endl;
        }
    }
    
    fout.close();
}
