#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <vector>
#include <algorithm>
#include <exception>
#include <ostream>
#include <fstream>

#include <cfdlib/real.hpp>

template <typename T>
class Field1D
{
public:
    Field1D(size_t n)
        : m_data(n), m_n(n)
    {}
    
    Field1D(size_t n, const T& fill_value)
        : m_data(n, fill_value), m_n(n)
    {}
    
    /// Returns the size of the field (Number of elements)
    size_t n() const
    {
        return m_n;
    }
    
    /// Sets all values of the field to fill_value
    void fill(const T& fill_value)
    {
        std::fill(std::begin(m_data), std::end(m_data), fill_value);
    }

    /// Read only access to the i-th value of the field
    /// Usage: std::cout << x(3);
    const T& operator() (size_t i) const
    {
        // row major order
        if (i < 0 || i >= m_n) {
            throw std::out_of_range("Index i out of range in Field1D(i)");
        }
        
        return m_data[i];
    }

    /// Read and write access to the i-th value of the field
    /// Usage: std::cout << x(3); x(2) = 5.;
    T& operator() (size_t i)
    {
        return const_cast<T&>(const_cast<const Field1D*>(this)->operator()(i));
    }
    
private:
    std::vector<T> m_data;
    size_t m_n;
};

template <typename T>
class Field2D
{
public:
    Field2D(size_t m, size_t n)
        : m_data(m*n), m_m(m), m_n(n)
    {}

    Field2D(size_t m, size_t n, const T& fill_value)
        : m_data(m*n, fill_value), m_m(m), m_n(n)
    {}

    /// Returns the number of rows of the field
    size_t m() const
    {
        return m_m;
    }

    /// Returns the number of columns of the field
    size_t n() const
    {
        return m_n;
    }
    
    /// Sets all values of the field to fill_value
    void fill(const T& fill_value)
    {
        std::fill(std::begin(m_data), std::end(m_data), fill_value);
    }
    
    /// Read only access to the (i,j)-th value of the field
    /// Usage: std::cout << x(3,4);
    const T& operator() (size_t i, size_t j) const
    {
        // row major order
        if (i < 0 || i >= m_m) {
            throw std::out_of_range("Index i out of range in Field2D(i,j)");
        }
        if (j < 0 || j >= m_n) {
            throw std::out_of_range("Index j out of range in Field2D(i,j)");
        }

        return m_data[i * m_n + j];
    }
    
    /// Read and write access to the (i,j)-th value of the field
    /// Usage: std::cout << x(3,4); x(2,1) = 5.;
    T& operator() (size_t i, size_t j)
    {
        return const_cast<T&>(const_cast<const Field2D*>(this)->operator()(i, j));
    }
    
    void reshape(size_t m, size_t n)
    {
        if (m * n != m_m * m_n) {
            throw std::runtime_error("Dimension error in Field2d::reshape");
        }
        
        m_m = m;
        m_n = n;
    }

private:
    // 1D storage of 2D field
    std::vector<T> m_data;
    size_t m_m, m_n;
};

template <typename T>
bool is_equal(const Field2D<T>& f1, const Field2D<T>& f2, double eps)
{
    if (f1.m() != f2.m() || f1.n() != f2.n()) {
        return false;
    }
    
    for (size_t i = 0; i < f1.m(); ++i) {
        for (size_t j = 0; j < f1.n(); ++j) {
            if (std::abs(f1(i,j) - f2(i,j)) > eps) {
                return false;
            }
         }
    }
    
    return true;
}

template <typename T>
bool is_equal(const Field1D<T>& f1, const Field1D<T>& f2, double eps)
{
    if (f1.n() != f2.n()) {
        return false;
    }
    
    for (size_t i = 0; i < f1.n(); ++i) {
        if (std::abs(f1(i) - f2(i)) >= eps) {
            return false;
        }
    }
    
    return true;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Field1D<T>& x)
{
    for (size_t i = 0; i < x.n(); ++i) {
        os << x(i) << std::endl;
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Field2D<T>& x)
{
    for (size_t i = 0; i < x.m(); ++i) {
        for (size_t j = 0; j < x.n(); ++j) {
            os << x(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

/// Stores the field to a file in binary form
/// Note: this is not portable. The result will 
/// be different on varying platforms.
/// 
/// The corresponding read function should be used only on the same system
template <typename T>
void write(const Field1D<T>& x, const std::string& filename)
{
    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    uint32_t vsize = static_cast<uint32_t>(x.n());
    fout.write((char*)&vsize, sizeof (uint32_t));
    fout.write((char*)&x(0), sizeof(T)*x.n());
    fout.close();
}

/// Stores the field to a file in binary form
/// Note: this is not portable. The result will 
/// be different on varying platforms.
/// 
/// The corresponding read function should be used only on the same system
template <typename T>
void write(const Field2D<T>& x, const std::string& filename)
{
    // NOTE: we actually should also write what data type we are writing!!!
    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    uint32_t m = static_cast<uint32_t>(x.m());
    uint32_t n = static_cast<uint32_t>(x.n());
    fout.write((char*)&m, sizeof (uint32_t));
    fout.write((char*)&n, sizeof (uint32_t));
    fout.write((char*)&x(0,0), sizeof(T)*x.n()*x.m());
    fout.close();
}

template <typename T>
Field1D<T> read_field1d(const std::string& filename)
{
    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t n;
    fin.read((char*)&n, sizeof(uint32_t));
    
    Field1D<T> result(n);
    fin.read((char*)&result(0), sizeof(T)*n);
    fin.close();
    return result;
}

template <typename T>
Field2D<T> read_field2d(const std::string& filename)
{
    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t n, m;
    fin.read((char*)&m, sizeof(uint32_t));
    fin.read((char*)&n, sizeof(uint32_t));
    
    Field2D<T> result(m, n);
    fin.read((char*)&result(0,0), sizeof(T)*n*m);
    fin.close();
    return result;
}

using DMatrix = Field2D<real>;
using DVector = Field1D<real>;

namespace matrix
{

    /// Returns the identity matrix
    inline DMatrix eye(size_t n)
    {
        DMatrix result(n, n, 0.);
        for (size_t i = 0; i < n; ++i) {
            result(i,i) = 1.;
        }

        return result;
    }
    
    inline DMatrix diag(size_t n, real value)
    {
        DMatrix result(n, n, 0.);
        for (size_t i = 0; i < n; ++i) {
            result(i,i) = value;
        }
        
        return result;
    }

    template <typename T>
    Field1D<T> flatten(const Field2D<T>& mat)
    {
        auto& mat_tmp = const_cast<Field2D<T>&>(mat);
        
        const size_t m_old = mat.m();
        const size_t n_old = mat.n();
        mat_tmp.reshape(1, m_old*n_old);
        
        Field1D<T> vec(m_old*n_old);
        for (size_t i = 0; i < m_old*n_old; ++i) {
            vec(i) = mat_tmp(0, i);
        }
        mat_tmp.reshape(m_old, n_old);
        
        return vec;
    }

    /// Block-view of a matrix
    class Block {
    public:
        Block(DMatrix& mat, size_t imin, size_t jmin, size_t m, size_t n)
            : m_mat(mat), m_imin(imin), m_jmin(jmin), m_m(m), m_n(n)
        {
            if (imin < 0 || imin >= mat.m() || m < 0 || imin + m > mat.m() ) {
                throw std::runtime_error("Illegal row indices of BlockMatrix");
            }
            if (jmin < 0 || jmin >= mat.n() || n < 0 || jmin + n > mat.n() ) {
                throw std::runtime_error("Illegal row indices of BlockMatrix");
            }
        }
        
        /// Read only access to the (i,j)-th value of the field
        /// Usage: std::cout << x(3,4);
        const real& operator() (size_t i, size_t j) const
        {
            // row major order
            if (i < 0 || i >= m_m) {
                throw std::out_of_range("Index i out of range in Block(i,j)");
            }
            if (j < 0 || j >= m_n) {
                throw std::out_of_range("Index j out of range in Block(i,j)");
            }
            
            return m_mat(i + m_imin, j + m_jmin);
        }
        
        /// Read and write access to the i-th value of the field
        /// Usage: std::cout << x(3,4); x(2, 3) = 5.;
        real& operator() (size_t i, size_t j)
        {
            return const_cast<real&>(const_cast<const Block*>(this)->operator()(i, j));
        }
        
        size_t m() const
        {
            return m_m;
        }
        
        size_t n() const
        {
            return m_n; 
        }
        
        /// Writes the matrix mat into the defined block
        Block& operator= (const DMatrix& mat)
        {
            if (mat.m() != m() || mat.n() != n()) {
                throw std::runtime_error("Cannot assign matrix to block: Dimensions do not match");
            }
            
            for (size_t i = 0; i < m(); ++i) {
                for (size_t j = 0; j < n(); ++j) {
                    this->operator()(i, j) = mat(i,j);
                }
            }

            return *this;
        }
        
        /// Converts the block into a matrix
        DMatrix matrix() const
        {
            DMatrix mat(m(), n());
            
            for (size_t i = 0; i < m(); ++i) {
                for (size_t j = 0; j < n(); ++j) {
                    mat(i,j) = this->operator()(i, j);
                }
            }
            return mat;
        }

     private:
         DMatrix& m_mat;
         size_t m_imin, m_jmin, m_m, m_n;
    };

}

#endif // FIELDS_HPP
