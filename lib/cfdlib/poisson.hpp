#ifndef POISSON_H
#define POISSON_H

#include <cfdlib/real.hpp>
#include <cfdlib/fields.hpp>

enum class SampleType
{
    OnCellCorner,
    OnCellCenter
};

class RegularGrid2d
{
public:
    /**
     * @brief Creates a regular grid where the dicretization points are
     * located on the cell corners
     * @param lenx X-length of the grid
     * @param leny X-length of the grid
     * @param imax Number of inner cell corners in x direction
     * @param jmax Number of inner cell corners in y direction
     */
    RegularGrid2d(real lenx, real leny, size_t imax, size_t jmax, SampleType type = SampleType::OnCellCorner)
        : m_lenx(lenx)
        , m_leny(leny)
        , m_imax(imax)
        , m_jmax(jmax)
        , m_evalOnCellCenter(type == SampleType::OnCellCenter)
    {}

    template <typename Function2P>
    DMatrix sampleFunction(Function2P func) const
    {
        real offset = m_evalOnCellCenter ?  0.5 : 1;
        
        DMatrix result(m_imax, m_jmax);
        for (size_t i = 0; i < m_imax; ++i) {
            real xp = x(i);
            for (size_t j = 0; j < m_jmax; ++j) {
                real yp = y(j);
                result(i, j) = func(xp, yp);
            }
        }

        return result;
    }

    real x(size_t idx) const
    {
       real offset = m_evalOnCellCenter ?  0.5 : 1;
       return static_cast<real>(idx + offset)*dx();
    }

    real y(size_t idx) const
    {
        real offset = m_evalOnCellCenter ?  0.5 : 1;
        return static_cast<real>(idx + offset)*dy();
    }

    real dx() const;
    real dy() const;
    
    size_t imax() const
    {
        return m_imax;
    }

    size_t jmax() const
    {
        return m_jmax;
    }
    
    bool matches(const DMatrix& m) const
    {
        return m_imax == m.m() && m_jmax == m.n();
    }

private:
    real m_lenx, m_leny;
    size_t m_imax, m_jmax;
    bool m_evalOnCellCenter;
};

/**
 * @brief Interface class to apply boundary conditions 
 * on a grid
 */
class BoundaryCondition
{
public:
    virtual void apply(DMatrix& m) const = 0;
};

class Poisson
{
public:

    /// Creates the linear system matrix of the poisson problem with Dirichlect
    /// boundary conditions p = 0.
    /// This can be used with an ordinary linear solver
    static DMatrix createSystemMatrix(const RegularGrid2d &grid);

    
    /** @brief Efficient sor solver for the possion problem
     *
     * This exploits the sparsity of the system matrix which should be
     * much faster than using the upper approach
     * 
     * @param grid       Grid parameters
     * @param rhs [in]   Right hand side of the poisson problem
     * @param p [in/out] The solution matrix
     * @param bc         Boundary conditions that are applied before each iteration
     * @param omega      SOR relaxation parameter
     * @param epsilon    L2-Error tolerance of the residuum
     * @param itermax    Maximum allowed iterations
     * 
     * @returns L2-Norm of the residuum of the solution
     */
    static real solve_sor(const RegularGrid2d& grid, const DMatrix& rhs, DMatrix &p,
                          const BoundaryCondition& bc,
                          real omega, real epsilon, unsigned int itermax);
    
    static DMatrix ghostCellsAdded(const DMatrix& m);
    static DMatrix ghostCellsRemoved(const DMatrix &m);
};


class HomogenousNeumannGridCenter : public BoundaryCondition
{
public:
    void apply(DMatrix& m) const override;
};

class HomogenousDirichletGridCorner : public BoundaryCondition
{
public:
    void apply(DMatrix& m) const override;
};

class HomogenousDirichletGridCenter : public BoundaryCondition
{
public:
    void apply(DMatrix& m) const override;
};


#endif // POISSON_H
