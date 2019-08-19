#include "poisson.hpp"

#include <cfdlib/blas.hpp>
#include <cfdlib/solver.hpp>

namespace  {
    inline real square(real v)
    {
        return v*v;
    }
}


DMatrix Poisson::createSystemMatrix(const RegularGrid2d& grid)
{
    // The resulting matrix is a block matrix of (imax * jmax) x (imax * jmax) elements
    // There are two characteristic blocks, each of the size jmax x jmax
    
    size_t jmax = grid.jmax();
    size_t imax = grid.imax();
    
    // First is the one on the diagonal
    DMatrix center_mat = matrix::diag(jmax, 2./square(grid.dx()) + 2./square(grid.dy()));

    real off_diag_value = - 1. / square(grid.dy());
    for (int i = 1; i < center_mat.m(); ++i) {
        center_mat(i-1, i) = off_diag_value;
        center_mat(i, i-1) = off_diag_value;
    }

    // Second are two next to the diaginal
    DMatrix off_center_mat = matrix::diag(jmax, -1./square(grid.dx()));

    // Assemble the final matrix
    DMatrix result = DMatrix(imax * jmax, imax * jmax);
    for (size_t i = 0; i < result.m(); i+=jmax)
    {
        matrix::Block(result, i, i, jmax, jmax) = center_mat;
        if (i + jmax < result.m()) {
            matrix::Block(result, i + jmax, i, jmax, jmax) = off_center_mat;
            matrix::Block(result, i, i + jmax, jmax, jmax) = off_center_mat;
        }
    }

    return result;
}

real Poisson::solve_sor(const RegularGrid2d &grid, const DMatrix &rhs, DMatrix& sol, const BoundaryCondition& bc, real omega, real epsilon, unsigned int itermax)
{
    if (!grid.matches(rhs)) {
        throw std::runtime_error("RHS matrix does not match grid"); 
    }
    
    if (!grid.matches(sol)) {
        throw std::runtime_error("Solution matrix does not match grid"); 
    }
    
    DMatrix p = ghostCellsAdded(sol);

    unsigned int iter = 0;
    real error = 0.;
    
    // coefficients of the poisson problem
    const real c1 = 1./square(grid.dx());
    const real c2 = 1./square(grid.dy());
    const real c3 = -2./square(grid.dx()) - 2./square(grid.dy());
    
    auto local_residuum = [&](const size_t& i, const size_t& j)
    {
        // make sure to use the same numerics as the generic sor solver (i.e. same addition order)
        const real lhs_ij = c1 * p(i-1, j)
                          + c2 * p(i, j-1) + c3 * p(i,j)  + c2 * p(i, j+1)
                          + c1 * p(i+1, j);
        
        // Residuum for the current grid cell
        return rhs(i - 1,j - 1) - lhs_ij;
    };
    
    Timer t;
    t.start();

    do
    {
        error = 0.;
        bc.apply(p);

        for (size_t i = 1; i < p.m() - 1; ++i) {
            for (size_t j = 1; j < p.n() - 1; ++j) {
                // Residuum for the current grid cell
                const real res_ij = local_residuum(i, j);

                p(i, j) += omega/c3 * res_ij;
                // compute l2 norm
                error += square(res_ij);
            }
        }
        error = sqrt(error);
        notify_iteration_finished(iter, error);
        iter++;
    } while(iter < itermax && error > epsilon);
    
    // Compute final residuum
    error = 0.;
    bc.apply(p);
    for (size_t i = 1; i < p.m() - 1; ++i) {
        for (size_t j = 1; j < p.n() - 1; ++j) {
            // Residuum for the current grid cell
            error += square(local_residuum(i, j));
        }
    }
    error = sqrt(error);

    t.stop();
    notify_finished(iter, error, t.milliseconds());

    // Remove ghost cells from solution
    sol = ghostCellsRemoved(p);
    return error;
}

real RegularGrid2d::dx() const
{
    if (!m_evalOnCellCenter) {
        return m_lenx / static_cast<real>(m_imax + 1);
    }
    else {
        return m_lenx/static_cast<real>(m_imax);
    }
}

real RegularGrid2d::dy() const
{
    if (!m_evalOnCellCenter) {
        return m_leny / static_cast<real>(m_jmax + 1);
    }
    else {
        return m_leny / static_cast<real>(m_jmax);
    }
}



void HomogenousNeumannGridCenter::apply(DMatrix &m) const
{
    // first and last row
    for (size_t i = 1; i < m.m() - 1; ++i) {
        m(i, 0) = m(i, 1);
        m(i, m.n() - 1) = m(i, m.n() - 2);
    }
    
    // first and last col
    for (size_t j = 1; j < m.n() -1; ++j) {
        m(0, j) = m(1, j);
        m(m.m() - 1, j) = m(m.m() - 2, j);
    }
}

void HomogenousDirichletGridCenter::apply(DMatrix &m) const
{
    // first and last row
    for (size_t i = 1; i < m.m() - 1; ++i) {
        m(i, 0) = -m(i, 1);
        m(i, m.n() - 1) = -m(i, m.n() - 2);
    }
    
    // first and last col
    for (size_t j = 1; j < m.n() -1; ++j) {
        m(0, j) = -m(1, j);
        m(m.m() - 1, j) = -m(m.m() - 2, j);
    }
}

void HomogenousDirichletGridCorner::apply(DMatrix &m) const
{
    // first and last row
    for (size_t i = 1; i < m.m() - 1; ++i) {
        m(i, 0) = 0.;
        m(i, m.n() - 1) = 0.;
    }
    
    // first and last col
    for (size_t j = 1; j < m.n() -1; ++j) {
        m(0, j) = 0.;
        m(m.m() - 1, j) = 0.;
    }
}

DMatrix Poisson::ghostCellsAdded(const DMatrix &m)
{
    DMatrix result(m.m() + 2, m.n() + 2, 0.);
    matrix::Block(result, 1, 1, m.m(), m.n()) = m;
    return result;
}

DMatrix Poisson::ghostCellsRemoved(const DMatrix &m)
{
    return matrix::Block(const_cast<DMatrix&>(m), 1, 1, m.m() - 2, m.n() - 2).matrix();
}
