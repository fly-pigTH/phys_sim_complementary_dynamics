#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_MATH
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_MATH

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {

const real ToReal(const float value);
const real ToReal(const double value);
const real ToReal(const integer value);
const double ToDouble(const real value);
const float ToFloat(const real value);

const real Pi();

// |a - b| <= abs_tol + rel_tol * |b|.
const bool IsClose(const real a, const real b, const real abs_tol,
    const real rel_tol);

// For all i, |ai - bi| <= abs_tol + rel_tol * |bi|.
const bool IsClose(const VectorXr& a, const VectorXr& b, const real abs_tol,
    const real rel_tol);

// Build a local frame.
const Matrix3r BuildFrameFromUnitNormal(const Vector3r& n);

const Matrix3r ToCrossProductMatrix(const Vector3r& a);

const Matrix3r ToRotationMatrix(const Vector3r& angle);

// Sparse matrices.
const SparseMatrixXr FromTriplet(const integer row_num, const integer col_num,
    const std::vector<Eigen::Triplet<real>>& nonzeros);
const std::vector<Eigen::Triplet<real>> ToTriplet(const SparseMatrixXr& mat);
const SparseMatrixXr FromDiagonal(const VectorXr& diagonal);
const SparseMatrixXr Transpose(const SparseMatrixXr& mat);

// Rotation matrix.
const bool IsRotationMatrix(const Matrix3r& R);

// Polar decomposition.
void PolarDecomposition(const Matrix3r& F, Matrix3r& R, Matrix3r& S);
const std::pair<Matrix3r, Matrix3r> PolarDecompositionDifferential(
    const Matrix3r& F, const Matrix3r& R, const Matrix3r& S,
    const Matrix3r& dF);
const std::pair<Matrix9r, Matrix9r> PolarDecompositionDifferential(
    const Matrix3r& F, const Matrix3r& R, const Matrix3r& S);

// Determinant gradient and Hessian.
const Matrix3r DeterminantGradient(const Matrix3r& A);
const Matrix9r DeterminantHessian(const Matrix3r& A);

// Check gradients.
const bool CheckGradient(const VectorXr& x0,
    const std::function<const real(const VectorXr&)>& func,
    const std::function<const VectorXr(const VectorXr&)>& func_grad,
    const Options& opt);

const bool CheckGradient(const VectorXr& x0,
    const std::function<const real(const VectorXr&)>& func,
    const std::function<const VectorXr(const VectorXr&)>& func_grad,
    const std::function<const bool(const integer)>& skip_x_dof,
    const Options& opt);

const bool CheckJacobian(const VectorXr& x0,
    const std::function<const VectorXr(const VectorXr&)>& func,
    const std::function<const MatrixXr(const VectorXr&)>& func_jac,
    const Options& opt);

const bool CheckJacobian(const VectorXr& x0,
    const std::function<const VectorXr(const VectorXr&)>& func,
    const std::function<const MatrixXr(const VectorXr&)>& func_jac,
    const std::function<const bool(const integer)>& skip_x_dof,
    const std::function<const bool(const integer)>& skip_func_dof,
    const Options& opt);

}

#endif