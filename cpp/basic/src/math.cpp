#include "basic/include/math.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {

const real ToReal(const float value) {
    return static_cast<real>(value);
}

const real ToReal(const double value) {
    return static_cast<real>(value);
}

const real ToReal(const integer value) {
    return static_cast<real>(value);
}

const double ToDouble(const real value) {
    return static_cast<double>(value);
}

const float ToFloat(const real value) {
    return static_cast<float>(value);
}

const real Pi() {
    return ToReal(3.141592653589793238462643383);
}

const bool IsClose(const real a, const real b, const real abs_tol,
    const real rel_tol) {
    return std::abs(a - b) <= std::abs(b) * rel_tol + abs_tol;
}

const bool IsClose(const VectorXr& a, const VectorXr& b, const real abs_tol,
    const real rel_tol) {
    CheckCondition(a.size() == b.size(), "basic::math::IsClose",
        "Sizes of two vectors are different.");
    return ((a - b).cwiseAbs().array() <= b.cwiseAbs().array() * rel_tol
        + abs_tol).all();
}

const Matrix3r BuildFrameFromUnitNormal(const Vector3r& n) {
    real max_len = -1;
    Matrix3r R; R.setZero();
    R.col(2) = n;
    for (integer i = 0; i < 3; ++i) {
        const Vector3r x = n.cross(Vector3r::Unit(i));
        const real x_len = x.norm();
        if (x_len > max_len) {
            max_len = x_len;
            R.col(0) = x / x_len;
        }
    }
    R.col(1) = n.cross(R.col(0));
    return R;
}

const Matrix3r ToCrossProductMatrix(const Vector3r& a) {
    Matrix3r A;
    A << 0, -a(2), a(1),
        a(2), 0, -a(0),
        -a(1), a(0), 0;
    return A;
}

const Matrix3r ToRotationMatrix(const Vector3r& angle) {
    // Define
    //     a_0 = cos(theta).
    //     a_1 = sin(theta) / theta.
    //     a_2 = (1 - cos(theta)) / theta^2.
    // Then, R = I + a_1(theta) A + a_2(theta) A^2. Here A = [axis_].
    // To see this, consider any vector v:
    // Rv = v + sin(theta) * unit_axis x v +
    //      (1 - cos(theta)) unit_axis x (unit_axis x v).
    const real theta = angle.stableNorm();
    const real theta_sqr = theta * theta;
    const real a0 = std::cos(theta);
    // 1 rad = 57 degrees.
    // 2e-3 rad = 0.1 degrees, which should be sufficiently small for motions
    // we care about.
    const real tol = ToReal(2e-3);
    real a1 = 0;
    real a2 = 0;
    if (theta < tol) {
        // sinx = x - x3 / 3! + x5 / 5! - ...
        // sinx / x = 1 - x2 / 3! + x4 / 5! - ...
        a1 = 1 - theta_sqr / 6;
        // sinx / x = a1 + O(theta^4).
        // cosx = 1 - x2 / 2! + x4 / 4! - x6 / 6! + ...
        // (1 - cosx) / x2 = 1 / 2! - x2 / 4! + x4 / 6! - ...
        a2 = ToReal(0.5) - theta_sqr / 24;
        // (1 - cosx) / x2 = a2 + O(theta^4).
    } else {
        a1 = std::sin(theta) / theta;
        a2 = (1 - std::cos(theta)) / theta_sqr;
    }
    const Matrix3r A = ToCrossProductMatrix(angle);
    return Matrix3r::Identity() + a1 * A + a2 * A * A;
}

// Sparse matrices.
const SparseMatrixXr FromTriplet(const integer row_num, const integer col_num,
    const std::vector<Eigen::Triplet<real>>& nonzeros) {
    SparseMatrixXr mat(row_num, col_num);
    mat.setFromTriplets(nonzeros.begin(), nonzeros.end());
    mat.makeCompressed();
    return mat;
}

const std::vector<Eigen::Triplet<real>> ToTriplet(const SparseMatrixXr& mat) {
    SparseMatrixXr mat_compressed = mat;
    mat_compressed.makeCompressed();
    std::vector<Eigen::Triplet<real>> nonzeros;
    const integer outer_size = static_cast<integer>(mat_compressed.outerSize());
    for (integer k = 0; k < outer_size; ++k)
        for (SparseMatrixXr::InnerIterator it(mat_compressed, k); it; ++it) {
            nonzeros.push_back(Eigen::Triplet<real>(it.row(), it.col(),
                it.value()));
        }
    return nonzeros;
}

const SparseMatrixXr FromDiagonal(const VectorXr& diagonal) {
    CheckCondition(diagonal.size() > 0, "basic::FromDiagonal",
        "Empty diagonal.");
    const integer matrix_size = static_cast<integer>(diagonal.size());
    std::vector<Eigen::Triplet<real>> nonzeros;
    for (integer i = 0; i < matrix_size; ++i) {
        nonzeros.push_back(Eigen::Triplet<real>(i, i, diagonal(i)));
    }
    return FromTriplet(matrix_size, matrix_size, nonzeros);
}

const SparseMatrixXr Transpose(const SparseMatrixXr& mat) {
    const std::vector<Eigen::Triplet<real>> nonzeros = ToTriplet(mat);
    std::vector<Eigen::Triplet<real>> nonzeros_t;
    for (const auto& tuple : nonzeros)
        nonzeros_t.push_back(Eigen::Triplet<real>(
            tuple.col(), tuple.row(), tuple.value()));
    return FromTriplet(static_cast<integer>(mat.cols()),
        static_cast<integer>(mat.rows()), nonzeros_t);
}

const bool IsRotationMatrix(const Matrix3r& R) {
    const Matrix3r I = Matrix3r::Identity();
    return IsClose(R.determinant(), 1, ToReal(1e-3), 0) &&
        IsClose((R * R.transpose() - I).squaredNorm() / (3 * 3), 0,
            ToReal(1e-6), 0);
}

void PolarDecomposition(const Matrix3r& F, Matrix3r& R, Matrix3r& S) {
    const Eigen::JacobiSVD<Matrix3r>
        svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3r Sig = svd.singularValues().asDiagonal();
    const Matrix3r U = svd.matrixU();
    const Matrix3r V = svd.matrixV();
    R = U * V.transpose();
    S = V * Sig * V.transpose();
}

const std::pair<Matrix3r, Matrix3r> PolarDecompositionDifferential(
    const Matrix3r& F, const Matrix3r& R, const Matrix3r& S,
    const Matrix3r& dF) {
    const Matrix3r lhs = R.transpose() * dF - dF.transpose() * R;
    Matrix3r A = Matrix3r::Zero();
    A(0, 0) = S(0, 0) + S(1, 1);
    A(1, 1) = S(0, 0) + S(2, 2);
    A(2, 2) = S(1, 1) + S(2, 2);
    A(0, 1) = A(1, 0) = S(1, 2);
    A(0, 2) = A(2, 0) = -S(0, 2);
    A(1, 2) = A(2, 1) = S(0, 1);
    const Matrix3r A_inv = A.inverse();
    const Vector3r b(lhs(0, 1), lhs(0, 2), lhs(1, 2));
    const Vector3r xyz = A_inv * b;
    const real x = xyz(0), y = xyz(1), z = xyz(2);
    Matrix3r W = Matrix3r::Zero();
    W(0, 0) = W(1, 1) = W(2, 2) = 0;
    W(0, 1) = x; W(0, 2) = y;
    W(1, 0) = -x; W(1, 2) = z;
    W(2, 0) = -y; W(2, 1) = -z;
    const Matrix3r dR = R * W;
    const Matrix3r dS = R.transpose() * (dF - dR * S);
    return std::make_pair(dR, dS);
}

const std::pair<Matrix9r, Matrix9r> PolarDecompositionDifferential(
    const Matrix3r& F, const Matrix3r& R, const Matrix3r& S) {
    // lhs01 = R.col(0).dot(dF.col(1)) - dF.col(0).dot(R.col(1)).
    Vector9r lhs01, lhs02, lhs12;
    lhs01 << -R.col(1), R.col(0), Vector3r::Zero();
    lhs02 << -R.col(2), Vector3r::Zero(), R.col(0);
    lhs12 << Vector3r::Zero(), -R.col(2), R.col(1);
    Matrix3r A = Matrix3r::Zero();
    A(0, 0) = S(0, 0) + S(1, 1);
    A(1, 1) = S(0, 0) + S(2, 2);
    A(2, 2) = S(1, 1) + S(2, 2);
    A(0, 1) = A(1, 0) = S(1, 2);
    A(0, 2) = A(2, 0) = -S(0, 2);
    A(1, 2) = A(2, 1) = S(0, 1);
    const Matrix3r A_inv = A.inverse();
    Matrix3Xr b(3, 9);
    b.row(0) = lhs01; b.row(1) = lhs02; b.row(2) = lhs12;
    const Matrix3Xr xyz = A_inv * b;
    const Vector9r x = xyz.row(0), y = xyz.row(1), z = xyz.row(2);
    Matrix3r W = Matrix3r::Zero();
    W(0, 0) = W(1, 1) = W(2, 2) = 0;
    // R01 * -x + R02 * -y
    // R11 * -x + R12 * -y
    // R21 * -x + R22 * -y
    // R00 * x + R02 * -z
    // R10 * x + R12 * -z
    // R20 * x + R22 * -z
    // R00 * y + R01 * z
    // R10 * y + R11 * z
    // R20 * y + R21 * z
    const Matrix9r dRdF = lhs01 * x.transpose() + lhs02 * y.transpose()
        + lhs12 * z.transpose();
    // F = RS.
    // S(F) = R(F).T * F.
    // dS = dR.T * F + R.T * dF.
    //    = (dRdF * dF).reshape(dim, dim).T * F + R.T * dF.
    Matrix9r dSdF = Matrix9r::Zero();
    for (integer i = 0; i < 3; ++i)
        for (integer j = 0; j < 3; ++j) {
            Matrix3r dF = Matrix3r::Zero();
            dF(i, j) = 1;
            dSdF.col(i + j * 3) +=
                (dRdF.col(i + j * 3).reshaped(3, 3).transpose() * F
                    + R.transpose() * dF).reshaped();
        }
    return std::make_pair(dRdF, dSdF);
}

const Matrix3r DeterminantGradient(const Matrix3r& A) {
    // dJ/dA = JA^-T
    // A = [ a0 | a1 | a2 ].
    // J = a0.dot(a1 x a2).
    // dJ/dA = [ a1 x a2 | a2 x a0 | a0 x a1 ]
    Matrix3r dJdA;
    dJdA.col(0) = A.col(1).cross(A.col(2));
    dJdA.col(1) = A.col(2).cross(A.col(0));
    dJdA.col(2) = A.col(0).cross(A.col(1));
    return dJdA;
}

const Matrix9r DeterminantHessian(const Matrix3r& A) {
    Matrix9r H = Matrix9r::Zero();
    const Matrix3r A0 = ToCrossProductMatrix(A.col(0));
    const Matrix3r A1 = ToCrossProductMatrix(A.col(1));
    const Matrix3r A2 = ToCrossProductMatrix(A.col(2));
    H.block<3, 3>(0, 3) += -A2;
    H.block<3, 3>(0, 6) += A1;
    H.block<3, 3>(3, 0) += A2;
    H.block<3, 3>(3, 6) += -A0;
    H.block<3, 3>(6, 0) += -A1;
    H.block<3, 3>(6, 3) += A0;
    return H;
}

const bool CheckGradient(const VectorXr& x0,
    const std::function<const real(const VectorXr&)>& func,
    const std::function<const VectorXr(const VectorXr&)>& func_grad,
    const Options& opt) {

    return CheckGradient(x0, func, func_grad,
        [](const integer i){ return false; }, opt);
}

const bool CheckGradient(const VectorXr& x0,
    const std::function<const real(const VectorXr&)>& func,
    const std::function<const VectorXr(const VectorXr&)>& func_grad,
    const std::function<const bool(const integer)>& skip_x_dof,
    const Options& opt) {

    const std::string error_location = "basic::CheckGradient";

    const integer verbose = opt.integer_option().at("verbose");
    const real grad_check_abs_tol = opt.real_option().at("grad_check_abs_tol");
    const real grad_check_rel_tol = opt.real_option().at("grad_check_rel_tol");

    const integer n = static_cast<integer>(x0.size());
    CheckCondition(n > 0, error_location, "Empty x0.");

    const VectorXr g0 = func_grad(x0);
    CheckCondition(g0.size() == x0.size(), error_location,
        "Dimensions of f and g do not agree.");

    // Check gradients.
    for (integer i = 0; i < n; ++i) {
        if (skip_x_dof(i)) continue;

        bool success = false;
        // Try a series of eps.
        for (const real eps : { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
                                1e-10, 1e-11, 1e-12 }) {
            VectorXr x_pos = x0;
            x_pos(i) += eps;
            VectorXr x_neg = x0;
            x_neg(i) -= eps;
            real f_pos, f_neg;
            try {
                f_pos = func(x_pos);
                f_neg = func(x_neg);
            } catch (std::runtime_error& ) {
                continue;
            }
            const real grad = (f_pos - f_neg) / (2 * eps);
            if (verbose > 0) {
                std::cout << "Active gradient at x" << i << " (eps = "
                    << eps << "): " << g0(i) << " (" << std::abs(grad - g0(i))
                    << ")" << std::endl;
            }
            if (IsClose(grad, g0(i), grad_check_abs_tol, grad_check_rel_tol)) {
                success = true;
                break;
            }
        }
        CheckCondition(success, error_location, "Gradient check failed.");
        if (!success) return false;
    }

    return true;
}

const bool CheckJacobian(const VectorXr& x0,
    const std::function<const VectorXr(const VectorXr&)>& func,
    const std::function<const MatrixXr(const VectorXr&)>& func_jac,
    const Options& opt) {

    return CheckJacobian(x0, func, func_jac,
        [](const integer i){ return false; },
        [](const integer i){ return false; },
        opt);
}

const bool CheckJacobian(const VectorXr& x0,
    const std::function<const VectorXr(const VectorXr&)>& func,
    const std::function<const MatrixXr(const VectorXr&)>& func_jac,
    const std::function<const bool(const integer)>& skip_x_dof,
    const std::function<const bool(const integer)>& skip_func_dof,
    const Options& opt) {

    const std::string error_location = "basic::CheckJacobian";

    const integer verbose = opt.integer_option().at("verbose");
    const real jac_check_abs_tol = opt.real_option().at("jac_check_abs_tol");
    const real jac_check_rel_tol = opt.real_option().at("jac_check_rel_tol");

    const integer n = static_cast<integer>(x0.size());
    CheckCondition(n > 0, error_location, "Empty x0.");

    const VectorXr f0 = func(x0);
    const integer m = static_cast<integer>(f0.size());
    const MatrixXr J0 = func_jac(x0);
    CheckCondition(J0.rows() == f0.size() && J0.cols() == x0.size(),
        error_location, "Dimensions of J, f, and x do not agree.");

    VectorXr active_f_dofs = VectorXr::Ones(m);
    for (integer i = 0; i < m; ++i) {
        active_f_dofs(i) = (active_f_dofs(i) ? 1 : 0);
    }

    // Check Jacobian.
    for (integer i = 0; i < n; ++i) {
        if (skip_x_dof(i)) continue;

        bool success = false;
        // Try a series of eps.
        for (const real eps : { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
                                1e-10, 1e-11, 1e-12 }) {
            VectorXr x_pos = x0;
            x_pos(i) += eps;
            VectorXr x_neg = x0;
            x_neg(i) -= eps;
            VectorXr f_pos, f_neg;
            try {
                f_pos = func(x_pos);
                f_neg = func(x_neg);
            } catch (std::runtime_error& ) {
                continue;
            }
            const VectorXr active_Jac =
                (f_pos - f_neg).cwiseProduct(active_f_dofs) / (2 * eps);
            const VectorXr active_J0 = J0.col(i).cwiseProduct(active_f_dofs);
            if (verbose > 0) {
                std::cout << "Active Jacobian at x" << i << " (eps = "
                    << eps << "): " << active_J0.norm() << " ("
                    << (active_Jac - active_J0).norm() << ")" << std::endl;
            }
            if (IsClose(active_Jac, active_J0, jac_check_abs_tol,
                        jac_check_rel_tol)) {
                success = true;
                break;
            }
        }
        CheckCondition(success, error_location, "Jacobian check failed.");
        if (!success) return false;
    }

    return true;
}

}