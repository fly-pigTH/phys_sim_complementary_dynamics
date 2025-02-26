#include "sim/include/simulator.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_complementary_dynamics {
namespace sim {

void Simulator::ClearForce(const integer vertex_index) {
    CheckVertexIndex(vertex_index, "ClearForce");
    forces_.col(vertex_index).setZero();
}

void Simulator::ApplyForce(const integer vertex_index, const Vector3r& force) {
    CheckVertexIndex(vertex_index, "ApplyForce");
    forces_.col(vertex_index) += force;
}

void Simulator::Step(const real time_step, const Options& opt) {
    const std::string error_location = "sim::Simulator::Step";
    CheckCondition(time_step > 0, error_location,
        "Time step must be positive.");
    const real h = time_step;

    const VectorXr& p = rig_parameters_;
    const Matrix3Xr ur = rig_->ur(p);
    const Matrix3Xr& X = rest_vertices_;

    if (opt.bool_option().at("grad_check")) {
        // Rig derivatives.
        const auto func = [&](const VectorXr& x) -> const VectorXr {
            return rig_->ur(x).reshaped();
        };
        const auto func_jac = [&](const VectorXr& x) -> const MatrixXr {
            return rig_->J(x).toDense();
        };
        CheckJacobian(rig_parameters_, func, func_jac, opt);

        // Elasticity --- this test is expensive!
        CheckElastics(positions_, opt);
    }

    if (opt.bool_option().at("skip_complement")) {
        velocities_ = (X + ur - positions_) / h;
        positions_ = X + ur;
        return;
    }

    const integer verbose = opt.integer_option().at("verbose");
    const Matrix3Xr M = Vector3r::Ones() * masses_.transpose();
    const SparseMatrixXr M_diag = FromDiagonal(M.reshaped());
    const SparseMatrixXr J = rig_->J(p);
    const integer n = vertex_num_;
    const integer m = rig_->m();
    // J is a tall and thin matrix.
    // Define u = ur + uc. Consider a perturbation at ur(p):
    // \delta ur = J * \delta p.
    // Complementary dynamics require that uc be limited to
    // \int_\Omega <\delta ur, uc> dX = 0 for any \delta ur. Assuming the soft
    // body is homogeneous, the above constraint becomes
    // J.T * M * uc = 0.        (Eq. 9 in Complementary Dynamics.)
    //
    // In summary, the dynamics of u is constrained by
    // J.T * M * (u - ur) = 0.
    // Let C = J.T * M. C is a short and fat matrix. C.T = M * J.
    const SparseMatrixXr Ct = M_diag * J;
    const SparseMatrixXr C = Transpose(Ct);
    // Now the constraint becomes C * (x - X - ur) = 0, or
    // C * x = C * (X + ur).
    const VectorXr b = C * (X + ur).reshaped();
    // Now consider the incremental potential:
    // min_x 1 / 2h^2 (x - y) M (x - y) + E(x),
    // where y = x0 + h * v0 + h2 M^-1 f_ext. To double check:
    // g = M / h2 (x - y) - f(x) = 0 =>
    // M / h2 (x - x0 - h v0) = f_ext + f(x).
    const Matrix3Xr y = positions_ + h * velocities_
        + h * h * forces_.cwiseQuotient(M);

    // With the kinematic constraints, the IP optimization becomes:
    // min_x 1 / 2h2 (x - y) M (x - y) + E(x),
    // s.t. C * x = b.
    //
    // Consider the Lagrangian multipler method:
    // min_x 1 / 2h2 (x - y) M (x - y) + E(x) + la.T * (C x - b).
    // The stationary point satisfies
    // M / h2 (x - y) - fe(x) + C.T * la = 0.
    // C x = b.
    // We apply Newton's method to find (x, la):
    // - Pick an initial guess x0 and la0.
    // - At each iteration, we compute (dx, dla) from the linearization at the
    //   current (xk, lak):
    //   M / h2 (xk + dx - y) - fe(xk) + K(xk) dx + C.T lak + C.T dla = 0.
    //   C xk + C dx = b.
    //   [ M / h2 + K, C.T ] [ dx  ] = [ -M / h2 (xk - y) + fe(xk) - C.T lak ]
    //   [ C,            0 ] [ dla ]   [ b - C xk ].
    // - After solving (dx, dla), run line search:
    //   xk = xk + alpha * dx,
    //   lak = lak + alpha * dla.
    //   We use the IP energy as the merit function. Note that the optimization
    //   problem has a constraint, so we additionally enforce xk always falls in
    //   the feasible set C x = b.
    //   If xk is feasible, i.e., b - C xk = 0, then the line search is always
    //   feasible. This is because C dx = 0, so any (xk + alpha * dx) satisfies
    //   C (xk + alpha * dx) = b. Therefore, it is sufficient to ensure x0 is
    //   feasible, after which we can use the IP energy as the merit function.

    const integer max_newton_iter = opt.integer_option().at("max_newton_iter");
    CheckCondition(max_newton_iter >= 0, error_location,
        "Negative max_newton_iter.");
    const integer max_ls_iter = opt.integer_option().at("max_ls_iter");
    CheckCondition(max_ls_iter >= 0, error_location, "Negative max_ls_iter.");

    const auto merit = [&](const Matrix3Xr& x) -> const real {
        const real ip = M.cwiseProduct(x - y).cwiseProduct(x - y).sum()
            / (2 * h * h) + E(x);
        return ip;
    };
    const auto merit_grad = [&](const Matrix3Xr& x) -> const Matrix3Xr {
        const Matrix3Xr grad = M.cwiseProduct(x - y) / (h * h) - fe(x);
        return grad;
    };
    const auto residual = [&](const Matrix3Xr& x, const VectorXr& la)
        -> const std::pair<VectorXr, VectorXr> {
        const VectorXr momentum_residual = merit_grad(x).reshaped() + Ct * la;
        const VectorXr constraint_residual = C * x.reshaped() - b;
        return { momentum_residual, constraint_residual };
    };
    // Functions useful for debugging purposes.
    const auto converged = [&](const std::pair<VectorXr, VectorXr>& residual)
        -> const bool {
        // What is a good threshold for the force?
        // One possibility is a certain percentage of gravity, say, 1%;
        // Also, for most of our scenes, the gravity of 10g mass should be
        // sufficiently small.
        // We will use the larger of the two as the force threshold.
        // We multiple it by h since we are doing impulse-based solve.
        const real force_tol = std::max(masses_.mean() * 9.81 * 0.01,
            1e-2 * 9.81);
        // For the displacement, we use 1mm as the tolerance.
        const real disp_tol = 1e-3;
        return residual.first.cwiseAbs().maxCoeff() < force_tol &&
            residual.second.cwiseAbs().maxCoeff() < disp_tol;
    };

    // Initial guess: note that this x0 satisfies C x0 = b.
    Matrix3Xr xk = X + ur;
    VectorXr lak = MatrixXr(C * Ct).colPivHouseholderQr().solve(
        -C * merit_grad(xk).reshaped());
    std::pair<VectorXr, VectorXr> rk = residual(xk, lak);
    real mk = merit(xk);
    if (verbose > 0) {
        std::cout << "------ Newton iter 0 ------" << std::endl;
        std::cout << "|x0| = " << xk.stableNorm() << ", |la0| = "
            << lak.stableNorm() << std::endl;
        std::cout << "|r0|_inf = " << rk.first.cwiseAbs().maxCoeff() << ", "
            << rk.second.cwiseAbs().maxCoeff() << std::endl;
        std::cout << "m0 = " << mk << std::endl;
    }
    if (converged(rk)) {
        if (verbose > 0) std::cout << "CONVERGED." << std::endl;
        // Converged.
        velocities_ = (xk - positions_) / h;
        positions_ = xk;
        return;
    }

    // The Newton loop.
    for (integer k = 0; k < max_newton_iter; ++k) {
        // Assemble the linear system of equations.
        // [ M / h2 + K, C.T ] [ dx  ] = [ -M / h2 (xk - y) + fe(xk) - C.T lak ]
        // [ C,            0 ] [ dla ]   [ b - C xk ].
        //
        // Apply Schur's complement:
        // [ A, C.T ] [ dx  ] = [ c ]
        // [ C,   0 ] [ dla ]   [ d ]
        // =>
        // [ C, C inv(A) C.T ] [ dx  ] = [ C inv(A) c ]
        // [ C,            0 ] [ dla ]   [ d ]
        // =>
        // C inv(A) C.T dla = C inv(A) c - d.
        const SparseMatrixXr A = M_diag / (h * h) + K(xk);
        const VectorXr c = -M.cwiseProduct(xk - y).reshaped() / (h * h)
            + fe(xk).reshaped() - Ct * lak;
        const VectorXr d = b - C * xk.reshaped();

        ////////////////////////////////////////////////////////////////////////
        // Task 3.2 (3 points).
        ////////////////////////////////////////////////////////////////////////
        //
        // Solve (dx, dla). The above comment describes a Schur-complement way,
        // but feel free to use other linear solvers in Eigen or even implement
        // your own solver from scratch.
        //
        // TODO.
        VectorXr dla = VectorXr::Zero(d.size());
        Matrix3Xr dx = Matrix3Xr::Zero(3, n);

        // Line search.
        const Matrix3Xr merit_gk = merit_grad(xk);
        const real ddk = merit_gk.cwiseProduct(dx).sum();
        real alpha = 1;
        real mk_next = 0;
        bool ls_succeed = false;
        for (integer j = 0; j < max_ls_iter; ++j) {
            mk_next = merit(xk + alpha * dx);
            if (verbose > 1) {
                std::cout << "Line search " << j << ", alpha = " << alpha
                    << ", mk(alpha) = " << mk_next << ", mk = "
                    << mk << std::endl;
            }
            if (mk_next < mk) {
                ls_succeed = true;
                xk = xk + alpha * dx;
                lak = lak + alpha * dla;
                mk = mk_next;
                break;
            }
            alpha /= 2;
        }
        if (!ls_succeed) {
            CheckCondition(ddk < 0, error_location, "Not a descent direction: "
                "dd = " + std::to_string(ddk) + ".");
            CheckCondition(false, error_location, "Line search failed. "
                "Consider increasing max ls iterations.");
        }

        // Convergence check.
        rk = residual(xk, lak);
        if (verbose > 0) {
            std::cout << "------ Newton iter " << k + 1 << " ------"
                << std::endl;
            std::cout << "|x" << k + 1 << "| = " << xk.stableNorm()
                << ", |la" << k + 1 << "| = " << lak.stableNorm() << std::endl;
            std::cout << "|r" << k + 1 << "|_inf = "
                << rk.first.cwiseAbs().maxCoeff() << ", "
                << rk.second.cwiseAbs().maxCoeff() << std::endl;
            std::cout << "m" << k + 1 << " = " << mk << std::endl;
        }
        if (converged(rk)) {
            if (verbose > 0) std::cout << "CONVERGED." << std::endl;
            // Converged.
            velocities_ = (xk - positions_) / h;
            positions_ = xk;
            return;
        }
    }

    CheckCondition(false, error_location, "Newton fails. Consider increasing "
        "the max iterations or relaxing the convergence tolerance.");
}

}
}