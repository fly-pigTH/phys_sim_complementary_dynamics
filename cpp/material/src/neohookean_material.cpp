#include "material/include/neohookean_material.hpp"
#include "basic/include/math.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

// Implemented the Stable neohookean material which is derived as follows.
// First, some notations:
// - J = |F|,
// - C = Ft * F,
// - Ic = tr(C).
// dIc/dF = 2F
// dJ/dF = JF^-T
//
// \Psi(F) = 0.5mu * (Ic - dim) + lambda / 2 * (J - alpha)^2
//           - 0.5mu * log(Ic + delta).
// P = mu * F + lambda * (J - alpha) * JF^-T - mu * F / (Ic + delta)
//   = (1 - 1 / (Ic + delta)) * mu * F + lambda * (J - alpha) * JF^{-T}.
//
// When F = R (J = 1, Ic = dim, dJ/dF = R), P should be zero:
// P = [1 - 1 / (dim + delta)] * mu * R + lambda * (1 - alpha) * R.
// So we obtain [1 - 1 / (dim + delta)] * mu + lambda * (1 - alpha) = 0.
// alpha = [1 - 1 / (dim + delta)] * mu / lambda + 1 <--------------- Eqn. (1).
// Note that the motivation requires delta > 0, so alpha is also positive.
//
// Furthermore, we want the Hessian to reach semi-definiteness at F = 0.
// This leads to a second constraint on alpha and delta, allowing us to solve
// both. See the functions below for details about the Hessian.

const real NeohookeanMaterial::Psi(const Matrix3r& F) const {
    const Matrix3r C = F.transpose() * F;
    const real J = F.determinant();
    const real Ic = C.trace();
    const real delta = 1;
    const real alpha = (1 - 1 / (3 + delta)) * mu() / lambda() + 1;
    return mu() / 2 * (Ic - 3) + lambda() / 2 * (J - alpha) * (J - alpha)
        - 0.5 * mu() * std::log(Ic + delta);
}

const Matrix3r NeohookeanMaterial::P(const Matrix3r& F) const {
    const Matrix3r C = F.transpose() * F;
    const real J = F.determinant();
    const real Ic = C.trace();
    const real delta = 1;
    const real alpha = (1 - 1 / (3 + delta)) * mu() / lambda() + 1;
    const Matrix3r dJdF = DeterminantGradient(F);
    return (1 - 1 / (Ic + delta)) * mu() * F + lambda() * (J - alpha) * dJdF;
}

const Matrix3r NeohookeanMaterial::dP(const Matrix3r& F,
    const Matrix3r& dF) const {
    const Matrix3r C = F.transpose() * F;
    const real J = F.determinant();
    const real Ic = C.trace();
    const real delta = 1;
    const real alpha = (1 - 1 / (3 + delta)) * mu() / lambda() + 1;
    // dJ/dF = JF^-T
    // F = [ f0 | f1 | f2 ].
    // J = f0.dot(f1 x f2).
    // dJ/dF = [ f1 x f2 | f2 x f0 | f0 x f1 ]
    const Matrix3r dJdF = DeterminantGradient(F);
    // P = (1 - 1 / (Ic + delta)) * mu * F + la * (J - alpha) * dJdF.
    // dIcdF = 2F.
    const real dIc = 2 * F.cwiseProduct(dF).sum();
    const real dJ = dJdF.cwiseProduct(dF).sum();
    const Matrix3r ddJdF = (DeterminantHessian(F)
        * dF.reshaped()).reshaped(3, 3);
    return (1 - 1 / (Ic + delta)) * mu() * dF
        + dIc / ((Ic + delta) * (Ic + delta)) * mu() * F
        + lambda() * (J - alpha) * ddJdF + lambda() * dJ * dJdF;
}

const Matrix9r NeohookeanMaterial::dPdF(const Matrix3r& F) const {
    const Matrix3r C = F.transpose() * F;
    const real J = F.determinant();
    const real Ic = C.trace();
    const real delta = 1;
    const real alpha = (1 - 1 / (3 + delta)) * mu() / lambda() + 1;
    // dJ/dF = JF^-T
    // F = [ f0 | f1 | f2 ].
    // J = f0.dot(f1 x f2).
    // dJ/dF = [ f1 x f2 | f2 x f0 | f0 x f1 ]
    const Matrix3r dJdF = DeterminantGradient(F);
    const Vector3r f0 = F.col(0);
    const Vector3r f1 = F.col(1);
    const Vector3r f2 = F.col(2);
    // P = (1 - 1 / (Ic + delta)) * mu * F + la * (J - alpha) * dJdF.
    // Part I:
    const Vector9r f = F.reshaped();
    Matrix9r dPdF = ((1 - 1 / (Ic + delta)) * mu()) * Matrix9r::Identity()
        + (2 * mu() / ((Ic + delta) * (Ic + delta))) * (f * f.transpose());
    // Part II: la * (J - alpha) * dJdF.
    const Vector9r djdf = dJdF.reshaped();
    dPdF += lambda() * djdf * djdf.transpose();
    // The trickiest part in part II:
    dPdF += lambda() * (J - alpha) * DeterminantHessian(F);
    return dPdF;
    // Regarding delta and alpha:
    // F = 0, J = 0, Ic = 0.
    // eig = mu * (1 - 1 / delta) = 0 => delta = 1.
}

}
}