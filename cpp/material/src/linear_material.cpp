#include "material/include/linear_material.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

const real LinearMaterial::Psi(const Matrix3r& F) const {
    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r e = (F + F.transpose()) / 2 - I;
    const real trace_e = e.trace();
    return mu() * e.array().square().sum() + lambda() / 2 * (trace_e * trace_e);
}

const Matrix3r LinearMaterial::P(const Matrix3r& F) const {
    const Matrix3r I = Matrix3r::Identity();
    return mu() * (F + F.transpose() - 2 * I) + lambda() * (F - I).trace() * I;
}

const Matrix3r LinearMaterial::dP(const Matrix3r& F, const Matrix3r& dF) const {
    const Matrix3r I = Matrix3r::Identity();
    return mu() * (dF + dF.transpose()) + lambda() * dF.trace() * I;
}

const Matrix9r LinearMaterial::dPdF(const Matrix3r& F) const {
    Matrix9r ret; ret.setZero();
    // mu * dF.
    for (integer i = 0; i < 9; ++i) ret(i, i) = mu();
    // mu * dF.transpose().
    for (integer i = 0; i < 3; ++i)
        for (integer j = 0; j < 3; ++j)
            ret(i * 3 + j, j * 3 + i) += mu();
    // la * dF.trace() * I.
    for (integer i = 0; i < 3; ++i)
        for (integer j = 0; j < 3; ++j)
            ret(i * 3 + i, j * 3 + j) += lambda();
    return ret;
}

}
}