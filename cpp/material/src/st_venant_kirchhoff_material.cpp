#include "material/include/st_venant_kirchhoff_material.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

const real StVenantKirchhoffMaterial::Psi(const Matrix3r& F) const {
    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r E = (F.transpose() * F - I) / 2;
    const real trace_E = E.trace();
    return mu() * E.array().square().sum() + lambda() / 2 * (trace_E * trace_E);
}

const Matrix3r StVenantKirchhoffMaterial::P(const Matrix3r& F) const {
    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r E = (F.transpose() * F - I) / 2;
    const real trace_E = E.trace();
    return F * (2 * mu() * E + lambda() * trace_E * I);
}

const Matrix3r StVenantKirchhoffMaterial::dP(const Matrix3r& F,
    const Matrix3r& dF) const {

    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r E = (F.transpose() * F - I) / 2;
    const Matrix3r dE = (dF.transpose() * F + F.transpose() * dF) / 2;
    const real trace_E = E.trace();
    const real trace_dE = dE.trace();
    return dF * (2 * mu() * E + lambda() * trace_E * I)
        + F * (2 * mu() * dE + lambda() * trace_dE * I);
}

const Matrix9r StVenantKirchhoffMaterial::dPdF(const Matrix3r& F) const {

    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r E = (F.transpose() * F - I) / 2;
    const real trace_E = E.trace();

    Matrix9r dPdF; dPdF.setZero();
    for (integer i = 0; i < 3; ++i)
        for (integer j = 0; j < 3; ++j) {
            Matrix3r dF; dF.setZero();
            dF(i, j) = 1;
            Matrix3r dE; dE.setZero();
            dE.col(j) += F.row(i).transpose();
            dE.row(j) += F.row(i);
            dE /= 2;
            const real trace_dE = F(i, j);
            dPdF.col(i + j * 3) += (dF * (2 * mu() * E + lambda() * trace_E * I)
                + F * (2 * mu() * dE + lambda() * trace_dE * I)).reshaped();
        }
    return dPdF;
}

}
}