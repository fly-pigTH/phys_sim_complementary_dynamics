#include "material/include/corotated_material.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

const real CorotatedMaterial::Psi(const Matrix3r& F) const {
    Matrix3r R, S;
    PolarDecomposition(F, R, S);

    ////////////////////////////////////////////////////////////////////////////
    // Task 1.1 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // TODO.
    return 0;
}

const Matrix3r CorotatedMaterial::P(const Matrix3r& F) const {
    ////////////////////////////////////////////////////////////////////////////
    // Task 1.2 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // TODO.
    return Matrix3r::Zero();
}

const Matrix3r CorotatedMaterial::dP(const Matrix3r& F,
    const Matrix3r& dF) const {
    Matrix3r R, S;
    PolarDecomposition(F, R, S);
    const auto dRdS = PolarDecompositionDifferential(F, R, S, dF);
    const Matrix3r dR = dRdS.first;
    const Matrix3r dS = dRdS.second;
    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r e_c = S - I;
    const Matrix3r de_c = dS;
    const real trace_e_c = e_c.trace();
    const real dtrace_e_c = de_c.trace();
    return dR * (2 * mu() * e_c + lambda() * trace_e_c * I)
        + R * (2 * mu() * de_c + lambda() * dtrace_e_c * I);
}

const Matrix9r CorotatedMaterial::dPdF(const Matrix3r& F) const {
    ////////////////////////////////////////////////////////////////////////////
    // Task 1.3 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // TODO.
    return Matrix9r::Zero();
}

}
}