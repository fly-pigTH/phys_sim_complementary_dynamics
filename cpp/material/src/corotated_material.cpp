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
    const Matrix3r I = Matrix3r::Identity();
    const Matrix3r e_c = S - I;
    const real trace_e_c = e_c.trace();
    // Task End
    return mu() * e_c.array().square().sum() + lambda() / 2 * (trace_e_c * trace_e_c);
}

const Matrix3r CorotatedMaterial::P(const Matrix3r& F) const {
    ////////////////////////////////////////////////////////////////////////////
    // Task 1.2 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // TODO.
    return Matrix3r::Zero();
     // Task Begin
     Matrix3r R, S;
     PolarDecomposition(F, R, S);
     const auto dXdF = PolarDecompositionDifferential(F, R, S);
     const Matrix9r dRdF= dXdF.first;
     const Matrix9r dSdF = dXdF.second;
     const Matrix3r I = Matrix3r::Identity();
     const Matrix3r e_c = S - I;
     const real trace_e_c = e_c.trace();
     const Vector9r dPsidS = (2 * mu() * e_c + lambda() * trace_e_c * I).transpose().reshaped();
     Matrix3r dPsidF = Matrix3r::Zero();
 
     // using paper result directly
     dPsidF = 2* mu() * (F-R) + lambda() * (S-I).trace() * R;
     // self calculation
     // dPsidF = (dSdF.transpose() * dPsidS).reshaped(3, 3).transpose();
     // Task End
     return dPsidF;
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
    // Task begin
    Matrix3r R, S;
    PolarDecomposition(F, R, S);
    const Matrix9r I9 = Matrix9r::Identity();

    const auto dXdF = PolarDecompositionDifferential(F, R, S);
    const Matrix9r dRdF= dXdF.first;
    const Matrix9r dSdF = dXdF.second;

    Eigen::Matrix<double, 9, 1> vecdRdF = Eigen::Map<Eigen::Matrix<double, 9, 1>>(R.data());
    Matrix9r RRT = vecdRdF * vecdRdF.transpose();

    Matrix9r dPdF = Matrix9r::Zero();
    dPdF = 2* mu() * (I9-dRdF) + lambda()*(R.transpose()*F).trace()*dRdF + lambda() * RRT;
    //Task end
    return dPdF;
}

}
}