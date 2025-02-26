#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_COROTATED_MATERIAL
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_COROTATED_MATERIAL

#include "material/include/material.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

class CorotatedMaterial : public Material {
public:
    CorotatedMaterial() : Material(MaterialType::kCorotated) {}

    const real Psi(const Matrix3r& F) const override;
    const Matrix3r P(const Matrix3r& F) const override;
    const Matrix3r dP(const Matrix3r& F, const Matrix3r& dF) const override;
    const Matrix9r dPdF(const Matrix3r& F) const override;
};

}
}

#endif