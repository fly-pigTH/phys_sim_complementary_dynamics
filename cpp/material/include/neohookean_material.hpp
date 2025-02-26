#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_NEOHOOKEAN_MATERIAL
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_NEOHOOKEAN_MATERIAL

#include "material/include/material.hpp"

// Implementing the stable Neohookean material described in
// "Stable Neo-Hookean Flesh Simulation"
// https://graphics.pixar.com/library/StableElasticity/paper.pdf
// We implemented Eqn. (14) for the 3D version and Eqn. (13) for the 2D version.
// Note that the paper from dl.acm.org contains a typo in alpha used by Eqn.
// (14). The Pixar version (the link above) fixed this typo and should be used.

namespace phys_sim_complementary_dynamics {
namespace material {

class NeohookeanMaterial : public Material {
public:
    NeohookeanMaterial() : Material(MaterialType::kNeohookean) {}

    const real Psi(const Matrix3r& F) const override;
    const Matrix3r P(const Matrix3r& F) const override;
    const Matrix3r dP(const Matrix3r& F, const Matrix3r& dF) const override;
    const Matrix9r dPdF(const Matrix3r& F) const override;
};

}
}

#endif
