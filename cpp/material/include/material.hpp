#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_MATERIAL
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_MATERIAL_MATERIAL

#include "basic/include/config.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

enum class MaterialType {
    kLinear = 0,
    kStVenantKirchhoff,
    kCorotated,
    kNeohookean,
    kTotalNum
};

class Material {
public:
    Material(const MaterialType type);
    virtual ~Material() {}

    void Initialize(const real density, const real youngs_modulus,
        const real poissons_ratio);

    const MaterialType type() const { return type_; }
    const real density() const { return density_; }
    const real rho() const { return density_; }
    const real youngs_modulus() const { return youngs_modulus_; }
    const real k() const { return youngs_modulus_; }
    const real poissons_ratio() const { return poissons_ratio_; }
    const real nu() const { return poissons_ratio_; }
    const real lambda() const { return lambda_; }
    const real mu() const { return mu_; }

    const real ComputeEnergyDensity(const Matrix3r& deformation_grad) const {
        return Psi(deformation_grad);
    }
    virtual const real Psi(const Matrix3r& F) const = 0;
    const Matrix3r ComputeStressTensor(const Matrix3r& deformation_grad) const {
        return P(deformation_grad);
    }
    virtual const Matrix3r P(const Matrix3r& F) const = 0;
    const Matrix3r ComputeStressTensorDifferential(
        const Matrix3r& deformation_grad,
        const Matrix3r& delta_deformation_grad) const {
        return dP(deformation_grad, delta_deformation_grad);
    }
    virtual const Matrix3r dP(const Matrix3r& F, const Matrix3r& dF) const = 0;
    const Matrix9r ComputeStressTensorDifferential(
        const Matrix3r& deformation_grad) const {
        return dPdF(deformation_grad);
    }
    virtual const Matrix9r dPdF(const Matrix3r& F) const = 0;

private:
    MaterialType type_;

    real density_;
    real youngs_modulus_;
    real poissons_ratio_;
    real lambda_;
    real mu_;
};

const std::string ToString(const MaterialType type);

}
}

#endif