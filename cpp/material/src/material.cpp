#include "material/include/material.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {
namespace material {

Material::Material(const MaterialType type) :
    type_(type), density_(0), youngs_modulus_(0), poissons_ratio_(0),
    lambda_(0), mu_(0) {}

void Material::Initialize(const real density, const real youngs_modulus,
    const real poissons_ratio) {

    CheckCondition(density > 0 && youngs_modulus > 0 &&
        0 < poissons_ratio && poissons_ratio < 0.5,
        "material::Material::Initialize", "Invalid material parameters.");

    density_ = density;
    youngs_modulus_ = youngs_modulus;
    poissons_ratio_ = poissons_ratio;
    const real k = youngs_modulus_;
    const real nu = poissons_ratio_;
    lambda_ = k * nu / ((1 + nu) * (1 - 2 * nu));
    mu_ = k / (2 * (1 + nu));
}

const std::string ToString(const MaterialType type) {
    switch (type) {
        case MaterialType::kLinear:
            return "linear";
        case MaterialType::kStVenantKirchhoff:
            return "st_venant_kirchhoff";
        case MaterialType::kCorotated:
            return "corotated";
        case MaterialType::kNeohookean:
            return "neohookean";
        default:
            return "unsupported type";
    }
}

}
}