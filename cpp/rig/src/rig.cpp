#include "rig/include/rig.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {
namespace rig {

Rig::Rig(const RigType type)
    : type_(type), parameter_num_(0),
    vertex_num_(0), vertices_(Matrix3Xr::Zero(3, 0)),
    element_num_(0), elements_(Matrix4Xi::Zero(4, 0)) {}

void Rig::Initialize(const integer parameter_num, const Matrix3Xr& vertices,
    const Matrix4Xi& elements, const Options& opt) {
    const std::string error_location = "rig::Rig::Initialize";

    CheckCondition(parameter_num > 0, error_location,
        "Zero or negative parameter number.");
    parameter_num_ = parameter_num;

    vertex_num_ = static_cast<integer>(vertices.cols());
    CheckCondition(vertex_num_ > 0, error_location, "Invalid vertices.");
    vertices_ = vertices;

    element_num_ = static_cast<integer>(elements.cols());
    CheckCondition(element_num_ > 0, error_location, "Invalid elements.");
    elements_ = elements;

    CheckCondition(elements_.minCoeff() >= 0 &&
        elements_.maxCoeff() < vertex_num_, error_location,
        "Out-of-bound vertex indices.");

    InitializeDerived(opt);
}

void Rig::CheckParameters(const VectorXr& parameters) const {
    CheckCondition(static_cast<integer>(parameters.size()) == parameter_num_,
        "rig::Rig::CheckParameters", "Incompatible parameter size.");
}

const std::string ToString(const RigType type) {
    switch (type) {
        case RigType::kLinearBlendSkinning:
            return "linear_blend_skinning";
        default:
            return "unsupported type";
    }
}

}
}