#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_RIG_LINEAR_BLEND_SKINNING_RIG
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_RIG_LINEAR_BLEND_SKINNING_RIG

#include "rig/include/rig.hpp"

namespace phys_sim_complementary_dynamics {
namespace rig {

class LinearBlendSkinningRig : public Rig {
public:
    LinearBlendSkinningRig() : Rig(RigType::kLinearBlendSkinning),
        transformation_number_(0), weights_(MatrixXr::Zero(0, 0)) {}

    const Matrix3Xr ur(const VectorXr& p) const override;
    const SparseMatrixXr J(const VectorXr& p) const override;
    const bool IsConstantJacobian() const override { return true; }

    const integer transformation_number() const {
        return transformation_number_;
    }
    const integer k() const { return transformation_number_; }
    const MatrixXr& weights() const { return weights_; }
    const MatrixXr& w() const { return weights_; }

private:
    void InitializeDerived(const Options& opt) override;

    // Number of local transformations.
    integer transformation_number_;
    MatrixXr weights_;
};

}
}

#endif