#include "rig/include/linear_blend_skinning_rig.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_complementary_dynamics {
namespace rig {

const Matrix3Xr LinearBlendSkinningRig::ur(const VectorXr& p) const {
    CheckParameters(p);

    ////////////////////////////////////////////////////////////////////////////
    // Task 3.1 (2 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Implement the displacements from linear-blend skinning. Please refer to
    // J and InitializeDerived below to figure out the meaning of p. Your
    // implementation should be compatible with J, which computes the Jacobian
    // of this function with respect to its input p.
    //
    // TODO.
    return Matrix3Xr::Zero(3, n());
}

const SparseMatrixXr LinearBlendSkinningRig::J(const VectorXr& p) const {
    CheckParameters(p);

    std::vector<Eigen::Triplet<real>> nonzeros;
    const Matrix3Xr& X = vertices();
    for (integer k = 0; k < transformation_number_; ++k) {
        const Matrix3r R = p.segment<9>(12 * k).reshaped(3, 3);
        const Vector3r t = p.segment<3>(12 * k + 9);
        for (integer i = 0; i < n(); ++i) {
            const Vector3r Xw = X.col(i) * weights_(k, i);
            for (integer j = 0; j < 3; ++j) {
                nonzeros.push_back(Eigen::Triplet<real>(
                    3 * i + j, 12 * k + j, Xw(0)));
                nonzeros.push_back(Eigen::Triplet<real>(
                    3 * i + j, 12 * k + 3 + j, Xw(1)));
                nonzeros.push_back(Eigen::Triplet<real>(
                    3 * i + j, 12 * k + 6 + j, Xw(2)));
                nonzeros.push_back(Eigen::Triplet<real>(
                    3 * i + j, 12 * k + 9 + j, weights_(k, i)));
            }
        }
    }

    return FromTriplet(3 * n(), m(), nonzeros);
}

void LinearBlendSkinningRig::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "rig::LinearBlendSkinningRig::InitializeDerived";

    // Initialize transformation_number_.
    transformation_number_ = m() / 12;
    CheckCondition(transformation_number_ > 0 &&
        transformation_number_ * 12 == m(), error_location,
        "Invalid k number.");

    // Initialize weights_.
    weights_ = MatrixXr::Zero(transformation_number_, n());
    // A simple weight scheme.
    CheckCondition(opt.GetMatrixOptionRows("samples") == 3 &&
        opt.GetMatrixOptionCols("samples") == transformation_number_,
        error_location, "Incompatible sample number.");
    const Matrix3Xr samples = opt.matrix_option().at("samples");
    const Matrix3Xr& X = vertices();
    const real bbox_size = (X.rowwise().maxCoeff()
        - X.rowwise().minCoeff()).maxCoeff();
    // TODO: make 0.3 an input parameter?
    const real sigma_sqr = 0.3 * bbox_size * bbox_size;
    for (integer k = 0; k < transformation_number_; ++k) {
        const RowVectorXr dist_sqr = (X.colwise() - samples.col(k)
            ).cwiseAbs2().colwise().sum();
        // A few reference numbers:
        // - exp(-1) \approx 0.37;
        // - exp(-2) \approx 0.14;
        // - exp(-3) \approx 0.05;
        weights_.row(k) = (-dist_sqr / sigma_sqr).array().exp();
    }

    // Normalize weights_.
    const RowVectorXr weights_sum = weights_.colwise().sum();
    weights_ = weights_.cwiseQuotient(
        VectorXr::Ones(transformation_number_) * weights_sum);
}

}
}