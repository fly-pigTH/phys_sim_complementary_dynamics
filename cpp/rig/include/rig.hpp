#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_RIG_RIG
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_RIG_RIG

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"

namespace phys_sim_complementary_dynamics {
namespace rig {

enum class RigType {
    kLinearBlendSkinning = 0,
    kTotalNum
};

class Rig {
public:
    Rig(const RigType type);
    virtual ~Rig() {}

    void Initialize(const integer parameter_num, const Matrix3Xr& vertices,
        const Matrix4Xi& elements, const Options& opt);

    const RigType type() const { return type_; }
    const integer parameter_num() const { return parameter_num_; }
    const integer m() const { return parameter_num_; }
    const integer vertex_num() const { return vertex_num_; }
    const integer n() const { return vertex_num_; }
    const Matrix3Xr& vertices() const { return vertices_; }
    const integer element_num() const { return element_num_; }
    const Matrix4Xi& elements() const { return elements_; }

    const Matrix3Xr ComputeDisplacements(const VectorXr& parameters) const {
        return ur(parameters);
    }
    virtual const Matrix3Xr ur(const VectorXr& p) const = 0;

    const SparseMatrixXr ComputeJacobian(const VectorXr& parameters) const {
        return J(parameters);
    }
    virtual const SparseMatrixXr J(const VectorXr& p) const = 0;
    virtual const bool IsConstantJacobian() const = 0;

protected:
    virtual void InitializeDerived(const Options& opt) = 0;

    void CheckParameters(const VectorXr& parameters) const;

private:
    RigType type_;

    // Rig parameters.
    integer parameter_num_;

    // A 3D tetrahedron mesh.
    integer vertex_num_;
    Matrix3Xr vertices_;

    integer element_num_;
    Matrix4Xi elements_;
};

const std::string ToString(const RigType type);

}
}

#endif
