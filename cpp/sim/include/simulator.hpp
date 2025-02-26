#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_SIM_SIMULATOR
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_SIM_SIMULATOR

#include "basic/include/options.hpp"
#include "material/include/material.hpp"
#include "rig/include/rig.hpp"

namespace phys_sim_complementary_dynamics {
namespace sim {

class Simulator {
public:
    Simulator();
    ~Simulator() {}

    void AddFiniteElements(const Matrix3Xr& rest_vertices,
        const Matrix4Xi& elements,
        const material::MaterialType material_type, const real density,
        const real youngs_modulus, const real poissons_ratio);

    void AddRig(const rig::RigType type, const integer parameter_num,
        const Options& opt);

    // Time stepping.
    void ClearForce(const integer vertex_index);
    void ClearForces() { forces_.setZero(); }
    void ApplyForce(const integer vertex_index, const Vector3r& force);

    // The core time stepping function.
    void Step(const real time_step, const Options& opt);

    const integer vertex_num() const { return vertex_num_; }
    const Matrix3Xr& rest_vertices() const { return rest_vertices_; }
    const integer element_num() const { return element_num_; }
    const Matrix4Xi& elements() const { return elements_; }
    const std::shared_ptr<material::Material> material() const {
        return material_;
    }
    const VectorXr& masses() const { return masses_; }
    const VectorXr& volumes() const { return volumes_; }
    const std::shared_ptr<rig::Rig> rig() const { return rig_; }

    const Matrix3Xr& positions() const { return positions_; }
    const Matrix3Xr& q() const { return positions_; }
    const Matrix3Xr& velocities() const { return velocities_; }
    const Matrix3Xr& q_dot() const { return velocities_; }
    const VectorXr& rig_parameters() const { return rig_parameters_; }
    const VectorXr& p() const { return rig_parameters_; }
    const Matrix3Xr& forces() const { return forces_; }
    const Matrix3Xr& f() const { return forces_; }

    void set_positions(const Matrix3Xr& positions) { set_q(positions); }
    void set_q(const Matrix3Xr& q);
    void set_velocities(const Matrix3Xr& velocities) { set_q_dot(velocities); }
    void set_q_dot(const Matrix3Xr& q_dot);
    void set_rig_parameters(const VectorXr& rig_parameters) {
        set_p(rig_parameters);
    }
    void set_p(const VectorXr& p);
    void set_forces(const Matrix3Xr& forces) { set_f(forces); }
    void set_f(const Matrix3Xr& f);

private:
    void CheckVertexIndex(const integer vertex_index,
        const std::string& error_location) const;
    void CheckPositions(const Matrix3Xr& positions,
        const std::string& error_location) const;
    void CheckVelocities(const Matrix3Xr& velocities,
        const std::string& error_location) const;

    const real GetElasticEnergy(const Matrix3Xr& positions) const {
        return E(positions);
    }
    const real E(const Matrix3Xr& q) const;
    const Matrix3Xr GetElasticForce(const Matrix3Xr& positions) const {
        return fe(positions);
    }
    const Matrix3Xr fe(const Matrix3Xr& q) const;
    const SparseMatrixXr GetStiffnessMatrix(
        const Matrix3Xr& positions) const {
        return K(positions);
    }
    const SparseMatrixXr K(const Matrix3Xr& q) const;
    // Check E, fe, and K.
    void CheckElastics(const Matrix3Xr& positions, const Options& opt) const;

    static const real ComputeFiniteElementVolume(
        const Eigen::Matrix<real, 3, 4>& vertices);

    // Finite element information.
    integer vertex_num_;
    Matrix3Xr rest_vertices_;
    integer element_num_;
    Matrix4Xi elements_;
    std::shared_ptr<material::Material> material_;
    VectorXr masses_;
    VectorXr volumes_;
    // The dim x dim matrix inverse needed to compute F.
    std::vector<Matrix3r> rest_inverse_;

    // Rig information.
    std::shared_ptr<rig::Rig> rig_;

    // State information.
    Matrix3Xr positions_;
    Matrix3Xr velocities_;
    VectorXr rig_parameters_;
    Matrix3Xr forces_;
};

const Matrix3Xi GetSurfaceElements(const Matrix4Xi& elements);

}
}

#endif