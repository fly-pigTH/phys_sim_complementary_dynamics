#include "sim/include/simulator.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "material/include/linear_material.hpp"
#include "material/include/st_venant_kirchhoff_material.hpp"
#include "material/include/corotated_material.hpp"
#include "material/include/neohookean_material.hpp"
#include "rig/include/linear_blend_skinning_rig.hpp"

namespace phys_sim_complementary_dynamics {
namespace sim {

Simulator::Simulator()
    : vertex_num_(0), rest_vertices_(3, 0), element_num_(0),
    elements_(4, 0), material_(nullptr), masses_(0), volumes_(0),
    rest_inverse_(0), rig_(nullptr),
    positions_(3, 0), velocities_(3, 0), rig_parameters_(0), forces_(3, 0) {}

void Simulator::AddFiniteElements(const Matrix3Xr& rest_vertices,
    const Matrix4Xi& elements, const material::MaterialType material_type,
    const real density, const real youngs_modulus, const real poissons_ratio) {

    const std::string error_location = "sim::Simulator::AddFiniteElements";

    vertex_num_ = static_cast<integer>(rest_vertices.cols());
    CheckCondition(vertex_num_ > 3, error_location, "Too few vertices.");
    rest_vertices_ = rest_vertices;

    element_num_ = static_cast<integer>(elements.cols());
    CheckCondition(element_num_ > 0, error_location, "Too few elements.");
    elements_ = elements;

    // Create material.
    if (material_type == material::MaterialType::kLinear) {
        material_ = std::make_shared<material::LinearMaterial>();
    } else if (material_type == material::MaterialType::kStVenantKirchhoff) {
        material_ = std::make_shared<material::StVenantKirchhoffMaterial>();
    } else if (material_type == material::MaterialType::kCorotated) {
        material_ = std::make_shared<material::CorotatedMaterial>();
    } else if (material_type == material::MaterialType::kNeohookean) {
        material_ = std::make_shared<material::NeohookeanMaterial>();
    } else {
        CheckCondition(false, error_location, "Unsupported material type: "
            + material::ToString(material_type) + ".");
    }
    material_->Initialize(density, youngs_modulus, poissons_ratio);

    // Compute masses and volumes.
    masses_ = VectorXr::Zero(vertex_num_);
    volumes_ = VectorXr::Zero(element_num_);
    rest_inverse_.resize(element_num_);
    for (integer j = 0; j < element_num_; ++j) {
        // Obtain the volume of the j-th element.
        const Eigen::Matrix<real, 3, 4> v =
            rest_vertices_(Eigen::all, elements_.col(j));
        rest_inverse_[j] = (v.rightCols(3).colwise() - v.col(0)).inverse();
        volumes_(j) = ComputeFiniteElementVolume(v);
        const real vol_per_vertex = volumes_(j) / 4;
        for (integer i = 0; i < 4; ++i) {
            masses_(elements_(i, j)) += vol_per_vertex;
        }
    }
    // Multiple density.
    masses_ *= material_->density();

    const Matrix3Xr zeros = Matrix3Xr::Zero(3, vertex_num_);
    positions_ = rest_vertices;
    velocities_ = zeros;
    forces_ = zeros;

}

void Simulator::AddRig(const rig::RigType type, const integer parameter_num,
    const Options& opt) {
    const std::string error_location = "sim::Simulator::AddRig";
    if (type == rig::RigType::kLinearBlendSkinning) {
        rig_ = std::make_shared<rig::LinearBlendSkinningRig>();
    } else {
        CheckCondition(false, error_location, "Unsupported rig type: "
            + rig::ToString(type) + ".");
    }
    rig_->Initialize(parameter_num, rest_vertices_, elements_, opt);
    rig_parameters_ = VectorXr::Zero(parameter_num);
}

void Simulator::set_q(const Matrix3Xr& q) {
    CheckPositions(q, "set_q");
    positions_ = q;
}

void Simulator::set_q_dot(const Matrix3Xr& q_dot) {

    std::stringstream ss;
    ss << "Incompatible velocity column number (expected): "
        << q_dot.cols() << " (" << vertex_num_ << ").";
    CheckCondition(static_cast<integer>(q_dot.cols()) == vertex_num_,
        "sim::Simulator::set_q_dot", ss.str());

    velocities_ = q_dot;
}

void Simulator::set_p(const VectorXr& p) {

    std::stringstream ss;
    ss << "Incompatible rig parameter size (expected): "
        << p.size() << " (" << rig_->m() << ").";
    CheckCondition(static_cast<integer>(p.size()) == rig_->m(),
        "sim::Simulator::set_p", ss.str());

    rig_parameters_ = p;
}

void Simulator::set_f(const Matrix3Xr& f) {

    std::stringstream ss;
    ss << "Incompatible force column number (expected): " << f.cols() << " ("
        << vertex_num_ << ").";
    CheckCondition(static_cast<integer>(f.cols()) == vertex_num_,
        "sim::Simulator::set_f", ss.str());

    forces_ = f;
}

void Simulator::CheckVertexIndex(const integer vertex_index,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Out-of-bound vertex_index (vertex_num_): " << vertex_index << " (" <<
        vertex_num_ << ").";
    CheckCondition(0 <= vertex_index && vertex_index < vertex_num_,
        "sim::Simulator::" + error_location, ss.str());
}

void Simulator::CheckPositions(const Matrix3Xr& positions,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Incompatible position column number (expected): " << positions.cols()
        << " (" << vertex_num_ << ").";
    CheckCondition(static_cast<integer>(positions.cols()) == vertex_num_,
        "sim::Simulator::" + error_location, ss.str());
}

void Simulator::CheckVelocities(const Matrix3Xr& velocities,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Incompatible velocity column number (expected): "
        << velocities.cols() << " (" << vertex_num_ << ").";
    CheckCondition(static_cast<integer>(velocities.cols()) == vertex_num_,
        "sim::Simulator::" + error_location, ss.str());   
}

const real Simulator::ComputeFiniteElementVolume(
    const Eigen::Matrix<real, 3, 4>& vertices) {
    const Vector3r v10 = vertices.col(1) - vertices.col(0);
    const Vector3r v20 = vertices.col(2) - vertices.col(0);
    const Vector3r v30 = vertices.col(3) - vertices.col(0);
    return std::abs(v10.cross(v20).dot(v30)) / 6;
}

const real Simulator::E(const Matrix3Xr& q) const {
    CheckPositions(q, "E");

    const auto E_ele = [&](const std::shared_ptr<material::Material>& material,
        const integer element_index) -> const real {

        ////////////////////////////////////////////////////////////////////////
        // Task 2.1 (2 points).
        ////////////////////////////////////////////////////////////////////////
        //
        // Compute the strain energy from the deformation described by vertex
        // positions q. We have provided the for loop that collects the energy
        // from each finite element and leave the computation inside each finite
        // element for you to implement.
        //
        // TODO.
        return 0.0;
    };

    real energy = 0;
    // Loop over all elements.
    for (integer e = 0; e < element_num_; ++e) {
        energy += E_ele(material_, e);
    }

    return energy;
}

const Matrix3Xr Simulator::fe(const Matrix3Xr& q) const {
    CheckPositions(q, "fe");

    ////////////////////////////////////////////////////////////////////////////
    // Task 2.2 (3 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Compute the elastic force from deformed vertex positions q.
    //
    // TODO.
    return Matrix3Xr::Zero(3, vertex_num_);
}

const SparseMatrixXr Simulator::K(const Matrix3Xr& q) const {
    CheckPositions(q, "K");

    std::vector<Eigen::Triplet<real>> nonzeros(3 * 4 * 3 * 4 * element_num_);

    ////////////////////////////////////////////////////////////////////////////
    // Task 2.3 (4 points).
    ////////////////////////////////////////////////////////////////////////////
    //
    // Compute the stiffness matrix from deformed vertex positions q. The rows
    // and columns of K use the following order of vertices:
    //
    // q(0, 0), q(1, 0), q(2, 0), q(0, 1), q(1, 1), q(2, 1), ...
    //
    // The first three rows/columns correspond to the x, y, and z degrees of
    // freedom (Dofs) of the first vertex, the next three rows/columns
    // correspond to the second vertex, and so on.
    //
    // TODO.
    return FromTriplet(3 * vertex_num_, 3 * vertex_num_, nonzeros);
}

void Simulator::CheckElastics(const Matrix3Xr& positions,
    const Options& opt) const {
    CheckPositions(positions, "CheckElastics");

    // The following gradient and Hessian checks passed.
    const auto func = [&](const VectorXr& x) -> const real {
        return E(x.reshaped(3, vertex_num_));
    };
    const auto func_grad = [&](const VectorXr& x) -> const VectorXr {
        return -fe(x.reshaped(3, vertex_num_)).reshaped();
    };
    const auto func_hess = [&](const VectorXr& x) -> const MatrixXr {
        return K(x.reshaped(3, vertex_num_)).toDense();
    };
    CheckGradient(positions.reshaped(), func, func_grad, opt);
    CheckJacobian(positions.reshaped(), func_grad, func_hess, opt);
}

const Matrix3Xi GetSurfaceElements(const Matrix4Xi& elements) {
    std::set<std::array<integer, 3>> surface_elements;
    const integer element_num = static_cast<integer>(elements.cols());
    for (integer e = 0; e < element_num; ++e) {
        for (integer i = 0; i < 4; ++i) {
            std::array<integer, 3> ele;
            for (integer d = 0; d < 3; ++d)
                ele[d] = elements((i + d) % 4, e);
            // Build key to this element.
            std::sort(ele.begin(), ele.end());
            // Check if this key has been visited before.
            if (surface_elements.find(ele) != surface_elements.end()) {
                surface_elements.erase(ele);
            } else {
                surface_elements.insert(ele);
            }
        }
    }

    const integer face_num = static_cast<integer>(surface_elements.size());
    Matrix3Xi faces = Matrix3Xi::Zero(3, face_num);
    integer cnt = 0;
    for (const auto& f : surface_elements) {
        for (integer i = 0; i < 3; ++i)
            faces(i, cnt) = f[i];
        ++cnt;
    }

    return faces;
}

}
}