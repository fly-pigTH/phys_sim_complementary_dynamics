#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "sim/include/simulator.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

namespace phys_sim_complementary_dynamics {

void TestFish() {
    // The fish data comes from the C++ codebase of Complementary Dynamics:
    // https://github.com/seungbaebang/complementary-dynamics-cpp.
    PrintInfo("phys_sim_complementary_dynamics::sim::TestFish",
        "TestFish begins...");

    // Load fish mesh.
    const Matrix3Xr fish_vertices =
        LoadMatrix<real>("asset/mesh/fish.vert").transpose();
    const Matrix4Xi fish_elements =
        LoadMatrix<integer>("asset/mesh/fish.ele").transpose();
    const Matrix3Xi fish_triangles = sim::GetSurfaceElements(fish_elements);
    const integer vertex_num = static_cast<integer>(fish_vertices.cols());
    const integer element_num = static_cast<integer>(fish_elements.cols());
    const integer triangle_num = static_cast<integer>(fish_triangles.cols());
    std::cout << "Mesh information: " << vertex_num << " vertices, "
        << element_num << " tetrahedrons, "
        << triangle_num << " triangles." << std::endl;
    const Vector3r bbox_min = fish_vertices.rowwise().minCoeff();
    const Vector3r bbox_max = fish_vertices.rowwise().maxCoeff();
    std::cout << "Bounding box: (" << bbox_min.transpose() << "), ("
        << bbox_max.transpose() << ")" << std::endl;

    // Initialize the simulator.
    const real density = 1e3;
    const real youngs_modulus = 5e4;
    const real poissons_ratio = 0.4;
    const real time_step = 4e-2;
    const real max_time = 5;
    Options sim_opt;
    sim_opt.integer_option()["max_newton_iter"] = 20;
    sim_opt.integer_option()["max_ls_iter"] = 10;
    sim_opt.bool_option()["grad_check"] = false;
    sim_opt.real_option()["grad_check_abs_tol"] = 1e-4;
    sim_opt.real_option()["grad_check_rel_tol"] = 1e-2;
    sim_opt.real_option()["jac_check_abs_tol"] = 1e-4;
    sim_opt.real_option()["jac_check_rel_tol"] = 1e-2;
    sim_opt.integer_option()["verbose"] = 2;
    sim::Simulator sim;
    sim.AddFiniteElements(fish_vertices, fish_elements,
        material::MaterialType::kCorotated, density, youngs_modulus,
        poissons_ratio);
    const integer bone_num = 2;
    Options rig_opt;
    Matrix3Xr bone_positions = Matrix3Xr::Zero(3, bone_num);
    bone_positions.col(0) = Vector3r(0.0, 0.294604370396, 1.21966209344);
    bone_positions.col(1) = Vector3r(0.0, 1.75584204756, -3.77682802848);
    rig_opt.matrix_option()["samples"] = bone_positions;
    sim.AddRig(rig::RigType::kLinearBlendSkinning, bone_num * 12, rig_opt);

    const auto p_at_t = [&](const real t) -> const VectorXr {
        // Let (Rs, ts) be the rotation and translation applied to a sample s.
        // x = Rs * (X - sample) + ts + sample;
        //   = Rs * X + ts + sample - Rs * sample.
        const Vector3r y = Vector3r::UnitY();
        const Matrix3r R0 = Eigen::AngleAxis<real>(
            std::sin(2 * Pi() * t) / 3, y).toRotationMatrix();
        const Vector3r t0(0, std::sin(Pi() * t) * 0.5 + 0.2 * t, 2 * t);
        const Matrix3r R1 = Eigen::AngleAxis<real>(
            std::sin(2 * Pi() * (t - 0.2)) / 3, y).toRotationMatrix();
        const Vector3r t1(0, 0, 2 * t);

        VectorXr p = VectorXr::Zero(2 * 12);
        p.segment<9>(0) = R0.reshaped();
        p.segment<3>(9) = t0;
        p.segment<9>(12) = R1.reshaped();
        p.segment<3>(12 + 9) = t1;
        return p;
    };

    // Visualize the mesh.
    polyscope::init();
    polyscope::view::setUpDir(polyscope::UpDir::YUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::ZFront);
    polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Free);
    polyscope::view::lookAt({ 10.0, 12.0, 6.0 }, { 0.0, 0.0, 3.0 });
    // For taking screenshots.
    polyscope::options::screenshotExtension = ".jpg";

    // Scene extent.
    polyscope::options::automaticallyComputeSceneExtents = true;
    polyscope::state::lengthScale = 1.;
    polyscope::state::boundingBox =
        std::tuple<glm::vec3, glm::vec3>{ { -12.0, -4.0, -10.0 },
            { 12.0, 5.5, 5.0 } };

    polyscope::removeAllStructures();
    // Register surface mesh structures.
    polyscope::registerSurfaceMesh("fish", fish_vertices.transpose(),
        fish_triangles.transpose());
    polyscope::getSurfaceMesh("fish")
        ->setSurfaceColor({ 0x4C / 255.0, 0xC9 / 255.0, 0xFE / 255.0 });
    polyscope::getSurfaceMesh("fish")->setEdgeWidth(1.0);

    // The following UI design is borrowed from the SIGGRAPH course:
    // https://github.com/siggraphcontact/rigidBodyTutorial.
    bool paused = true;
    bool step_once = false;
    bool skip_complement = false;
    real time_so_far = 0;
    const auto callback = [&]() -> void {
        ImGui::PushItemWidth(100);

        if (ImGui::Button("Reset")) {
            sim.ClearForces();
            sim.set_q(fish_vertices);
            sim.set_q_dot(Matrix3Xr::Zero(3, sim.vertex_num()));
            sim.set_p(p_at_t(0));
            time_so_far = 0;
            paused = true;
            step_once = false;
        }
        ImGui::SameLine();
        if (ImGui::Button("Step once")) step_once = true;

        ImGui::Checkbox("Pause", &paused);
        ImGui::SameLine();
        ImGui::Checkbox("Skip complement", &skip_complement);
        ImGui::Text("Time: %.2f seconds.", time_so_far);
        if ((!paused || step_once) && time_so_far <= max_time) {
            Tic();
            sim.ClearForces();
            sim.set_p(p_at_t(time_so_far));
            sim_opt.bool_option()["skip_complement"] = skip_complement;
            sim.Step(time_step, sim_opt);
            time_so_far += time_step;

            step_once = false;
        }
        // Visualization.
        polyscope::getSurfaceMesh("fish")->updateVertexPositions(
            sim.q().transpose());

        ImGui::PopItemWidth();

        // Call the lines below to take a screenshot.
        /*
        {
            // The code below is modified from polyscope::screenshot.
            char buff[50];
            snprintf(buff, 50, "video/screenshot_%06zu%s",
                polyscope::state::screenshotInd,
                polyscope::options::screenshotExtension.c_str());
            polyscope::screenshot(std::string(buff), true);
            ++polyscope::state::screenshotInd;
        }
        */
    };

    polyscope::state::userCallback = callback;
    polyscope::show();

    PrintInfo("phys_sim_complementary_dynamics::sim::TestFish",
        "TestFish ends...");
}

}

int main() {
    phys_sim_complementary_dynamics::TestFish();
    return 0;
}