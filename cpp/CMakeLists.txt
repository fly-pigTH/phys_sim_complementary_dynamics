cmake_minimum_required(VERSION 3.20)

project(phys_sim_complementary_dynamics)

# Project-wide setup.
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

# Build external libraries.
add_subdirectory("${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/external/polyscope"
    "${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/build/polyscope")

# Configure include and lib paths.
include_directories("${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/cpp")
include_directories("${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/external")
include_directories("${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/external/eigen")
include_directories("${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}/external/polyscope/include")

# Compile each library.
# Dependency map:
# basic <- empty;
add_subdirectory("basic")
# material <- basic;
add_subdirectory("material")
# rig <- basic;
add_subdirectory("rig")
# sim <- rig, material, basic;
add_subdirectory("sim")

# Compile the main program.
add_executable(main "main.cpp")
target_link_libraries(main sim polyscope)