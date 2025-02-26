export PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT=$PWD
export PHYS_SIM_COMPLEMENTARY_DYNAMICS_BUILD_ROOT=$PWD/build/phys_sim_complementary_dynamics
mkdir -p $PHYS_SIM_COMPLEMENTARY_DYNAMICS_BUILD_ROOT

cd $PHYS_SIM_COMPLEMENTARY_DYNAMICS_BUILD_ROOT
CMAKE_OPTIONS="-DPHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT=${PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT}"
cmake ${CMAKE_OPTIONS} ../../cpp
make -j$(nproc)

cd $PHYS_SIM_COMPLEMENTARY_DYNAMICS_ROOT
./build/phys_sim_complementary_dynamics/main
