file(GLOB_RECURSE sim_cpp_files "./src/*.cpp")

add_library(sim SHARED ${sim_cpp_files})
target_link_libraries(sim rig material basic)

# Add a test project.
add_executable(test_sim "./main.cpp")
target_link_libraries(test_sim sim)