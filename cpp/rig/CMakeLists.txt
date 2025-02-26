file(GLOB_RECURSE rig_cpp_files "./src/*.cpp")

add_library(rig SHARED ${rig_cpp_files})
target_link_libraries(rig basic)

# Add a test project.
add_executable(test_rig "./main.cpp")
target_link_libraries(test_rig rig)