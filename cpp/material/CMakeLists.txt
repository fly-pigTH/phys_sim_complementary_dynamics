file(GLOB_RECURSE material_cpp_files "./src/*.cpp")

add_library(material SHARED ${material_cpp_files})
target_link_libraries(material basic)

# Add a test project.
add_executable(test_material "./main.cpp")
target_link_libraries(test_material material)