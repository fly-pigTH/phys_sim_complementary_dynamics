file(GLOB_RECURSE basic_cpp_files "./src/*.cpp")

add_library(basic SHARED ${basic_cpp_files})

# Add a test project.
add_executable(test_basic "./main.cpp")
target_link_libraries(test_basic basic)