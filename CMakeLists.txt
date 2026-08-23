cmake_minimum_required(VERSION 3.21)

project(
    light_matter_heatpipe_tube
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------
# C++QED
# ------------------------------------------------------------

add_subdirectory(
    qed-submodule/cppqed
    cppqed-build
)

# ------------------------------------------------------------
# threeL_level_master_test
# ------------------------------------------------------------

find_package(Eigen3 REQUIRED)

add_executable(
    three_level_master_test
    src/three_level_master_test.cc
)

target_link_libraries(
    three_level_master_test
    PRIVATE
        Eigen3::Eigen
)

# ------------------------------------------------------------
# Smoke test
# ------------------------------------------------------------

add_executable(
    cppqed_test
    tests/cppqed_test.cc
)

target_link_libraries(
    cppqed_test
    PRIVATE
        CPPQED::CPPQEDelements
        CPPQED::CPPQEDcore
        CPPQED::CPPQEDutils
)