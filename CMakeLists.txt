cmake_minimum_required(VERSION 3.24)

project(light-matter-heatpipe-tube
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersGeneral
    FiltersModeling
    FiltersSources
    IOGeometry
)

add_executable(domain_model_test
    src/domain_model.cc
)

target_link_libraries(domain_model_test
    PRIVATE
        ${VTK_LIBRARIES}
)

vtk_module_autoinit(
    TARGETS domain_model_test
    MODULES ${VTK_LIBRARIES}
)