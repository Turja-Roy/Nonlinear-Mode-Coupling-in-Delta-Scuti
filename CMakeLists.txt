cmake_minimum_required(VERSION 3.16)
project(stellar_coupling CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra)

find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(OpenMP)

# CLI11 and HighFive are vendored under third_party/ -- both header-only, and
# pinning them keeps a cluster build from needing anything installed.
add_library(coupling STATIC
    src/numeric.cpp src/angular.cpp src/model.cpp src/kappa.cpp
    src/triplets.cpp src/amplitude.cpp src/stability.cpp src/csv.cpp)
# Vendored, so the build does not depend on the cluster's Eigen or Boost.
target_include_directories(coupling PUBLIC src third_party third_party/eigen3 ${HDF5_INCLUDE_DIRS})
target_link_libraries(coupling PUBLIC GSL::gsl ${HDF5_C_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(coupling PUBLIC OpenMP::OpenMP_CXX)
endif()

foreach(app three_mode_search channels four_mode_search mixed_network make_inlists mw25)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE coupling)
endforeach()

enable_testing()
foreach(t numeric kappa amplitude)
    add_executable(test_${t} tests/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE coupling)
    # the checks are assert()s, and Release defines NDEBUG -- without this they
    # compile out and every test passes without testing anything
    target_compile_options(test_${t} PRIVATE -UNDEBUG)
    add_test(NAME ${t} COMMAND test_${t})
endforeach()
