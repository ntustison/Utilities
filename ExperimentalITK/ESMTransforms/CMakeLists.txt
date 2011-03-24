project(itkESMRegistration)

cmake_minimum_required(VERSION 2.6)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(ESMImageRegistration ESMImageRegistration.cxx)
target_link_libraries(ESMImageRegistration ITKIO ITKNumerics)
