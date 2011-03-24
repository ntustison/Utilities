PROJECT(TendonRegistration)

# Find ITK
INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF(USE_ITK_FILE)

SET(DeformableRegistration1
  DeformableRegistration1.cxx
)

# Header and template source files for project
# we can only include them if compiling on windows
IF(NOT UNIX)
IF(WIN32)
SET(DeformableRegistration1 ${DeformableRegistration1}
)
ENDIF(WIN32)
ENDIF(NOT UNIX)

LINK_DIRECTORIES (
  ${ITK_BINARY_DIR}/Code/Numerics/FEM
  ${ITK_BINARY_DIR}/Code/Numerics
  ${ITK_BINARY_DIR}/Code/Numerics/vxl
  ${ITK_BINARY_DIR}/Code/IO
  ${ITK_BINARY_DIR}/Utilities/png
  ${ITK_BINARY_DIR}/Utilities/zlib
  ${ITK_BINARY_DIR}/Utilities/MetaIO
  ${ITK_BINARY_DIR}/Code/BasicFilters
  ${ITK_BINARY_DIR}/Code/Common
)

ADD_EXECUTABLE(DeformableRegistration1 ${DeformableRegistration1})

TARGET_LINK_LIBRARIES(DeformableRegistration1 ITKCommon ITKIO ITKNumerics ITKAlgorithms ITKFEM ITKBasicFilters)


