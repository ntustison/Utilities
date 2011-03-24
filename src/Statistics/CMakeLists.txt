PROJECT(PermuationTest)

#  
# Find ITK 
# 
INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 
 

LINK_DIRECTORIES(
  ${Local_Libs}
  ${ITK_LIBRARY_PATH}
)

ADD_EXECUTABLE(PermutationTests PermutationTests.cxx )

TARGET_LINK_LIBRARIES(PermutationTests ITKNumerics ITKIO ITKCommon ITKStatistics)
