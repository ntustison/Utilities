PROJECT(ttest)

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

ADD_EXECUTABLE(ttest ttest.cxx )

TARGET_LINK_LIBRARIES(ttest ITKNumerics ITKIO ITKCommon ITKStatistics)
