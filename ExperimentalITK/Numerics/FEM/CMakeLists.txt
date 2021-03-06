# Source files for the FEM project

SET(FEM_SRCS 

  itkFEMElementBase.cxx

  itkFEMElement2DC0LinearLine.cxx
  itkFEMElement2DC0LinearLineStress.cxx
  itkFEMElement2DC1Beam.cxx
  itkFEMElement2DC0LinearQuadrilateral.cxx
  itkFEMElement2DC0LinearQuadrilateralLinearElasticity.cxx
  itkFEMElement2DC0LinearQuadrilateralMembrane.cxx
  itkFEMElement2DC0LinearQuadrilateralStress.cxx
  itkFEMElement2DC0LinearQuadrilateralStrain.cxx
  itkFEMElement2DC0LinearTriangular.cxx
  itkFEMElement2DC0LinearTriangularMembrane.cxx
  itkFEMElement2DC0LinearTriangularStress.cxx
  itkFEMElement2DC0LinearTriangularStrain.cxx
  itkFEMElement2DC0QuadraticTriangular.cxx
  itkFEMElement2DC0QuadraticTriangularStress.cxx
  itkFEMElement2DC0QuadraticTriangularStrain.cxx
  itkFEMElement3DC0LinearHexahedron.cxx
  itkFEMElement3DC0LinearHexahedronLinearElasticity.cxx
  itkFEMElement3DC0LinearHexahedronMembrane.cxx
  itkFEMElement3DC0LinearHexahedronStrain.cxx
  itkFEMElement3DC0LinearTetrahedron.cxx
  itkFEMElement3DC0LinearTetrahedronMembrane.cxx
  itkFEMElement3DC0LinearTetrahedronStrain.cxx
  
  itkFEMBSplineShapeFunctions.cxx
  itkFEMElement2DBSplinePatch.cxx
  itkFEMElement3DBSplinePatch.cxx
  itkFEMElement2DBSplinePatchLinearElasticity.cxx
  itkFEMElement2DBSplinePatchMembrane.cxx
  itkFEMElement2DBSplinePatchStrain.cxx
  itkFEMElement2DBSplinePatchStress.cxx
  itkFEMElement3DBSplinePatchLinearElasticity.cxx
  itkFEMElement3DBSplinePatchMembrane.cxx
  itkFEMElement3DBSplinePatchStrain.cxx

  itkFEMMaterialBase.cxx
  itkFEMMaterialLinearElasticity.cxx

  itkFEMLoadBase.cxx
  itkFEMLoadBC.cxx
  itkFEMLoadBCMFC.cxx
  itkFEMLoadNode.cxx
  itkFEMLoadElementBase.cxx
  itkFEMLoadGrav.cxx
  itkFEMLoadEdge.cxx
  itkFEMLoadPoint.cxx
  itkFEMLoadLandmark.cxx
 
  itkFEMLoadImplementationsRegister.cxx
  itkFEMLoadImplementationGenericBodyLoad.cxx
  itkFEMLoadImplementationGenericLandmarkLoad.cxx

  itkFEMSolver.cxx
  itkFEMSolverHyperbolic.cxx
  itkFEMSolverCrankNicolson.cxx

  itkFEMGenerateMesh.cxx
  itkFEMGenerateBSplineMesh.cxx

  itkFEMLinearSystemWrapper.cxx
  itkFEMLinearSystemWrapperVNL.cxx
  itkFEMLinearSystemWrapperDenseVNL.cxx
  itkFEMLinearSystemWrapperItpack.cxx
  itkFEMItpackSparseMatrix.cxx

  itkFEMLightObject.cxx
  itkFEMException.cxx
  itkFEMUtility.cxx
  itkFEMInitialization.cxx

)



# Header and template source files for FEM project
# we can only include them if compiling on windows
IF(NOT UNIX)
IF(WIN32)
SET(FEM_SRCS ${FEM_SRCS}

  itkFEMElements.h
  itkFEMElementBase.h
  itkFEMElementStd.h

  itkFEMElement1DStress.h
  itkFEMElement2DLinearElasticity.h
  itkFEMElement2DStress.h
  itkFEMElement2DStrain.h
  itkFEMElement2DMembrane.h
  itkFEMElement3DLinearElasticity.h
  itkFEMElement3DMembrane.h
  itkFEMElement3DStrain.h
  itkFEMElement2DC0LinearLine.h
  itkFEMElement2DC0LinearLineStress.h
  itkFEMElement2DC1Beam.h
  itkFEMElement2DC0LinearQuadrilateral.h
  itkFEMElement2DC0LinearQuadrilateralLinearElasticity.h
  itkFEMElement2DC0LinearQuadrilateralMembrane.h
  itkFEMElement2DC0LinearQuadrilateralStress.h
  itkFEMElement2DC0LinearQuadrilateralStrain.h
  itkFEMElement2DC0LinearTriangular.h
  itkFEMElement2DC0LinearTriangularMembrane.h
  itkFEMElement2DC0LinearTriangularStress.h
  itkFEMElement2DC0QuadraticTriangular.h
  itkFEMElement2DC0QuadraticTriangularStress.h
  itkFEMElement3DC0LinearHexahedron.h
  itkFEMElement3DC0LinearHexahedronLinearElasticity.h
  itkFEMElement3DC0LinearHexahedronMembrane.h
  itkFEMElement3DC0LinearHexahedronStrain.h  
  itkFEMElement3DC0LinearTetrahedron.h
  itkFEMElement3DC0LinearTetrahedronMembrane.h
  itkFEMElement3DC0LinearTetrahedronStrain.h  

  itkFEMBSplineShapeFunctions.h
  itkFEMElement2DBSplinePatch.h
  itkFEMElement2DBSplinePatchLinearElasticity.h
  itkFEMElement2DBSplinePatchMembrane.h
  itkFEMElement2DBSplinePatchStrain.h
  itkFEMElement2DBSplinePatchStress.h
  itkFEMElement3DBSplinePatch.h
  itkFEMElement3DBSplinePatchLinearElasticity.h
  itkFEMElement3DBSplinePatchMembrane.h
  itkFEMElement3DBSplinePatchStrain.h

  itkFEMElementStd.hxx
  itkFEMElement1DStress.hxx
  itkFEMELement2DLinearElasticity.hxx
  itkFEMElement2DMembrane.hxx
  itkFEMElement2DStrain.hxx
  itkFEMElement2DStress.hxx
  itkFEMELement3DLinearElasticity.hxx
  itkFEMElement3DMembrane.hxx
  itkFEMElement3DStrain.hxx

  itkFEMMaterials.h
  itkFEMMaterialBase.h
  itkFEMMaterialLinearElasticity.h

  itkFEMLoads.h
  itkFEMLoadBase.h
  itkFEMLoadBC.h
  itkFEMLoadBCMFC.h
  itkFEMLoadNode.h
  itkFEMLoadElementBase.h
  itkFEMLoadGrav.h
  itkFEMLoadEdge.h
  itkFEMLoadPoint.h
  itkFEMLoadLandmark.h
  itkFEMImageMetricLoad.h

  itkFEMLoadImplementationGenericBodyLoad.h
  itkFEMLoadImplementationGenericLandmarkLoad.h
  itkFEMLoadImplementationTest.h
  itkFEMImageMetricLoadImplementation.h

  itkFEMLoadTest.h

  itkFEMSolver.h
  itkFEMSolverHyperbolic.h
  itkFEMSolverCrankNicolson.h

  itkFEMGenerateMesh.h
  itkFEMGenerateBSplineMesh.h

  itkFEMSolution.h

  itkFEMLinearSystemWrappers.h
  itkFEMLinearSystemWrapper.h
  itkFEMLinearSystemWrapperVNL.h
  itkFEMLinearSystemWrapperDenseVNL.h
  itkFEMLinearSystemWrapperItpack.h
  itkFEMItpackSparseMatrix.h
  
  itkFEM.h
  itkFEMLightObject.h
  itkFEMP.h
  itkFEMPArray.h
  itkFEMObjectFactory.h
  itkFEMMacro.h
  itkFEMUtility.h
  itkFEMException.h
  itkFEMInitialization.h

  itkVisitorDispatcher.h

  itkFEMImageMetricLoad.hxx
  
  ../../Algorithms/itkFEMRegistrationFilter.h
  ../../Algorithms/itkFEMRegistrationFilter.hxx

)
ENDIF(WIN32)
ENDIF(NOT UNIX)

IF(CMAKE_BUILD_TOOL MATCHES ".*dev")
  SET(FEM_SRCS ${FEM_SRCS};ReadMe.txt;ToDo.txt)
ENDIF(CMAKE_BUILD_TOOL MATCHES ".*dev")

  
# Define source groups that look nice in DSP files
SOURCE_GROUP("Elements" REGULAR_EXPRESSION itkFEMElement*)
SOURCE_GROUP("Nodes" REGULAR_EXPRESSION itkFEMNode*)
SOURCE_GROUP("Loads" REGULAR_EXPRESSION itkFEMLoad*)
SOURCE_GROUP("Materials" REGULAR_EXPRESSION itkFEMMaterial*)
SOURCE_GROUP("Linear Systems" REGULAR_EXPRESSION itkFEMLinearSystem*)
SOURCE_GROUP("Source Files Template" REGULAR_EXPRESSION ./*.hxx)

INCLUDE_DIRECTORIES(${VXL_NETLIB_INCLUDE_DIR}) # Some .cxx files need f2c.h.

# Build the library
ADD_LIBRARY(ITKFEM ${FEM_SRCS})
TARGET_LINK_LIBRARIES(ITKFEM 
                      ITKNumerics
                      ITKCommon
                      ITKBasicFilters
                      ITKIO)
IF(UNIX)
  TARGET_LINK_LIBRARIES (ITKFEM -lm )
ENDIF(UNIX)

INSTALL_TARGETS(/lib/InsightToolkit ITKFEM)
INSTALL_FILES(/include/InsightToolkit/Numerics/FEM "(\\.h|\\.hxx)$")
