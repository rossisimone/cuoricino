
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
 
TRIBITS_ADD_EXECUTABLE(
  IntegratedHeart
  NAME IntegratedHeart
  SOURCES main.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(IH_mesh  
  SOURCE_FILES leftventricle.mesh biventricle2.mesh
  DEST_FILES Mesh/leftventricle.mesh Mesh/biventricle2.mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/mesh/inria
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(IH_data
  CREATE_SYMLINK
  SOURCE_FILES Run_IntegratedHeart.dat biventricle2fibers.dat AlievPanfilovParameters.xml MonodomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/multiscale/IntegratedHeart
)

#TRIBITS_COPY_FILES_TO_BINARY_DIR(IH_fibers
#  CREATE_SYMLINK
#  SOURCE_FILES biventricle2fibers.dat
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/multiscale/IntegratedHeart
#)

TRIBITS_COPY_FILES_TO_BINARY_DIR(IH_multiscale
  CREATE_SYMLINK
  SOURCE_FILES MultiscaleDatabase FunctionsDatabase
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/multiscale/IntegratedHeart
)


#TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_mesh
#  CREATE_SYMLINK
#  SOURCE_FILES inria
#  DEST_FILES Mesh
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data/mesh
#)
#
#TRIBITS_COPY_FILES_TO_BINARY_DIR(MS_functions_database
#  CREATE_SYMLINK
#  SOURCE_FILES FunctionsDatabase MultiscaleDatabase
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
#)

