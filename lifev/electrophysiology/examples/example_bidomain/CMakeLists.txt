
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_bidomain
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_test_bidomain_data
  CREATE_SYMLINK
  SOURCE_FILES BidomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(NiedererBM11.mesh_test_bidomain
  SOURCE_FILES benchmark_05mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/inria/
)
