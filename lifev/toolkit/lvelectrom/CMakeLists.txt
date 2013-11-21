
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  lvelectrom
  SOURCES lvelectrom.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_defaults_lvelectrom
  CREATE_SYMLINK
  SOURCE_FILES MonodomainSolverParamList.xml MonodomainSolverPreconditioner
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../defaults/
)

