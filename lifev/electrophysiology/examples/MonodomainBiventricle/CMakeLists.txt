
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  MonodomainBiventricle
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_MonodomainBiventricle
  CREATE_SYMLINK
  SOURCE_FILES data MonodomainSolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

