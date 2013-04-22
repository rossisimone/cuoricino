
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  test_monodomain
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_test_monodomain_data
  CREATE_SYMLINK
  SOURCE_FILES data MonodomainSolverParamList.xml NegroniLascano96Parameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

