
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  IonicFitzHughNagumoModel
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_IonicFitzHughNagumoModel_data
  CREATE_SYMLINK
  SOURCE_FILES data MonodomainSolverParamList.xml FitzHughNagumoParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

