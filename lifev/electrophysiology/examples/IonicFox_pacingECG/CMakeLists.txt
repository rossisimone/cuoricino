
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  IonicFoxModel_pacingECG
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_IonicFoxModelpacingECG_data
  CREATE_SYMLINK
  SOURCE_FILES data MonodomainSolverParamList.xml FoxParameters.xml SolverParamList2.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

