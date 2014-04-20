
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  example_10TusscherwStimPro
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_0DTTStimPro_data
  CREATE_SYMLINK
  SOURCE_FILES MonodomainSolverParamList.xml StimulationParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

