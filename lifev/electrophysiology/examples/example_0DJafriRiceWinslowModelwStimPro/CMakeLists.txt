
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  test_0DJafriRiceWinslowModelwStimPro
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DJafriRiceWinslowModelwStimPro_data
  CREATE_SYMLINK
  SOURCE_FILES JafriRiceWinslowParameters.xml StimulationParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


