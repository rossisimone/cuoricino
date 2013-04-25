
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  0DJafriRiceWinslowwStimPro
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_0DJafriRiceWinslowModelwStimPro_data
  CREATE_SYMLINK
  SOURCE_FILES JafriRiceWinslowParameters.xml StimulationParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

