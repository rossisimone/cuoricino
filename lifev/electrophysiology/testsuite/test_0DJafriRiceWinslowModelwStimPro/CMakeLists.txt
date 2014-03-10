
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_0DJafriRiceWinslowModelwStimPro
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DJafriRiceWinslowModelwStimPro_data
  CREATE_SYMLINK
  SOURCE_FILES JafriRiceWinslowParameters.xml StimulationParameters.xml data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

