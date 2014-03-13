
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_0DJafriRiceWinslow
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DJafriRiceWinslowModel_data
  CREATE_SYMLINK
  SOURCE_FILES JafriRiceWinslowParameters.xml 
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

