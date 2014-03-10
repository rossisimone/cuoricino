
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  0DNegroniLascanoModelJafriRiceWinslow
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_0DNegroniLascanoModelJafriRiceWinslow_data
  CREATE_SYMLINK
  SOURCE_FILES NegroniLascano96Parameters.xml JafriRiceWinslowParameters.xml StimulationParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


