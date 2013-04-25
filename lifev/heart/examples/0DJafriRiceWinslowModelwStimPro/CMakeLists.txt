
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  0DTenTusscher
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_0DTenTusscherModel_data
  CREATE_SYMLINK
  SOURCE_FILES TenTusscherParameters.xml StimulationParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

