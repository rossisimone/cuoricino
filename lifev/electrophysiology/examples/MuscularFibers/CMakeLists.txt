
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  MuscularFibers
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_MuscularFibers_data
  CREATE_SYMLINK
  SOURCE_FILES data ParamList.xml 
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

