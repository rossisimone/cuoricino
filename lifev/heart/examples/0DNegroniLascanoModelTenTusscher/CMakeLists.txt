
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  0DNegroniLascanoModelTenTusscher
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_0DNegroniLascanoModelTenTusscher_data
  CREATE_SYMLINK
  SOURCE_FILES NegroniLascano96Parameters.xml TenTusscherParameters.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

