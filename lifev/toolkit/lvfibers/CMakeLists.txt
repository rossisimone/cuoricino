
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  lvfibers
  SOURCES lvfibers.cpp  
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_defaults_lvfibers
  CREATE_SYMLINK
  SOURCE_FILES FiberGenerationParamList.xml FiberGenerationPreconditioner
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../defaults/
)

