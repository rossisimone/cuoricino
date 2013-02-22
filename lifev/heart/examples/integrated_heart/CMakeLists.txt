
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  IntegratedHeart
  SOURCES integratedHeart.cpp flowConditions.cpp ud_functions.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_integrated_heart_data
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_heart_integrated_heart_mesh
  CREATE_SYMLINK
  SOURCE_FILES integratedHeart_prototype1_fluid.mesh integratedHeart_prototype1_structure.mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../data/mesh/inria
)

