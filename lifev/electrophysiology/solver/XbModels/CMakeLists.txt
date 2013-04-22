SET( solver_XbModels_HEADERS
  solver/XbModels/ElectroXbModel.hpp
  solver/XbModels/XbNegroniLascano96.hpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_XbModels_HEADERS})
  STRING(REPLACE "solver/XbModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/solver/XbModels"
  COMPONENT ${PACKAGE_NAME}
)
