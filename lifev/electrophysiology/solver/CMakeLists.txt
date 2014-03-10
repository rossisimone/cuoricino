ADD_SUBDIRECTORIES(
	XbModels
	IonicModels
)
SET(solver_HEADERS
  solver/ElectroIonicData.hpp
  solver/ElectroIonicSolver.hpp
  solver/ElectroFunctors.hpp
  solver/ElectroETAMonodomainSolver.hpp
  solver/ElectroETABidomainSolver.hpp
  solver/StimulationProtocol.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/ElectroIonicData.cpp
  solver/ElectroFunctors.cpp
  solver/StimulationProtocol.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/solver"
  COMPONENT ${PACKAGE_NAME}
)
