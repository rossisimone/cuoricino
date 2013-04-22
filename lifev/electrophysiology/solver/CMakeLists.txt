ADD_SUBDIRECTORIES(
	XbModels
	IonicModels
)
SET(solver_HEADERS
#  solver/ElectroAlievPanfilov.hpp
  solver/ElectroMonodomainData.hpp
  solver/ElectroStiffnessFibers.hpp
  solver/ElectroIonicData.hpp
  solver/ElectroBidomainData.hpp
  solver/ElectroIonicSolver.hpp
  solver/ElectroLuoRudy.hpp
  solver/ElectroBidomainSolver.hpp
  solver/ElectroFunctors.hpp
  solver/ElectroMonodomainSolver.hpp
  solver/ElectroETAMonodomainSolver.hpp
  solver/ElectroUpdatedMonodomainSolver.hpp
#  solver/ElectroODEModel.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/ElectroIonicData.cpp
  solver/ElectroMonodomainData.cpp
  solver/ElectroBidomainData.cpp
  solver/ElectroFunctors.cpp
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
