SET( solver_IonicModels_HEADERS
  solver/IonicModels/IonicAlievPanfilov.hpp
  solver/IonicModels/IonicFitzHughNagumo.hpp
  solver/IonicModels/IonicMitchellSchaeffer.hpp
  solver/IonicModels/IonicJafriRiceWinslow.hpp
  solver/IonicModels/IonicTenTusscher.hpp
  solver/IonicModels/IonicMinimalModel.hpp
  solver/IonicModels/ElectroIonicModel.hpp
  solver/IonicModels/IonicFox.hpp
CACHE INTERNAL "")

SET( solver_IonicModels_SOURCES
  solver/IonicModels/IonicAlievPanfilov.cpp
  solver/IonicModels/IonicFitzHughNagumo.cpp
  solver/IonicModels/IonicMitchellSchaeffer.cpp
  solver/IonicModels/IonicJafriRiceWinslow.cpp
  solver/IonicModels/IonicTenTusscher.cpp
  solver/IonicModels/IonicMinimalModel.cpp
  solver/IonicModels/IonicFox.cpp
  solver/IonicModels/ElectroIonicModel.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_IonicModels_HEADERS})
  STRING(REPLACE "solver/IonicModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/solver/IonicModels"
  COMPONENT ${PACKAGE_NAME}
)
