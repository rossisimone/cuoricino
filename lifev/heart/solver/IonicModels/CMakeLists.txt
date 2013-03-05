SET( IonicModels_HEADERS
  IonicModels/IonicAlievPanfilov.hpp
  IonicModels/IonicFitzHughNagumo.hpp
  IonicModels/IonicMinimalModel.hpp
  IonicModels/HeartIonicModel.hpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${IonicModels_HEADERS})
  STRING(REPLACE "IonicModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/heart/solver/IonicModels"
  COMPONENT ${PACKAGE_NAME}
)
