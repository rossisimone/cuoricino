SET(util_HEADERS
  util/ElectrophysiologyUtility.hpp
  util/CardiacStimulus.hpp  
  util/CardiacStimulusPMJ.hpp
  util/CardiacStimulusSingleSource.hpp
  util/CardiacStimulusPacingProtocol.hpp
CACHE INTERNAL "")

SET(util_SOURCES
  util/CardiacStimulus.cpp  
  util/CardiacStimulusPMJ.cpp
  util/CardiacStimulusSingleSource.cpp
  util/CardiacStimulusPacingProtocol.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${util_HEADERS})
  STRING(REPLACE "util/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/util"
  COMPONENT ${PACKAGE_NAME}
)
