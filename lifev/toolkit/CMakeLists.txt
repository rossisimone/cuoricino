INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)
INCLUDE(AddSubdirectories)
INCLUDE(TribitsAddExecutable)

#
# A) Define the package
#

TRIBITS_PACKAGE(Toolkit)
SET(${PACKAGE_NAME}_VERSION "3.1")

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

SET(HEADERS "")
SET(SOURCES "")

ADD_SUBDIRECTORIES(lvelectrom)
ADD_SUBDIRECTORIES(lvfibers)

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevtoolkit
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

#TRIBITS_ADD_TEST_DIRECTORIES(testsuite)
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
