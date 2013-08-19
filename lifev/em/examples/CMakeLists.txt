
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  basicExample
  basicExampleHO
  basicExampleGHO
  exampleGHO_withFiberRotation
  exampleHO_withFiberRotation
  exampleStatic_GAHO
   exampleHO_ellipsoid
  example_ActiveStrain
  example_CalciumDependentActivation
  example_laplacian
  example_sheet
  linearizedGHO
  linearizedGHOrotation
  linearizedGHOventricle
)
