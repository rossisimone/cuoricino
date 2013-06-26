
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  basicExample
  basicExampleHO
  basicExampleGHO
  exampleGHO_withFiberRotation
  exampleHO_withFiberRotation
   exampleHO_ellipsoid
  example_ActiveStrain
  example_CalciumDependentActivation
  example_laplacian
)
