
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  benchmarkPurkinje
  FiberInterpolation
  IonicRogersMcCulloch_pacingECG
  IonicMinimalModel_pacingECG
  IonicFox_pacingECG
)
