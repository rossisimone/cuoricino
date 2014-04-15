
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DNegroniLascanoModel
  #0DNegroniLascanoModelJafriRiceWinslow
  benchmarkPurkinje
  FiberInterpolation
  XbNegroniLascano96
  example_0DJafriRiceWinslowModel
  example_0DJafriRiceWinslowModelwStimPro
  #IonicRogersMcCulloch_pacingECG
  #IonicRogersMcCullochSpiral
  IonicMinimalModel_pacingECG
  IonicFox
  IonicFox_pacingECG
)
