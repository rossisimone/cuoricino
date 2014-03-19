
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DNegroniLascanoModel
  0DNegroniLascanoModelJafriRiceWinslow
  ActivationWith0DTenTusscher06Model
  FiberInterpolation
  XbNegroniLascano96
  example_0DJafriRiceWinslowModelwStimPro
  IonicRogersMcCulloch_pacingECG
  IonicRogersMcCullochSpiral
  IonicMinimalModel
  IonicMinimalModel_pacingECG
  IonicFox
  IonicFox_pacingECG
  Restart
  IonicHodgkinHuxley
)
