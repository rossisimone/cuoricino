
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DAlievPanfilovModel
  0DLuoRudyIModel
  0DMitchellSchaefferModel
  0DJafriRiceWinslowModel
  0DJafriRiceWinslowModelwStimPro
  0DTenTusscher06Model
  0DMinimalModel
  0DNegroniLascanoModel
  0DNegroniLascanoModelJafriRiceWinslow
  0DFitzHughNagumoModel
  0DFoxModel
  ActivationWith0DTenTusscher06Model
  FiberInterpolation
  ETA_MonodomainSolver
  XbNegroniLascano96
  IonicRogersMcCulloch_pacingECG
  IonicRogersMcCullochSpiral
  IonicMinimalModel
  IonicMinimalModel_pacingECG
  IonicFox
  IonicFox_pacingECG
  Restart
)
