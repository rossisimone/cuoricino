
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DAlievPanfilovModel
  0DMitchellSchaefferModel
  0DJafriRiceWinslowModel
  0DJafriRiceWinslowModelwStimPro
  0DTenTusscherModel
  0DTenTusscherModelwStimPro
  0DMinimalModel
  0DNegroniLascanoModel
  0DNegroniLascanoModelTenTusscher
  0DNegroniLascanoModelJafriRiceWinslow
  0DFitzHughNagumoModel
  0DFoxModel
#  CreateFiberEllipsoid
  FiberInterpolation
#  ActiveStrainODE
#  AlievPanfilov
  basic_example
  ETA_MonodomainSolver
  XbNegroniLascano96
  IonicAlievPanfilov
  IonicFitzHughNagumo
  IonicRogersMcCulloch_pacingECG
  IonicRogersMcCullochSpiral
  IonicMinimalModel
  IonicMinimalModel_pacingECG
  IonicFox
  IonicFox_pacingECG
  IonicTenTusscherModel
  Restart
  RosenbrockTest
#  RiceModel
)
