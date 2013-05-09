
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DAlievPanfilovModel
  0DMitchellSchaefferModel
  0DJafriRiceWinslowModel
  0DTenTusscherModel
  0DFoxModel
  0DMinimalModel
  0DNegroniLascanoModel
  0DFitzHughNagumoModel
#  CreateFiberEllipsoid
  FiberInterpolation
#  ActiveStrainODE
#  AlievPanfilov
  basic_example
  ETA_MonodomainSolver
#  oneway_EM
  XbNegroniLascano96
  IonicAlievPanfilov
  IonicFitzHughNagumo
  IonicFitzHughNagumo_test
  IonicRogersMcCulloch
  IonicMinimalModel
  IonicMinimalModel_pacingECG
  IonicFox
  IonicFox_pacingECG
  Restart
  RosenbrockTest  
#  RiceModel
)
