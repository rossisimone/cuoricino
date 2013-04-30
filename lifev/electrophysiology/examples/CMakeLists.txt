
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
  IonicRogersMcCulloch
  IonicMinimalModel
  IonicFox
  Restart
  RosenbrockTest  
#  RiceModel
)
