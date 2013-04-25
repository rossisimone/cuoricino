
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
  FiberInterpolation
#  ActiveStrainODE
#  AlievPanfilov
  basic_example
  ETA_MonodomainSolver
#  IntegratedHeart
  oneway_EM
  XbNegroniLascano96
  IonicAlievPanfilov
  IonicFitzHughNagumo
  IonicMinimalModel
  Restart
  RosenbrockTest  
#  RiceModel
)
