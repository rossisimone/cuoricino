
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  0DAlievPanfilovModel
  0DMitchellSchaefferModel
  0DJafriRiceWinslowModel
  0DMinimalModel
  0DNegroniLascanoModel
  0DFitzHughNagumoModel
  FiberInterpolation
#  ActiveStrainODE
#  AlievPanfilov
  basic_example
  ETA_MonodomainSolver
  integrated_heart
  oneway_EM
  XbNegroniLascano96
  IonicAlievPanfilov
  IonicFitzHughNagumo
  IonicMinimalModel
  Restart
#  RiceModel
)
