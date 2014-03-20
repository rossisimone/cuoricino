INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
	test_benchmark
	benchmarkPurkinje	 
    test_bidomain
	test_restart
	test_0DAlievPanfilovModel
	test_0DLuoRudyIModel
	test_0DJafriRiceWinslowModel
	test_0DMinimalModel
	test_0DFitzHughNagumoModel
	test_0DFoxModel
	test_0DTenTusscher06Model
	test_0DHodgkinHuxley
	test_0DNoblePurkinje
)
