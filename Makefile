CC = g++
CXX = g++
CXXFLAGS = -march=native -Wall -ffast-math -pthread -Winline -O3 -g -std=c++11 -Iinclude
#CXXFLAGS = -march=native -Wall -pthread -Winline -O3 -g -std=c++11 -Iinclude
H_CXXFLAGS = -march=native -Wall -pthread  -O3 -g -std=c++11 -Iinclude

#OPTIONS = -msse2 -mfpmath=sse -Wall -ffast-math -pthread -Winline  -O3 -fno-omit-frame-pointer -g -std=c++11



LDFLAGS := -lmpfr -lgmp -lgmpxx -pthread -g
ifdef PROFILE_BUILD
    LDFLAGS := $(LDFLAGS) -Wl,--no-as-needed -lprofiler
endif

HOSTNAME = $(shell hostname)
ifeq ($(HOSTNAME),lmfdb5.maths.bris.ac.uk)
    LDFLAGS := -L/data/local/lib $(LDFLAGS)
endif

THETA_SUM_OBJECTS = \
		    theta_sums/cache.o \
		    theta_sums/derivative_computations.o \
		    theta_sums/direct_evaluation.o \
		    theta_sums/exp_sum_euler_maclaurin.o \
		    theta_sums/G_functions.o \
		    theta_sums/H_and_J_integrals.o \
		    theta_sums/H_functions.o \
		    theta_sums/ICn.o \
		    theta_sums/misc.o \
		    theta_sums/theta_algorithm.o \
		    theta_sums/theta_sums.o \
		    theta_sums/w_coefficient.o

MAIN_SUM_OBJECTS = \
		   main_sum/main_sum.o \
		   main_sum/ms_misc.o \
		   main_sum/stage1.o \
		   main_sum/stage2.o \
		   main_sum/stage3.o

OTHER_OBJECTS = \
		log/log.o \
		misc/pow.o

OBJECTS = $(MAIN_SUM_OBJECTS) \
	  $(THETA_SUM_OBJECTS) \
	  $(OTHER_OBJECTS) \
	  test.o \
	  zetacalc.o \
	  tests/mainsum_tests.o

EXECUTABLES = zetacalc \
	      test

TESTS =       tests/mainsum_tests \
	      tests/Htest \
	      tests/thetatime \
	      tests/thetatest \
	      tests/write_testfile \
	      tests/write_testfile2 \
	      tests/run_testfile \

zetacalc: $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) zetacalc.o
	$(CXX) -o zetacalc zetacalc.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

$(THETA_SUM_OBJECTS): include/theta_sums.h include/log.h include/misc.h theta_sums/precomputed_tables.h
$(MAIN_SUM_OBJECTS): include/theta_sums.h include/log.h include/main_sum.h include/misc.h theta_sums/precomputed_tables.h
log/log.o: include/log.h
$(EXECUTABLES): include/theta_sums.h include/log.h include/main_sum.h include/misc.h

theta_sums/H_functions.o: theta_sums/H_functions.cc
	$(CXX) -c theta_sums/H_functions.cc $(H_CXXFLAGS) -o theta_sums/H_functions.o

test: $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) test.o
	$(CXX) -o test test.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#tests: tests/mainsum_tests tests/Htest tests/thetatime tests/thetatest tests/write_testfile
tests: $(TESTS)

tests/%: tests/%.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
	$(CXX) -o $@ $< $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#tests/mainsum_tests:  tests/mainsum_tests.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
#	$(CXX) -o tests/mainsum_tests tests/mainsum_tests.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#tests/thetatime: tests/thetatime.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
#	$(CXX) -o tests/thetatime tests/thetatime.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#tests/thetatest: tests/thetatest.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
#	$(CXX) -o tests/thetatest tests/thetatest.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#tests/write_testfile: tests/write_testfile.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
#	$(CXX) -o tests/write_testfile tests/write_testfile.o $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

#build/mainsum_tests.o: tests/mainsum_tests.cc include/main_sum.h include/theta_sums.h
#	$(CXX) -c tests/mainsum_tests.cc $(OPTIONS) $(INCLUDEDIR) -o build/mainsum_tests.o

.PHONY: clean
clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
	-rm $(TESTS)


build/blfi.o: blfi.cc include/main_sum.h
	$(CXX) -c blfi.cc $(OPTIONS) $(INCLUDEDIR) -o build/blfi.o

blfi:  build/theta_sums.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/log.o include/main_sum.h include/log.h build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/w_coefficient.o build/cache.o build/zetacalc.o build/main_sum.o build/stage1.o build/stage2.o build/stage3.o build/log.o build/ms_misc.o build/blfi.o gmpfrxx/gmpfrxx.o
	$(CXX) -O3 -o blfi build/theta_sums.o build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/ms_misc.o build/stage1.o build/stage2.o build/stage3.o build/main_sum.o build/log.o build/w_coefficient.o build/cache.o build/blfi.o gmpfrxx/gmpfrxx.o $(LIBS)
