CC = g++
CXX = g++
CXXFLAGS = -march=native -Wall -pthread -Winline -O3 -ffast-math -g -std=c++11 -Iinclude
H_CXXFLAGS = -march=native -Wall -pthread  -O3 -g -std=c++11 -Iinclude

LDFLAGS := -L$(HOME)/lib -std=c++11
LDLIBS := -lmpfr -lgmp -lgmpxx -pthread -g
ifdef PROFILE_BUILD
    LDFLAGS := $(LDFLAGS) -Wl,--no-as-needed
    LDLIBS := $(LDLIBS) -lprofiler
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

EXECUTABLES = zetacalc \
	      test \
	      blfi

TESTS =       tests/mainsum_tests \
	      tests/Htest \
	      tests/thetatime \
	      tests/thetatime2 \
	      tests/thetatest \
	      tests/write_testfile \
	      tests/write_testfile2 \
	      tests/run_testfile \

zetacalc: $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) zetacalc.o

blfi: $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) blfi.o gmpfrxx/gmpfrxx.o

theta_sums.a: $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
	ar rs theta_sums.a $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)

$(THETA_SUM_OBJECTS): include/theta_sums.h include/log.h include/misc.h theta_sums/precomputed_tables.h
$(MAIN_SUM_OBJECTS): include/theta_sums.h include/log.h include/main_sum.h include/misc.h theta_sums/precomputed_tables.h
log/log.o: include/log.h

theta_sums/H_functions.o: theta_sums/H_functions.cc
	$(CXX) -c theta_sums/H_functions.cc $(H_CXXFLAGS) -o theta_sums/H_functions.o

test: $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) test.o
	$(CXX) -o test test.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDFLAGS)

tests: $(TESTS)

tests/%: tests/%.o $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $< $(MAIN_SUM_OBJECTS) $(THETA_SUM_OBJECTS) $(OTHER_OBJECTS) $(LDLIBS)

.PHONY: clean
clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
	-rm $(TESTS)
