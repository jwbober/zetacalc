CXX = clang++

#OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range -g
OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -ffast-math -g -pthread -Winline
H_OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -g -pthread

#OPTIONS = -msse2 -mfpmath=sse -Wall -g
#H_OPTIONS = -msse2 -mfpmath=sse -Wall -g

#OPTIONS = -O3 -Wall -ffast-math -g
LIBS = -lmpfr -lgmp -lgmpxx -pthread -g
#LIBS = -lmpfr -lgmp -lgmpxx
INCLUDEDIR = -Iinclude
#INCLUDEDIR = -I/usr/local/sage/local/include

zetacalc:  build/theta_sums.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/log.o include/rs_sum.h include/log.h build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/w_coefficient.o build/cache.o build/zetacalc.o build/riemann_siegel_sum.o build/stage1.o build/stage2.o build/stage3.o build/log.o build/rs_misc.o
	$(CXX) -O3 -o zetacalc build/theta_sums.o build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/rs_misc.o build/stage1.o build/stage2.o build/stage3.o build/riemann_siegel_sum.o build/log.o build/w_coefficient.o build/cache.o build/zetacalc.o $(LIBS)

blfi:  build/theta_sums.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/log.o include/rs_sum.h include/log.h build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/w_coefficient.o build/cache.o build/zetacalc.o build/riemann_siegel_sum.o build/stage1.o build/stage2.o build/stage3.o build/log.o build/rs_misc.o build/blfi.o gmpfrxx/gmpfrxx.o
	$(CXX) -O3 -o blfi build/theta_sums.o build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/rs_misc.o build/stage1.o build/stage2.o build/stage3.o build/riemann_siegel_sum.o build/log.o build/w_coefficient.o build/cache.o build/blfi.o gmpfrxx/gmpfrxx.o $(LIBS)


test:  build/test.o build/theta_sums.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/log.o include/rs_sum.h include/log.h build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/w_coefficient.o build/cache.o build/riemann_siegel_sum.o build/stage1.o build/stage2.o build/stage3.o build/log.o build/rs_misc.o
	$(CXX) -O3 -o test build/theta_sums.o build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/rs_misc.o build/stage1.o build/stage2.o build/stage3.o build/riemann_siegel_sum.o build/log.o build/w_coefficient.o build/cache.o build/test.o $(LIBS)

build/test.o: test.cc include/theta_sums.h include/log.h
	$(CXX) -c test.cc $(OPTIONS) $(INCLUDEDIR) -o build/test.o

build/blfi.o: blfi.cc include/rs_sum.h
	$(CXX) -c blfi.cc $(OPTIONS) $(INCLUDEDIR) -o build/blfi.o

build/theta_sums.o: theta_sums/theta_sums.cc include/theta_sums.h
	$(CXX) -c theta_sums/theta_sums.cc $(OPTIONS) $(INCLUDEDIR) -o build/theta_sums.o

build/direct_evaluation.o: theta_sums/direct_evaluation.cc include/theta_sums.h
	$(CXX) -c theta_sums/direct_evaluation.cc $(OPTIONS) $(INCLUDEDIR) -o build/direct_evaluation.o

build/exp_sum_euler_maclaurin.o: theta_sums/exp_sum_euler_maclaurin.cc include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/exp_sum_euler_maclaurin.cc $(OPTIONS) $(INCLUDEDIR) -o build/exp_sum_euler_maclaurin.o

build/theta_algorithm.o: theta_sums/theta_algorithm.cc include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/theta_algorithm.cc $(OPTIONS) $(INCLUDEDIR) -o build/theta_algorithm.o

build/G_functions.o: theta_sums/G_functions.cc include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/G_functions.cc $(OPTIONS) $(INCLUDEDIR) -o build/G_functions.o

build/H_functions.o: theta_sums/H_functions.cc include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/H_functions.cc $(H_OPTIONS) $(INCLUDEDIR) -o build/H_functions.o

build/ICn.o: theta_sums/ICn.cc include/theta_sums.h theta_sums/precomputed_tables.h include/log.h
	$(CXX) -c theta_sums/ICn.cc $(OPTIONS) $(INCLUDEDIR) -o build/ICn.o

build/H_and_J_integrals.o: theta_sums/H_and_J_integrals.cc include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/H_and_J_integrals.cc $(OPTIONS) $(INCLUDEDIR) -o build/H_and_J_integrals.o

build/cache.o: theta_sums/cache.cc include/theta_sums.h
	$(CXX) -c theta_sums/cache.cc $(OPTIONS) $(INCLUDEDIR) -o build/cache.o

build/derivative_computations.o: theta_sums/derivative_computations.cc include/theta_sums.h
	$(CXX) -c theta_sums/derivative_computations.cc $(OPTIONS) $(INCLUDEDIR) -o build/derivative_computations.o

build/stage1.o: rs_sum/stage1.cc include/rs_sum.h
	$(CXX) -c rs_sum/stage1.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage1.o

build/stage2.o: rs_sum/stage2.cc include/rs_sum.h
	$(CXX) -c rs_sum/stage2.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage2.o

build/stage3.o: rs_sum/stage3.cc include/rs_sum.h
	$(CXX) -c rs_sum/stage3.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage3.o

build/rs_misc.o: rs_sum/rs_misc.cc include/rs_sum.h
	$(CXX) -c rs_sum/rs_misc.cc $(OPTIONS) $(INCLUDEDIR) -o build/rs_misc.o

build/riemann_siegel_sum.o: rs_sum/riemann_siegel_sum.cc include/rs_sum.h include/theta_sums.h
	$(CXX) -c rs_sum/riemann_siegel_sum.cc $(OPTIONS) $(INCLUDEDIR) -o build/riemann_siegel_sum.o

build/zetacalc.o: include/rs_sum.h zetacalc.cc
	$(CXX) -c zetacalc.cc $(OPTIONS) $(INCLUDEDIR) -o build/zetacalc.o

build/misc.o: theta_sums/misc.cc theta_sums/precomputed_tables.h include/misc.h include/log.h
	$(CXX) -c theta_sums/misc.cc $(OPTIONS) $(INCLUDEDIR) -o build/misc.o

build/w_coefficient.o: theta_sums/w_coefficient.cc theta_sums/w_coefficient.h include/misc.h include/theta_sums.h theta_sums/precomputed_tables.h
	$(CXX) -c theta_sums/w_coefficient.cc $(OPTIONS) $(INCLUDEDIR) -o build/w_coefficient.o

build/log.o: log/log.cc include/log.h
	$(CXX) -c log/log.cc $(OPTIONS) $(INCLUDEDIR) -o build/log.o

clean:
	rm zetacalc
	rm build/*
