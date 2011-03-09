#OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range -g
OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -ffast-math -g -pthread
H_OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -g -pthread

#OPTIONS = -msse2 -mfpmath=sse -Wall -g
#H_OPTIONS = -msse2 -mfpmath=sse -Wall -g

#OPTIONS = -O3 -Wall -ffast-math -g
LIBS = -lmpfr -lgmp -lgmpxx -pthread -g
#LIBS = -lmpfr -lgmp -lgmpxx
INCLUDEDIR = -Iinclude
#INCLUDEDIR = -I/usr/local/sage/local/include

zetacalc:  build/theta_sums.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/log.o include/rs_sum.h include/log.h build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/stats.o build/w_coefficient.o build/cache.o build/compute_zeta.o build/riemann_siegel_sum.o build/stage1.o build/stage2.o build/stage3.o build/log.o build/rs_misc.o
	g++ -O3 -o zetacalc build/theta_sums.o build/direct_evaluation.o build/exp_sum_euler_maclaurin.o build/theta_algorithm.o build/G_functions.o build/H_functions.o build/ICn.o build/H_and_J_integrals.o build/derivative_computations.o build/misc.o build/rs_misc.o build/stage1.o build/stage2.o build/stage3.o build/riemann_siegel_sum.o build/log.o build/stats.o build/w_coefficient.o build/cache.o build/compute_zeta.o $(LIBS)

build/test.o: test.cc include/theta_sums.h include/log.h
	g++ -c test.cc $(OPTIONS) $(INCLUDEDIR) -o build/test.o

build/main.o: main.cc include/theta_sums.h
	g++ -c main.cc $(OPTIONS) $(INCLUDEDIR) -o build/main.o

build/stats.o: stats.cc include/theta_sums.h include/rs_sum.h
	g++ -c stats.cc $(OPTIONS) $(INCLUDEDIR) -o build/stats.o

build/theta_sums.o: theta_sums.cc include/theta_sums.h
	g++ -c theta_sums.cc $(OPTIONS) $(INCLUDEDIR) -o build/theta_sums.o

build/direct_evaluation.o: direct_evaluation.cc include/theta_sums.h
	g++ -c direct_evaluation.cc $(OPTIONS) $(INCLUDEDIR) -o build/direct_evaluation.o

build/exp_sum_euler_maclaurin.o: exp_sum_euler_maclaurin.cc include/theta_sums.h precomputed_tables.h
	g++ -c exp_sum_euler_maclaurin.cc $(OPTIONS) $(INCLUDEDIR) -o build/exp_sum_euler_maclaurin.o

build/theta_algorithm.o: theta_algorithm.cc include/theta_sums.h precomputed_tables.h
	g++ -c theta_algorithm.cc $(OPTIONS) $(INCLUDEDIR) -o build/theta_algorithm.o

build/G_functions.o: G_functions.cc include/theta_sums.h precomputed_tables.h
	g++ -c G_functions.cc $(OPTIONS) $(INCLUDEDIR) -o build/G_functions.o

build/H_functions.o: H_functions.cc include/theta_sums.h precomputed_tables.h
	g++ -c H_functions.cc $(H_OPTIONS) $(INCLUDEDIR) -o build/H_functions.o

build/ICn.o: ICn.cc include/theta_sums.h precomputed_tables.h include/log.h
	g++ -c ICn.cc $(OPTIONS) $(INCLUDEDIR) -o build/ICn.o

build/H_and_J_integrals.o: H_and_J_integrals.cc include/theta_sums.h precomputed_tables.h
	g++ -c H_and_J_integrals.cc $(OPTIONS) $(INCLUDEDIR) -o build/H_and_J_integrals.o

build/cache.o: cache.cc include/theta_sums.h
	g++ -c cache.cc $(OPTIONS) $(INCLUDEDIR) -o build/cache.o

build/derivative_computations.o: derivative_computations.cc include/theta_sums.h
	g++ -c derivative_computations.cc $(OPTIONS) $(INCLUDEDIR) -o build/derivative_computations.o

build/stage1.o: rs_sum/stage1.cc include/rs_sum.h
	g++ -c rs_sum/stage1.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage1.o

build/stage2.o: rs_sum/stage2.cc include/rs_sum.h
	g++ -c rs_sum/stage2.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage2.o

build/stage3.o: rs_sum/stage3.cc include/rs_sum.h
	g++ -c rs_sum/stage3.cc $(OPTIONS) $(INCLUDEDIR) -o build/stage3.o

build/rs_misc.o: rs_sum/rs_misc.cc include/rs_sum.h
	g++ -c rs_sum/rs_misc.cc $(OPTIONS) $(INCLUDEDIR) -o build/rs_misc.o

build/riemann_siegel_sum.o: rs_sum/riemann_siegel_sum.cc include/rs_sum.h include/theta_sums.h
	g++ -c rs_sum/riemann_siegel_sum.cc $(OPTIONS) $(INCLUDEDIR) -o build/riemann_siegel_sum.o

build/compute_zeta.o: include/rs_sum.h compute_zeta.cc
	g++ -c compute_zeta.cc $(OPTIONS) $(INCLUDEDIR) -o build/compute_zeta.o

build/misc.o: misc.cc precomputed_tables.h include/misc.h include/log.h
	g++ -c misc.cc $(OPTIONS) $(INCLUDEDIR) -o build/misc.o

build/w_coefficient.o: w_coefficient.cc w_coefficient.h include/misc.h include/theta_sums.h precomputed_tables.h
	g++ -c w_coefficient.cc $(OPTIONS) $(INCLUDEDIR) -o build/w_coefficient.o

build/log.o: log/log.cc include/log.h
	g++ -c log/log.cc $(OPTIONS) $(INCLUDEDIR) -o build/log.o

clean:
	rm zetacalc
	rm build/*
