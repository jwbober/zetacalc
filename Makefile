#OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range -g
OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -ffast-math -g
H_OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -g
#OPTIONS = -O3 -Wall -ffast-math -g
LIBS = -lmpfr -lgmp -lprofiler -lgmpxx
#LIBS = -lmpfr -lgmp -lgmpxx
INCLUDEDIR = 
#INCLUDEDIR = -I/usr/local/sage/local/include

a.out: theta_sums.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o zeta.h log.h direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o stats.o w_coefficient.o cache.o
	g++ -O3 theta_sums.o direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o stats.o w_coefficient.o cache.o $(LIBS)

test:  theta_sums.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o misc.o zeta.o log.o zeta.h log.h direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o test.o stats.o w_coefficient.o cache.o
	g++ -O3 -o test theta_sums.o direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o test.o misc.o zeta.o log.o stats.o w_coefficient.o cache.o $(LIBS)

test.o: test.cc theta_sums.h log.h
	g++ -c test.cc $(OPTIONS) $(INCLUDEDIR)

main.o: main.cc theta_sums.h
	g++ -c main.cc $(OPTIONS) $(INCLUDEDIR)

stats.o: stats.cc theta_sums.h
	g++ -c stats.cc $(OPTIONS) $(INCLUDEDIR)

theta_sums.o: theta_sums.cc theta_sums.h
	g++ -c theta_sums.cc $(OPTIONS) $(INCLUDEDIR)

direct_evaluation.o: direct_evaluation.cc theta_sums.h
	g++ -c direct_evaluation.cc $(OPTIONS) $(INCLUDEDIR)

exp_sum_euler_maclaurin.o: exp_sum_euler_maclaurin.cc theta_sums.h precomputed_tables.h
	g++ -c exp_sum_euler_maclaurin.cc $(OPTIONS) $(INCLUDEDIR)

theta_algorithm.o: theta_algorithm.cc theta_sums.h
	g++ -c theta_algorithm.cc $(OPTIONS) $(INCLUDEDIR)

G_functions.o: G_functions.cc theta_sums.h precomputed_tables.h
	g++ -c G_functions.cc $(OPTIONS) $(INCLUDEDIR)

H_functions.o: H_functions.cc theta_sums.h precomputed_tables.h
	g++ -c H_functions.cc $(H_OPTIONS) $(INCLUDEDIR)

ICn.o: ICn.cc theta_sums.h precomputed_tables.h log.h
	g++ -c ICn.cc $(OPTIONS) $(INCLUDEDIR)

H_and_J_integrals.o: H_and_J_integrals.cc theta_sums.h precomputed_tables.h
	g++ -c H_and_J_integrals.cc $(OPTIONS) $(INCLUDEDIR)

cache.o: cache.cc theta_sums.h
	g++ -c cache.cc $(OPTIONS) $(INCLUDEDIR)

derivative_computations.o: derivative_computations.cc theta_sums.h
	g++ -c derivative_computations.cc $(OPTIONS) $(INCLUDEDIR)

zeta.o: zeta.cc zeta.h
	g++ -c zeta.cc $(OPTIONS) $(INCLUDEDIR)

misc.o: misc.cc precomputed_tables.h misc.h log.h
	g++ -c misc.cc $(OPTIONS) $(INCLUDEDIR)

w_coefficient.o: w_coefficient.cc w_coefficient.h misc.h
	g++ -c w_coefficient.cc $(OPTIONS) $(INCLUDEDIR)

log.o: log.cc log.h
	g++ -c log.cc $(OPTIONS) $(INCLUDEDIR)

clean:
	rm theta_sums.o direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o test.o test stats.o w_coefficient.o
