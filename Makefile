#OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range
OPTIONS = -O3 -msse2 -mfpmath=sse -Wall -ffast-math

a.out: theta_sums.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o zeta.h log.h direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o
	g++ -O3 theta_sums.o direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o -lmpfr -lgmp -msse2 -mfpmath=sse -lprofiler -lgmpxx

main.o: main.cc theta_sums.h
	g++ -c main.cc $(OPTIONS)

theta_sums.o: theta_sums.cc theta_sums.h
	g++ -c theta_sums.cc $(OPTIONS)

direct_evaluation.o: direct_evaluation.cc theta_sums.h
	g++ -c direct_evaluation.cc $(OPTIONS)

exp_sum_euler_maclaurin.o: exp_sum_euler_maclaurin.cc theta_sums.h precomputed_tables.h
	g++ -c exp_sum_euler_maclaurin.cc $(OPTIONS)

theta_algorithm.o: theta_algorithm.cc theta_sums.h
	g++ -c theta_algorithm.cc $(OPTIONS)

G_functions.o: G_functions.cc theta_sums.h precomputed_tables.h
	g++ -c G_functions.cc $(OPTIONS)

H_functions.o: H_functions.cc theta_sums.h precomputed_tables.h
	g++ -c H_functions.cc $(OPTIONS)

ICn.o: ICn.cc theta_sums.h precomputed_tables.h log.h
	g++ -c ICn.cc $(OPTIONS)

H_and_J_integrals.o: H_and_J_integrals.cc theta_sums.h precomputed_tables.h
	g++ -c H_and_J_integrals.cc $(OPTIONS)

derivative_computations.o: derivative_computations.cc theta_sums.h
	g++ -c derivative_computations.cc $(OPTIONS)

zeta.o: zeta.cc zeta.h
	g++ -c zeta.cc $(OPTIONS)

misc.o: misc.cc precomputed_tables.h w_coefficient.h misc.h
	g++ -c misc.cc $(OPTIONS)

log.o: log.cc log.h
	g++ -c log.cc $(OPTIONS)

clean:
	rm theta_sums.o direct_evaluation.o exp_sum_euler_maclaurin.o theta_algorithm.o G_functions.o H_functions.o ICn.o H_and_J_integrals.o derivative_computations.o main.o misc.o zeta.o log.o 
