/*
 * Functions for computing
 *   SUM(v, K, t) = sum_{n=v}^{v + K - 1} n^{-(1/2 + it)}
 *
 * These "stage 3" functions use Ghaith Hiary's "theta algorithm" and
 * are very fast when v is large. They should only be used when
 * v >> t^(1/3) (our current choice is v > 1100 t^(1/3), but the specific
 * choice of 1100 depends on the quality of the theta algorithm
 * implementation and may also be machine dependent, so it is subject
 * to change.)
 *
 */

#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "rs_sum.h"
#include "log.h"

using namespace std;

void compute_taylor_coefficients(mpfr_t t, Complex Z[30]);

const int MAX_THREADS = 30;
extern string NUM_THREADS_FILE;
extern bool use_num_threads_file;
extern int default_number_of_threads;

unsigned int stage_3_block_size(Double v, Double t) {
    // With this choice, the block size will hit N when
    // v = N/.9 * t^(1/3)
    //
    // So if we want to start stage 3 with a block size of 100, then we should
    // set stage_2_bound to 112 * t^(1/3)
    //unsigned int block_size = (unsigned int)(  min(.9 * v * pow(t, -.3333333333333333333333), pow(v, .75)/sqrt(500)) );
    unsigned int block_size = (unsigned int)(  min(.9 * v * pow(t, -.3333333333333333333333), v/500.0 * 500.0) );
    //unsigned int block_size = (unsigned int)(2 * v * pow(t, -.3333333333333333333333));
    return block_size;
}

Complex zeta_block_stage3_basic(mpz_t v, unsigned int *K, mpfr_t t, Complex ZZ[30], Double epsilon, int Kmin) {
    if(*K == 0) {
        return 0.0;
    }
    if(*K == 1) {
        if(0) {
            cout << v << " " << v << endl;
        }
        return exp_itlogn(v)/sqrt(mpz_get_d(v));
    }
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    unsigned int block_size = min(stage_3_block_size(vv, tt), *K);
    if(block_size < 50 && *K >= 50) {
        cout << "Error: in stage 3, computed a block size of " << block_size << ", which is smaller than 50. Shouldn't reach stage two until we can use a larger block size." << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;
    if(0) {
        mpz_t n;
        mpz_init(n);
        mpz_add_ui(n, v, block_size - 1);
        cout << v << " " << n << endl;
        mpz_clear(n);
    }




    Double w = (block_size-1)/vv;
    Double w_power = 1;

    Complex Z[30];

    for(int l = 0; l < 19; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    int j = 18;

    // Compute Z[l]
 
    mpfr_t a, b, x;

    int precision = mpfr_get_prec(t);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(x, precision);

    mpfr_const_pi(x, GMP_RNDN);             // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);         // x = 2 pi
    mpfr_mul_z(x, x, v, GMP_RNDN);            // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN);            // a = t / (2 pi v)

//    mpfr_mul_si(x, x, -2, GMP_RNDN);        // x = -4 pi v
//    mpfr_mul(x, x, v, GMP_RNDN);            // x = -4 pi v^2
//    mpfr_div(b, t, x, GMP_RNDN);            // b = -t/ (4 pi v^2)

//    mpfr_mul_si(b, v, -2, GMP_RNDN);
    mpfr_div_z(b, a, v, GMP_RNDN);
    mpfr_div_si(b, b, -2, GMP_RNDN);

//    cout << mpfr_get_d(b, GMP_RNDN) << endl;

    Complex S = compute_exponential_sums(a, b, j, block_size-1, Z, epsilon, Kmin, 0);

    // we don't need the previous values of a and b anymore, so
    // we can erase them.

    Complex z = exp_itlogn(v);
    z = z / sqrt(mpz_get_d(v));
    S = S * z;

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(x);


//    cout << *K << endl;

    return S;
}


class stage3_data_t {
public:
    mpz_t v;
    mpfr_t t;
    unsigned int blocksize;
    Double delta;
    int M;
    Complex * S;
    
    queue<int> * thread_queue;
    pthread_mutex_t * queue_mutex;
    pthread_cond_t * queue_nonempty_signaler;
    int id;
    int Kmin;
    Complex * Z;


    stage3_data_t() {
        mpz_init(v);
        mpfr_init2(t, 150);
    }

    void set(   mpz_t _v,
                unsigned int _blocksize,
                mpfr_t _t,
                Double _delta,
                int _M, Complex * _S,
                queue<int> * _thread_queue,
                pthread_mutex_t * _queue_mutex,
                pthread_cond_t * _queue_nonempty_signaler,
                int _id, 
                int _Kmin,
                Complex * _Z) {

        mpz_set(v, _v);
        mpfr_set(t, _t, GMP_RNDN);

        blocksize = _blocksize;
        delta = _delta;
        M = _M;
        S = _S;

        thread_queue = _thread_queue;
        queue_mutex = _queue_mutex;
        queue_nonempty_signaler = _queue_nonempty_signaler;

        id = _id;
        Kmin = _Kmin;
        Z = _Z;
    }

    ~stage3_data_t() {
        mpz_clear(v);
        mpfr_clear(t);
    }

};

void * zeta_block_stage3(void * thread_data) {
    stage3_data_t * data = (stage3_data_t * )(thread_data);

    zeta_block_stage3(data->v, data->blocksize, data->t, data->Z, data->delta, data->M, data->S, data->Kmin);

    // when we are done we need to put our thread number in the queue, and if the queue
    // was empty, we need to signal that it is no longer empty

    pthread_mutex_lock(data->queue_mutex);

    bool queue_was_empty = data->thread_queue->empty();

    data->thread_queue->push(data->id);

    if(queue_was_empty) {
        pthread_cond_signal(data->queue_nonempty_signaler);
    }

    pthread_mutex_unlock(data->queue_mutex);

    pthread_exit(NULL);
}


Complex zeta_block_stage3(mpz_t n, unsigned int N, mpfr_t t, Complex Z[30], Double delta, int M, Complex * S, int Kmin) {
    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
    }

    if(N == 0) {
        return 0;
    }


    mpz_t v;
    mpz_init(v);

    mpz_set(v, n);
    unsigned int K = N;
    while(N > 0) {
        Complex current_term = zeta_block_stage3_basic(v, &K, t, Z, exp(-20), Kmin);
        Complex multiplier = exp(I * delta * log(mpz_get_d(v)));
        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }
        N = N - K;
        mpz_add_ui(v, v, K);
        K = N;
    }

    mpz_clear(v);

    return S[0];
}


Complex zeta_sum_stage3(mpz_t n, mpz_t N, mpfr_t t, Double delta, int M, Complex * S, int Kmin, int verbose) {
    //
    // Compute and return the sum
    //
    // \sum_{k=n}^{n + N - 1} n^{.5 + it}
    //
    // We do this by repeatedly calling zeta_block_stage3 with a block size of 1000000
    // and add up all the terms. We could just call it directly, but we anticipate
    // shortly changing this routine to use multiple threads.
    //

    time_t start_wall_time = time(NULL);
    clock_t last_cpu_time = clock();
    double total_cpu_time = 0;

    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
    }


    Complex S2[M];
    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
        S2[l] = 0.0;
    }


    const unsigned int block_size = 100000000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, N, block_size);

    Complex Z[30];
    compute_taylor_coefficients(t, Z);

    int num_threads = default_number_of_threads;
    if(use_num_threads_file) {
        ifstream num_threads_file;
        num_threads_file.open(NUM_THREADS_FILE.c_str());
        if(num_threads_file) {
            num_threads_file >> num_threads;
            num_threads_file.close();
        }
    }

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set(v, n);
    mpz_sub_ui(v, v, block_size);
    if(num_threads == 1 || mpz_cmp_si(number_of_blocks, num_threads) < 0 ) {
        for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
            mpz_add_ui(v, v, block_size);
            zeta_block_stage3(v, block_size, t, Z, delta, M, S2, Kmin);
            for(int l = 0; l < M; l++) {
                S[l] += S2[l];
            }
            //S = S + zeta_block_mpfr(v, block_size, t);
            if(mpz_divisible_ui_p(k, 20u)) {
                time_t current_wall_time = time(NULL);
                clock_t current_cpu_time = clock();
                time_t elapsed_wall_time = current_wall_time - start_wall_time;
                double elapsed_cpu_time = ((double)current_cpu_time - (double)last_cpu_time)/CLOCKS_PER_SEC;
                total_cpu_time += elapsed_cpu_time;
                if(verbose) {
                    cout << "In stage3, completed " << k << " large blocks out of " << number_of_blocks << "." << endl;
                    cout << "        In stage3 thus far: " << elapsed_wall_time << " real seconds; " << total_cpu_time << " cpu seconds; " << elapsed_cpu_time << "cpu seconds this chunk. " << endl;
                }
                last_cpu_time = current_cpu_time;
                int current_blocksize = stage_3_block_size(mpz_get_d(v), mpfr_get_d(t, GMP_RNDN));
                if(verbose)
                    cout << "        Current blocksize ~= " << current_blocksize << endl;
            }
        }
        mpz_add_ui(v, v, block_size);

        zeta_block_stage3(v, remainder, t, Z, delta, M, S2, Kmin);
        for(int l = 0; l < M; l++) {
            S[l] += S2[l];
        }
        //S = S + zeta_block_mpfr(v, remainder, t);
    }
    else {
        pthread_t threads[MAX_THREADS];
        bool unjoined[MAX_THREADS];
        for(int l = 0; l < MAX_THREADS; l++) {
            unjoined[l] = false;
        }

        queue<int> thread_queue;
        pthread_mutex_t queue_mutex;
        pthread_cond_t queue_nonempty_signaler;

        pthread_mutex_init(&queue_mutex, NULL);
        pthread_cond_init(&queue_nonempty_signaler, NULL);
        
        stage3_data_t thread_data[MAX_THREADS];
        Complex S2[MAX_THREADS][M];

        // start by spawning a bunch of threads
        //
        // we lock the queue in case any of these threads would finish before
        // we start waiting for them to finish.

        pthread_mutex_lock(&queue_mutex);
        mpz_set_ui(k, 0u);
        for(int n = 0; n < num_threads; n++, mpz_add_ui(k, k, 1u)) {
            mpz_add_ui(v, v, block_size);
            thread_data[n].set(v, block_size, t, delta, M, S2[n], &thread_queue, &queue_mutex, &queue_nonempty_signaler, n, Kmin, Z);
            pthread_create(&threads[n], NULL, zeta_block_stage3, (void *)(&thread_data[n]));
            unjoined[n] = true;
        }

        // at this point the queue should be empty, and we want to wait for
        // it to be nonempty, so we wait
        
        pthread_cond_wait(&queue_nonempty_signaler, &queue_mutex);



        for(; mpz_cmp(k, number_of_blocks) < 0;) {
            //if(next_thread == num_threads) {
            //    for(int n = 0; n < num_threads; n++) {
            //        void * status;
            //        pthread_join(threads[n], &status);
            //        for(int l = 0; l < M; l++) {
            //            S[l] += S2[n][l];
            //        }
            //    }
            //    next_thread = 0;
            //}

            // at this point it should be guaranteed that the thread queue is
            // nonempty, so there is a thread that has finished and is ready to do
            // some more work. also, we have the queue locked.

            int next_thread = thread_queue.front();
            thread_queue.pop();


            void * status;
            if(unjoined[next_thread]) {
                pthread_join(threads[next_thread], &status);  // join this thread so that it can be destroyed.
                unjoined[next_thread] = false;
            }


            // record the data from this thread
            for(int l = 0; l < M; l++) {
                S[l] += S2[next_thread][l];
            }
        
            // now we create a new thread

            if(next_thread < num_threads) {
                mpz_add_ui(k, k, 1u);
                mpz_add_ui(v, v, block_size);
                thread_data[next_thread].set(v, block_size, t, delta, M, S2[next_thread], &thread_queue, &queue_mutex, &queue_nonempty_signaler, next_thread, Kmin, Z);
                pthread_create(&threads[next_thread], NULL, zeta_block_stage3, (void *)(&thread_data[next_thread]));
                unjoined[next_thread] = true;
            }

            // now we have created a new thread to do some work.
            //
            // if the thread queue is empty, we wait for a signal
            // on the condition variable queue_nonempty_signaler. otherwise, we just
            // go on to the next iteration of the loop.

            if(thread_queue.empty()) {
                pthread_cond_wait(&queue_nonempty_signaler, &queue_mutex);
            }

            //S = S + zeta_block_mpfr(v, block_size, t);
            if(mpz_divisible_ui_p(k, 20u)) {
                time_t current_wall_time = time(NULL);
                clock_t current_cpu_time = clock();
                time_t elapsed_wall_time = current_wall_time - start_wall_time;
                double elapsed_cpu_time = ((double)current_cpu_time - (double)last_cpu_time)/CLOCKS_PER_SEC;
                total_cpu_time += elapsed_cpu_time;
                if(verbose) {
                    cout << "In stage3, completed " << k << " large blocks out of " << number_of_blocks << "." << endl;
                    cout << "        In stage3 thus far: " << elapsed_wall_time << " real seconds; " << total_cpu_time << " cpu seconds; " << elapsed_cpu_time << "cpu seconds this chunk. " << endl;
                }
                last_cpu_time = current_cpu_time;
                int current_blocksize = stage_3_block_size(mpz_get_d(v), mpfr_get_d(t, GMP_RNDN));
                if(verbose) {
                    cout << "        Current blocksize ~= " << current_blocksize << endl;
                    cout << "        Sum so far: " << S[0] << endl;
                }
                if(use_num_threads_file) {
                    ifstream num_threads_file;
                    int new_num_threads = num_threads;
                    num_threads_file.open(NUM_THREADS_FILE.c_str());
                    if(!num_threads_file) {
                        new_num_threads = 2;
                    }
                    else {
                        num_threads_file >> new_num_threads;
                        num_threads_file.close();
                    }
                    if(new_num_threads > num_threads) {
                        if(new_num_threads > MAX_THREADS)
                            new_num_threads = MAX_THREADS;
                        if(verbose)
                            cout << "Increasing the number of threads to " << new_num_threads << endl;
                        // We actually need to spawn new threads now. First we check to see how many blocks
                        // are remaining. If it is less than the number of new threads we should
                        // spawn, we actually do nothing. (Assuming that MAX_THREADS isn't very large
                        // this can never be too wasteful.)
                        
                        mpz_sub(number_of_blocks, number_of_blocks, k);
                        if(mpz_cmp_si(number_of_blocks, num_threads - new_num_threads) >= 0) {
                            mpz_add(number_of_blocks, number_of_blocks, k);
                            mpz_add_ui(k, k, new_num_threads - num_threads);
                            for(int k = num_threads; k < new_num_threads; k++) {
                                mpz_add_ui(v, v, block_size);
                                thread_data[k].set(v, block_size, t, delta, M, S2[k], &thread_queue, &queue_mutex, &queue_nonempty_signaler, k, Kmin, Z);
                                pthread_create(&threads[k], NULL, zeta_block_stage3, (void *)(&thread_data[k]));
                                unjoined[k] = true;
                            }
                        }
                        else {
                            mpz_add(number_of_blocks, number_of_blocks, k);
                            // Do nothing. This is the case where the number of blocks
                            // remaining is very small.
                        }

                        num_threads = new_num_threads;   
                    }
                    else if(new_num_threads < num_threads) {
                        if(new_num_threads < 1)
                            new_num_threads = 1;
                        if(verbose)
                            cout << "Decreasing the number of threads to " << new_num_threads << endl;
                        pthread_mutex_unlock(&queue_mutex);
                        for(int k = new_num_threads; k < num_threads;k++) {
                            void * status;
                            if(verbose) {
                                cout << "Attempting to join thread " << k << endl;
                                cout.flush();
                            }
                            if(unjoined[k]) {
                                pthread_join(threads[k], &status);
                                unjoined[k] = false;
                            }
                        }
                        if(verbose)
                            cout << "Reacquiring lock" << endl;
                        pthread_mutex_lock(&queue_mutex);
                        if(verbose) {
                            cout << "Lock acquired." << endl;
                            cout.flush();
                        }
                        num_threads = new_num_threads;
                    }
                }
            }
        }

        // at this point there are still some threads working which need to finish.
        // we are guaranteed that the thread queue is nonempty, so
        // will don't have to worry about that, but we do need to wait for all
        // threads to finish.

        // there is a good chance that below we will join threads that
        // have already been joined. i don't know what happens in this case.
        // hopefully pthread_join just returns an error and continues.

        // we have the queue locked, and there may be threads waiting for it,
        // so we have to unlock it here.

        pthread_mutex_unlock(&queue_mutex);

        for(int n = 0; n < num_threads; n++) {
            void * status;
            if(unjoined[n]) {
                pthread_join(threads[n], &status);
                unjoined[n] = false;
            }
                //for(int l = 0; l < M; l++) {
                //    S[l] += S2[n][l];
                //}
        }

        // now all the threads are done and we need to record their data

        // at this point there should be no other threads running, so
        // we don't have to worry about any locks or conditions.

        while(!thread_queue.empty()) {
            int next_thread = thread_queue.front();
            thread_queue.pop();

            // record the data from this thread
            for(int l = 0; l < M; l++) {
                S[l] += S2[next_thread][l];
            }
        }


        // finally, we still have one more block to compute.
        mpz_add_ui(v, v, block_size);

        zeta_block_stage3(v, remainder, t, Z, delta, M, S2[0], Kmin);
        for(int l = 0; l < M; l++) {
            S[l] += S2[0][l];
        }
        
        pthread_mutex_destroy(&queue_mutex);
        pthread_cond_destroy(&queue_nonempty_signaler);


    }



    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    time_t end_wall_time = time(NULL);
    time_t elapsed_wall_time = end_wall_time - start_wall_time;
    if(verbose)
        cout << "Spent " << elapsed_wall_time << " seconds in stage 3." << endl;


    return S[0];

}


void compute_taylor_coefficients(mpfr_t t, Complex Z[30]) {
    Double tt1 = mpfr_get_d(t, GMP_RNDN);
    Double tt2 = tt1 * tt1;
    Double tt3 = tt2 * tt1;
    Double tt4 = tt3 * tt1;
    Double tt5 = tt4 * tt1;
    Double tt6 = tt5 * tt1;
    Double tt7 = tt6 * tt1;

    Z[0] = 1.0;
    Z[1] = -.5;
    Z[2] =  .375;
    Z[3] = I * tt1 / 3.0;
    Z[4] = -I * tt1 * 5.0/12.0;
    Z[5] = I * tt1 * 9.0/20.0;
    Z[6] = - tt2 / 18.0;
    Z[7] =  tt2 / 9.0;
    Z[8] = - tt2 * 77.0/480.0;
    Z[9] = -I * tt3 / 162.0;
    Z[10] = I  * tt3 * 11.0/648.0;
    Z[11] = -I * tt3 * 133.0/4320.0;
    Z[12] = tt4 / 1944.0;
    Z[13] = 0.0;
    Z[14] = 0.0;
    Z[15] = tt5 / 29160.0;
    Z[16] = 0.0;
    Z[17] = 0.0;
    Z[18] = tt6 / 524880.0;
    Z[19] = 0;
    Z[20] = 0;
    Z[21] = tt7 / 11022480.0;
    Z[22] = 0;

}
