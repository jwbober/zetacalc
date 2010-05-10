                                                                                //-----------------------------------------------------------
                                                                                //
                                                                                // Derivatives of the function
                                                                                //
                                                                                // g(t) = exp(2 pi i alpha t + 2 pi i b t^2)
                                                                                //
void initialize_power_arrays(int n, Complex alpha, Complex b);                  // at 0, 1, and K. Only supports derivatives up
Complex g_derivative_at_1(int n);                                               // to the 21st, and initialize_power_arrays() must
Complex g_derivative_at_0(int n);                                               // be called before the other functions, with n
Complex g_derivative_at_K_without_exponential_factor(int n, int K);             // at least as large as any derivatives to be
                                                                                // computed. For the function ...without_exponential_factor()
                                                                                // the answer needs to be multiplied by
                                                                                //
                                                                                //  exp(2 pi i a K + 2 pi i b K^2)
                                                                                //
                                                                                // to get the actual derivative.
                                                                                // ----------------------------------------------------------



