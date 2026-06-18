#include <cmath>
#include <stdexcept>

extern "C" double besselK_c(double nu, double x) {
    // 1. Intercept huge values. K_nu(x) asymptotically approaches 0.
    // The C++ library aborts around x > 700 due to float underflow.
    if (x > 700.0) {
        return 0.0; 
    }
    
    // 2. Wrap the function in a try-catch block. 
    // If ANY math exception happens, return 0.0 instead of crashing the whole C program.
    try {
        return std::cyl_bessel_k(nu, x);
    } catch (...) {
        return 0.0;
    }
}