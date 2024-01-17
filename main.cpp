#include "Jacobi.h"
#include "Domain.h"

int main(){
    // Coefficients
    const double alpha = 1.e-3;
    const double nu = 1.0;
    const double a = 6.28, b = 6.28;

    // Spatial domain
    const double ax = 0.0, bx = 1.0;
    const double ay = 0.0, by = 1.0;
    const uint64_t M = 10, N = 10;
    const double hx = (bx-ax)/M, hy = (by-ay)/N;
    
    // travail
    int iter_Jmax = 10*M*N;
    double Tmax = 1000;
    
    // Tolerance
    double tol_j = 1.e-5, tol_t = 1.e-8;

    Jacobi j(alpha, nu, a, b, ax, bx, ay, by, M, N, hx, hy, iter_Jmax, Tmax, tol_j, tol_t);

    return 0;
}