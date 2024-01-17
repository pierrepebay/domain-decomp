#include "Jacobi.h"

void Jacobi::run() {
    // Coefficients
    double ht = 1.e-3;
    double alpha = 1/ht, nu = 1.0;
    double a = 3.14, b = 3.14;
    double lambda;
    // Domaine spatial
    double ax = 0.0, bx = 1.0;
    double ay = 0.0, by = 1.0;
    const uint64_t M = 10, N = 10;
    double hx = (bx-ax)/(double)M, hy = (by-ay)/(double)N;

    // travail
    int iter_Jmax = 10*M*N;
    int i, j, ij, imj, ipj, ijm, ijp;
    int iter_t, iter_j;
    double s, t;
    double hx2, hy2, w0, x, y;
    double Tmax = 1000;

    // Tolerance
    double tol_j = 1.e-5, tol_t = 1.e-8;
    double err_j, err_t;


    // Allocation tableaux
    std::array<double, (M+2)*(N+2)> u = {};
    std::array<double, (M+2)*(N+2)> uu = {};
    std::array<double, (M+2)*(N+2)> v = {};
    std::array<double, (M+2)*(N+2)> vv = {};
    std::array<double, (M+2)*(N+2)> f = {};
    std::array<double, (M+2)*(N+2)> ff = {};

    // Initialisation u et f
    for (i = 0; i < M+2; i++)
        for (j = 0; j < N+2; j++){
            x = ax + (double)i*hx;
            y = ay +(double)j*hy;
            w0 = sin(a*x)*cos(b*y);
            ij = i*(N+2) + j;
            f[ij] = (a*a+b*b)*w0;
            u[ij] = 0.0;
            if ( i == 0 || i == M+1 || j == 0 || j == N+1) u[ij] = w0;
            uu[ij] = 0.0;
        }
    lambda = alpha+2.0/hx/hx+2.0/hy/hy;
    hx2 = nu/hx/hx;
    hy2 = nu/hy/hy;

    // Boucle en temps
    iter_t = 0;
    t = 0.0;
    err_t = 1.0;
    while (err_t > tol_t && t < Tmax){
        iter_t++;
        t = (double)iter_t*ht;

        // nouveau 2nd membre (Euler implicite)

        for (i = 0; i < M+2; i++)
            for (j = 0; j < N+2; j++){
                ij = i*(N+2)+j;
                ff[ij] = f[ij] + u[ij]/ht;
                vv[ij] = u[ij];
                v[ij] = u[ij];
            }

        iter_j = 0;
        err_j = 1.0;
        // Boucle de Jacobi
        while (err_j > tol_j && iter_j<iter_Jmax){
            iter_j++;

            // mise à jour Jacobi
            for (i = 1; i < M+1; i++)
                for (j = 1; j < N+1; j++){
                    ij = i*(N+2)+j;
                    imj = (i-1)*(N+2)+j; ipj = (i+1)*(N+2)+j;
                    ijm = i*(N+2)+j-1; ijp = i*(N+2)+j+1;
                    v[ij] =(ff[ij] + hx2*(vv[imj]+vv[ipj]) + hy2*(vv[ijm]+vv[ijp]))/lambda;
                }

            // test d'arret Jacobi & vv <- v
            err_j=0;
            s = 0.0;
            for (i = 1; i < M+1; i++)
                for (j = 1; j < N+1; j++){
                    ij = i*(N+2)+j;
                    err_j += (v[ij]-vv[ij])*(v[ij]-vv[ij]);
                    s += v[ij]*v[ij];
                    vv[ij] = v[ij];
                }
            err_j = sqrt(err_j/s);
            std::cout << "----Jacobi iter=" << iter_j << " err=" << err_j << std::endl;
        }

        // récupération de la solution de Jacobi
        for (i = 1; i < M+1; i++)
            for (j = 1; j < N+1; j++){
                ij = i*(N+2)+j;
                u[ij] = v[ij];
            }

        // test solution stationnaire et uu <- u
        err_t = 0.0;
        s = 0.0;
        for (i = 1; i < M+1; i++)
            for (j = 1; j < N+1; j++){
                ij = i*(N+2)+j;
                err_t += (u[ij]-uu[ij])*(u[ij]-uu[ij]);
                s += u[ij]*u[ij];
                uu[ij] = u[ij];
            }
        err_t = sqrt(err_t/s);

        std::cout << "----iter=" << iter_t << " err=" << err_t << std::endl;
    }
}
