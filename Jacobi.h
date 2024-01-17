#ifndef DOMAIN_DECOMP_JACOBI_H
#define DOMAIN_DECOMP_JACOBI_H

#include <cmath>
#include <iostream>
#include <array>

struct Jacobi {
private:
    double alpha, nu, a, b, ax, bx, ay, by, hx, hy, tol_j, tol_t;
    uint64_t M, N, iter_Jmax;
    double Tmax;
public:
    Jacobi(
        double alpha,
        double nu,
        double a,
        double b,
        double ax,
        double bx,
        double ay,
        double by,
        uint64_t M,
        uint64_t N,
        double hx,
        double hy,
        uint64_t iterJmax,
        double Tmax,
        double tolJ,
        double tolT
    ) : alpha(alpha)
    , nu(nu)
    , a(a)
    , b(b)
    , ax(ax)
    , bx(bx)
    , ay(ay)
    , by(by)
    , hx(hx)
    , hy(hy)
    , tol_j(tolJ)
    , tol_t(tolT)
    , M(M)
    , N(N)
    , iter_Jmax(iterJmax)
    , Tmax(Tmax) {}

    void run();

};


#endif //DOMAIN_DECOMP_JACOBI_H
