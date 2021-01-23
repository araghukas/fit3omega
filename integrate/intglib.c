#include "intglib.h"


// -----------------------------------------------------------------------------
//             IMPORTANT NOTE: i = 1, 2, 3, ... STRICTLY POSITIVE NONZERO!!
//                       this is true everywhere, i indicates layer number.
//
//                       There is no 0th layer; layers 1 is the uppermost one!
// -----------------------------------------------------------------------------


/* A and B coefficients */
double complex B(int i, double L, double omega, const double kxs[], const double kys[], const double Cvs[]){
    // i is the LAYER index so shift down for array index
    return csqrt((kxs[i-1] / kys[i-1])*L*L + I*2.*omega*Cvs[i-1]/kys[i-1]);
}


// semi-infinite case
double complex A_s(int i, int n_layers, double L, double omega, const double kxs[], const double kys[],  const double Cvs[],  const double ds[]){

    double complex A_n = -1. + 0.*I;

    double complex A_i, A_ii, B_i, B_ii;
    double complex num, den, kB_factor, tanh_term;

    if (i == n_layers) {
        A_i = A_n;
    }
    else {
        B_i = B(i, L, omega, kxs, kys, Cvs);
        B_ii = B(i+1, L, omega, kxs, kys, Cvs);
        A_ii = A_s(i+1, n_layers, L, omega, kxs, kys, Cvs, ds);
        kB_factor = kys[i] * B_ii / kys[i-1] / B_i;
        tanh_term = ctanh(B_i * ds[i-1]);

        num = A_ii * kB_factor - tanh_term;
        den = 1. - A_ii * kB_factor * tanh_term;

        A_i = num / den;
    }

    return A_i;
}


// adiabatic case
double complex A_a(int i, int n_layers, double L, double omega, const double kxs[], const double kys[], const double Cvs[], const double ds[]){

    double complex B_n = B(n_layers, L, omega, kxs, kys, Cvs);
    double complex A_n = -ctanh(B_n * ds[n_layers - 1]);

    double complex A_i, A_ii, B_i, B_ii;
    double complex num, den, kB_factor, tanh_term;

    if (i == n_layers) {
        A_i = A_n;
    }
    else {
        B_i = B(i, L, omega, kxs, kys, Cvs);
        B_ii = B(i+1, L, omega, kxs, kys, Cvs);
        A_ii = A_a(i+1, n_layers, L, omega, kxs, kys, Cvs, ds);
        kB_factor = kys[i] * B_ii / kys[i-1] / B_i;
        tanh_term = ctanh(B_i * ds[i-1]);

        num = A_ii * kB_factor - tanh_term;
        den = 1. - A_ii * kB_factor * tanh_term;

        A_i = num / den;
    }

    return A_i;
}


// isothermal case
double complex A_o(int i, int n_layers, double L, double omega, const double kxs[], const double kys[], const double Cvs[], const double ds[]){

    double complex B_n = B(n_layers, L, omega, kxs, kys, Cvs);
    double complex A_n = -1. / ctanh(B_n * ds[n_layers - 1]);

    double complex A_i, A_ii, B_i, B_ii;
    double complex num, den, kB_factor, tanh_term;

    if (i == n_layers) {
        A_i = A_n;
    }
    else {
        B_i = B(i, L, omega, kxs, kys, Cvs);
        B_ii = B(i+1, L, omega, kxs, kys, Cvs);
        A_ii = A_o(i+1, n_layers, L, omega, kxs, kys, Cvs, ds);
        kB_factor = kys[i] * B_ii / kys[i-1] / B_i;
        tanh_term = ctanh(B_i * ds[i-1]);

        num = A_ii * kB_factor - tanh_term;
        den = 1. - A_ii * kB_factor * tanh_term;

        A_i = num / den;
    }

    return A_i;
}






/* INTEGRAND FUNCTIONS */
// semi-infinite case
double f_s(int n_layers, double L, double b, double omega, const double kxs[], const double kys[], const double Cvs[], const double ds[]){
    double complex B_1 = B(1, L, omega, kxs, kys, Cvs);
    double complex A_1 = A_s(1, n_layers, L, omega, kxs, kys, Cvs, ds);
    double complex result =  1. / (A_1 * B_1) * (sin(b * L) / (b * L)) * (sin(b * L) / (b * L));
    return creal(result);
}


// adiabatic case
double f_a(int n_layers, double L, double b, double omega, const double kxs[], const double kys[], const double Cvs[], const double ds[]){
    double complex B_1 = B(1, L, omega, kxs, kys, Cvs);
    double complex A_1 = A_a(1, n_layers, L, omega, kxs, kys, Cvs, ds);
    double complex result = 1. / (A_1 * B_1) * (sin(b * L) / (b * L)) * (sin(b * L) / (b * L));
    return creal(result);
}


// isothermal case
double f_o(int n_layers, double L, double b, double omega, const double kxs[], const double kys[], const double Cvs[], const double ds[]){
    double complex B_1 = B(1, L, omega, kxs, kys, Cvs);
    double complex A_1 = A_o(1, n_layers, L, omega, kxs, kys, Cvs, ds);
    double complex result = 1. / (A_1 * B_1) * (sin(b * L) / (b * L)) * (sin(b * L) / (b * L));
    return creal(result);
}






/* INTEGRAL FUNCTIONS */
// semi-infinite case
double integrate_f_s(int n_layers, double xi, double xf, int N, double b, double omega, const double ds[], const double kxs[], const double kys[], const double Cvs[]){

    double h = (xf - xi) / N;

    double result = .5 * h * (f_s(n_layers, xi, b, omega, kxs, kys, Cvs, ds) + f_s(n_layers, xf, b, omega, kxs, kys, Cvs, ds));

    for (int k = 1; k < N; k++) {
        result += h * f_s(n_layers, xi + k*h, b, omega, kxs, kys, Cvs, ds);
    }

    return result;
}


// adiabatic case
double integrate_f_a(int n_layers, double xi, double xf, int N, double b, double omega, const double ds[], const double kxs[], const double kys[], const double Cvs[]){

    double h = (xf - xi) / N;

    double result = .5 * h * (f_a(n_layers, xi, b, omega, kxs, kys, Cvs, ds) + f_a(n_layers, xf, b, omega, kxs, kys, Cvs, ds));

    for (int k = 1; k < N; k++) {
        result += h * f_a(n_layers, xi + k*h, b, omega, kxs, kys, Cvs, ds);
    }

    return result;
}


// isothermal case
double integrate_f_o(int n_layers, double xi, double xf, int N, double b, double omega, const double ds[], const double kxs[], const double kys[], const double Cvs[]){

    double h = (xf - xi) / N;

    // add endpoints, `k == 0` and `k == N` cases
    double result = .5 * h * (f_o(n_layers, xi, b, omega, kxs, kys, Cvs, ds) + f_o(n_layers, xf, b, omega, kxs, kys, Cvs, ds));

    for (int k = 1; k < N; k++) {
        result += h * f_o(n_layers, xi + k*h, b, omega, kxs, kys, Cvs, ds);
    }

    return result;
}
