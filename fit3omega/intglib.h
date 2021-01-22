#include <stdio.h>
#include <complex.h>

double complex B(int i, double L, double omega,
  const double kxs[], const double kys[], const double Cvs[]);


// semi-infinite BC's function variants
double complex A_s(int i, int n_layers, double L, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double f_s(int n_layers, double L, double b, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double integrate_f_s(int n_layers, double xi, double xf, int N, double b, double omega,
  const double ds[], const double kxs[], const double kys[], const double Cvs[]);


// adiabatic BC's function variants
double complex A_a(int i, int n_layers, double L, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double f_a(int n_layers, double L, double b, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double integrate_f_a(int n_layers, double xi, double xf, int N, double b, double omega,
  const double ds[], const double kxs[], const double kys[], const double Cvs[]);


// isothermal BC's function variants
double complex A_o(int i, int n_layers, double L, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double f_o(int n_layers, double L, double b, double omega,
  const double kxs[], const double kys[], const double Cvs[], const double ds[]);

double integrate_f_o(int n_layers, double xi, double xf,
  int N, double b, double omega,
  const double ds[], const double kxs[], const double kys[], const double Cvs[]);
