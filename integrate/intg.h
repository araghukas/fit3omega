#include <complex>
#include <vector>
using namespace std;


class IntegralTermBT_EQ1 {
  const int nf;
  const vector<double> omegas;
  const double b;
  const double lambda_i;
  const double lambda_f;
  const int N;

  int nl;
  double* ds;
  double* kxs;
  double* kys;
  double* Cvs;
  char b_type;

  complex<double>* result;
  complex<double>* fB_omegas_;
  complex<double>* phi_omegas_;
  complex<double>* fA_omegas_;
  complex<double>* _AB_next;
  complex<double>* _kkB;
  complex<double>* _tanh_term;
  complex<double>* _integrand;

  complex<double>* fB(int i_layer, double lambda);
  complex<double>* phi(int i_layer, double lambda);
  complex<double>* fA(int i_layer, double lambda);
  complex<double>* integrand(double lambda);
  complex<double> sinc_sq(double x);

  void init_memory()
  {
    result      = new complex<double>[nf];
    fB_omegas_  = new complex<double>[nf];
    phi_omegas_ = new complex<double>[nf];
    fA_omegas_  = new complex<double>[nf];

    _AB_next    = new complex<double>[nf];
    _kkB        = new complex<double>[nf];
    _tanh_term  = new complex<double>[nf];
    _integrand  = new complex<double>[nf];
  }

public:
  IntegralTermBT_EQ1(vector<double> omegas,
                     double b,
                     double lambda_i,
                     double lambda_f,
                     int N)

    : nf{static_cast<int>(omegas.size())},
      omegas{omegas},
      b{b},
      lambda_i{lambda_i},
      lambda_f{lambda_f},
      N{N}
  {
    init_memory();
    nl = 0;
  }

  ~IntegralTermBT_EQ1()
  {
    delete[] result;
    delete[] fB_omegas_;
    delete[] phi_omegas_;
    delete[] fA_omegas_;

    delete[] _AB_next;
    delete[] _kkB;
    delete[] _tanh_term;
    delete[] _integrand;
  }

  void debug(int i);
  int get_nf() { return nf; }
  complex<double>* integrand(int nl,
                             double lambda,
                             double* d_,
                             double* kx_,
                             double* ky_,
                             double* Cv_,
                             char b_type);

  complex<double>* integral(int nl,
                            double* d_,
                            double* kx_,
                            double* ky_,
                            double* Cv_,
                            char b_type);
};
