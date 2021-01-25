#include <complex>
#include <vector>
using namespace std;


class IntegralTermBT_EQ1 {
  const int nf;
  const int nl;
  const vector<double> omegas;
  const vector<double> ds;
  const vector<double> kxs;
  const vector<double> kys;
  const vector<double> Cvs;
  const double b;
  const double lambda_i;
  const double lambda_f;
  const int N;
  const char b_type;

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
  complex<double> sinc_sq(double x);  

public:
  //  void debug();  
  complex<double>* integrand(double lambda);
  complex<double>* integral();
  
  IntegralTermBT_EQ1(vector<double> omegas,
                     vector<double> ds,
                     vector<double> kxs,
                     vector<double> kys,
                     vector<double> Cvs,
                     double b,
                     double lambda_i,
                     double lambda_f,
                     int N,
                     char b_type)
    : nf{static_cast<int>(omegas.size())},
      nl{static_cast<int>(ds.size())},
      omegas{omegas},
      ds{ds},
      kxs{kxs},
      kys{kys},
      Cvs{Cvs},
      b{b},
      lambda_i{lambda_i},
      lambda_f{lambda_f},
      N{N},
      b_type{b_type}
  {
    result = new complex<double>[nf];
    fB_omegas_ = new complex<double>[nf];
    phi_omegas_ = new complex<double>[nf];
    fA_omegas_ = new complex<double>[nf];

    _AB_next = new complex<double>[nf];
    _kkB = new complex<double>[nf];
    _tanh_term = new complex<double>[nf];
    _integrand = new complex<double>[nf];
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

  complex<double> integrate(int lambda_i, int lambda_f, int N);
};
