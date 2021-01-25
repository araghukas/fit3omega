#include <iostream>
#include "intg.h"


void print(double* arr, int n, string s)
{
  cout << s << "==========================" << endl;
  for (int i = 0; i < n; i++)
    cout << arr[i] << endl;
  cout << endl;
}


void print(vector<double>* arr, int n, string s)
{
  cout << s << "==========================" << endl;
  for (int i = 0; i < n; i++)
    cout << (*arr)[i] << endl;
  cout << endl;
}


void print(vector<double> arr, int n, string s)
{
  cout << s << "==========================" << endl;
  for (int i = 0; i < n; i++)
    cout << arr[i] << endl;
  cout << endl;
}


void print(complex<double>* arr, int n, string s)
{
  cout << s << "==========================" << endl;
  for (int i = 0; i < n; i++)
    cout << arr[i] << endl;
  cout << endl;
}


void IntegralTermBT_EQ1::debug(int i)
{
  // cout << "DEBUG::::: " << i << endl;
  // cout << "address: " << result << endl;
  // print(omegas, nf, "omegas");
  // print(result, nf, "result");
  // print(kxs, nl, "kxs");
  // print(Cvs, nl, "Cvs");
  cout << "result_memory: " << result;
}


int main()
{
  vector<double> omegas {3141.592654, 3453.252493, 3795.830362,
                         4172.393466, 4586.313290, 5041.295785,
                         5541.414549, 6091.147298, 6695.415960,
                         7359.630737, 8089.738547, 8892.276271,
                         9774.429274, 10744.095744, 11809.957402,
                         12981.557234, 14269.384935, 15684.970819,
                         17240.989064, 18951.371178, 22898.000490,
                         25169.583073, 27666.516662, 30411.157069,
                         33428.077902, 36744.290581, 40389.486176,
                         44396.301242, 48800.610024, 53641.845651,
                         58963.353193, 64812.777741, 71242.490988,
                         78310.060131, 86078.763286, 94618.156040,
                         104004.694196, 114322.418313, 125663.706144};
  double b = 3.4285e-05;
  double lambda_i = 1e-3;
  double lambda_f = 1e7;
  int N = 200;

  double d_[] {5e-08, 3e-06, 0.0003};
  double kx_[] {1, 10, 100};
  double ky_[] {1, 10, 100};
  double Cv_[] {1706100.0, 2124260.0, 1630300.0};
  char b_type = 's';
  int nl = 3;

  IntegralTermBT_EQ1 INTG {omegas, b, lambda_i, lambda_f, N};

  for (int k = 0; k < 1000; k++) {
    complex<double>* res = INTG.integral(nl,d_,kx_,ky_,Cv_,b_type);
    for (int i = 0; i < 8; i++)
      cout << res[i].real() << ", ";
    cout << endl;
  }
}
