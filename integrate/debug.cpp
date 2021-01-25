#include "intg.h"


void print(vector<double>* arr, int n, string s)
{
  cout << s << "==========================" << endl;
  for (int i = 0; i < n; i++)
    cout << (*arr)[i] << endl;
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
  cout << "DEBUG::::: " << i << endl;
  cout << "address: " << result << endl;
  print(result, nf, "result");
  print(fB_omegas_, nf, "fB");
  print(fA_omegas_, nf, "fA");
}


int main()
{
  //
}
