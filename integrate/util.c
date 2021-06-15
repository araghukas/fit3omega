#include "integrate.h"


// writes into array points linear in log-space, between min and max
void make_logspace(double *arr, double min, double max, int size)
{
  for (int k = 0; k < size; k++)
    arr[k] = min * pow(max / min, ((double) k) / (size - 1));
}


// a function that computes (sinc(x))^2
double sinc_sq(double x)
{
  double sinc = sin(x) / x;
  return sinc*sinc;
}
