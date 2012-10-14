
/*------------------------------------------------------------------------------
 * FILE: ornuhl.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for sampling Ornstein-Uhlenbeck deviates
 *
 * REFERENCES: http://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process
 *
 *------------------------------------------------------------------------------
 */

#ifndef __OrnsteinUhlenbeckProcess_HEADER__
#define __OrnsteinUhlenbeckProcess_HEADER__

#include <complex>
#include <cmath>
#include <fstream>
#include "random.hpp"


class OrnsteinUhlenbeckProcess
/* 
 * The time correlation of O-U process is <x(t)*x(t+s)> = exp(-|s|)
 *
 * The variance is given by (sigma^2) / (2 * theta)
 *
*/
{
private:
  RandomNumberStream stream;
  double theta, sigma;

public:
  OrnsteinUhlenbeckProcess(double theta, double sigma, int seed)
    : stream(seed),
      theta(theta),
      sigma(sigma) { }

  double Deviate(double x, double dt)
  {
    const double dW = SampleNormal(0, dt);
    return -theta*x*dt + sigma*dW;
  }
  std::complex<double> Deviate(std::complex<double> z, double dt)
  {
    const std::complex<double> dW(SampleNormal(0, dt), SampleNormal(0, dt));
    return -theta*z*dt + sigma*dW;
  }

private:
  double SampleNormal(double mu, double sig2)
  {
    const double C = stream.RandomDouble(-1,1);
    return sqrt(2*sig2) * erfinv(C) + mu;
  }
  double erfinv(double y) const
  {
    double x = 0.5*sqrt(M_PI)*(y + M_PI/12.0*pow(y,3));
    double f,g;
    do {
      f  = erf(x) - y;
      g  = 2.0/sqrt(M_PI) * exp(-x*x);
      x -= f/g;
    } while (fabs(f/g) > 1e-12);
    return x;
  }
} ;

#endif // __OrnsteinUhlenbeckProcess_HEADER__
