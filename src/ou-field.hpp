
/*------------------------------------------------------------------------------
 * FILE: ou-field.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for generating vector fields from an
 *      Ornstein–Uhlenbeck process
 *
 * REFERENCES:
 *
 *   http://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process
 *   Wolfram Schmidt, Wolfgang Hillebrandt, Jens C. Niemeyerb (2005)
 *
 *------------------------------------------------------------------------------
 */

#ifndef __StochasticVectorField_HEADER__
#define __StochasticVectorField_HEADER__

#include <valarray>
#include <vector>
#include <complex>
#include <fstream>
#include "ornuhl.hpp"

class StochasticVectorField2d
{
private:
  int k1, numk;
  double P0, zeta;
  std::vector<OrnsteinUhlenbeckProcess> OuProcesses;
  std::valarray<std::complex<double> > Fx, Fy;
  std::valarray<std::complex<double> > Kx, Ky;

public:
  void Initialize(double P0, double zeta, int k1, int seed);
  void DeSerialize(std::string serial);
  void AdvanceField(double dt);
  void SampleField(double x, double y, double z, double *F);
  std::string Serialize();
  double GetZeta() {return zeta;}
  double GetP0() {return P0;}
  int GetK1() {return k1;}
} ;

class StochasticVectorField3d
{
private:
  int k1, numk;
  double P0, zeta;
  std::vector<OrnsteinUhlenbeckProcess> OuProcesses;
  std::valarray<std::complex<double> > Fx, Fy, Fz;
  std::valarray<std::complex<double> > Kx, Ky, Kz;

public:
  void Initialize(double P0, double zeta, int k1, int seed);
  void DeSerialize(std::string serial);
  void AdvanceField(double dt);
  void SampleField(double x, double y, double z, double *F);
  std::string Serialize();
  double GetZeta() {return zeta;}
  double GetP0() {return P0;}
  int GetK1() {return k1;}
} ;

#endif // __StochasticVectorField_HEADER__
