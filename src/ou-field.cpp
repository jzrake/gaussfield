
/*------------------------------------------------------------------------------
 * FILE: ou-field.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for generating vector fields from an OU-process
 *
 * REFERENCES:
 *
 *   http://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process
 *   Wolfram Schmidt, Wolfgang Hillebrandt, Jens C. Niemeyerb (2005)
 *
 *------------------------------------------------------------------------------
 */

#include "ou-field.hpp"


static double StoneProfile(double k)
{
  const double KvecPeak = 2.0;
  return pow(k,6) * exp(-8*k/(2*M_PI*KvecPeak));
}

// Constructors
// -----------------------------------------------------------------------------

void StochasticVectorField2d::Initialize(double P0_, double zeta_, int k1_,
					 int seed)
{
  k1 = k1_;
  numk = (2*k1+1)*(2*k1+1);
  P0 = P0_;
  zeta = zeta_;
  Fx.resize(numk);
  Fy.resize(numk);
  Kx.resize(numk);
  Ky.resize(numk);

  const int si = (2*k1+1);
  const int sj = 1;
  double totpower = 0.0;

  for (int i=-k1; i<=k1; ++i) {
    for (int j=-k1; j<=k1; ++j) {

      const int m = (i+k1)*si + (j+k1)*sj;
      Kx[m] = 2*M_PI*i;
      Ky[m] = 2*M_PI*j;

      const double K[2] = { Kx[m].real(), Ky[m].real() };
      const double k = sqrt(K[0]*K[0] + K[1]*K[1]);

      totpower += StoneProfile(k);
    }
  }

  for (int m=0; m<numk; ++m) {
    const double K[3] = { Kx[m].real(), Ky[m].real() };
    const double k = sqrt(K[0]*K[0] + K[1]*K[1]);
    const double Pk = StoneProfile(k);

    OrnsteinUhlenbeckProcess process(1.0, sqrt(2*P0*Pk/totpower), seed+m);
    OuProcesses.push_back(process);
  }
}

void StochasticVectorField3d::Initialize(double P0_, double zeta_, int k1_,
					 int seed)
{
  k1 = k1_;
  numk = (2*k1+1)*(2*k1+1)*(2*k1+1);
  P0 = P0_;
  zeta = zeta_;
  Fx.resize(numk);
  Fy.resize(numk);
  Fz.resize(numk);
  Kx.resize(numk);
  Ky.resize(numk);
  Kz.resize(numk);

  const int si = (2*k1+1)*(2*k1+1);
  const int sj = (2*k1+1);
  const int sk = 1;
  double totpower = 0.0;

  for (int i=-k1; i<=k1; ++i) {
    for (int j=-k1; j<=k1; ++j) {
      for (int k=-k1; k<=k1; ++k) {

        const int m = (i+k1)*si + (j+k1)*sj + (k+k1)*sk;
        Kx[m] = 2*M_PI*i;
        Ky[m] = 2*M_PI*j;
        Kz[m] = 2*M_PI*k;

        const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
        const double k = sqrt(K[0]*K[0] + K[1]*K[1] + K[2]*K[2]);

	totpower += StoneProfile(k);
      }
    }
  }

  for (int m=0; m<numk; ++m) {
    const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
    const double k = sqrt(K[0]*K[0] + K[1]*K[1] + K[2]*K[2]);
    const double Pk = StoneProfile(k);

    OrnsteinUhlenbeckProcess process(1.0, sqrt(2*P0*Pk/totpower), seed+m);
    OuProcesses.push_back(process);
  }
}

void StochasticVectorField2d::DeSerialize(std::string serial)
{
  std::stringstream stream;
  stream << serial;

  stream.read((char*)&k1, sizeof(int));
  numk = (2*k1+1)*(2*k1+1);

  Kx.resize(numk);
  Ky.resize(numk);

  Fx.resize(numk);
  Fy.resize(numk);

  OrnsteinUhlenbeckProcess dummy(0,0,0);
  OuProcesses.assign(numk, dummy);

  size_t bytes = numk*sizeof(std::complex<double>);

  stream.read((char*)&P0, sizeof(double));
  stream.read((char*)&zeta, sizeof(double));

  stream.read((char*)&Kx[0], bytes);
  stream.read((char*)&Ky[0], bytes);

  stream.read((char*)&Fx[0], bytes);
  stream.read((char*)&Fy[0], bytes);

  stream.read((char*)&OuProcesses[0],
              sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}

void StochasticVectorField3d::DeSerialize(std::string serial)
{
  std::stringstream stream;
  stream << serial;

  stream.read((char*)&k1, sizeof(int));
  numk = (2*k1+1)*(2*k1+1)*(2*k1+1);

  Kx.resize(numk);
  Ky.resize(numk);
  Kz.resize(numk);

  Fx.resize(numk);
  Fy.resize(numk);
  Fz.resize(numk);

  OrnsteinUhlenbeckProcess dummy(0,0,0);
  OuProcesses.assign(numk, dummy);

  size_t bytes = numk*sizeof(std::complex<double>);

  stream.read((char*)&P0, sizeof(double));
  stream.read((char*)&zeta, sizeof(double));

  stream.read((char*)&Kx[0], bytes);
  stream.read((char*)&Ky[0], bytes);
  stream.read((char*)&Kz[0], bytes);

  stream.read((char*)&Fx[0], bytes);
  stream.read((char*)&Fy[0], bytes);
  stream.read((char*)&Fz[0], bytes);

  stream.read((char*)&OuProcesses[0],
              sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());
}

// AdvanceField
// -----------------------------------------------------------------------------
void StochasticVectorField2d::AdvanceField(double dt)
{
  for (int m=0; m<numk; ++m) {
    std::complex<double> dW[2], dV[2];

    const double K[2] = { Kx[m].real(), Ky[m].real() };
    const double kdotk = K[0]*K[0] + K[1]*K[1];

    if (kdotk < 1e-8) continue; // don't advance the zero mode

    dW[0] = OuProcesses[m].Deviate(Fx[m], dt);
    dW[1] = OuProcesses[m].Deviate(Fy[m], dt);

    for (int p=0; p<2; ++p) {
      dV[p] = 0.0;
      for (int q=0; q<2; ++q) {
        const double P = (1 - 2*zeta)*K[p]*K[q]/kdotk + zeta*(p==q);
        dV[p] += P * dW[q];
      }
    }

    Fx[m] += dV[0];
    Fy[m] += dV[1];
  }
}
void StochasticVectorField3d::AdvanceField(double dt)
{
  for (int m=0; m<numk; ++m) {
    std::complex<double> dW[3], dV[3];

    const double K[3] = { Kx[m].real(), Ky[m].real(), Kz[m].real() };
    const double kdotk = K[0]*K[0] + K[1]*K[1] + K[2]*K[2];

    if (kdotk < 1e-8) continue; // don't advance the zero mode

    dW[0] = OuProcesses[m].Deviate(Fx[m], dt);
    dW[1] = OuProcesses[m].Deviate(Fy[m], dt);
    dW[2] = OuProcesses[m].Deviate(Fz[m], dt);

    for (int p=0; p<3; ++p) {
      dV[p] = 0.0;
      for (int q=0; q<3; ++q) {
        const double P = (1 - 2*zeta)*K[p]*K[q]/kdotk + zeta*(p==q);
        dV[p] += P * dW[q] / sqrt(1 - 1.25*zeta + 3*zeta*zeta);
      }
    }

    Fx[m] += dV[0];
    Fy[m] += dV[1];
    Fz[m] += dV[2];
  }
}

// SampleField
// -----------------------------------------------------------------------------
void StochasticVectorField2d::SampleField(double x, double y, double z, double *F)
{
  const std::complex<double> X(0,x);
  const std::complex<double> Y(0,y);
  std::valarray<std::complex<double> > expfac = exp(Kx*X + Ky*Y);

  F[0] = (Fx * expfac).sum().real();
  F[1] = (Fy * expfac).sum().real();
}
void StochasticVectorField3d::SampleField(double x, double y, double z, double *F)
{
  const std::complex<double> X(0,x);
  const std::complex<double> Y(0,y);
  const std::complex<double> Z(0,z);
  std::valarray<std::complex<double> > expfac = exp(Kx*X + Ky*Y + Kz*Z);

  F[0] = (Fx * expfac).sum().real();
  F[1] = (Fy * expfac).sum().real();
  F[2] = (Fz * expfac).sum().real();
}



// Serialize field and random number states to bit array.
// -----------------------------------------------------------------------------
std::string StochasticVectorField2d::Serialize()
{
  std::stringstream stream;
  size_t bytes = numk*sizeof(std::complex<double>);

  stream.write((char*)&k1, sizeof(int));
  stream.write((char*)&P0, sizeof(double));
  stream.write((char*)&zeta, sizeof(double));

  stream.write((char*)&Kx[0], bytes);
  stream.write((char*)&Ky[0], bytes);

  stream.write((char*)&Fx[0], bytes);
  stream.write((char*)&Fy[0], bytes);

  stream.write((char*)&OuProcesses[0],
               sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());

  return stream.str();
}
std::string StochasticVectorField3d::Serialize()
{
  std::stringstream stream;
  size_t bytes = numk*sizeof(std::complex<double>);

  stream.write((char*)&k1, sizeof(int));
  stream.write((char*)&P0, sizeof(double));
  stream.write((char*)&zeta, sizeof(double));

  stream.write((char*)&Kx[0], bytes);
  stream.write((char*)&Ky[0], bytes);
  stream.write((char*)&Kz[0], bytes);

  stream.write((char*)&Fx[0], bytes);
  stream.write((char*)&Fy[0], bytes);
  stream.write((char*)&Fz[0], bytes);

  stream.write((char*)&OuProcesses[0],
               sizeof(OrnsteinUhlenbeckProcess)*OuProcesses.size());

  return stream.str();
}
