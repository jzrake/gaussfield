
/*------------------------------------------------------------------------------
 * FILE: random.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Contains a class for random number streams
 *
 * REFERENCES: Numerical Recipes, Third Edition: Section 7.1
 *
 *
 *------------------------------------------------------------------------------
 */

#include "random.hpp"

RandomNumberStream::RandomNumberStream(int n)
{
  v = 4101842887655102017LL;
  w = 1;

  u = n^v; int64();
  v =   u; int64();
  w =   v; int64();
}
double RandomNumberStream::RandomDouble()
{
  return int64() * 5.42101086242752217E-20f;
}
double RandomNumberStream::RandomDouble(double a, double b)
{
  return int64() * 5.42101086242752217E-20f * (b-a) + a;
}
unsigned long long RandomNumberStream::RandomInteger()
{
  return int64();
}
unsigned long long RandomNumberStream::RandomInteger(int a, int b)
{
  return a + int64() % (b-a);
}
unsigned long long RandomNumberStream::int64()
{
  u  = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v>>17; v ^= v<<31; v ^= v>>8;
  w  = 4294957665U*(w & 0xffffffff) + (w>>32);
  unsigned long long x = u^(u<<21); x ^= x>>35; x ^= x<<4;
  return (x+v)^w;
}
