
/*------------------------------------------------------------------------------
 * FILE: random.hpp
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

#ifndef __RandomNumberStream_HEADER__
#define __RandomNumberStream_HEADER__

class RandomNumberStream
{
private:
  unsigned long long u,v,w;

public:
  RandomNumberStream(int n=1);
  double RandomDouble();
  double RandomDouble(double a, double b);
  unsigned long long RandomInteger();
  unsigned long long RandomInteger(int a, int b);
private:
  unsigned long long int64();
} ;

#endif // __RandomNumberStream_HEADER__
