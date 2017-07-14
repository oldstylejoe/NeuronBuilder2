/***************************************************************************
                          rand3.h  -  description
                             -------------------
    begin                : Tue Apr 18 2000
    copyright            : (C) 2000 by joe snider
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/*random number generator class
 * 
 *I haven't read it, but NR cites
 * Knuth, D.E. 1981, Seminumerical Algorithms, 2nd ed., vol2 of The Art of 
 *     Computer Programming (Reading, MA: Addison-Wessley), chpt 3.2-3.3
 * */

//random number generator yanked from Num Rec ch 7.1
//Uses a subtractive algorithm by Knuth

#ifndef RAND3
#define RAND3

#include <math.h>

#define MBIG 1000000000
#define MZ 0
#define FAC (1.0/MBIG)
#define MSEED 161803398
#define SEED_BIG 1000000

//According to Knuth, any large MBIG, and any smaller (but still large)
//MSEED can be substituted for the above values.

class CRand3{
 public:
   //nothing special about 7
   CRand3(unsigned long iseed = 7) {
      long mj,mk;
      int i,ii,k;
   
      iseed %= SEED_BIG;             //make sure the seed is large enough
      mj=MSEED - long(iseed);
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
     }
      for (k=1;k<=4;k++) {
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
      }
      inext=0;
      inextp=31; //The constant 31 is special see Knuth   

      //Default power (1) and range (1., 100.) for the power law generator.
      powerExponentHelper = 1./(1.+1.);
      powerHelper1 = pow(100., 1.+1.) - pow(1., 1.+1.);
      powerHelper2 = pow(1., 1.+1.);
   }


 public:
   //returns a random deviate between 0 and 1
   double next() {
      if (++inext == 56) inext=1;
      if (++inextp == 56) inextp=1;
      long mj = ma[inext]-ma[inextp];
      if (mj < MZ) mj += MBIG;
      ma[inext]=mj;
      return mj*FAC;
   }
   double operator() () {
      return next();
   }

   //returns an integer in the range [0,n)
   //appropriate for use in std::random_shuffle
   int operator() (const int& inN) {
      return int(double(inN)*(*this)());
   }

   //returns an Gaussian random variable using the Box-Muller algorithm.
   //inStdDev is the standard deviation.
   double Gaussian(const double& inStdDev) {
      double u, v, s;
      do {
         u = 2.*next()-1.;
         v = 2.*next()-1.;
         s = u*u + v*v;
      } while ( s >= 1. || s == 0.);
      return inStdDev * u * sqrt(-2.*log(s)/s);
      //return inStdDev * v * sqrt(-2.*log(s)/s);  //This is another random, but ignore it.
   }

   double GetRand() {return (*this)();}

   //Set the parameters for power law generation.
   //Default are power 1, low = 1., high = 100.
   //Note, the power is the actual power; no implicit minus.
   //No checking, but must satisfy 0 < inLow < inHigh to give usable results.
   //Power = -1 has to be handled seperately (or use something like -0.999)
   void SetPowerParameters(const double& inLow, const double& inHigh, const double& inPower) {
      powerExponentHelper = 1./(inPower+1.);
      powerHelper1 = pow(inHigh, inPower+1.) - pow(inLow, inPower+1.);
      powerHelper2 = pow(inLow, inPower+1.);
   }

   double Power() {
      return pow(powerHelper1*GetRand() + powerHelper2, powerExponentHelper);
   }
   
 private:
   int inext,inextp;
   // the value 56 is special and should not be changed
   long ma[56];
   double powerHelper1, powerHelper2, powerExponentHelper;
};

#endif
