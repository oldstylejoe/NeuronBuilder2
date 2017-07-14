//Joe Snider
//1/05
//
//Store a point.  Just a sized vector<double> (may add on the ability to add, etc...)

#include <iostream>
#include <vector>
#include <math.h>

#ifndef POINT__
#define POINT__

using namespace std;

const double EQUAL_DISTANCE = 1.e-12; //used in testing equality of two points

class CPoint: public vector<double>{
public:
   //typedefs
   typedef vector<double>::iterator iterator;
   typedef vector<double>::const_iterator const_iterator;

public:
   //constructors
   CPoint(): vector<double>() {}
   CPoint(const CPoint& inPoint): vector<double>() {
      *this = inPoint;
   }

   CPoint& operator=(const CPoint& inPoint) {
      clear();
      const_iterator iter = inPoint.begin();
      const_iterator iter_end = inPoint.end();
      for(; iter != iter_end; ++iter) {
         push_back(*iter);
      }
      return *this;
   }

   virtual ~CPoint() {}

public:
   //gets and sets
   int GetDimension() const  {return (int)size();}
   void SetDimension(const int& inSize) {resize(inSize);}

public:
   //interface

   //Return the length (distance from 0).
   double Length() const {
      return sqrt(Dot(*this));
   }

   //Return the length squared (distance from 0).
   //Saves a sqrt if possible instead of Length().
   double LengthSquared() const {
      return Dot(*this);
   }

   //Return true if this vector has exactly zero length.
   //Assumes that the points are real, that is it just tests that the square is non-zero.
   bool IsZero() const {
      return Dot(*this) == 0.0;
   }

   //different dimensions means not equal
   bool operator==(const CPoint& inPoint) const {
      if( inPoint.GetDimension() != GetDimension() ) { 
         return false;
      }
      for(int i = 0; i < GetDimension(); ++i) {
         if( fabs( inPoint[i] - at(i) ) > EQUAL_DISTANCE ) {
            return false;
         }
      }
      return true;
      /*const_iterator iterR = begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         if( fabs(*iterR - *iterA) > EQUAL_DISTANCE ) {
            return false;
         }
      }
      return true;*/
   }

   //different dimensions means not equal
   bool operator!=(const CPoint& inPoint) const {
      return !(*this == inPoint);
   }

   CPoint operator+(const CPoint& inPoint) const {
      CPoint pntReturn = *this;
      iterator iterR = pntReturn.begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         *iterR += *iterA;
      }
      return pntReturn;
   }

   CPoint operator-(const CPoint& inPoint) const {
      CPoint pntReturn = *this;
      iterator iterR = pntReturn.begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         *iterR -= *iterA;
      }
      return pntReturn;
   }

   CPoint& operator+=(const CPoint& inPoint) {
      iterator iterR = begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         *iterR += *iterA;
      }
      return *this;
   }

   CPoint& operator-=(const CPoint& inPoint) {
      iterator iterR = begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         *iterR -= *iterA;
      }
      return *this;
   }

   CPoint& operator*=(const double& inScale) {
      iterator iter = begin();
      iterator iter_end = end();
      for(; iter != iter_end; ++iter) {
         *iter *= inScale;
      }
      return *this;
   }

   CPoint operator*(const double& inScale) const {
      CPoint pntReturn = *this;
      iterator iter = pntReturn.begin();
      iterator iter_end = pntReturn.end();
      for(; iter != iter_end; ++iter) {
         *iter *= inScale;
      }
      return pntReturn;
   }

   double operator*(const CPoint& inPoint) const {
      return Dot(inPoint);
   }

   double Dot(const CPoint& inPoint) const {
      double dReturn = 0.;
      const_iterator iterR = begin();
      const_iterator iterA = inPoint.begin();
      const_iterator iterA_end = inPoint.end();
      for(; iterA != iterA_end; ++iterA, ++iterR) {
         dReturn += (*iterR)*(*iterA);
      }
      return dReturn;
   }

   //Scale in each direction seperately.
   //inScale may be longer than necessary (i.e. have dimension 3 when this has dim 2).
   CPoint& Scale(const CPoint& inScale) {
      iterator iterThis = begin();
      iterator iterThis_end = end();
      const_iterator iterScaler = inScale.begin();
      for(; iterThis != iterThis_end; ++iterScaler, ++iterThis) {
         *iterThis *= *iterScaler;
      }
      return *this;
   }

   //Apply the input matrix.
   //Type T must have operator[][] access to a double (it's a matrix).
   template<class T>
   CPoint& LinearTransform(const T& inM) {
      CPoint pTemp = *this;
      for(unsigned i = 0; i < size(); ++i) {
         at(i) = 0.;
         for(unsigned j = 0; j < size(); ++j) {
            const double x = (pTemp[j]) * (inM[i][j]);
            at(i) += x;
         }
      }
      return *this;
   }

private:
   //helpers

private:
   //data members

};

//allow streaming
ostream& operator<<(ostream& inOut, const CPoint& inPoint) {
   if(inPoint.size() > 0) {
      CPoint::const_iterator iter = inPoint.begin();
      CPoint::const_iterator iter_end = inPoint.end();
      --iter_end;
      for(; iter != iter_end; ++iter) {
         inOut << *iter << " ";
      }
      inOut << *iter_end;
   }
   return inOut;
}

//have to handle double*CPoint externally
CPoint operator*(const double& inScale, const CPoint& inPoint) {
   CPoint pntReturn = inPoint;
   CPoint::iterator iter = pntReturn.begin();
   CPoint::iterator iter_end = pntReturn.end();
   for(; iter != iter_end; ++iter) {
      *iter *= inScale;
   }
   return pntReturn;
}

#endif //POINT__
