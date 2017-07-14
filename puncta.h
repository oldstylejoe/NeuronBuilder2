//Joe Snider
//2/05
//
//Puncta allow neuron branches to be stabilized.
//They're just derived from vector and have an occupied member.
//
//3/06
//Changed to include a strength member.
//Branches are retracted past 'weak' puncta.
//
//12/08
//Strength was apparently a bad idea because I never did anything with it.
//It's used now to label puncta as stabilizers (str=1.) or removed (str=0)

#include <iostream>
#include <vector>

#include "point.h"

#ifndef PUNCTA__
#define PUNCTA__

using namespace std;

class CPuncta: public CPoint{
public:
   //typedefs
   //typedef vector<double>::iterator iterator;
   //typedef vector<double>::const_iterator const_iterator;

public:
   //constructors
   CPuncta(): m_bOccupied(false), m_dStrength(0.), CPoint() {}
   CPuncta(const CPoint& inPoint): m_bOccupied(false), m_dStrength(0.), CPoint(inPoint) {}
   CPuncta(const CPuncta& inPuncta): CPoint() {
      *this = inPuncta;
   }

   CPuncta& operator=(const CPuncta& inPuncta) {
      const_iterator iter = inPuncta.begin();
      const_iterator iter_end = inPuncta.end();
      for(; iter != iter_end; ++iter) {
         push_back(*iter);
      }
      m_bOccupied = inPuncta.GetOccupied();
      m_dStrength = inPuncta.GetStrength();
      return *this;
   }

   virtual ~CPuncta() {}

public:
   //gets and sets
   bool GetOccupied() const {return m_bOccupied;}
   void SetOccupied(const bool& inOccupied) {m_bOccupied = inOccupied;}

   double GetStrength() const {return m_dStrength;}
   void SetStrength(const double& inStrength) {m_dStrength = inStrength;}

public:
   //interface
   bool IsOccupied() const {
      return m_bOccupied;
   }

private:
   //helpers

private:
   //data members
   bool m_bOccupied;
   double m_dStrength;

};

#endif //PUNCTA__

