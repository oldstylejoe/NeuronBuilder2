//Joe Snider
//5/05
//
//Contains a branch. A multiply forward linked list.

#include <iostream>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <gsl/gsl_integration.h>

#include "point.h"
#include "line.h"
#include "histogram.h"

#ifndef BRANCH
#define BRANCH

using namespace std;
using namespace boost;

//Evaluate the line integral in 3D.
//Used in CBranch::ProductMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction3D(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((bx-ax)*x+ax)*((by-ay)*x+ay)*((bz-az)*x+az), n);
   return dRet;
}

//Evaluate the line integral in 2D.
//Used in CBranch::ProductMoment by the gsl integration routine.
//Takes 5 parameters: power, ax, bx, ay, by.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction2D(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double dRet = pow( ((bx-ax)*x+ax)*((by-ay)*x+ay), n);
   return dRet;
}

class CBranch: public enable_shared_from_this<CBranch> {
public:
   //typedefs
   typedef shared_ptr<CBranch> BranchPointerType;
   typedef vector<BranchPointerType> DaughtersHolderType;
   typedef CLine LineType;

public:
   //constructors
   CBranch(): m_uStabilizerCount(0) {
      m_uLevel = 0;
   }

   CBranch(BranchPointerType inParent): m_brParent(inParent), m_Line(), m_uStabilizerCount(0) {
      if( inParent != NULL ) {
         m_uLevel = inParent->GetLevel() + 1;
         inParent->AddDaughter(shared_from_this());
      } else {
         m_uLevel = 0;
      }
   }

public:
   //gets and sets
   unsigned GetLevel() const {return m_uLevel;}

   unsigned GetStabilizerCount() const {return m_uStabilizerCount;}
   void SetStabilizerCount(const unsigned& inStabilizerCount) {m_uStabilizerCount = inStabilizerCount;}

   const DaughtersHolderType& GetDaughters() const {return m_Daughters;}

   const CLine& GetLine() const {return m_Line;}

   BranchPointerType GetParent() const {return m_brParent;}
   void SetParent(BranchPointerType inParent) {
      m_brParent = inParent;
      if( inParent != NULL ) {
         inParent->AddDaughter(shared_from_this());
      } else {
         m_uLevel = 0;
      }
   }

   const CPoint& GetMinimum() const {return m_pntMinimum;}
   const CPoint& GetMaximum() const {return m_pntMaximum;}

public:
   //interface
    
    //Recursivesly scale the arbor by the input factors.
    //This can potentially create dangling parts 
    //   (daughters that don't touch their parents).
    void Scale(const double& x, const double& y, const double& z) {
        CPoint pntScale;
        pntScale.push_back(x);
        pntScale.push_back(y);
        pntScale.push_back(z);
        m_Line.Scale(pntScale);

        DaughtersHolderType::iterator iterD = m_Daughters.begin();
        DaughtersHolderType::iterator iterD_end = m_Daughters.end();
        for(; iterD != iterD_end; ++iterD) {
            (*iterD)->Scale(x,y,z);
        }
    }

   //Update the levels of this branch and all daughters.
   void UpdateLevels(int inLevel) {
      m_uLevel = inLevel;
      DaughtersHolderType::iterator iterD = m_Daughters.begin();
      DaughtersHolderType::iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->UpdateLevels(inLevel+1);
      }
   }

   //Return the product moment. Uses the gsl integrator and the public function
   // to get the moments.
   //Specialized for 2 or 3 dimensions.
   double ProductMoment(const int& inN) const {
      double dRet = 0.;
      if(GetLine().front().size() == 2) {
         DoProductMoment2D(inN, dRet);
      } else if(GetLine().front().size() == 3) {
         DoProductMoment3D(inN, dRet);
      } else {
         cerr << "Error in CBranch::ProductMoment: dimension must be 2 or 3...continuing\n" << flush;
      }
      return dRet;
   }

   //Remove the last daughter (used by CArbor::TrimDaughter)
   void RemoveLastDaughter() {
      m_Daughters.pop_back();
   }

   //Save the branch in a Neurolucida like format (swc).
   //Good enough to fool my NL reader, but not tested on NL itself.
   //inPoints has a list of the currently added points.
   //This is not terribly efficient, but shouldn't matter much.
   //  If it's a bottle neck, then make the search for a parent faster.
   void RecursiveDisplay(ostream& inOut) const {
      vector<CPoint> vecPoints;
      //Starts at the root (by definition)
      LineType::const_iterator iterL = GetLine().begin();
      LineType::const_iterator iterL_end = GetLine().end();
      vecPoints.push_back(*iterL);
      inOut << "1 1 " << *iterL << ((iterL->size() == 2)?" 1 1 -1\n":" 1 -1\n");
      ++iterL;
      for(; iterL != iterL_end; ++iterL) {
         vecPoints.push_back(*iterL);
         inOut << vecPoints.size() << " 3 "
			 << *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 ") 
            << vecPoints.size()-1 << "\n";
      }

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         //call the daughter version of the display
         (*iterD)->RecursiveDisplay(inOut, vecPoints);
      }
   }

   //Undivide the line.
   void UndivideLine() {
      m_Line.UnDivide();
   }

   //Trim this branch down to the last daughter.
   //i.e. send          to
   //      _________               ________
   //         /   \       =>           /   \
   //        /     \                  /     \
   //
   //Does nothing if there are no daughters.
   bool TrimToDaughter() {
      if( DaughterCount() > 0 ) {
         //cerr << "gh1 " << m_Line << "\n" << flush;
         //for(int i = 0; i < DaughterCount(); ++i) {
         //  cerr << "gh2 " << m_Daughters[i]->FirstPoint() << "\n" << flush;
         //}
         //The last daughter is the last one on the list of daughers.
         //maintain the number of line divisions, just in case.
         if(m_Line.Shrink(m_Daughters.back()->FirstPoint())) {
            //cout << "|" << flush;
         } else {
            cerr << "\nWarning in CBranch::TrimToDaughter: unable to trim.\n" << flush;
            cerr << m_Daughters.back()->FirstPoint() << " "
                  << m_Line << "\n" << flush;
            //cerr << m_Line.GetTurningPoints() << "\n" << flush;
            return false;
         }
      }
      return true;
   }

   //Recursively find the number of stabilizers at each level.
   //Modifies inL_to_S.
   void LevelStabilizers(map<int, vector<unsigned> >& inL_to_S) const {
      vector<unsigned> vecBlank;
      inL_to_S.insert( map<int, vector<unsigned> >::value_type(GetLevel(), vecBlank) ).
         first->second.push_back(GetStabilizerCount());
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->LevelStabilizers(inL_to_S);
      }
   }

   //Recursively find the number of stabilizers of the branches.
   //Modifies inS_to_N.
   void Stabilizers(map<int, unsigned>& inS_to_N) const {
      inS_to_N.insert( map<int, unsigned>::value_type(GetStabilizerCount(), 0) ).first->second += 1;
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->Stabilizers(inS_to_N);
      }
   }

   //Recursively find the angle formed between this branch and all daughters.
   //Modifies inHisto.
   void Angles(CHistogram<double>& inHisto) const {
      CPoint pntCurrent = m_Line.StartToEnd();
      double dCurrentLength = pntCurrent.Length();
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if(dCurrentLength > 1.e-12) {
            CPoint pntDaughter = (*iter)->GetLine().StartToEnd();
            double dDaughterLength = pntDaughter.Length();
            if(dDaughterLength > 1.e-12) {
               inHisto.increment(acos(pntCurrent.Dot(pntDaughter)/(dCurrentLength*dDaughterLength)));
            }
         }
         (*iter)->Angles(inHisto);
      }
   }

   //Recursively find the angle formed between this branch and all daughters.
   //Dumps angles to the input stream.
   //Changed to measure the angle locally (3/07 Joe Snider).
   void RawAngles(ostream& inOut) const {
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if( GetLine().size() > 1 ) {
            CPoint pntCurrent;
            CLine::const_iterator iterL = GetLine().begin();
            CLine::const_iterator iterL_end = GetLine().end();
            for(; iterL != iterL_end; ++iterL) {
               if( (*iterL - (*iter)->GetLine().at(0)).LengthSquared() < 1.e-12) {
                  pntCurrent = *iterL;
                  ++iterL;
                  if(iterL != iterL_end) {
                     pntCurrent = *iterL - pntCurrent;
                  } else {
                     --iterL;
                  }
               }
            }
            double dCurrentLength = pntCurrent.Length();
            if(dCurrentLength > 1.e-12 && (*iter)->GetLine().size() > 1) {
               CPoint pntDaughter = (*iter)->GetLine().at(1) - (*iter)->GetLine().at(0);
               double dDaughterLength = pntDaughter.Length();
               if(dDaughterLength > 1.e-12) {
                  inOut << acos(pntCurrent.Dot(pntDaughter)/(dCurrentLength*dDaughterLength)) << "\n";
               }
            }
         }
         (*iter)->RawAngles(inOut);
      }
   }

   //Recursively find the lengths of this branch and all daughters.
   //Changed (3/06) to define the length as the distance between branches.
   //This is what's commonly measured, so it is better for comparison to experiment.
   //Modifies inHisto.
   void Lengths(CHistogram<double>& inHisto) const {
      CPoint pntStart = FirstPoint();

      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if(!IsRoot()) {
            inHisto.increment( m_Line.DistanceToStart((*iter)->FirstPoint()) - m_Line.DistanceToStart(pntStart) );
            pntStart = (*iter)->FirstPoint();
         }
         (*iter)->Lengths(inHisto);
      }

      //include the end
      if(DaughterCount() == 0) {
         //this is a leaf just add the length of this branch
         inHisto.increment( m_Line.StartToEnd().Length() );
      } else if (!IsRoot()) {
         //add the last daughter to the end
         inHisto.increment( m_Line.DistanceToStart(LastPoint()) - m_Line.DistanceToStart(m_Daughters.back()->FirstPoint()) );
      }
   }

   //Recursively find the lengths of this branch and all daughters.
   //Dump them to the input stream.
   void RawLengths(ostream& inOut) const {
      CPoint pntStart = FirstPoint();

      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if(!IsRoot()) {
            inOut << m_Line.DistanceToStart((*iter)->FirstPoint()) - m_Line.DistanceToStart(pntStart) << "\n";
            pntStart = (*iter)->FirstPoint();
         }
         (*iter)->RawLengths(inOut);
      }

      //include the end
      if(DaughterCount() == 0) {
         //this is a leaf just add the length of this branch
         inOut << m_Line.Length() << "\n";
      } else if (!IsRoot()) {
         //add the last daughter to the end
         inOut << m_Line.DistanceToStart(LastPoint()) - m_Line.DistanceToStart(m_Daughters.back()->FirstPoint()) << "\n";
      }
   }

   //Find the max level (searches all daughters).
   int MaxLevel() const {
      int iLevel = GetLevel();
      DoMaxLevel(iLevel);
      return iLevel;
   }

   //Find the center of mass of this and all daughter branches.
   //Note: it's worth doing this alone, because many other measurements depend on it.
   CPoint TotalMean() const {
      unsigned uCount = 0;
      CPoint pCOM;
      pCOM.assign(m_pntMinimum.GetDimension(), 0.);
      DoTotalMean(pCOM, uCount);
      return pCOM*(1./double(uCount));
   }

   //--------Joe Snider, 6/06. Removed. This causes all sorts of problems...
   //Loop over the Daughters for access to the lines.
   ////Find all the segments in this and all daughter branches.
   ////Modifies inReturn.
   //void AllSegments(vector<CLine>& inReturn) const {
   //   inReturn.push_back(m_Line);
   //   DaughtersHolderType::const_iterator iter = m_Daughters.begin();
   //   DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
   //   for(; iter != iter_end; ++iter) {
   //      (*iter)->AllSegments(inReturn);
   //   }
   //}

   //Return the number of segments in the branch.
   int NumberSegments() const {
      return int(m_Line.size() - 1);
   }

   //Divide the current branch into inNum smaller branches.
   //Allows insertion along the branch rather than just at the end.
   void Divide(const int& inNum) {
      m_Line.Divide(inNum);
   }

   //Get the inN'th point of the line.
   //Returns the last point if the line is too short.
   CPoint NthPoint(const int& inN) const {
      if(int(m_Line.size()) > inN) {
         return m_Line[inN];
      }
      return LastPoint();
   }

   //Get the inN'th turning point of the line.
   //Returns the last point if the line is too short.
   CPoint NthTurningPoint(const int& inN) const {
      if(int(m_Line.GetTurningPoints().size()) > inN) {
         return m_Line.GetTurningPoints().at(inN);
      }
      return LastPoint();
   }

   //Get the first point on the line.
   //Returns an empty point if the line is too short.
   CPoint FirstPoint() const {
	   if(m_Line.GetTurningPoints().size() > 0) {
		  return m_Line.GetTurningPoints().at(0);
      }
      CPoint pRet;
      return pRet;
   }

   //Get the last point on the line (makes insertion easier)
   //Returns an empty point if the line is too short.
   CPoint LastPoint() const {
      if(m_Line.size() > 0) {
         return m_Line[m_Line.size()-1];
      }
      CPoint pRet;
      return pRet;
   }

   //Check if this is a root branch (no parent)
   bool IsRoot() const {
      return (m_uLevel == 0);
   }

   //Add a daughter
   //Sorts the daughters in the order they appear along the branch.
   void AddDaughter(BranchPointerType inDaughter) {
      inDaughter->m_brParent = shared_from_this();
      inDaughter->m_uLevel = m_uLevel+1;
      m_Daughters.push_back(inDaughter);

      //Do the sorting.
      multimap<double, BranchPointerType> mapSorter;
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         mapSorter.insert(
            multimap<double, BranchPointerType>::
            value_type(m_Line.DistanceToStart((*iter)->FirstPoint()), *iter) );
      }

      m_Daughters.clear();
      multimap<double, BranchPointerType>::const_iterator iterM = mapSorter.begin();
      multimap<double, BranchPointerType>::const_iterator iterM_end = mapSorter.end();
      for(; iterM != iterM_end; ++iterM) {
         m_Daughters.push_back(iterM->second);
      }
   }

   //Remove the daughter pointed at by the input pointer
   void RemoveDaughter(BranchPointerType inD) {
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         //note don't worry about the order, it can't be affected
         if( (*iter)->GetLine() == inD->GetLine() ) {
            m_Daughters.erase(iter);
            return;
         }
      }
      cerr << "Warning CBranch::RemoveDaughter: requested daughter not found...ignoring\n" << flush;
   }

   //Get the number of daughters.
   int DaughterCount() const {
      return (int)m_Daughters.size();
   }

   //Add a point to the line.
   void AddPoint(const CPoint& inPoint) {
      m_Line.push_back(inPoint);
      if(m_Line.size() == 1) {
         m_pntMinimum = inPoint;
         m_pntMaximum = inPoint;
      } else {
         for(int i = 0; i < inPoint.GetDimension(); ++i) {
            m_pntMinimum[i] = min(m_pntMinimum[i], inPoint[i]);
            m_pntMaximum[i] = max(m_pntMaximum[i], inPoint[i]);
         }
      }
      UpdateMinimumMaximum(m_pntMinimum, m_pntMaximum);
   }

   //Add a set of points to the line.
   //The type must satisfy the requirements for a forward stl iterator to CPoints.
   //This has the advantage of only updating the max/mins of the parent once all the new points have been added.
   template <class T>
   void AddPoint(const T& inBegin, const T& inEnd) {
      if(m_Line.size() == 0) {
         m_pntMinimum = *inBegin;
         m_pntMaximum = *inBegin;
      } else {
         for(T iter = inBegin; iter < inEnd; ++iter) {
            m_Line.push_back(*iter);
            for(int i = 0; i < inBegin->GetDimension(); ++i) {
               m_pntMinimum[i] = min(m_pntMinimum[i], iter->operator[](i));
               m_pntMaximum[i] = max(m_pntMaximum[i], iter->operator[](i));
            }
         }
      }
      UpdateMinimumMaximum(m_pntMinimum, m_pntMaximum);
   }

   //Recursively check this branch and all daughter branches for intersection with the input sphere,
   //  centered at inCenter with radius inRadius.
   //Note: the max/min values are used to possibly stop the search early.
   bool DoesSphereIntersect(const CPoint& inCenter, const double& inRadius) const {
      //first check that the input sphere is within the bounding box
      //(extend the bounding box by the radius and just check the point)
      for(int i = 0; i < inCenter.GetDimension(); ++i) {
         if( (inCenter[i]+inRadius < m_pntMinimum[i]) ||
            (inCenter[i] > inRadius+m_pntMaximum[i]) ) {
               return false;
         }
      }

      //next check if it intersects with the current line
      if( m_Line.DoesSphereIntersect(inCenter, inRadius) ) {
         return true;
      }

      //finally, recurse and check all the daughter branches
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if((*iter)->DoesSphereIntersect(inCenter, inRadius)) {
            return true;
         }
      }

      //everything failed, so no intersection
      //hopefully we won't get here very often
      return false;
   }

   //Recursively check this branch and all daughter branches for intersection in a cylinder
   //  with the input sphere, centered at inCenter with radius inRadius.
   //See Cline::DoesSphereHover for more detail, but this checks that inCenter is within
   //  inRadius of at least one segment and is directly above it.
   //Note: the max/min values are used to possibly stop the search early.
   bool DoesSphereHover(const CPoint& inCenter, const double& inRadius) const {
      //First check that the input sphere is within the bounding box.
      //(extend the bounding box by the radius and just check the point)
      for(int i = 0; i < inCenter.GetDimension(); ++i) {
         if( (inCenter[i]+inRadius < m_pntMinimum[i]) ||
            (inCenter[i] > inRadius+m_pntMaximum[i]) ) {
               return false;
         }
      }

      //Next check if it intersects with the current line.
      if( m_Line.DoesSphereHover(inCenter, inRadius) ) {
         return true;
      }

      //Finally, recurse and check all the daughter branches.
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if((*iter)->DoesSphereHover(inCenter, inRadius)) {
            return true;
         }
      }

      //Everything failed, so no intersection.
      return false;
   }

   //Clear the branch.
   //Note, smart pointerification makes this safe.
   void Clear() {
      m_brParent.reset();
      m_uLevel = 0;
      m_Daughters.clear();
      m_Line.clear();
   }

private:
   //helpers

   //Do the product moment recursion in 2d.
   void DoProductMoment2D(const int& inN, double& inRet) const {
      if(GetLine().size() > 1) {
         //Create some things needed for the gsl integrator.
         gsl_integration_workspace* w =
            gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

         double dResult, dError;
         double dAlpha[5];
         dAlpha[0] = (double)inN;

         CLine::const_iterator iterL = m_Line.begin();
         CLine::const_iterator iterR = m_Line.begin();
         ++iterR;
         CLine::const_iterator iter_end = m_Line.end();
         for(; iterR != iter_end; ++iterL, ++iterR) {
            dAlpha[1] = iterL->at(0);
            dAlpha[2] = iterR->at(0);
            dAlpha[3] = iterL->at(1);
            dAlpha[4] = iterR->at(1);

            gsl_function LF_2D;
            LF_2D.function = &LineFunction2D;
            LF_2D.params = &dAlpha;
            gsl_integration_qag(&LF_2D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
            inRet += (*iterR-*iterL).Length()*dResult;
         }

         //clean up the workspace
         gsl_integration_workspace_free(w);
      }

      //Do the recursion on the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoProductMoment2D(inN, inRet);
      }
   }

   //Do the product moment recursion in 3d.
   void DoProductMoment3D(const int& inN, double& inRet) const {
      if(GetLine().size() > 1) {
         //Create some things needed for the gsl integrator.
         gsl_integration_workspace* w =
            gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

         double dResult, dError;
         double dAlpha[7];
         dAlpha[0] = (double)inN;

         CLine::const_iterator iterL = m_Line.begin();
         CLine::const_iterator iterR = m_Line.begin();
         ++iterR;
         CLine::const_iterator iter_end = m_Line.end();
         for(; iterR != iter_end; ++iterL, ++iterR) {
            dAlpha[1] = iterL->at(0);
            dAlpha[2] = iterR->at(0);
            dAlpha[3] = iterL->at(1);
            dAlpha[4] = iterR->at(1);
            dAlpha[5] = iterL->at(2);
            dAlpha[6] = iterR->at(2);

            gsl_function LF_3D;
            LF_3D.function = &LineFunction3D;
            LF_3D.params = &dAlpha;

            gsl_integration_qag(&LF_3D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
            inRet += (*iterR-*iterL).Length()*dResult;
         }

         //clean up the workspace
         gsl_integration_workspace_free(w);
      }

      //Do the recursion on the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoProductMoment3D(inN, inRet);
      }
   }

   //Do the total mean recursion.
   //Modifies inPoint (sum of all points) and inCount (=number of points).
   void DoTotalMean(CPoint& inPoint, unsigned& inCount) const {
      CLine::const_iterator iterL = m_Line.begin();
      CLine::const_iterator iterL_end = m_Line.end();
      for(; iterL != iterL_end; ++iterL) {
         inPoint += *iterL;
         ++inCount;
      }
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoTotalMean(inPoint, inCount);
      }
   }

   //Do the calculation to find the max level
   void DoMaxLevel(int& inLevel) const {
      inLevel = max((int)inLevel, (int)GetLevel());
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->DoMaxLevel(inLevel);
      }
   }

   //Recursively update the max and min of the branch.
   //The min (max) is guarunteed to cover all daughter branches, so it can be used for optimization later.
   void UpdateMinimumMaximum(const CPoint& inMinimum, const CPoint& inMaximum) {
      for(int i = 0; i < inMinimum.GetDimension(); ++i) {
         m_pntMinimum[i] = min(m_pntMinimum[i], inMinimum[i]);
         m_pntMaximum[i] = max(m_pntMaximum[i], inMaximum[i]);
      }

      //update the parent (unless this is the root)
      if(!IsRoot()) {
         m_brParent->UpdateMinimumMaximum(m_pntMinimum, m_pntMaximum);
      }
   }

   //Do recursion
   void RecursiveDisplay(ostream& inOut, vector<CPoint>& inPoints) const {
      //cout << "branch.h gh1 " << GetLine().size() << "\n" << flush;
      //it is possible to have a blank daughter
      if(GetLine().size() > 1) {
         //find the parent id (-1 if it's not found)
         int iParent = -1;
         vector<CPoint>::const_iterator iterP 
            = find(inPoints.begin(), inPoints.end(), GetLine().at(0));
         if(iterP != inPoints.end()) {
            iParent = (int)distance((vector<CPoint>::const_iterator)inPoints.begin(), iterP) + 1;
         } else {
            cerr << "Failed to find a parent in CBranch::RecursiveDisplay. This is probably a bug...continuing\n" << flush;
            //////////////////////////////testing only/////////////////////////////////////////
            cerr << "CBranch.h gh1 " << m_uLevel << " " << GetLine().at(0) << "\n" << flush;
            double dMinDistTest = 1.e12;
            int iClosestID = -1;
            for(int l = 0; l < inPoints.size(); ++l) {
               if(dMinDistTest > (inPoints[l]-GetLine().at(0)).LengthSquared()) {
                  dMinDistTest = (inPoints[l]-GetLine().at(0)).LengthSquared();
                  iClosestID = l;
               }
            }
            cerr << "CBranch.h gh2 " << inPoints[iClosestID] << "\n" << flush;
            /////////////////////////////end testing only//////////////////////////////////////
         }

         LineType::const_iterator iterL = GetLine().begin();
         LineType::const_iterator iterL_end = GetLine().end();
         //the first point on the line is a duplicate
         ++iterL;

         //the second point has a special parent
         inPoints.push_back(*iterL);
         inOut << inPoints.size() << " 3 "
			 << *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 " )
            << iParent << "\n";
         ++iterL;

         //the rest of the points go in order
         for(; iterL != iterL_end; ++iterL) {
            inPoints.push_back(*iterL);
            inOut << inPoints.size() << " 3 "
				<< *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 ") 
               << inPoints.size()-1 << "\n";
         }
      }

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         //call the daughter version of the display
         (*iterD)->RecursiveDisplay(inOut, inPoints);
      }
   }

private:
   //data members
   unsigned m_uLevel;
   BranchPointerType m_brParent;
   DaughtersHolderType m_Daughters;
   CPoint m_pntMinimum;
   CPoint m_pntMaximum;
   CLine m_Line;
   unsigned m_uStabilizerCount;
};

#endif //BRANCH
