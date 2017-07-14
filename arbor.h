//Joe Snider
//1/05
//
//Simulate a neuronal arbor. A sort of tree containing branches which can contain daughters.
//
//TODO: this should really be derived from CBranch. That is a serious change, though.

#include <iostream>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <time.h>
#include <boost/shared_ptr.hpp>

#include "branch.h"
#include "puncta.h"
#include "KDTree.h"
#include "histogram.h"
#include "tnt.h"
#include "jama_eig.h"

#ifndef ARBOR
#define ARBOR

const double PUNCTA_REMOVE_PROB = 1.;
const double PUNCTA_REMOVE_RADIUS_SIGMA = 0.1;

#define HOVER  //Set whether to use hovering puncta or all withing a range (DoesSphereIntersect or DoesSphereHover).

const int MAX_DAUGHTERS = 1;
const int MAX_DAUGHTERS_ROOT = 3;
//const int SEGMENTS_PER_BRANCH = 5; set at stabilizers now

const double POISSON_TEST_VALUE = 0.05;

using namespace std;
using namespace boost;
using namespace TNT;
using namespace JAMA;

class CArbor {
public:
   //typedefs
   typedef shared_ptr<CBranch> BranchPointerType;
   typedef vector<BranchPointerType> DaughterHolderType;
   typedef list<CBranch*> AvailableBranchType; //This is not allowed to delete, so just use the pointer.

public:
   //constructors
   CArbor(): m_dProbExtendLinearly(0.5) {}

   ~CArbor() {
      //Safe pointers, so nothing to do.
   }

public:
   //gets and sets

   const CBranch& GetBranch() const {return *m_brBranch;}
   
   const double GetProbExtendLinearly() const {return m_dProbExtendLinearly;}
   void SetProbExtendLinearly(const double& inPEL) {m_dProbExtendLinearly = inPEL;}

   const AvailableBranchType& GetAvailableBranches() const {return m_listAvailableBranches;}

   const DaughterHolderType& GetAllBranches() const {return m_vecAllBranches;}

public:
   //interface
    
    //scale the arbor
    void Scale(const double& x, const double& y, const double& z) {
        m_brBranch->Scale(x,y,z);
    }
   
   //Return the branch with the max level
   int MaxLevel() const {
      return m_brBranch->MaxLevel();
   }

   //Return the product moment.
   //Just passes the calculation on to the branch.
   double ProductMoment(const int& inN) const {
      return m_brBranch->ProductMoment(inN);
   }

   //Write the data in a form that NeuroLeucidia can read.
   void Write(ostream& inOut) const {
      //Put some identifiers at the top
      time_t timeCurrent;
      time( &timeCurrent );
      //char cstrBuffer[128];
      //ctime_s(cstrBuffer, 26, timeCurrent);
      inOut << "# Created by NeuronBuilder2\n";
      inOut << "# Date: " << ctime(&timeCurrent);
      m_brBranch->RecursiveDisplay(inOut);
   }

   //Return the starting point.
   //Return a blank point if there is no starting point.
   CPoint Start() const {
      if(m_brBranch->GetLine().size() > 0) {
         return m_brBranch->GetLine().front();
      }

      //empty arbor
      CPoint pRet;
      return pRet;
   }

   //Update the level structure since it can be broken in various places.
   void UpdateLevels() {
      m_brBranch->UpdateLevels(m_brBranch->GetLevel());
   }

   //Count the number of stabilizers at each level.
   //Modifies inL_to_S.
   void LevelStabilizers(map<int, vector<unsigned> >& inL_to_S) const {
      m_brBranch->LevelStabilizers(inL_to_S);
   }

   //Count the stabilizers of this arbor.
   //Modifies inS_to_N.
   void Stabilizers(map<int, unsigned>& inS_to_N) const {
      m_brBranch->Stabilizers(inS_to_N);
   }

   //Find the angles formed by this branch with the daughter branches.
   //Modifies inHisto.
   void Angles(CHistogram<double>& inHisto) const {
      m_brBranch->Angles(inHisto);
   }

   //Dump the raw lengths to the input stream.
   void RawAngles(ostream& inOut) const {
      m_brBranch->RawAngles(inOut);
   }

   //Find the lengths of this branch.
   //Modifies inHisto.
   void Lengths(CHistogram<double>& inHisto) const {
      m_brBranch->Lengths(inHisto);
   }

   //Dump the raw lengths to the input stream.
   void RawLengths(ostream& inOut) const {
      m_brBranch->RawLengths(inOut);
   }

   //Calculate the probability density funciton of the rescaled neuron.
   void InertialHistogram(CHistogram<double>& inMass, CHistogram<double>& inCount) const {
      /*//get all the segments
      vector<CLine> vecSegments;
      m_brBranch->AllSegments(vecSegments);

      //find the center of mass to use as the origin (shift theorm may be faster, but I'd be more likely to screw it up)
      CPoint pntCOM = m_brBranch->TotalMean();

      //check for too large of a dimension
      int iDimension = pntCOM.GetDimension();
      if( (iDimension != 2) && (iDimension != 3) ) {
         return;
      }

      //Shift the segment copies.
      vector<CLine>::iterator iter = vecSegments.begin();
      vector<CLine>::iterator iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         iter->Shift(-1.*pntCOM);
      }

      //Create and initialize the moment of inertia tensor.
      //Uses TNT library from NIST.
      Array2D<double> a2dInertiaTensor(3,3);
      a2dInertiaTensor = 0.0;

      //do the sum to find the moment of inertia
      //note: some bad variable names (x, y, z are the position of the point, m its mass (length))
      double m, x, y, z;
      CPoint pV;
      double dM = 0.;
      iter = vecSegments.begin();
      iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         if(iter->size() > 1) {
            CLine::const_iterator iterLL = iter->begin();
            CLine::const_iterator iterLR = iter->begin();
            ++iterLR;
            CLine::const_iterator iterL_end = iter->end();
            for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
               //minimize copying (the compiler should help significantly)
               pV = *iterLL;
               pV += *iterLR;
               pV *= 0.5;
               x = pV[0];
               y = pV[1];
               z = (iDimension==2)?0.:pV[2];
               m = (*iterLL - *iterLR).Length();
               dM += m;
               a2dInertiaTensor[0][0] += m*(y*y+z*z);
               a2dInertiaTensor[0][1] -= m*x*y;
               a2dInertiaTensor[0][2] -= m*x*z;
               a2dInertiaTensor[1][1] += m*(x*x+z*z);
               a2dInertiaTensor[1][0] -= m*y*x;
               a2dInertiaTensor[1][2] -= m*y*z;
               a2dInertiaTensor[2][2] += m*(x*x+y*y);
               a2dInertiaTensor[2][0] -= m*z*x;
               a2dInertiaTensor[2][1] -= m*z*y;
            }
         }
      }

      //Do the eigenvalue problem (note: moment of inertia is hermitian).
      //Not the best variable names, but there'll be alot of V[1][2]*x+..., so it's best.
      Eigenvalue<double> E1(a2dInertiaTensor);
      Array2D<double> D(3,3);
      Array2D<double> V(3,3);
      E1.getD(D);
      E1.getV(V);

      //Rotate the segment copies, and find the max/min for scaling.
      double dXMax = -1.e100;
      double dYMax = -1.e100;
      double dZMax = -1.e100;
      double dXMin = 1.e100;
      double dYMin = 1.e100;
      double dZMin = 1.e100;
      iter = vecSegments.begin();
      iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         iter->Rotate(V);
         dXMax = max( dXMax, iter->Max(0) );
         dYMax = max( dYMax, iter->Max(1) );
         dZMax = (iDimension==2) ? 1. : max( dZMax, iter->Max(2) );
         dXMin = min( dXMin, iter->Min(0) );
         dYMin = min( dYMin, iter->Min(1) );
         dZMin = (iDimension==2) ? 0. : min( dZMin, iter->Min(2) );
      }

      //Scale the segments (don't have to worry about the dimension).
      CPoint pScale;
      pScale.push_back( 2./(dXMax-dXMin) );
      pScale.push_back( 2./(dYMax-dYMin) );
      pScale.push_back( 2./(dZMax-dZMin) );
      iter = vecSegments.begin();
      iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         iter->Scale(pScale);
      }

      //Do the incrementing.
      iter = vecSegments.begin();
      iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         if(iter->size() > 1) {
            CLine::const_iterator iterLL = iter->begin();
            CLine::const_iterator iterLR = iter->begin();
            ++iterLR;
            CLine::const_iterator iterL_end = iter->end();
            for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
               pV = *iterLL;
               pV += *iterLR;
               pV *= 0.5;
               x = pV[0];
               y = pV[1];
               z = (iDimension==2)?0.:pV[2];
               m = (*iterLL - *iterLR).Length();
               inMass.increment( sqrt(x*x+y*y+z*z), m );
               inCount.increment( sqrt(x*x+y*y+z*z), 1. );
            }
         }
      }*/
   }

   //Return the inertial volume of the neuron.
   //Other results have shown that this is linearly related to the box-counting volume, and it's better defined.
   //Note: this can only be defined (at least by me) in 2 or 3 dimensions, and any others return 0.
   //This calculates the total length (it's called the total mass dM here), so it may as well be returned too.
   //Modifies inInertialVolume and inTotalLength
   void InertialVolume(double& inInertialVolume, double& inTotalLength) const {
      //find the center of mass to use as the origin (shift theorm may be faster, but I'd be more likely to screw it up)
      const CPoint pntCOM = m_brBranch->TotalMean();

      //check for too large of a dimension
      if( (pntCOM.GetDimension() != 2) && (pntCOM.GetDimension() != 3) ) {
         inInertialVolume = 0;
         inTotalLength = 0;
         return;
      }

      //create and initialize the moment of inertia tensor (just a vector<vector<double> >)
      vector<vector<double> > vecInertia;
      for(int i = 0; i < 3; ++i) {
         vector<double> vecInsert(3, 0.);
         vecInertia.push_back(vecInsert);
      }

      //do the sum to find the moment of inertia
      //note: some bad variable names (x, y, z are the position of the point, m its mass (length))
      //double m, x, y, z;
      CPoint V;
      double dM = 0.;
      DaughterHolderType::const_iterator iter = m_vecAllBranches.begin();
      DaughterHolderType::const_iterator iter_end = m_vecAllBranches.end();
      if(pntCOM.GetDimension() == 2) {
         for(; iter != iter_end; ++iter) {
            //z = 0 by definition
            if( (*iter)->GetLine().size() > 1) {
               CLine::const_iterator iterLL = (*iter)->GetLine().begin();
               CLine::const_iterator iterLR = (*iter)->GetLine().begin();
               ++iterLR;
               CLine::const_iterator iterL_end = (*iter)->GetLine().end();
               for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
                  //minimize copying (the compiler should help significantly)
                  V = *iterLL;
                  V += *iterLR;
                  V *= 0.5;
                  V -= pntCOM;
                  const double x = V[0];
                  const double y = V[1];
                  const double m = (*iterLL - *iterLR).Length();
                  dM += m;
                  vecInertia[0][0] += m*(y*y);
                  vecInertia[0][1] -= m*x*y;
                  vecInertia[1][1] += m*(x*x);
                  vecInertia[1][0] -= m*y*x;
               }
            }
         }
         inInertialVolume = sqrt((vecInertia[0][0]*vecInertia[1][1]-vecInertia[0][1]*vecInertia[1][0])/(dM*dM));
         inTotalLength = dM;
      } else {
         //z is non-zero
         for(; iter != iter_end; ++iter) {
            //z = 0 by definition
            if( (*iter)->GetLine().size() > 1) {
               CLine::const_iterator iterLL = (*iter)->GetLine().begin();
               CLine::const_iterator iterLR = (*iter)->GetLine().begin();
               ++iterLR;
               CLine::const_iterator iterL_end = (*iter)->GetLine().end();
               for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
                  //minimize copying (the compiler should help significantly)
                  V = *iterLL;
                  V += *iterLR;
                  V *= 0.5;
                  V -= pntCOM;
                  const double x = V[0];
                  const double y = V[1];
                  const double z = V[2];
                  const double m = (*iterLL - *iterLR).Length();
                  dM += m;
                  vecInertia[0][0] += m*(y*y+z*z);
                  vecInertia[0][1] -= m*x*y;
                  vecInertia[0][2] -= m*x*z;
                  vecInertia[1][1] += m*(x*x+z*z);
                  vecInertia[1][0] -= m*y*x;
                  vecInertia[1][2] -= m*y*z;
                  vecInertia[2][2] += m*(x*x+y*y);
                  vecInertia[2][0] -= m*z*x;
                  vecInertia[2][1] -= m*z*y;
               }
            }
         }
         inInertialVolume = sqrt(Determinant(vecInertia)/(dM*dM*dM));
         inTotalLength = dM;
      }
   }

   //Set the starting point.
   //Also clears the branch and deletes any daughters.
   void CreateStart(const CPoint& inCenter) {
      //reset
      m_vecAllBranches.clear();
      BranchPointerType brTemp( new CBranch() );
      m_brBranch = brTemp;
      m_listAvailableBranches.clear();

      //add the input point to the main branch, and the list of available branches
      m_brBranch->AddPoint(inCenter);
      m_listAvailableBranches.push_back(m_brBranch.get());
      m_vecAllBranches.push_back(m_brBranch);
   }

   //Do a single step of neuron growth using the input puncta.
   //May modify the occupation flag of the inPuncta.
   //The templated type T is the random number generator ( operator() returns a double in [0,1) )
   //Returns true if an insertion occured.
   template<class T>
   bool Step(const double& inLength, const double& inPunctaRange, const double& inPunctaRemoveRadius, T& inRand, CKDTree<CPuncta>& inPuncta) {
      //Used by the kd-tree (I really don't remember why I made the search return a deque, but it shouldn't hurt anything).
      deque<CPuncta*> dqPunctaInRange;

      if(m_listAvailableBranches.size() > 0) {
         //pick a random end-point at which to insert the new branch
         //TODO: need a better way to get a random element
         int iID = int(double(m_listAvailableBranches.size())*inRand());
         AvailableBranchType::iterator iterID = m_listAvailableBranches.begin();
         for(int iSkip = 0; iSkip < iID; ++iSkip, ++iterID);

         //First try to extend the branch linearly.
         int iNumSegments = (*iterID)->NumberSegments();
         if(inRand() < m_dProbExtendLinearly) {
            //Have to make sure this isn't the starting point.
            if( iNumSegments > 0 ) {
               CPoint pExtendL = (*iterID)->LastPoint();
               CPoint pExtendR = pExtendL;
               pExtendR -= (*iterID)->NthPoint((int)(*iterID)->GetLine().size()-2);
               pExtendR *= inRand() * inLength / pExtendR.Length();
               pExtendR += pExtendL;
               FindStabilizers(pExtendL, pExtendR, inPunctaRange, inPuncta, dqPunctaInRange);
               if(dqPunctaInRange.size() > 0) {
                  //Accept the new branch.
                  //Add turning points at the stabilizers (projected onto the new branch).
                  //Use a temporary line to find the correct points
                  CLine lineTemp;
                  lineTemp.push_back(pExtendL);
                  lineTemp.push_back(pExtendR);

                  //Maintain order with a map
                  multimap<double, CPoint> mapNewTurningPoints;
                  deque<CPuncta*>::const_iterator iter = dqPunctaInRange.begin();
                  deque<CPuncta*>::const_iterator iter_end = dqPunctaInRange.end();
                  for(; iter != iter_end; ++iter) {
                     //Find the new turning point.
                     const CPoint pTemp = lineTemp.FindClosestPoint(**iter);
                     const double dDist = (pTemp-pExtendL).LengthSquared();
                     //stl magic
                     mapNewTurningPoints.insert(multimap<double, CPoint>::value_type(dDist, pTemp));
                  }
                  
                  //Insert the new turning points
                  multimap<double, CPoint>::const_iterator iterTP = mapNewTurningPoints.begin();
                  multimap<double, CPoint>::const_iterator iterTP_end = mapNewTurningPoints.end();
                  for(; iterTP != iterTP_end; ++iterTP) {
                     (*iterID)->AddPoint(iterTP->second);
                  }

                  //Remove puncta near the stabilizers.
                  //TrimRandomPuncta(inPunctaRemoveRadius, dqPunctaInRange, PUNCTA_REMOVE_PROB, inRand, inPuncta);
                  TrimRandomPuncta(inPunctaRemoveRadius/HackRR(pExtendR.Length()), dqPunctaInRange, PUNCTA_REMOVE_PROB, inRand, inPuncta);
                  (*iterID)->SetStabilizerCount( (*iterID)->GetStabilizerCount() + dqPunctaInRange.size() );
                  return true;
               }
            }
         }

         //Couldn't extend linearly, so try to create a new branch randomly along the current branch.
         //CPoint pTrialL = (*iterID)->NthTurningPoint(inRand(iNumSegments));
         //int iGaussianTrial = 0;
         //if(iNumSegments > 0) {
         //   do{
         //      iGaussianTrial = int(0.5*double(iNumSegments) + inRand.Gaussian( 0.5*double(iNumSegments) ) + 0.5);
         //   } while(iGaussianTrial < 0 || iGaussianTrial >= iNumSegments);
         //}
         //CPoint pTrialL = (*iterID)->NthTurningPoint(iGaussianTrial);
         CPoint pTrialL = (*iterID)->NthTurningPoint(iNumSegments-int(double(inRand(iNumSegments)+inRand(iNumSegments)+inRand(iNumSegments))/3.+.5));
         /*int iTest = 0;
         for(int k = 0; inRand() > POISSON_TEST_VALUE && k < iNumSegments; ++k) {
            iTest = k;
         }
         CPoint pTrialL = (*iterID)->NthTurningPoint(iTest);*/
         double dLengthTest;
         CPoint pTrialR;
         do {
            pTrialR = pTrialL;
            for(int i = 0; i < pTrialR.GetDimension(); ++i) {
               //pTrialR[i] += inLength*(2.*inRand()-1.);
               pTrialR[i] += inLength*inRand.Gaussian(1./2.);
            }
            dLengthTest = (pTrialR-pTrialL).LengthSquared();
         } while ( dLengthTest > 4.*inLength*inLength );//&& dLengthTest < 0.25*inLength*inLength);
         //pTrialR -= pTrialL;
         //pTrialR *= inLength / pTrialR.Length();
         //pTrialR += pTrialL;
         FindStabilizers(pTrialL, pTrialR, inPunctaRange, inPuncta, dqPunctaInRange);

         //Add the current branch if it was accepted.
         if(dqPunctaInRange.size() > 0) {
            BranchPointerType brTemp(new CBranch());
            //Add turning points at the stabilizers (projected onto the new branch).
            //Use a temporary line to find the correct points
            CLine lineTemp;
            lineTemp.push_back(pTrialL);
            lineTemp.push_back(pTrialR);

            //Maintain order with a map
            multimap<double, CPoint> mapNewTurningPoints;
            deque<CPuncta*>::const_iterator iter = dqPunctaInRange.begin();
            deque<CPuncta*>::const_iterator iter_end = dqPunctaInRange.end();
            for(; iter != iter_end; ++iter) {
               //Find the new turning point.
               const CPoint pTemp = lineTemp.FindClosestPoint(**iter);
               const double dDist = (pTemp-pTrialL).LengthSquared();
               //stl magic
               mapNewTurningPoints.insert(multimap<double, CPoint>::value_type(dDist, pTemp));
            }

            //Insert the new turning points
            brTemp->AddPoint(pTrialL);
            multimap<double, CPoint>::const_iterator iterTP = mapNewTurningPoints.begin();
            multimap<double, CPoint>::const_iterator iterTP_end = mapNewTurningPoints.end();
            for(; iterTP != iterTP_end; ++iterTP) {
               brTemp->AddPoint(iterTP->second);
            }

            brTemp->SetStabilizerCount( (int)dqPunctaInRange.size() );
            m_listAvailableBranches.push_back(brTemp.get());
            m_vecAllBranches.push_back(brTemp);
            (*iterID)->AddDaughter(brTemp);

            if( !(*iterID)->IsRoot() ) {
               if( (*iterID)->DaughterCount() > MaxDaughters(*iterID) ) {
                  m_listAvailableBranches.erase(iterID);
               }
            }

            //TrimRandomPuncta(inPunctaRemoveRadius, dqPunctaInRange, PUNCTA_REMOVE_PROB, inRand, inPuncta);
            TrimRandomPuncta(inPunctaRemoveRadius/HackRR(pTrialR.Length()), dqPunctaInRange, PUNCTA_REMOVE_PROB, inRand, inPuncta);
            return true;
         }
      }

      //Couldn't add any of the trials.
      return false;
   }

   //Trim the ends of this arbor such that only ends with enough stabilizers survive.
   //Returns the number of branches modified.
   int TrimEnds(const double& inPunctaRange, const CKDTree<CPuncta>& inPuncta, const int& inMinStabilizers = 1) {
      int iCount = 0;
      DaughterHolderType::iterator iter = m_vecAllBranches.begin();
      DaughterHolderType::iterator iter_end = m_vecAllBranches.end();
      for(; iter != iter_end; ++iter) {
         int iNumStabilizers;
         if( (*iter)->DaughterCount() > 0 ) {
            //A branch with a potentially small nub at the end
            CPoint pStart = (*iter)->GetDaughters().back()->FirstPoint();
            CPoint pEnd = (*iter)->LastPoint();
            if( (pEnd-pStart).LengthSquared() < 1.e-12 ) {
               //skip if the branch ends at a daughter
               iNumStabilizers = inMinStabilizers + 1;
            } else {
               iNumStabilizers = CountStabilizers(pStart, pEnd, inPunctaRange, inPuncta);
            }
         } else {
            //A leaf that may be small
            iNumStabilizers = (*iter)->GetStabilizerCount();
         }
         //Find the stabilizers of the end.
         if(iNumStabilizers < inMinStabilizers) {
            //trim the branch.
            ++iCount;
            if( (*iter)->DaughterCount() > 0 ) {
               //just make a line shorter
               if(! ((*iter)->TrimToDaughter()) ) {
                  cerr << "Warning in CArbor::TrimEnds: daughter trim failed...stopping trim\n" << flush;
                  return -1;
               }
            } else {
               //remove an entire leaf (unless it's the root)
               //(have to restart the search because branches change)
               if( (*iter)->IsRoot() ) {
                  cerr << "Warning in CArbor::TrimEnds: attempted to remove the root \
                        (arbor is probably too sparse)...stopping trim\n" << flush;
                  return -2;
               } else {
                  (*iter)->GetParent()->RemoveDaughter(*iter);
                  m_vecAllBranches.erase(iter);
                  //cout << ";" << flush;
                  iter = m_vecAllBranches.begin();
                  iter_end = m_vecAllBranches.end();
               }
            }
         }
      }
      return iCount;
   }

   //Merge daughters uphill.
   //Takes branches that would have ended at a single daughter and turns them into one branch.
   //This is a memory leak for now. I think the daughter could just be deleted, but it's dangerous.
   //Pointers are smart now, so no memory leak (Joe 3/07).
   //Be sure to update the levels elsewhere.
   bool TrimDaughter() {
      DaughterHolderType::iterator iter = m_vecAllBranches.begin();
      DaughterHolderType::iterator iter_end = m_vecAllBranches.end();
      for(; iter != iter_end; ++iter) {
         if( (*iter)->DaughterCount() > 0 && (*iter)->GetLine().size() > 1) {
            if( (*iter)->LastPoint() == (*iter)->GetDaughters().back()->FirstPoint() ) {
               //Find the daughter branch in the list of all branches.
               DaughterHolderType::iterator iterDaughter = m_vecAllBranches.end();
               DaughterHolderType::iterator iterD = m_vecAllBranches.begin();
               DaughterHolderType::iterator iterD_end = m_vecAllBranches.end();
               while( iterDaughter == iterD_end ) {
                  if(iterD == iterD_end) {
                     cerr << "Error: unknown daughter error in CArbor::TrimDaughter...stopping trim\n" << flush;
                     return false;
                  }
                  if( (*iterD)->GetLine() == (*iter)->GetDaughters().back()->GetLine() ) {
                     iterDaughter = iterD;
                  }
                  ++iterD;
               }

               //Remove the daughter from the list of available branches (if it's there).
               //Add this branch to the available branches (if it's not there).

               //Record the number of divisions, just in case, then undivide.
               int iDivisions = int((*iter)->GetLine().size() + (*iterDaughter)->GetLine().size());
               (*iter)->UndivideLine();
               (*iterDaughter)->UndivideLine();

               //Merge the daughter line onto this branch's line.
               for(unsigned i = 1; i < (*iterDaughter)->GetLine().GetTurningPoints().size(); ++i) {
                  (*iter)->AddPoint( (*iterDaughter)->GetLine().GetTurningPoints().at(i) );
               }

               //Remove the daughter from the parent (this)
               (*iter)->RemoveLastDaughter();

               //Add the daughter's daughters to this branch.
               for(unsigned j = 0; j < (*iterDaughter)->GetDaughters().size(); ++j) {
                  (*iter)->AddDaughter( (*iterDaughter)->GetDaughters().at(j) );
               }

               //Remove the daughter from the list of all branches.
               m_vecAllBranches.erase(iterDaughter);

               //Have to restart the search.
               iter = m_vecAllBranches.begin();
               iter_end = m_vecAllBranches.end();

               //cout << ">" << flush;
            }
         }
      }
      return true;
   }

   //Check if it is possible to extend any of the available endpoints by the input length to the input kd-tree.
   bool Extendable(CKDTree<CPuncta>& inPuncta, const double& inRange) const {
      deque<CPuncta*> dqPunctaInRange;
      //cout << "-" << flush;
      int iPossible = 0;
      AvailableBranchType::const_iterator iter = m_listAvailableBranches.begin();
      AvailableBranchType::const_iterator iter_end = m_listAvailableBranches.end();
      for(; (iter != iter_end) && (iPossible == 0); ++iter) {
         dqPunctaInRange.clear();
         CPuncta pcTemp( (*iter)->LastPoint() );
         inPuncta.FindPointsInRange( pcTemp, inRange, dqPunctaInRange );
         deque<CPuncta*>::iterator iterPIR = dqPunctaInRange.begin();
         deque<CPuncta*>::iterator iterPIR_end = dqPunctaInRange.end();
         for(; iterPIR != iterPIR_end; ++iterPIR) {
            if(!(*iterPIR)->IsOccupied()) {
               return true;
            }
         }
      }
      return false;
   }

private:
   //helpers

   //Return the maximum allowed daughters.
   int MaxDaughters(const CBranch* inBranch) const {
      if(m_brBranch->IsRoot()){
         return MAX_DAUGHTERS_ROOT;
      }
      return MAX_DAUGHTERS;
   }

   //Mark any puncta in range of the input puncta as occupied.
   //May modify the occupied flags of the input kdTree.
   void TrimPuncta(const double& inPunctaRemoveRadius, const deque<CPuncta*>& inStabilized, CKDTree<CPuncta>& inPuncta) const {
      deque<CPuncta*> dqRemove;
      deque<CPuncta*>::const_iterator iter = inStabilized.begin();
      deque<CPuncta*>::const_iterator iter_end = inStabilized.end();
      for(; iter != iter_end; ++iter) {
         dqRemove.clear();
         inPuncta.FindPointsInRange( **iter, inPunctaRemoveRadius, dqRemove );
         deque<CPuncta*>::iterator iterDoRemove = dqRemove.begin();
         deque<CPuncta*>::iterator iterDoRemove_end = dqRemove.end();
         for(; iterDoRemove != iterDoRemove_end; ++iterDoRemove) {
            (*iterDoRemove)->SetOccupied(true);
         }
      }
   }

   //Mark some puncta in range of the input puncta as occupied.
   //May modify the occupied flags of the input kdTree.
   //The type T is a random number generator.
   template<class T>
   void TrimRandomPuncta(const double& inPunctaRemoveRadius, 
      const deque<CPuncta*>& inStabilized, 
      const double& inRemoveProbability,
      T& inRand,
      CKDTree<CPuncta>& inPuncta) const {
      deque<CPuncta*> dqRemove;
      deque<CPuncta*>::const_iterator iter = inStabilized.begin();
      deque<CPuncta*>::const_iterator iter_end = inStabilized.end();
      for(; iter != iter_end; ++iter) {
         //const double dRemoveRadius = max(0., inPunctaRemoveRadius + inRand.Gaussian(PUNCTA_REMOVE_RADIUS_SIGMA));
         const double dRemoveRadius = inPunctaRemoveRadius;
         dqRemove.clear();
         inPuncta.FindPointsInRange( **iter, dRemoveRadius, dqRemove );
         deque<CPuncta*>::iterator iterDoRemove = dqRemove.begin();
         deque<CPuncta*>::iterator iterDoRemove_end = dqRemove.end();
         for(; iterDoRemove != iterDoRemove_end; ++iterDoRemove) {
            if(inRand() < inRemoveProbability) {
               (*iterDoRemove)->SetOccupied(true);
            }
         }
      }
   }

   //Find the puncta that stabilize the input branch (if any).
   //Modifies the input deque, and may change the occupied flags of the input kdTree.
   void FindStabilizers(const CPoint& inLeft, 
                        const CPoint& inRight, 
                        const double& inPunctaRange, 
                        CKDTree<CPuncta>& inPuncta, 
                        deque<CPuncta*>& inStabilizers) const {
      //First, do a logarithmic search for "nearby" points in the kd-tree.
      CPuncta pcCenterBranch( (inLeft+inRight)*0.5 );
      double dNearbyRange = inPunctaRange + 0.5*((inLeft-inRight).Length());
      deque<CPuncta*> dqNearbyPuncta;
      inPuncta.FindPointsInRange( pcCenterBranch, dNearbyRange, dqNearbyPuncta );

      //Now check all nearby puncta for intersection, and label them appropriately.
      if(dqNearbyPuncta.size() > 0) {
         CLine lBranch;
         lBranch.push_back(inLeft);
         lBranch.push_back(inRight);
         deque<CPuncta*>::iterator iter = dqNearbyPuncta.begin();
         deque<CPuncta*>::iterator iter_end = dqNearbyPuncta.end();
         for(; iter != iter_end; ++iter) {
            if(!((*iter)->GetOccupied())) {
#ifdef HOVER
               if(lBranch.DoesSphereHover(**iter, inPunctaRange)) {
#else
               if(lBranch.DoesSphereIntersect(**iter, inPunctaRange)) {
#endif
                  (*iter)->SetOccupied(true);
                  (*iter)->SetStrength(1.);
                  inStabilizers.push_back(*iter);
               }
            }
         }
      }
   }

   //Return true if the input line (inLeft to inRight) is stabilized.
   bool IsStabilized(const CPoint& inLeft, 
                     const CPoint& inRight, 
                     const double& inPunctaRange, 
                     const CKDTree<CPuncta>& inPuncta) const {
      //First, do a logarithmic search for "nearby" points in the kd-tree.
      CPuncta pcCenterBranch( (inLeft+inRight)*0.5 );
      double dNearbyRange = inPunctaRange + 0.5*((inLeft-inRight).Length());
      deque<CPuncta> dqNearbyPuncta;
      inPuncta.FindPointsInRange( pcCenterBranch, dNearbyRange, dqNearbyPuncta );

      //Now check all nearby puncta for intersection, and label them appropriately.
      if(dqNearbyPuncta.size() > 0) {
         CLine lBranch;
         lBranch.push_back(inLeft);
         lBranch.push_back(inRight);
         deque<CPuncta>::const_iterator iter = dqNearbyPuncta.begin();
         deque<CPuncta>::const_iterator iter_end = dqNearbyPuncta.end();
         for(; iter != iter_end; ++iter) {
#ifdef HOVER
            if(lBranch.DoesSphereHover(*iter, inPunctaRange)) {
#else
            if(lBranch.DoesSphereIntersect(*iter, inPunctaRange)) {
#endif
               return true;
            }
         }
      }

      //did not find any stabilizers
      return false;
   }

   //Return the number of stabilizers of the input line (left to right)
   //Note: it's faster to call IsStabilized if all you need to know is whether
   //  or not the input line is stabailized
   int CountStabilizers(const CPoint& inLeft, 
                        const CPoint& inRight, 
                        const double& inPunctaRange, 
                        const CKDTree<CPuncta>& inPuncta) const {
      //First, do a logarithmic search for "nearby" points in the kd-tree.
      CPuncta pcCenterBranch( (inLeft+inRight)*0.5 );
      double dNearbyRange = inPunctaRange + 0.5*((inLeft-inRight).Length());
      deque<CPuncta> dqNearbyPuncta;
      inPuncta.FindPointsInRange( pcCenterBranch, dNearbyRange, dqNearbyPuncta );

      //ofstream testFile("test.txt");
      //testFile << inLeft << "\n" << inRight << "\n\n\n";
      //for(int i = 0; i < dqNearbyPuncta.size(); ++i) {
      //   testFile << dqNearbyPuncta[i] << "\n";
      //}
      //testFile << flush;
      //testFile.close();

      //Now check all nearby puncta for intersection
      int iCount = 0;
      if(dqNearbyPuncta.size() > 0) {
         CLine lBranch;
         lBranch.push_back(inLeft);
         lBranch.push_back(inRight);
         deque<CPuncta>::const_iterator iter = dqNearbyPuncta.begin();
         deque<CPuncta>::const_iterator iter_end = dqNearbyPuncta.end();
         for(; iter != iter_end; ++iter) {
#ifdef HOVER
            if(lBranch.DoesSphereHover(*iter, inPunctaRange)) {
#else
            if(lBranch.DoesSphereIntersect(*iter, inPunctaRange)) {
#endif
               ++iCount;
            }
         }
      }

      //Return the count.
      return iCount;
   }

   //return the determinant of the input matrix (specialized for 3x3)
   double Determinant(const vector<vector<double> >& inM) const {
      return inM[0][0]*inM[1][1]*inM[2][2] + 
         inM[0][1]*inM[1][2]*inM[2][0] +
         inM[0][2]*inM[1][0]*inM[2][1] -
         inM[2][0]*inM[1][1]*inM[0][2] -
         inM[2][1]*inM[1][2]*inM[0][0] -
         inM[2][2]*inM[1][0]*inM[0][1];
   }
   
   //hack the removal radius to break scaling
   //remmoval radius is divided by whatever is returned
   double HackRR(const double& x) const {
       //return (x<20.)?1.:5.;
      return 1.;
   }

private:
   //data members
   BranchPointerType m_brBranch;
   DaughterHolderType m_vecAllBranches;
   AvailableBranchType m_listAvailableBranches;
   double m_dProbExtendLinearly;
};

//allow streaming
//format is x1 y1 (z1...) x2 y2 (z2...) to make it easy for root to draw the set of lines.
ostream& operator<<(ostream& inOut, const CArbor& inArbor) {
   CArbor::DaughterHolderType::const_iterator iter = inArbor.GetAllBranches().begin();
   CArbor::DaughterHolderType::const_iterator iter_end = inArbor.GetAllBranches().end();
   for(; iter != iter_end; ++iter) {
      if( (*iter)->GetLine().size() > 1 ) {
         CLine::const_iterator iterLL = (*iter)->GetLine().begin();
         CLine::const_iterator iterLR = (*iter)->GetLine().begin();
         ++iterLR;
         CLine::const_iterator iterL_end = (*iter)->GetLine().end();
         for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
            for(int i = 0; i < iterLL->GetDimension(); ++i) {
               inOut << (*iterLL)[i] << " ";
            }
            for(int i = 0; i < iterLR->GetDimension(); ++i) {
               inOut << (*iterLR)[i] << " ";
            }
            inOut << "\n";
         }
      }
   }
   return inOut;
}

#endif //ARBOR
