//Joe Snider
//2/05
//
//Model of a brain; just a bunch of arbors and some puncta.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "arbor.h"
#include "KDTree.h"
//#include "rand3.h"
#include "histogram.h"

#ifndef BRAIN
#define BRAIN

#define RANDOM_SPACING

using namespace std;

const int MAX_INSERTION_WITHOUT_PLACEMENT = 100;
const int MAX_STEPS_WITHOUT_PLACEMENT = 1000;

class CBrain{
public:
   //templates

public:
   //constructors (T is a random number generator, must have operator() return a double in [0,1))
   template<class T>
   CBrain(const int& inSize,
          const int& inNumberPuncta,
          const CPoint& inLow,
          const CPoint& inHigh,
          const double& inProbExtendLinearly,
          T& inRand) {
      m_pntLow = inLow;
      m_pntHigh = inHigh;
      m_iNumberPuncta = inNumberPuncta;
      m_dProbExtendLinearly = inProbExtendLinearly;

      BuildRandomPuncta(inRand);
      InitializeArbors(inSize, inRand);
      //BuildTestPuncta(inRand);
      //InitializeTestArbors(inSize, inRand);
   }

public:
   //gets and sets

   int GetCount() const {return (int)m_vecArbors.size();}

   const vector<CArbor>& GetArbors() const {return m_vecArbors;}

   CPoint GetLow() const {return m_pntLow;}
   CPoint GetHigh() const {return m_pntHigh;}
   int GetNumberPuncta() const {return m_iNumberPuncta;}
   
   double GetProbExtendLinearly() const {return m_dProbExtendLinearly;}
   void SetProbExtendLinearly(const double& inPEL) {
      m_dProbExtendLinearly = inPEL;
      //update all the arbors
      vector<CArbor>::iterator iter = m_vecArbors.begin();
      vector<CArbor>::iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->SetProbExtendLinearly(inPEL);
      }
   }

public:
   //interface
    
    //scale the arbors
    void Scale(const double& x, const double& y, const double& z) {
        vector<CArbor>::iterator iter = m_vecArbors.begin();
        vector<CArbor>::iterator iter_end = m_vecArbors.end();
        for(; iter != iter_end; ++iter) {
            iter->Scale(x,y,z);
        }
    }

   //Dump the product moments to the stream.
   //inMinLevel optionally ignores any without a level at least that (enough branchings).
   void ProductMoment(ostream& inOut, const int& inLow, const int& inHigh, const int& inStep, int inMinLevel = 0) const {
      int i = 0;
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter, ++i) {
         if(iter->MaxLevel() > inMinLevel) {
            inOut << "simulated_" << i << " ";
            for(int j = inLow; j <= inHigh; j += inStep) {
               inOut << iter->ProductMoment(j);
               if(j < inHigh) {
                  inOut << " ";
               } else {
                  inOut << "\n" << flush;
               }
            }
         }
      }
   }

   //Dump the arbor to stdout in .swc format.
   void DumpSWC() const {
      for(unsigned i = 0; i < m_vecArbors.size(); ++i) {
         m_vecArbors[i].Write(cout);
      }
   }

   //Save all the arbors to a format the Neuroleucida can recognize.
   //Creates a file for each arbor.
   void Save(const string& inFileBase = "simulated") const {
      char cstrBuffer[100]; //100 characters is plenty, just stores integers
      for(unsigned i = 0; i < m_vecArbors.size(); ++i) {
         sprintf(cstrBuffer, "%d", i);
         string strFileName = inFileBase;
         strFileName.append("_");
         strFileName.append(cstrBuffer);
         strFileName.append(".swc");
         ofstream ofData(strFileName.c_str());
         if(ofData.good()) {
            m_vecArbors[i].Write(ofData);
         } else {
            cerr << "Unable to create file "
               << strFileName << " ... "
               << "stopping write\n" << flush;
            return;
         }
      }
   }

   //Return the number of occupied puncta.
   //Unfortunately, this gets a bit messy with the kd-tree, but it's worth it to do the point look up fast.
   int CountOccupiedPuncta() const {
      int iRet = 0;
      //Loop over all the leaves of the kd-tree and count the occupied.
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL = m_kdPuncta.GetLeaves().begin();
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL_end = m_kdPuncta.GetLeaves().end();
      for(; iterL != iterL_end; ++iterL) {
         CKDNode<CPuncta>::DataType::const_iterator iterD = (*iterL)->GetData().begin();
         CKDNode<CPuncta>::DataType::const_iterator iterD_end = (*iterL)->GetData().end();
         for(; iterD != iterD_end; ++iterD) {
            iRet += iterD->IsOccupied()?1:0;
         }
      }
      return iRet;
   }

   //Count the number of stabilizers at each level of all branches.
   //Modifies inL_to_S to contain the level, the average number of stabilizers, and the standard error.
   void LevelStabilizers(map<int, pair<double, double> >& inL_to_S) const {
      //Get the data in vectors.
      map<int, vector<unsigned> > mapL_S;
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->LevelStabilizers(mapL_S);
      }

      //Do the averaging.
      map<int, vector<unsigned> >::const_iterator iterL = mapL_S.begin();
      map<int, vector<unsigned> >::const_iterator iterL_end = mapL_S.end();
      for(; iterL != iterL_end; ++iterL) {
         double dMean = 0.;
         vector<unsigned>::const_iterator iterS = iterL->second.begin();
         vector<unsigned>::const_iterator iterS_end = iterL->second.end();
         for(; iterS != iterS_end; ++iterS) {
            dMean += double(*iterS);
         }
         dMean /= double(iterL->second.size());
         
         double dStdErr = 0.;
         iterS = iterL->second.begin();
         iterS_end = iterL->second.end();
         for(; iterS != iterS_end; ++iterS) {
            dStdErr += (double(*iterS) - dMean)*(double(*iterS) - dMean);
         }
         dStdErr = sqrt(dStdErr)/double(iterL->second.size());

         inL_to_S[iterL->first].first = dMean;
         inL_to_S[iterL->first].second = dStdErr;
      }
   }

   //Count the stabilizers of all branches.
   //Modifies inS_to_N.
   void Stabilizers(map<int, unsigned>& inS_to_N) const {
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->Stabilizers(inS_to_N);
      }
   }

   //Find the angles formed by all branches.
   //Modifies inHisto.
   void Angles(CHistogram<double>& inHisto) const {
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->Angles(inHisto);
      }
   }

   //Dump the raw lengths to the input stream.
   void RawAngles(ostream& inOut) const {
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->RawAngles(inOut);
      }
   }

   //Find the lengths of all branches.
   //Modifies inHisto.
   void Lengths(CHistogram<double>& inHisto) const {
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->Lengths(inHisto);
      }
   }

   //Dump the raw lengths to the input stream.
   void RawLengths(ostream& inOut) const {
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->RawLengths(inOut);
      }
   }

   //Trim the ends of the branches back to their daughters.
   //May modify the branches, but won't change any of the stabilizers.
   void Trim(const double& inPunctaRange, const int& inMinStabilizers = 1) {
      vector<CArbor>::iterator iter = m_vecArbors.begin();
      vector<CArbor>::iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->TrimDaughter();
         iter->UpdateLevels();
         iter->TrimEnds(inPunctaRange, m_kdPuncta, inMinStabilizers);
         while( iter->TrimEnds(inPunctaRange, m_kdPuncta, inMinStabilizers) > 0 ) {
            if(!iter->TrimDaughter()) {
               cerr << "Warning in CBrain::Trim: TrimDaughter failed...stopping trim\n" << flush;
               return;
            }
            iter->UpdateLevels();
         }
      }
   }

   //Grow the neurons.
   //The type T must have operator() return a double in [0,1), e.g. a random number generator.
   //Returns the number of puncta used.
   //Changed 4/06 - inMaxCount now specifies the max number of steps allowed without a placement.
   template <class T>
   int Grow(const int& inMaxCount, 
            const double& inLength, 
            const double& inPunctaRange, 
            const double& inPunctaRemoveRadius, 
            T& inRand) {
      int iNumPlaced = 0;
      int iStepsWithoutPlacement = 0;
      int iNumArbors = (int)m_vecArbors.size();
      int iMaxTrials = MAX_STEPS_WITHOUT_PLACEMENT*iNumArbors;
      int iArborToExtend;

      //The main loop where placement occurs.
      while( iStepsWithoutPlacement < inMaxCount ) {
         //Choose a random arbor.
         iArborToExtend = inRand(iNumArbors);
         //iArborToExtend = SmallestArbor();
         CArbor& arborTemp = m_vecArbors[iArborToExtend];

         if( arborTemp.Step(RandomValue(inLength, inRand), inPunctaRange, inPunctaRemoveRadius, inRand, m_kdPuncta) ) {
            ++iNumPlaced;
            iStepsWithoutPlacement = 0;
            //Check if there are available puncta (return if not).
            //TODO: (this could be optimized by letting CArbor::Step return the number of stabilizers).
            if( !ArePunctaAvailable() ) {
               return m_iNumberPuncta;
            }
         } else {
            ++iStepsWithoutPlacement;
         }

         //Give up eventually.
         if( iStepsWithoutPlacement > iMaxTrials ) {
               return CountOccupiedPuncta();
         }
      }

      return CountOccupiedPuncta();
   }

   //dump the puncta.
   //Unfortunately, this gets a bit messy with the kd-tree, but it's worth it to do the point look up fast.
   void DumpPuncta(ostream& inOut) const {
      //Loop over all the leaves of the kd-tree and dump their data (awkward).
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL = m_kdPuncta.GetLeaves().begin();
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL_end = m_kdPuncta.GetLeaves().end();
      for(; iterL != iterL_end; ++iterL) {
         CKDNode<CPuncta>::DataType::const_iterator iterD = (*iterL)->GetData().begin();
         CKDNode<CPuncta>::DataType::const_iterator iterD_end = (*iterL)->GetData().end();
         for(; iterD != iterD_end; ++iterD) {
            for(int i = 0; i < iterD->GetDimension(); ++i) {
               inOut << (*iterD)[i] << " ";
            }
            inOut << ((iterD->IsOccupied())?1:0) << " "
                    << iterD->GetStrength() << "\n";
         }
      }
   }

   //Dump the segments.
   void DumpSegments(ostream& inOut, int iID = 0, const string& header = "") const {
      inOut << header;
      int i;
      vector<CArbor>::const_iterator iterA = m_vecArbors.begin();
      vector<CArbor>::const_iterator iterA_end = m_vecArbors.end();
      for(; iterA != iterA_end; ++iterA, ++iID) {
         CArbor::DaughterHolderType::const_iterator iter = iterA->GetAllBranches().begin();
         CArbor::DaughterHolderType::const_iterator iter_end = iterA->GetAllBranches().end();
         for(; iter != iter_end; ++iter) {
            if( (*iter)->GetLine().size() > 1 ) {
               CLine::const_iterator iterLL = (*iter)->GetLine().begin();
               CLine::const_iterator iterLR = (*iter)->GetLine().begin();
               ++iterLR;
               CLine::const_iterator iterL_end = (*iter)->GetLine().end();
               for(; iterLR != iterL_end; ++iterLL, ++iterLR) {
                  for(i = 0; i < iterLL->GetDimension(); ++i) {
                     inOut << (*iterLL)[i] << " ";
                  }
                  for(i = 0; i < iterLR->GetDimension(); ++i) {
                     inOut << (*iterLR)[i] << " ";
                  }
                  inOut << " " << iID << " ";
                  inOut << " " << (*iter)->GetLevel() << "\n";
               }
            }
         }
      }
   }

   //Dump the inertial volumes and total lengths of the arbors.
   void DumpInertialVolume(ostream& inOut) const {
      double dIV, dTL;
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->InertialVolume(dIV, dTL);
         if(dIV > 1.e-6 && dTL > 0.) {
            inOut << dTL << " " << dIV << "\n";
         }
      }
   }

   //Find the profile with inertial histogram normalization.
   //Aligns and scales to an uniform neuron (rotate by inertial and scale to [-1,1]).
   //Then finds the average length at a given distance.
   void InertialHistogram(ostream& inOut) const {
      //Initialize and get the values.
      CHistogram<double> histoMass(0., 2., 20);
      histoMass.zero();
      CHistogram<double> histoCount(0., 2., 20);
      histoCount.zero();
      vector<CArbor>::const_iterator iter = m_vecArbors.begin();
      vector<CArbor>::const_iterator iter_end = m_vecArbors.end();
      for(; iter != iter_end; ++iter) {
         iter->InertialHistogram(histoMass, histoCount);
      }

      //build the histogram (and normalize).
      double a, b;
      double dNormalization = 0.;
      int iDimension = m_pntLow.GetDimension();
      CHistogram<double>::iterator iterH = histoMass.bin_begin();
      CHistogram<double>::iterator iterHC = histoCount.bin_begin();
      CHistogram<double>::iterator iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         if(*iterHC > 0.5) {
            a = histoMass.GetBinLeftEdge(iterH);
            b = histoMass.GetBinRightEdge(iterH);
            *iterH /= (iDimension==2)?acos(-1.)*(b*b-a*a):4./3.*acos(-1.)*(b*b*b-a*a*a);
            dNormalization += *iterH;
         }
      }

      //dump the histogram
      dNormalization = 1./dNormalization;
      iterH = histoMass.bin_begin();
      iterHC = histoCount.bin_begin();
      iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         if(*iterHC > 0.5) {
            inOut << histoMass.GetBinCenter(iterH) << " "
               << *iterH * dNormalization << " "
               << *iterH * dNormalization / sqrt(*iterHC) << "\n";
         }
      }
   }

private:
   //helpers

   //Return true if there is at least one available puncta.
   int ArePunctaAvailable() const {
      int iRet = 0;
      //Loop over all the leaves of the kd-tree and count the occupied.
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL = m_kdPuncta.GetLeaves().begin();
      CKDTree<CPuncta>::LeafHolderType::const_iterator iterL_end = m_kdPuncta.GetLeaves().end();
      for(; iterL != iterL_end; ++iterL) {
         CKDNode<CPuncta>::DataType::const_iterator iterD = (*iterL)->GetData().begin();
         CKDNode<CPuncta>::DataType::const_iterator iterD_end = (*iterL)->GetData().end();
         for(; iterD != iterD_end; ++iterD) {
            if( !iterD->IsOccupied() ) {
               return true;
            }
         }
      }
      return false;
   }

   //Insert the puncta.
   //Randomly puts inNumber of them down on the rectangle (or prism in 3d, or whatever in 4d...)
   //   defined by inLow, inHigh.
   //The type T must have operator() overloaded to return a number between [0,1).
   template <class T>
   void BuildRandomPuncta(T& inRand) {
      vector<CPuncta> vecPuncta;
      CPuncta pcInsert;
      for(int i = 0; i < m_iNumberPuncta; ++i) {
         pcInsert.clear();
         for(int j = 0; j < m_pntLow.GetDimension(); ++j) {
             pcInsert.push_back(m_pntLow[j] + (m_pntHigh[j]-m_pntLow[j])*inRand());
         }
         vecPuncta.push_back(pcInsert);
      }
      m_kdPuncta.Clear();
      m_kdPuncta.InsertData(vecPuncta.begin(), vecPuncta.end(), (int)vecPuncta.size(), m_pntLow.GetDimension());
   }

   //Insert the puncta in a test configuration.
   //Makes a line of them through the middle (m_pntLow to m_pntHigh).
   //The type T must have operator() overloaded to return a number between [0,1).
   template <class T>
   void BuildTestPuncta(T& inRand) {
      vector<CPuncta> vecPuncta;
      for(int i = 0; i < m_iNumberPuncta; ++i) {
         CPuncta pcInsert;
         pcInsert = m_pntLow + (m_pntHigh-m_pntLow)*inRand();
         vecPuncta.push_back(pcInsert);
      }
      m_kdPuncta.Clear();
      m_kdPuncta.InsertData(vecPuncta.begin(), vecPuncta.end(), vecPuncta.size(), m_pntLow.GetDimension());
   }

   //Initialize the arbors.
   //Randomly places inSize of them on a rectangle with some self avoidance.
   //Or if inSize == 1, then just put one in the center.
   //Also, weakly avoid the boundary (arbitrarily pad by 10% for now)
   //The type T must have operator() overloaded to return a number between [0,1).
   template<class T>
   void InitializeArbors(const int& inSize, T& inRand) {
      m_vecArbors.resize(inSize);

      //If there's only one arbor, then put it in the center.
      if(inSize == 1) {
         CPoint pntCenter = 0.5*(m_pntLow+m_pntHigh);
         m_vecArbors[0].CreateStart(pntCenter);
         m_vecArbors[0].SetProbExtendLinearly(m_dProbExtendLinearly);
         return;
      }
#ifdef RANDOM_SPACING
      for(int i = 0; i < inSize; ++i) {
         //Just put them down randomly 
         CPoint pntCenter;
         for(int j = 0; j < m_pntLow.GetDimension(); ++j) {
            pntCenter.push_back(m_pntLow[j] + (m_pntHigh[j]-m_pntLow[j])*inRand());
         }

         //Initialize the new arbor.
         m_vecArbors[i].CreateStart(pntCenter);
         m_vecArbors[i].SetProbExtendLinearly(m_dProbExtendLinearly);
      }
#else //RANDOM_SPACING
      //Start with a spacing that assigns equal area to each arbor.
      double dSpacing = 1.;
      for(int j = 0; j < m_pntLow.GetDimension(); ++j) {
         dSpacing *= m_pntHigh[j]-m_pntLow[j];
      }
      dSpacing = pow( dSpacing / double(inSize), 1./double(m_pntLow.GetDimension()) );

      //Do the insertion.
      for(int i = 0; i < inSize; ++i) {
         int iTrial = 0;
         double dMinDistance;
         CPoint pntCenter;
         do {
            //Choose a new random point.
            pntCenter.clear();
            for(int j = 0; j < m_pntLow.GetDimension(); ++j) {
               pntCenter.push_back(m_pntLow[j] + (m_pntHigh[j]-m_pntLow[j])*(0.8*inRand()+0.1));
               //pntCenter.push_back(m_pntLow[j] + (m_pntHigh[j]-m_pntLow[j])*inRand());
            }
            //if( m_pntLow.GetDimension() == 3 ) {
            //   pntCenter[2] = 0.;
            //}

            //Find the minimum distance to other points.
            dMinDistance = 1.e100;
            for(int k = 0; k < i; ++k) {
               double dDist = (m_vecArbors[k].Start() - pntCenter).Length();
               dMinDistance = min( dDist, dMinDistance );
            }

            //Decrease the desired spacing if this is taking too long.
            ++iTrial;
            if( iTrial > MAX_INSERTION_WITHOUT_PLACEMENT) {
               dSpacing *= 0.9;
               iTrial = 0;
            }
         } while ( dMinDistance < dSpacing );

         //Initialize the new arbor.
         m_vecArbors[i].CreateStart(pntCenter);
         m_vecArbors[i].SetProbExtendLinearly(m_dProbExtendLinearly);
      }
#endif //RANDOM_SPACING
   }

   //Initialize the arbors.
   //Place them all along the center line.
   //The type T must have operator() overloaded to return a number between [0,1).
   template<class T>
   void InitializeTestArbors(const int& inSize, T& inRand) {
      m_vecArbors.resize(inSize);

      for(int i = 0; i < inSize; ++i) {
         CPoint pntCenter = m_pntLow + (m_pntHigh-m_pntLow)*inRand();
         m_vecArbors[i].CreateStart(pntCenter);
         m_vecArbors[i].SetProbExtendLinearly(m_dProbExtendLinearly);
      }
   }

   //Return a value with some distribution.
   //The type T must have operator() overloaded to return a number between [0,1).
   template<class T>
   double RandomValue(const double& inX, T& inRand) const {
      /*//A Gaussian in the range [0, 2*inX].
      double dGauss = 0.;
      int iNumIterations = 4;
      for(int i = 0; i < iNumIterations; ++i) {
         dGauss += inRand();
      }
      return inX*dGauss*2./double(iNumIterations);*/
      //uniform in the range inX*[0.5, 1.5]
      return inX*(0.5+inRand());
   }

   //Return the id of the smallest arbor.
   template<class T>
   int SmallestArbor() const {
      int iRet = -1;
      double dMin = 1.e100;
      int i = 0;
      for(int i = 0; i < m_vecArbors.size(); ++i) {
         double dTest = double(m_vecArbors[i].GetAllBranches().size());
         if(dMin > dTest) {
            dMin = dTest;
            iRet = i;
         }
      }
      return iRet;
   }

private:
   //data members
   vector<CArbor> m_vecArbors;
   CKDTree<CPuncta> m_kdPuncta;
   CPoint m_pntLow;
   CPoint m_pntHigh;
   int m_iNumberPuncta;
   double m_dProbExtendLinearly;

};

#endif //BRAIN
