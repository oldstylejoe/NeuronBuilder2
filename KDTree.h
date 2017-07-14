//Joe Snider
//6/02
//
//updated 2/04 - now really a kd-tree
//
//Impelment a kd-tree container
//
//TODO: allow insertion and deletion of single points.

#ifndef KDTREE
#define KDTREE

#include <vector>
#include <deque>
#include <math.h>

#include "Median.h"
#include "KDNode.h"

using namespace std;

ostream& operator<<(ostream& osOut, const vector<double>& x) {
   vector<double>::const_iterator iter = x.begin();
   vector<double>::const_iterator iter_end = x.end();
   if(iter != iter_end) {
      osOut << "(" << *iter;
      ++iter;
      for(; iter != iter_end; ++iter) {
         osOut << "," << *iter;
      }
      osOut << ")";
   }

   return osOut;
}

//number of data points allowed in a leaf
const unsigned LEAF_SIZE = 10;

//Comparison class for sorting.
//Warning: be sure that m_iElement is set.
//The type T must have operator[] defined with return type sortable by operator<.
template <class T>
class CComp{
public:
   bool operator()(const T& x, const T& y) const {
      return x[m_iElement] < y[m_iElement];
   }

   int m_iElement;   //element on which to sort; bad hiding, but fast
};

//Type PointType must be a valid type for class CComp (above).
//Also must have size() overloaded which returns the dimension.
template <class PointType>
class CKDTree{
public:
   //typedefs
   typedef vector<CKDNode<PointType>*> NodeHolderType;
   typedef vector<CKDNode<PointType>*> LeafHolderType;

public:
   //constructor

   CKDTree() {
      m_Compare.m_iElement = 0;
   }

   ~CKDTree() {
      Clear();
   }

public:
   //interface

   //clear the structure for new data
   void Clear() {
      //delete all nodes (gets leaves too)
      typename NodeHolderType::iterator iterNode = m_NodeHolder.begin();
      typename NodeHolderType::iterator iterNode_end = m_NodeHolder.end();
      for(; iterNode != iterNode_end; ++iterNode) {
         delete *iterNode;
      }
      m_NodeHolder.clear();
      m_LeafHolder.clear();

      m_Compare.m_iElement = 0;
   }

   //get leaf data
   const LeafHolderType& GetLeaves() const {
      return m_LeafHolder;
   }

   //insert the data
   //IterInsert must be an STL style iterator to a
   // set of PointTypes
   template <class IterInsert>
   void InsertData(IterInsert inBegin,
                   IterInsert inEnd,
                   const int& inSize,
                   const int& inDimension) {
      if( m_NodeHolder.size() > 0 ) {
         Clear();
      }
      m_Compare.m_iElement = 0;
      m_iDimension = inDimension;

      //build the root
      CKDNode<PointType> * nodeRoot = new CKDNode<PointType>(0);
      nodeRoot->SetParent(NULL);
      m_NodeHolder.push_back(nodeRoot);

      //call the local helper that can be called recursively
      DoInsert(inBegin, inEnd, 0, inSize, nodeRoot);
   }

   //display the tree structure
   void Display(ostream& inOut) {
      DoDisplay(inOut, *(m_NodeHolder.begin()));
   }

   //Display the points without any formatting.
   void DisplayPoints(ostream& inOut) const {
      typename LeafHolderType::const_iterator iter = m_LeafHolder.begin();
      typename LeafHolderType::const_iterator iter_end = m_LeafHolder.end();
      for(; iter != iter_end; ++iter) {
         typename CKDNode<PointType>::DataType::const_iterator iterD = (*iter)->GetData().begin();
         typename CKDNode<PointType>::DataType::const_iterator iterD_end = (*iter)->GetData().end();
         for(;iterD != iterD_end; ++iterD) {
            for(int i = 0; i < m_iDimension; ++i) {
               inOut << iterD->at(i) << " ";
            }
            inOut << "\n";
         }
      }
   }

   //Search for all points within a distance of the input point
   //modifies inReturn
   //returns the number of points examined (for testing the time)
   int FindPointsInRange(const PointType& inPoint,
                          const double& inRange,
                          deque<PointType>& inReturn) const {
      deque<PointType> dqTemp;

      ////check dimensionality of input point
      //if(inPoint.size() != m_iDimension) {
      //   cerr << "Error: dimensions do not match...cancelling search\n" << flush;
      //   return 0;
      //}

      //Find all nodes in all possible leaves in the range (O(area*log(n)) elimination)
      //Recursive search modifying vecReturn.
      FindLeavesInRange(inPoint,
                        inRange,
                        *(m_NodeHolder.begin()),
                        dqTemp);

      //eliminate nodes not contained by the circle (O(n) search of left overs)
      double dRangeSquared = inRange*inRange;
      typename deque<PointType>::const_iterator iter = dqTemp.begin();
      typename deque<PointType>::const_iterator iter_end = dqTemp.end();
      for(; iter != iter_end; ++iter) {
         if(DistanceSquared(inPoint, *iter) < dRangeSquared) {
            inReturn.push_back(*iter);
         }
      }

      return dqTemp.size();
   }

   //Search for all points within a distance of the input point
   //modifies inReturn, and inReturn contains references to the actual kd-tree data.
   //Note: changing the kd-tree coordinates leads to unknown results (points won't be found, etc...),
   // but this allows a more complex object for the point (it can be a class, and we can easily modify it).
   //returns the number of points examined (for testing the time)
   int FindPointsInRange(const PointType& inPoint,
                         const double& inRange,
                         deque<PointType*>& inReturn) {
      deque<typename CKDNode<PointType>::DataType::iterator> dqTemp;

      ////check dimensionality of input point
      //if(inPoint.size() != m_iDimension) {
      //   cerr << "Error: dimensions do not match...cancelling search\n" << flush;
      //   return 0;
      //}

      //Find all nodes in all possible leaves in the range (O(area*log(n)) elimination)
      //Recursive search modifying vecReturn.
      FindLeavesInRange(inPoint,
                        inRange,
                        *(m_NodeHolder.begin()),
                        dqTemp);

      //eliminate nodes not contained by the circle (O(n) search of left overs)
      double dRangeSquared = inRange*inRange;
      typename deque<typename CKDNode<PointType>::DataType::iterator>::iterator iter = dqTemp.begin();
      typename deque<typename CKDNode<PointType>::DataType::iterator>::iterator iter_end = dqTemp.end();
      for(; iter != iter_end; ++iter) {
         if(DistanceSquared(inPoint, **iter) < dRangeSquared) {
            inReturn.push_back(&(**iter));
         }
      }

      return (int)dqTemp.size();
   }

   //Find the closest point to the input point.
   //Just checks all the points.
   //Testing only (use FindClosestPoint).
   PointType DumbFindClosestPoint(const PointType& inPoint) const {
      double dDistSquared = 1.e100;
      PointType pRet;
      typename LeafHolderType::const_iterator iter = m_LeafHolder.begin();
      typename LeafHolderType::const_iterator iter_end = m_LeafHolder.end();
      for(; iter != iter_end; ++iter) {
         typename CKDNode<PointType>::DataType::const_iterator iterD = (*iter)->GetData().begin();
         typename CKDNode<PointType>::DataType::const_iterator iterD_end = (*iter)->GetData().end();
         for(;iterD != iterD_end; ++iterD) {
            if(DistanceSquared( inPoint, *iterD ) < dDistSquared) {
               dDistSquared = DistanceSquared( inPoint, *iterD );
               pRet = *iterD;
            }
         }
      }
      return pRet;
   }

   //Find the closest point to the input point.
   //Must be at least 4 occupied leaves.
   PointType FindClosestPoint(const PointType& inPoint) const {
      if( m_LeafHolder.size() < 4 ) {
         return inPoint;
      }

      //Find the distance to the first hyperplane.
      double dR = fabs(inPoint[0] - m_NodeHolder[0]->GetHyperPlane()[0]);

      deque<PointType> dqTemp;

      //Find all nodes in all possible leaves in the range (O(area*log(n)) elimination)
      //Recursive search modifying vecReturn.
      FindLeavesInRange(inPoint, dR, *(m_NodeHolder.begin()), dqTemp);
      //cout << "KDTree.h gh0 " << inPoint << " " << dqTemp.size() << " " << dR << "\n" << flush;

      //Decrease (increase) the distance until there are at least 1 and at most 4*LEAF_SIZE+1
      // points to test.
      //Should be some sort of adaptive ranges here. It's possible this won't converge
      while(dqTemp.size() < 1 || dqTemp.size() > 4*LEAF_SIZE+1) {
         if(dqTemp.size() < 1) {
            dR *= 1.1;
         } else if(dqTemp.size() > 4*LEAF_SIZE+1) {
            dR *= .9;
         }
         dqTemp.clear();
         FindLeavesInRange(inPoint, dR, *(m_NodeHolder.begin()), dqTemp);
         //cout << "KDTree.h gh1 " << dqTemp.size() << " " << dR << "\n" << flush;
      }

      double dDistSquared = DistanceSquared(inPoint, dqTemp[0]);
      PointType pRet = dqTemp[0];
      typename deque<PointType>::iterator iter = dqTemp.begin();
      ++iter;
      typename deque<PointType>::iterator iter_end = dqTemp.end();
      for(; iter != iter_end; ++iter) {
         if(DistanceSquared(inPoint, *iter) < dDistSquared) {
            dDistSquared = DistanceSquared(inPoint, *iter);
            pRet = *iter;
         }
      }

      return pRet;
   }

private:

   //recurse the tree and display
   void DoDisplay(ostream& inOut, CKDNode<PointType>* inNode) const {
      inOut << inNode->GetDepth() << ": ";
      if( inNode->IsLeaf() ) {
         inOut << "Leaf ";
         typename CKDNode<PointType>::DataType::const_iterator iterData =
            inNode->GetData().begin();
         typename CKDNode<PointType>::DataType::const_iterator iterData_end =
            inNode->GetData().end();
         for(;iterData != iterData_end; ++iterData) {
            inOut << *iterData << " ";
         }
         inOut << "\n";
      } else if ( inNode->GetDepth() == 0 ) {
         inOut << "Root ";
         inOut << inNode->GetHyperPlane() << "\n";
         DoDisplay(inOut, inNode->GetLeftDaughter());
         DoDisplay(inOut, inNode->GetRightDaughter());
      } else {
         inOut << "Node ";
         inOut << inNode->GetHyperPlane() << "\n";
         DoDisplay(inOut, inNode->GetLeftDaughter());
         DoDisplay(inOut, inNode->GetRightDaughter());
      }
   }

   //do the inserting
   template <class IterInsert>
   void DoInsert(IterInsert inBegin,
                 IterInsert inEnd,
                 const int& inLevel,
                 const unsigned& inSize,
                 CKDNode<PointType>* inParent) {
      if(inSize < LEAF_SIZE) {
         //check if its a leaf and stop if so
         inParent->SetData(inBegin, inEnd);
         inParent->SetLeftDaughter(NULL);
         inParent->SetRightDaughter(NULL);
         m_LeafHolder.push_back(inParent);
      } else {
         //set the coordinate on which to sort
         m_Compare.m_iElement = (inLevel % m_iDimension);

         //find the median at the current level
         PointType ptMedian;
         HoareMedian(inBegin, inEnd, inSize, ptMedian, m_Compare);
         inParent->SetHyperPlane(ptMedian);

         //split data into left and right
         //Put data equal to the pivot into the smaller of left and right (default right).
         //   Maintains balance and prevents all elements from ending up in one bin.
         //   this gives an extra if statement, could be optimized?
         deque<PointType> dqLeft;
         deque<PointType> dqRight;
         IterInsert iter = inBegin;
         for(; iter != inEnd; ++iter) {
            if( m_Compare(*iter, ptMedian) ) {
               dqLeft.push_back(*iter);
            } else if( m_Compare(ptMedian, *iter) ) {
               dqRight.push_back(*iter);
            } else {
               if( dqLeft.size() < dqRight.size() ) {
                  dqLeft.push_back(*iter);
               } else {
                  dqRight.push_back(*iter);
               }
            }
         }

         //-Removed for production runs---Joe
         //check that daughters are smaller than the parent
         //this is fairly fast, so I keep it in even though it does
         //  slow things down slightly
         /*if( (dqLeft.size() >= inSize) || (dqRight.size() >= inSize) ) {
            // warn and mark as a leaf
            inParent->SetData(inBegin, inEnd);
            inParent->SetLeftDaughter(NULL);
            inParent->SetRightDaughter(NULL);
            m_LeafHolder.push_back(inParent);
            cerr << "Warning: daughter >= parent at level " << flush;
            cerr << inParent->GetDepth();
            cerr << ".\nLeft = "
                 << dqLeft.size() << " Right = "
                 << dqRight.size() << " Parent = "
                 << inSize << "\n"
                 << "pivot = " << ptMedian << " ...Continuing\n" << flush;
         } else {
            //build daughter nodes and recurse
            //left
            CKDNode<PointType>* nodeLeftInsert = new CKDNode<PointType>(inLevel+1);
            m_NodeHolder.push_back(nodeLeftInsert);
            nodeLeftInsert->SetParent(inParent);
            inParent->SetLeftDaughter(nodeLeftInsert);
            DoInsert(dqLeft.begin(),
                     dqLeft.end(),
                     inLevel+1,
                     dqLeft.size(),
                     nodeLeftInsert);
            //right
            CKDNode<PointType>* nodeRightInsert = new CKDNode<PointType>(inLevel+1);
            m_NodeHolder.push_back(nodeRightInsert);
            nodeRightInsert->SetParent(inParent);
            inParent->SetRightDaughter(nodeRightInsert);
            DoInsert(dqRight.begin(),
                     dqRight.end(),
                     inLevel+1,
                     dqRight.size(),
                     nodeRightInsert);
         }*/
         ////end removed---replacement follows
         //build daughter nodes and recurse
         //left
         CKDNode<PointType>* nodeLeftInsert = new CKDNode<PointType>(inLevel+1);
         m_NodeHolder.push_back(nodeLeftInsert);
         nodeLeftInsert->SetParent(inParent);
         inParent->SetLeftDaughter(nodeLeftInsert);
         DoInsert(dqLeft.begin(),
                  dqLeft.end(),
                  inLevel+1,
                  (int)dqLeft.size(),
                  nodeLeftInsert);
         //right
         CKDNode<PointType>* nodeRightInsert = new CKDNode<PointType>(inLevel+1);
         m_NodeHolder.push_back(nodeRightInsert);
         nodeRightInsert->SetParent(inParent);
         inParent->SetRightDaughter(nodeRightInsert);
         DoInsert(dqRight.begin(),
                  dqRight.end(),
                  inLevel+1,
                  (int)dqRight.size(),
                  nodeRightInsert);
         ////end replacement
      }
   }

   //Recurse the tree and search for leaves in a range.
   //modifies inNodes
   void FindLeavesInRange(const PointType& inPoint,
                          const double& inR,
                          CKDNode<PointType>* inNode,
                          deque<PointType>& inNodes) const {
      //check if this is a leaf
      if(inNode->IsLeaf()) {
         //insert the nodes of this leaf into the list and stop recursing
         typename CKDNode<PointType>::DataType::const_iterator iter = 
            inNode->GetData().begin();
         typename CKDNode<PointType>::DataType::const_iterator iter_end =
            inNode->GetData().end();
         for(; iter != iter_end; ++iter) {
            inNodes.push_back(*iter);
         }
      } else {
         //find the distance to the current hyperplane
         int iDepth = inNode->GetDepth();
         int iDimension = (iDepth % m_iDimension);
         double dDist = inPoint[iDimension] - (inNode->GetHyperPlane())[iDimension];

         if(fabs(dDist) > inR) {
            //the distance is greater than the range so only check one daughter
            if(dDist < 0.) {
               FindLeavesInRange(inPoint, inR, inNode->GetLeftDaughter(), inNodes);
            } else {
               FindLeavesInRange(inPoint, inR, inNode->GetRightDaughter(), inNodes);
            }
         } else {
            //the distance is less that the range so must check both daughters
            FindLeavesInRange(inPoint, inR, inNode->GetLeftDaughter(), inNodes);
            FindLeavesInRange(inPoint, inR, inNode->GetRightDaughter(), inNodes);
         }
      }
   }

   //Recurse the tree and search for leaves in a range.
   //modifies inNodes
   void FindLeavesInRange(const PointType& inPoint,
                          const double& inR,
                          CKDNode<PointType>* inNode,
                          deque<typename CKDNode<PointType>::DataType::iterator>& inNodes) {
      //check if this is a leaf
      if(inNode->IsLeaf()) {
         //insert the nodes of this leaf into the list and stop recursing
         typename CKDNode<PointType>::DataType::iterator iter = inNode->GetData().begin();
         typename CKDNode<PointType>::DataType::iterator iter_end = inNode->GetData().end();
         for(; iter != iter_end; ++iter) {
            inNodes.push_back(iter);
         }
      } else {
         //find the distance to the current hyperplane
         int iDepth = inNode->GetDepth();
         int iDimension = (iDepth % m_iDimension);
         double dDist = inPoint[iDimension] - (inNode->GetHyperPlane())[iDimension];

         if(fabs(dDist) > inR) {
            //the distance is greater than the range so only check one daughter
            if(dDist < 0.) {
               FindLeavesInRange(inPoint, inR, inNode->GetLeftDaughter(), inNodes);
            } else {
               FindLeavesInRange(inPoint, inR, inNode->GetRightDaughter(), inNodes);
            }
         } else {
            //the distance is less that the range so must check both daughters
            FindLeavesInRange(inPoint, inR, inNode->GetLeftDaughter(), inNodes);
            FindLeavesInRange(inPoint, inR, inNode->GetRightDaughter(), inNodes);
         }
      }
   }

   //return the Euclidean distance squared between two points
   double DistanceSquared(const PointType& inX, const PointType& inY) const {
      double dReturn = 0.;
      for(unsigned i = 0; i < inX.size(); ++i) {
         dReturn += (inX[i]-inY[i])*(inX[i]-inY[i]);
      }
      return dReturn;
   }

private:
   //data
   NodeHolderType m_NodeHolder;
   LeafHolderType m_LeafHolder;
   int m_iDimension;

   CComp<PointType> m_Compare;

};

#endif //KDTREE
