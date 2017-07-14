//Joe Snider
//6/02
//
//a node of a kd-tree

#ifndef KDNODE
#define KDNODE

#include <deque>
#include <vector>

using namespace std;

//see KDTree.h for requirements on PointType.
template<class PointType>
class CKDNode{
public:
   //typedefs
   typedef deque<PointType> DataType;

public:
   //constructor
   //depth specifies coordinate being sorted
   //hyperplane determines left and right sets
   //  empty if its a leaf
   //data contains the data in the node
   //  empty unless the node is a leaf
   CKDNode(const int& inDepth): m_iDepth(inDepth)
   {}

public:
   //interface

   //returns true if this node is a leaf
   bool IsLeaf() const {return m_dqData.size() > 0;}

public:
   //get and sets

   //parent
   CKDNode* GetParent() {return m_pParent;}
   void SetParent(CKDNode* inParent) {m_pParent = inParent;}

   //left daughter
   CKDNode* GetLeftDaughter() {return m_pLeftDaughter;}
   void SetLeftDaughter(CKDNode* inLeftDaughter)
      {m_pLeftDaughter = inLeftDaughter;}

   //right daughter
   CKDNode* GetRightDaughter() {return m_pRightDaughter;}
   void SetRightDaughter(CKDNode* inRightDaughter)
      {m_pRightDaughter = inRightDaughter;}

   //depth
   int GetDepth() const {return m_iDepth;}
   void SetDepth(const int& inDepth) {m_iDepth = inDepth;}

   //data
   template<class ForwardIterator>
   void SetData(ForwardIterator inBegin, ForwardIterator inEnd) {
      for(; inBegin != inEnd; ++inBegin) {
         m_dqData.push_back(*inBegin);
      }
   }
   const DataType& GetData() const {return m_dqData;}
   DataType& GetData() {return m_dqData;}

   //Hyper Plane
   const PointType& GetHyperPlane() const {
      return m_HyperPlane;
   }
   void SetHyperPlane(const PointType& inHyperPlane) {
      m_HyperPlane = inHyperPlane;
   }

private:
   //data
   CKDNode* m_pLeftDaughter;
   CKDNode* m_pRightDaughter;
   CKDNode* m_pParent;
   int m_iDepth;
   PointType m_HyperPlane;
   DataType m_dqData;
};

#endif //KDNODE
