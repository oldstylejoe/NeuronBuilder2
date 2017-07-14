//Joe Snider
//5/02
//
//uses Hoare's algorithm to find a median
//
//just a function
//
//generic stl style interface

#ifndef MEDIAN
#define MEDIAN

#include <deque>
#include <algorithm>
#include <vector>

using namespace std;

//IteratorType is an stl style iterator (forward)
//ValueType is the type of the data
//CompareType is a comparison operator with operator() defined (like std::map)
//modifies input variable
//
//finds the median of a set
//order nlog(n) algorithm so only use on small sets
//see HoareMedian below for an average O(n) algo
//
//for an even set it will return the left of the center 2
//  ie input of 1,2,3,4,5,6 returns 3
template <class IteratorType, class ValueType, class CompareType>
void Median(IteratorType inBegin,
            IteratorType inEnd,
            ValueType& inReturn,
            CompareType& inCompare) {
   vector<ValueType> vectorData;
   //copy the data
   for(; inBegin != inEnd; ++inBegin) {
      vectorData.push_back(*inBegin);
   }

   //sort the copy
   sort(vectorData.begin(), vectorData.end(), inCompare);

   //return the appropriate element (integer arithmatic is intentional)
   inReturn = vectorData[vectorData.size()/2];
}

//calls HoareElementFind of size/2
template <class IteratorType, class ValueType, class CompareType>
void HoareMedian(IteratorType inBegin,
                 IteratorType inEnd,
                 ValueType& inReturn,
                 CompareType& inCompare) {
   //have to count the elements
   IteratorType iter = inBegin;
   int iCount = 0;
   for(; iter != inEnd; ++iter, ++iCount);
   HoareElementFind(inBegin, inEnd, iCount/2, inReturn, inCompare);
}

//calls HoareElementFind of inSize/2
//input of inSize allowed for optimization
template <class IteratorType, class ValueType, class CompareType>
void HoareMedian(IteratorType inBegin,
                 IteratorType inEnd,
                 const unsigned& inSize,
                 ValueType& inReturn,
                 CompareType& inCompare) {
   //have to count the elements
   HoareElementFind(inBegin, inEnd, inSize/2, inReturn, inCompare);
}

//find the inElement smallest of the input set
template <class IteratorType, class ValueType, class CompareType>
void ElementFind(IteratorType inBegin,
                 IteratorType inEnd, 
                 size_t inElement,
                 ValueType& inReturn,
                 CompareType& inCompare) {
   vector<ValueType> vectorData;
   for(; inBegin != inEnd; ++inBegin) {
      vectorData.push_back(*inBegin);
   }
   sort(vectorData.begin(), vectorData.end());
   inReturn = vectorData[inElement];
}


//IteratorType is an stl style iterator (forward)
//ValueType is the type of the data, must have operator< defined
//modifies input variable
//
//finds the inElement smallest of a set using Hoare's 
//average O(n) algorithm
//
//for an even set it will return the left of the center 2
//  ie input of 1,2,3,4,5,6 returns 3
template <class IteratorType, class ValueType, class CompareType>
void HoareElementFind(IteratorType inBegin,
                      IteratorType inEnd,
                      size_t inElement,
                      ValueType& inReturn,
                      CompareType& inCompare) {
   //check for finished
   IteratorType itTemp = inBegin;
   //Joe Snider
   //6/06
   //Changed the following line.
   //It was unsafe (unsecure) to increment the iterator past the end.
   //if( ++(++(++(++(++(++itTemp))))) >= inEnd ) {
   if( (inEnd-inBegin) <= 6 ) {
	  //finished, set median and stop
      ElementFind(inBegin, inEnd, inElement, inReturn, inCompare);
   } else {
      //build left and right sets
      deque<ValueType> dqLeft;
      deque<ValueType> dqRight;
      ValueType vtCenter = *inBegin;
      ++inBegin;
      for(; inBegin != inEnd; ++inBegin) {
//Joe Snider
//4/04
//Fixed bug: default placement of median is now the left set.
//If the median appeared multiple times, it would be put in the wrong place.
//         if( inCompare(*inBegin, vtCenter) ) {
//            dqLeft.push_back(*inBegin);
//         } else {
//            dqRight.push_back(*inBegin);
//         }
         if( inCompare(vtCenter, *inBegin) ) {
            dqRight.push_back(*inBegin);
         } else {
            dqLeft.push_back(*inBegin);
         }
      }

      //now check if we got lucky with the center
      //or recurse on the larger deque
      size_t stLeftSize = dqLeft.size();
      //unsigned uRightSize = dqRight.size();
      if(inElement == stLeftSize) {
         //got lucky, done
         inReturn = vtCenter;
      } else if (inElement < stLeftSize) {
         //recurse appropriate direction
         HoareElementFind(dqLeft.begin(), 
                          dqLeft.end(), 
                          inElement, 
                          inReturn,
                          inCompare);
      } else {
         HoareElementFind(dqRight.begin(), 
                          dqRight.end(), 
                          inElement - (stLeftSize + 1), 
                          inReturn,
                          inCompare);
      }
   }
}

#endif //MEDIAN
