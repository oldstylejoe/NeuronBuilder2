//Joe Snider
//7/17
//
//Interface for generating neurons using the growth model.
//
//Not super general, but call like 
//    <executable> <removal radius> 
//Dumps out an swc file to stdout that you can pipe where you like.
//Units of removal readius are in terms of the trial length (1).
//
//Change other stuff by recompiling.

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "brain.h"
#include "arbor.h"
#include "rand3.h"
#include "histogram.h"
#include "ForEachUtils.h"

//a length scale for the brain, scales the boundary, puncta, and number of arbors
double LENGTH_SCALE = 20.;

//set the bounding box
double X_LOW = -1*LENGTH_SCALE;
double Y_LOW = -1*LENGTH_SCALE;
const double Z_LOW = -1*LENGTH_SCALE;//0.;
double X_HIGH = LENGTH_SCALE;
double Y_HIGH = LENGTH_SCALE;
const double Z_HIGH = LENGTH_SCALE;//3.;

double BRANCH_LENGTH = 1.;          //length of new branches
const double PUNCTA_RANGE = .05;           //range over which to consider a puncta intersecting a branch
const double PUNCTA_REMOVE_RADIUS = 1.3;   //remove all puncta within this range of all puncta that intersect a branch
const double PUNCTA_PER_BRANCH = 1.0;      //number of puncta per branch

const int PUNCTA = int(4.*LENGTH_SCALE*LENGTH_SCALE/(BRANCH_LENGTH*PUNCTA_RANGE)*PUNCTA_PER_BRANCH+0.5);
//const int PUNCTA = int(4.*(Z_HIGH-Z_LOW)*LENGTH_SCALE*LENGTH_SCALE/(acos(-1.)*BRANCH_LENGTH*PUNCTA_RANGE*PUNCTA_RANGE)*PUNCTA_PER_BRANCH+0.5);

const int MAIN_MAX_STEPS_WITHOUT_PLACEMENT = 1000000;

const int NUM_ARBORS = 1;//int(LENGTH_SCALE*LENGTH_SCALE/10.+0.5);            //the number of arbors to place

const double STEP_FORWARD_PROB = 1.; //the probablilty the a trial growth step is straight ahead

const int MIN_REPEATS = 0;
const int STEP_REPEATS = 1;
const int MAX_REPEATS = 100;

const int REQUIRED_STABILIZERS = 1;

using namespace std;

int main(int argc, char* argv[]) {
   unsigned uStartTime = (unsigned)time(NULL);

   srand( (unsigned) time(NULL) );
   CRand3 rand3(rand());

   if(argc != 2) {
      //cout << "usage: ./neuronbuilder2 <removal radius (double)>\nThe pipe the resulting output to a .swc file\n" << flush;
      cout << "usage: ./neuronbuilder2 <removal radius (double)>\nThen pipe the resulting output list of segments to a .txt file\n" << flush;
      return -1;
   }

   double removeRadius = atof(argv[1]);
   cerr << "Using removal radius: " << removeRadius << "\n";

   CPoint pntLow;
   pntLow.push_back(X_LOW);
   pntLow.push_back(Y_LOW);
   pntLow.push_back(Z_LOW);
   CPoint pntHigh;
   pntHigh.push_back(X_HIGH);
   pntHigh.push_back(Y_HIGH);
   pntHigh.push_back(Z_HIGH);
   //int iPunctaCount = int(4.*LENGTH_SCALE*LENGTH_SCALE/(BRANCH_LENGTH*PUNCTA_RANGE)*PUNCTA_PER_BRANCH+0.5);
   int iPunctaCount = int(4.*(Z_HIGH-Z_LOW)*LENGTH_SCALE*LENGTH_SCALE/(acos(-1.)*BRANCH_LENGTH*PUNCTA_RANGE*PUNCTA_RANGE)*PUNCTA_PER_BRANCH+0.5);
   cerr << "   Placing " << iPunctaCount << " puncta ... " << flush;
   CBrain B1(NUM_ARBORS, iPunctaCount, pntLow, pntHigh, STEP_FORWARD_PROB, rand3);
   cerr << "done\n" << flush;
   cerr << "Growing arbor ... " << flush;
   int iPlaced = B1.Grow(MAIN_MAX_STEPS_WITHOUT_PLACEMENT, BRANCH_LENGTH, PUNCTA_RANGE, removeRadius, rand3);
   cerr << "done. Consumed " << 100.0*float(iPlaced)/float(iPunctaCount) << "\% of puncta.\n" << flush;
   cerr << "Trimming danglers ... " << flush;
   B1.Trim(PUNCTA_RANGE, REQUIRED_STABILIZERS);
   cerr << "done\n" << flush;
   //B1.DumpSWC();
   B1.DumpSegments(cout, 0, "x1 y1 z1 x2 y2 z2 arbor depth\n");

   //dump the puncta for testing
   ofstream ofPuncta("puncta.txt");
   B1.DumpPuncta(ofPuncta);
   ofPuncta.close();

   cerr << "This took " << time(NULL)-uStartTime << " seconds\n" << flush;

   return 0;
}

