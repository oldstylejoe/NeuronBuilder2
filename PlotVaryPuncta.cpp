//Joe Snider
//3/07
//
//Fit some lines. Not general.
//Assumes pipes (cygwin)

#include <iostream>
#include <fstream>
#include <math.h>

#include "linear_regression.h"

int main(int inFileCount, char** inFiles) {
   double dX, dY, dM, dME, dB, dBE;

   if(inFileCount > 1) {
      for(int i = 1; i < inFileCount; ++i) {
         vector<double> vecX, vecY, vecYE;
         ifstream ifIn(inFiles[i]);
         if(ifIn.good()) {
            cerr << "Running " << inFiles[i] << "..." << flush;
            while(ifIn >> dX >> dY) {
               if(dX > 0. && dY > 0.) {
                  vecX.push_back(log(dX));
                  vecY.push_back(log(dY));
                  vecYE.push_back(0.);
               }
               LinearRegression(vecX, vecY, vecYE, dM, dME, dB, dBE);
               cout << dM-dME << " ";
               LinearRegression(vecY, vecX, vecYE, dM, dME, dB, dBE);
               cout << 1./(dM-dME) << "\n";
            }
            cerr << "done\n" << flush;
         } else {
            cerr << "Failed to read " << inFiles[i] << "...continuing\n" << flush;
         }
         ifIn.close();
      }
   }
   return 0;
}
