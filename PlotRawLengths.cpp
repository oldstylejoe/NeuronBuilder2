{
#include <fstream.h>

   TH1D * h1dData = new TH1D("h1dData", "", 200, 0., 2.);

   double dX;
   ifstream ifIn("raw_length_br50000_trial0.txt");
   while(ifIn >> dX) {
      if(dX > 1.e-12) {
         h1dData->Fill(log(dX));
      }
   }
   ifIn.close();

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   h1dData->Draw();

   pCanvas->Update();
   pCanvas->Draw();
}
