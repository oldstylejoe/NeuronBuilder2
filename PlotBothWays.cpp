{
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);

#include <fstream.h>
#include <string>

   TGraph * grXY = new TGraph();
   grXY->SetName("grXY");
   grXY->SetTitle("");

   TGraph * grYX = new TGraph();
   grYX->SetName("grYX");

   string strJunk;
   double dX, dY, dJunk;
   ifstream ifIn("total_lengths.txt");
   while(ifIn >> strJunk >> dY >> dJunk >> dX) {
      if(dX > 1. && dY > 1.) {
      //if(strJunk[strJunk.size()-1] == 'l') {
         grXY->SetPoint(grXY->GetN(), log10(dX), log10(dY));
         grYX->SetPoint(grYX->GetN(), log10(dY), log10(dX));
      //}
      }
   }
   ifIn.close();

   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   TF1 * fFitXY = new TF1("fFitXY", "1++x");
   grXY->Fit(fFitXY, "n");
   TF1 * fFitYX = new TF1("fFitYX", "1++x");
   grYX->Fit(fFitYX, "n");

   cout << "Slope is in the range [" 
      << fFitXY->GetParameter(1) - fFitXY->GetParError(1) << ", "
      << 1./(fFitYX->GetParameter(1) - fFitYX->GetParError(1)) << "]\n" << flush;
   double dBestM = sqrt( fFitXY->GetParameter(1) / fFitYX->GetParameter(1) );

   TF1 * fFitBest = new TF1("fFitBest", "x++1");
   fFitBest->FixParameter(0, dBestM);
   grXY->Fit(fFitBest, "+");

   grXY->SetMarkerStyle(20);
   grXY->Draw("ap");
   grXY->GetXaxis()->SetTitle("Area #left(log MC units^{2}#right)");
   grXY->GetXaxis()->CenterTitle();
   grXY->GetXaxis()->SetTitleSize(0.05);
   grXY->GetYaxis()->SetTitle("Length #left(log MC units#right)");
   grXY->GetYaxis()->CenterTitle();
   grXY->GetYaxis()->SetTitleSize(0.05);

   pCanvas->SetTickx();
   pCanvas->SetTicky();

   pCanvas->Update();
   pCanvas->Draw();

   cout << "Best estimate is " 
      << dBestM << " +- "
      << fFitXY->GetParError(1) << "\n" << flush;
   cout << "There are " << grXY->GetN() << " points.\n" << flush;
}
