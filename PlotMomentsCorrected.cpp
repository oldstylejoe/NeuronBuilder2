{
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);

#include <fstream>
#include <string>
#include <vector>

   int ipColor[] = {2,3,4,5,6,7};
   int ipMarker[] = {20,21,22,23,24,25};

   TMultiGraph * mgrData = new TMultiGraph("mgrData", "");

   TGraphErrors * greSlopes = new TGraphErrors();
   greSlopes->SetName("greSlopes");

   TF1 * fFit = new TF1("fFit", "x++1", -100., 100.);

   double dpMY[22];  //add more space if needed
   double dpOY[22];

   double dL, dX, dY, dOldValue, dOldError;
   string strJunk;
   ofstream ofOut("out.txt");
   for(int i = 4; i <= 20; i += 2) {
      TGraph * grXY = new TGraph();
      TString strID;
      strID += "grXY";
      strID += i;
      grXY->SetName(strID);
      grXY->SetMarkerColor(ipColor[i%5]);
      grXY->SetMarkerStyle(ipMarker[i%5]);
      fFit->SetLineColor(ipColor[i%5]);
      TGraph * grYX = new TGraph();
      TString strID;
      strID = "grYX";
      strID += i;
      grYX->SetName(strID);
      grYX->SetMarkerColor(ipColor[i%5]);
      grYX->SetMarkerStyle(ipMarker[i%5]);
      fFit->SetLineColor(ipColor[i%5]);

      ifstream ifPM("product_moment_merge.txt");
      while(ifPM >> strJunk >> dL >> dX >> dX) {
         for(int j = 3; j <= 20; ++j) {
            ifPM >> dY;
            if(i == j) {
               grXY->SetPoint(grXY->GetN(), log10(dX)-log10(dL), log10(dY));
               grYX->SetPoint(grYX->GetN(), log10(dY), log10(dX)-log10(dL));
            }
         }
      }
      ifPM.close();
      
      fFit->SetParameters(double(i), 0);
      grXY->Fit(fFit, "r e +");
      double dMX = fFit->GetParameter(0);
      double dDMX = fFit->GetParError(0);
      grYX->Fit(fFit, "r e +");
      double dMY = fFit->GetParameter(0);
      double dDMY = fFit->GetParError(0);
      mgrData->Add(grXY);
      dpMY[i] = dMY;
      dpOY[i] = fFit->GetParameter(1);
      //mgrData->Add(grYX);

      //double dSlopeLow = dMX-dDMX;
      //double dSlopeHigh = 1./(dMY-dDMY);
      double dSlopeLow = sqrt( (dMX-dDMX) / (dMY+dDMY) );
      double dSlopeHigh = sqrt( (dMX+dDMX) / (dMY-dDMY) );
      double dV = sqrt(dMX/dMY);//0.5*(dSlopeLow+dSlopeHigh);
      double dVE = dDMX;//0.5*(dSlopeHigh-dSlopeLow);

      greSlopes->SetPoint(greSlopes->GetN(), i, dV-0.5*double(i));
      greSlopes->SetPointError(greSlopes->GetN()-1, 0, dVE);

   }
   ofOut.close();
   
   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "Moments", 1);

   mgrData->Draw("ap");

   TF1 * fDrawer = new TF1("fDrawer", "[0]*x+[1]", 0., 1000.);
   for(int i = 4; i <= 20; i += 2) {
      fDrawer->SetLineColor(ipColor[i%5]);
      fDrawer->SetParameters(1./dpMY[i], -1.*dpOY[i]/dpMY[i]);
      fDrawer->DrawCopy("same");
   }

   mgrData->GetXaxis()->SetTitle("#frac{#LT2#GT}{#LT0#GT} #left(log #mu^{4}#right)");
   mgrData->GetXaxis()->SetTitleSize(0.05);
   mgrData->GetXaxis()->CenterTitle();
   mgrData->GetYaxis()->SetTitle("#LTn#GT #left(log #mum^{2n+3}#right)");
   mgrData->GetYaxis()->SetTitleSize(0.05);
   mgrData->GetYaxis()->CenterTitle();

   pCanvas->SetTickx();
   pCanvas->SetTicky();

   pCanvas->Update();
   pCanvas->Draw();

   TCanvas * pCanvas2 = new TCanvas("pCanvas2", "Slopes", 1);

   greSlopes->SetMarkerStyle(8);
   TF1 * fFit2 = new TF1("fFit2", "[0]");
   greSlopes->Fit("fFit2", "+");
   greSlopes->Draw("ap");

   greSlopes->GetXaxis()->SetTitle("Moment number (unitless)");
   greSlopes->GetXaxis()->SetTitleSize(0.05);
   greSlopes->GetXaxis()->CenterTitle();
   greSlopes->GetYaxis()->SetTitle("Slope (unitless)");
   greSlopes->GetYaxis()->SetTitleSize(0.05);
   greSlopes->GetYaxis()->CenterTitle();

   pCanvas2->SetTickx();
   pCanvas2->SetTicky();

   pCanvas2->Update();
   pCanvas2->Draw();

   cout << "The corrected (unbiased if there's no cutoff) scaling exponent is :\n     "
      << 4.*fFit2->GetParameter(0)-2. << " +- "
      << 4.*fFit2->GetParError(0) << "\n" << flush;

}
