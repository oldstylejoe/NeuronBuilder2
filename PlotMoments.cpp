{
#include <fstream>
#include <string>
#include <vector>

   int ipColor[] = {2,3,4,5,6,7};
   int ipMarker[] = {20,21,22,23,24,25};

   TMultiGraph * mgrData = new TMultiGraph("mgrData", "");

   TGraphErrors * greSlopes = new TGraphErrors();
   greSlopes->SetName("greSlopes");

   TGraphErrors * greRatios = new TGraphErrors();
   greRatios->SetName("greRatios");

   TF1 * fFit = new TF1("fFit", "[0]*x+[1]", 0., 100.);

   double dX, dY, dOldValue, dOldError;
   string strJunk;
   ofstream ofOut("out.txt");
   for(int i = 2; i <= 20; i += 2) {
      TGraph * grTemp = new TGraph();
      TString strID;
      strID += "grTemp";
      strID += i;
      grTemp->SetName(strID);
      grTemp->SetMarkerColor(ipColor[i%5]);
      grTemp->SetMarkerStyle(ipMarker[i%5]);
      fFit->SetLineColor(ipColor[i%5]);

      ifstream ifPM("product_moment_vary0.9.txt");
      while(ifPM >> strJunk >> dX) {
         for(int j = 1; j <= 20; ++j) {
            ifPM >> dY;
            if(i == j) {
               //if(strJunk[strJunk.size()-1] == 'l') {
                  grTemp->SetPoint(grTemp->GetN(), log(dX), log(dY));
                  //cout << "gh1 " << " " << dX << " " << dY << "\n" << flush;
               //}
            }
         }
      }
      ifPM.close();
      
      fFit->SetParameters(double(i), 0);
      grTemp->Fit(fFit, "q r e +");
      mgrData->Add(grTemp);

      greSlopes->SetPoint(greSlopes->GetN(), i, fFit->GetParameter(0));
      greSlopes->SetPointError(greSlopes->GetN()-1, 0, fFit->GetParError(0));

      //ofOut << i << " " << fFit->GetParameter(0) << " " << fFit->GetParError(0) << "\n";
      if(i > 3) {
         double dRatio = fFit->GetParameter(0)/dOldValue;
         double dRatioError = dRatio*sqrt(pow(dOldError/dOldValue,2)+
            pow(fFit->GetParError(0)/fFit->GetParameter(0),2));
         //greRatios->SetPoint(greRatios->GetN(), i-2, 2./(dRatio-1.)-(i-2.)-1.);
         //greRatios->SetPointError(greRatios->GetN()-1, 0, fabs(dRatioError/dRatio/dRatio));
         greRatios->SetPoint(greRatios->GetN(), i-2, dRatio);
         greRatios->SetPointError(greRatios->GetN()-1, 0, dRatioError);
         ofOut << i << " " << 2./(dRatio-1.)-(i-2.)-1. << " " << fabs(dRatioError/dRatio/dRatio) << "\n";
         //ofOut << i << " " << dRatio << " " << dRatioError << "\n";
         if(i == 4) {
            cout << "K4/K2 estimate = " 
               << 2./(dRatio-1.)-(i-2.)-1. << " +- "
               << fabs(dRatioError/dRatio/dRatio) << "\n" << flush;
         }
      }
      dOldValue = fFit->GetParameter(0);
      dOldError = fFit->GetParError(0);
   }
   ofOut.close();
   
   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "Moments", 1);

   mgrData->Draw("ap");

   pCanvas->Update();
   pCanvas->Draw();

   TCanvas * pCanvas2 = new TCanvas("pCanvas2", "Slopes", 1);

   greSlopes->SetMarkerStyle(8);
   TF1 * fFit2 = new TF1("fFit2", "[0]*x+1");
   greSlopes->Fit("fFit2", "q+");
   greSlopes->Draw("ap");

   cout << "Raw estimate (incorrect) "
      << 0.5*(1./(fFit2->GetParameter(0)-fFit2->GetParError(0))-1. + 
              1./(fFit2->GetParameter(0)+fFit2->GetParError(0))-1.) << " +- "
      << 0.5*fabs(1./(fFit2->GetParameter(0)-fFit2->GetParError(0)) -
              1./(fFit2->GetParameter(0)+fFit2->GetParError(0))) << "\n" << flush;

   pCanvas2->Update();
   pCanvas2->Draw();

   TCanvas * pCanvas2 = new TCanvas("pCanvas3", "Ratios", 1);

   greRatios->SetMarkerStyle(8);
   greRatios->Draw("ap");

   TF1 * fFit3 = new TF1("fFit3", "1.+2./([0]+x+1.)");
   fFit3->SetParameter(0,-1);
   greRatios->Fit(fFit3,"qe");

   pCanvas2->Update();
   pCanvas2->Draw();
}
