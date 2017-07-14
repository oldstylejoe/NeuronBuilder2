{
   gROOT->Reset();
   //gStyle->SetOptStat(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(2);
   gStyle->SetPadColor(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetTitleColor(0);
   //gStyle->SetStatColor(0);

   #include <iostream.h>
   #include <fstream.h>
   #include <vector.h>

   double dX, dY, dZ;
   int iID;
   int iCountEmpty = 0;
   int iCountOccupied = 0;

   TGraph2D * grPunctaOccupied = new TGraph2D();
   grPunctaOccupied->SetName("grPunctaOccupied");
   TGraph2D * grPunctaEmpty = new TGraph2D();
   grPunctaEmpty->SetName("grPunctaEmpty");
   ifstream ifPuncta("puncta.txt");
   while(ifPuncta >> dX >> dY >> dZ >> iID) {
      if(iID == 0) {
         grPunctaEmpty->SetPoint(iCountEmpty, dX, dY, dZ);
         ++iCountEmpty;
      }
      if(iID == 1) {
         grPunctaOccupied->SetPoint(iCountOccupied, dX, dY, dZ);
         ++iCountOccupied;
      }
   }
   ifPuncta.close();

   grPunctaEmpty->SetMarkerColor(2);
   grPunctaEmpty->SetMarkerStyle(8);
   grPunctaOccupied->SetMarkerColor(3);
   grPunctaOccupied->SetMarkerStyle(8);

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 2);

   grPunctaEmpty->Draw("p");
   grPunctaOccupied->Draw("psame");

   int iColor;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;

   ifstream ifSegments("segments.txt");
   TPolyLine3D plTemp(2);
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor) {
      plTemp.SetPoint(0, dX1, dY1, dZ1);
      plTemp.SetPoint(1, dX2, dY2, dZ2);
      plTemp.SetLineColor(2+(iColor%6));
      plTemp.Draw();
   }

   pCanvas->Update();
   pCanvas->Draw();

}