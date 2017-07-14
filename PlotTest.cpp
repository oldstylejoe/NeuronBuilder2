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

   TMultiGraph * mgrData = new TMultiGraph("mgrData", "");

   TGraph * grPunctaOccupied = new TGraph(100);
   grPunctaOccupied->SetName("grPunctaOccupied");
   TGraph * grPunctaEmpty = new TGraph(100);
   grPunctaEmpty->SetName("grPunctaEmpty");
   ifstream ifPuncta("puncta.txt");
   double dX, dY;
   int iID;
   int iCountOccupied = 0;
   int iCountEmpty = 0;
   while(ifPuncta >> dX >> dY >> iID) {
      if (grPunctaOccupied->GetN() <= iCountOccupied) {grPunctaOccupied->Set(2*iCountOccupied);}
      if (grPunctaEmpty->GetN() <= iCountEmpty) {grPunctaEmpty->Set(2*iCountEmpty);}
      if(iID == 0) {
         grPunctaEmpty->SetPoint(iCountEmpty, dX, dY);
         ++iCountEmpty;
      }
      if(iID == 1) {
         grPunctaOccupied->SetPoint(iCountOccupied, dX, dY);
         ++iCountOccupied;
      }
   }
   grPunctaEmpty->Set(iCountEmpty);
   grPunctaOccupied->Set(iCountOccupied);
   ifPuncta.close();

   grPunctaEmpty->SetMarkerColor(2);
   grPunctaEmpty->SetMarkerStyle(8);
   grPunctaOccupied->SetMarkerColor(3);
   grPunctaOccupied->SetMarkerStyle(8);

   mgrData->Add(grPunctaEmpty, "p");
   mgrData->Add(grPunctaOccupied, "p");

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 2);

   mgrData->Draw("a");

   TLine lineDrawer;
   ifstream ifArbor("segments.txt");
   double dX1, dY1, dX2, dY2;
   while( ifArbor >> dX1 >> dY1 >> dX2 >> dY2 ) {
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
   }
   ifArbor.close();

   pCanvas->Update();
   pCanvas->Draw();

}