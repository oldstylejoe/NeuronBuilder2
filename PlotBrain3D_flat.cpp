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

   TList * listSegments = new TList();
   TPolyMarker * pmCellBodies = new TPolyMarker();
   pmCellBodies->SetMarkerStyle(8);

   int iColor, iLevel;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;
   double dXMin = 1.e10;
   double dYMin = 1.e10;
   double dZMin = 1.e10;
   double dXMax = -1.e10;
   double dYMax = -1.e10;
   double dZMax = -1.e10;

   ifstream ifSegments("segments_vary1.25.txt");
   int iOldColor = -1;
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor >> iLevel) {
if(iColor == 0) {
      if(iOldColor != iColor) {
         iOldColor = iColor;
         pmCellBodies->SetPoint(iColor, dX1, dY1);
      }
      TLine * lineTemp = new TLine(dX1, dY1, dX2, dY2);
      lineTemp->SetLineColor(2+(iColor%6));
      if(iColor == 0) {
         lineTemp->SetLineColor(1);
         lineTemp->SetLineWidth(3);
      }
      listSegments->Add(lineTemp);

      dXMin = TMath::Min(dX1, TMath::Min(dX1, dXMin));
      dYMin = TMath::Min(dY1, TMath::Min(dY1, dYMin));
      dZMin = TMath::Min(dZ1, TMath::Min(dZ1, dZMin));
      dXMax = TMath::Max(dX1, TMath::Max(dX1, dXMax));
      dYMax = TMath::Max(dY1, TMath::Max(dY1, dYMax));
      dZMax = TMath::Max(dZ1, TMath::Max(dZ1, dZMax));
}
   }

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 800, 800);

   //TH2D * h2dAxis = new TH2D("h2dAxis", "", 10, 
   //   dXMin - 0.01*TMath::Abs(dXMin),
   //   dXMax + 0.01*TMath::Abs(dXMax),
   //   10,
   //   dYMin - 0.01*TMath::Abs(dYMin),
   //   dYMax + 0.01*TMath::Abs(dYMax));
   TH2D * h2dAxis = new TH2D("h2dAxis", "", 10, -6, 6, 10, -6, 6);
   h2dAxis->Draw();

   listSegments->Draw();
   pmCellBodies->Draw("same");

   pCanvas->Update();
   pCanvas->Draw();

}