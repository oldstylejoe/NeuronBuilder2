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

/*   TGraph2D * grPunctaOccupied = new TGraph2D();
   grPunctaOccupied->SetName("grPunctaOccupied");
   TGraph2D * grPunctaEmpty = new TGraph2D();
   grPunctaEmpty->SetName("grPunctaEmpty");
   ifstream ifPuncta("puncta.txt");
   double dX, dY, dZ;
   int iID;
   int iCountEmpty = 0;
   int iCountOccupied = 0;
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
   //grPunctaEmpty->SetMarkerStyle(8);
   grPunctaOccupied->SetMarkerColor(3);
   //grPunctaOccupied->SetMarkerStyle(8);
*/
   TList * listSegments = new TList();
   TPolyMarker3D * pm3CellBodies = new TPolyMarker3D();
   pm3CellBodies->SetMarkerStyle(8);

   int iColor, iLevel;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;
   double dXMin = 1.e10;
   double dYMin = 1.e10;
   double dZMin = 1.e10;
   double dXMax = -1.e10;
   double dYMax = -1.e10;
   double dZMax = -1.e10;

   ifstream ifSegments("segments_vary1.txt");
   int iOldColor = -1;
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor >> iLevel) {
//if((iColor%70) == 0) {
      if(iOldColor != iColor) {
         iOldColor = iColor;
         pm3CellBodies->SetPoint(iColor, dX1, dY1, dZ1);
      }
      TPolyLine3D * plTemp = new TPolyLine3D(2);
      plTemp->SetPoint(0, dX1, dY1, dZ1);
      plTemp->SetPoint(1, dX2, dY2, dZ2);
      plTemp->SetLineColor(2+(iColor%6));
      if(iColor == 2) {
         plTemp->SetLineColor(1);
         plTemp->SetLineWidth(3);
      }
      listSegments->Add(plTemp);

      dXMin = TMath::Min(dX1, TMath::Min(dX1, dXMin));
      dYMin = TMath::Min(dY1, TMath::Min(dY1, dYMin));
      dZMin = TMath::Min(dZ1, TMath::Min(dZ1, dZMin));
      dXMax = TMath::Max(dX1, TMath::Max(dX1, dXMax));
      dYMax = TMath::Max(dY1, TMath::Max(dY1, dYMax));
      dZMax = TMath::Max(dZ1, TMath::Max(dZ1, dZMax));
//}
   }

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 800, 800);

   TView3D * viewCanvas = new TView3D();
   viewCanvas->SetRange(dXMin - 0.01*TMath::Abs(dXMin),
      dYMin - 0.01*TMath::Abs(dYMin),
      dZMin - 0.01*TMath::Abs(dZMin),
      dXMax + 0.01*TMath::Abs(dXMax),
      dYMax + 0.01*TMath::Abs(dYMax),
      dZMax + 0.01*TMath::Abs(dZMax));
   //viewCanvas->SetOutlineToCube();

   listSegments->Draw();
   pm3CellBodies->Draw("same");

/*   grPunctaEmpty->Draw("psame");
   grPunctaOccupied->Draw("psame");
*/
//   pCanvas->Update();
//   pCanvas->Draw();

}