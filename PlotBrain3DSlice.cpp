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

   int iColor;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;
   double dXMin = 1.e10;
   double dYMin = 1.e10;
   double dZMin = 1.e10;
   double dXMax = -1.e10;
   double dYMax = -1.e10;
   double dZMax = -1.e10;

   ifstream ifSegments("segments.txt");
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor) {
      if(dZ1*dZ2 < 0.) {
         TPolyLine3D * plTemp = new TPolyLine3D(2);
         plTemp->SetPoint(0, dX1, dY1, dZ1);
         plTemp->SetPoint(1, dX2, dY2, dZ2);
         plTemp->SetLineColor(1);
         plTemp->SetLineWidth(4);
         listSegments->Add(plTemp);

         dXMin = TMath::Min(dX1, TMath::Min(dX1, dXMin));
         dYMin = TMath::Min(dY1, TMath::Min(dY1, dYMin));
         dZMin = -.1;//TMath::Min(dZ1, TMath::Min(dZ1, dZMin));
         dXMax = TMath::Max(dX1, TMath::Max(dX1, dXMax));
         dYMax = TMath::Max(dY1, TMath::Max(dY1, dYMax));
         dZMax = .1;//TMath::Max(dZ1, TMath::Max(dZ1, dZMax));
      }
   }

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 2);

   TView * viewCanvas = new TView(1);
   viewCanvas->SetRange(dXMin - 0.01*TMath::Abs(dXMin),
      dYMin - 0.01*TMath::Abs(dYMin),
      dZMin - 0.01*TMath::Abs(dZMin),
      dXMax + 0.01*TMath::Abs(dXMax),
      dYMax + 0.01*TMath::Abs(dYMax),
      dZMax + 0.01*TMath::Abs(dZMax));
   viewCanvas->SetOutlineToCube();

   listSegments->Draw();

   pCanvas->Update();
   pCanvas->Draw();

}