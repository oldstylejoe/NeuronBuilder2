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

   int colorChoice[] = {1,2,3,4,5,6,7,8};
   TMultiGraph * mgrData = new TMultiGraph("mgrData", "");

   TGraph * grPunctaOccupied = new TGraph(100);
   grPunctaOccupied->SetName("grPunctaOccupied");
   TGraph * grPunctaRemoved = new TGraph(100);
   grPunctaRemoved->SetName("grPunctaRemoved");
   TGraph * grPunctaEmpty = new TGraph(100);
   grPunctaEmpty->SetName("grPunctaEmpty");
   ifstream ifPuncta("puncta.txt");
   double dX, dY, str;
   int iID;
   int iCountEmpty = 0;
   int iCountOccupied = 0;
   int iCountRemoved = 0;
   while(ifPuncta >> dX >> dY >> iID >> str) {
      if(iID == 0) {
         if(iCountEmpty >= grPunctaEmpty->GetN()) {grPunctaEmpty->Set(2*iCountEmpty);}
         grPunctaEmpty->SetPoint(iCountEmpty, dX, dY);
         ++iCountEmpty;
      }
      if(iID == 1) {
         if(str > 0.5) {
            if(iCountOccupied >= grPunctaOccupied->GetN()) {grPunctaOccupied->Set(2*iCountOccupied);}
            grPunctaOccupied->SetPoint(iCountOccupied, dX, dY);
            ++iCountOccupied;
         } else {
            if(iCountRemoved >= grPunctaRemoved->GetN()) {grPunctaRemoved->Set(2*iCountRemoved);}
            grPunctaRemoved->SetPoint(iCountRemoved, dX, dY);
            ++iCountRemoved;
         }
      }
   }
   ifPuncta.close();
   grPunctaEmpty->Set(iCountEmpty);
   grPunctaOccupied->Set(iCountOccupied);
   grPunctaRemoved->Set(iCountRemoved);

   grPunctaEmpty->SetMarkerColor(19);
   grPunctaEmpty->SetMarkerStyle(8);
   grPunctaOccupied->SetMarkerColor(2);
   grPunctaOccupied->SetMarkerStyle(21);
   grPunctaRemoved->SetMarkerColor(20);
   grPunctaRemoved->SetMarkerStyle(8);

   mgrData->Add(grPunctaEmpty, "p");
   mgrData->Add(grPunctaOccupied, "p");
   mgrData->Add(grPunctaRemoved, "p");

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 800, 800);

   mgrData->Draw("a");

   TLine lineDrawer;
   lineDrawer.SetLineWidth(2);
   TEllipse ellipseDrawer;
   int iColor, iSize;
   double dX1, dY1, dX2, dY2;
   ifstream ifArborTest("segments.txt");
   while( ifArborTest >> dX1 >> dY1 >> dX2 >> dY2 >> iColor >> iSize) {
      //lineDrawer.SetLineWidth(max(8-iSize, 1));
      lineDrawer.SetLineColor(1);
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
   }
   ifArborTest.close();
   
   ifstream ifArbor("segments_vary1.txt");
   int iOldColor = -1;
   while( ifArbor >> dX1 >> dY1 >> dX2 >> dY2 >> iColor >> iSize) {
      //lineDrawer.SetLineWidth(max(8-iSize, 1));
      lineDrawer.SetLineColor(colorChoice[iSize%8]);
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
      if( iOldColor != iColor ) {
         iOldColor = iColor;
         ellipseDrawer.SetFillColor(2+((iColor+6)%6));
         ellipseDrawer.DrawEllipse(dX1, dY1, 1., 1., 0., 360., 0.);
      }
   }
   ifArbor.close();

   pCanvas->Update();
   pCanvas->Draw();

}