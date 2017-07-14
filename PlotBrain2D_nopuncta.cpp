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

   TCanvas * pCanvas = new TCanvas("pCanvas", "title", 800, 800);

   char fileName[] = "segments_vary1.5.txt";
   double lengthScale = 40.;
   TH2D * axis = new TH2D("axis", fileName, 10, -1.7*lengthScale, 1.7*lengthScale,
                          10, -1.7*lengthScale, 1.7*lengthScale);
   axis->Draw();
   double xOffset[9];
   double yOffset[9];
   int count = 0;
   for(int i = -1; i <= 1; ++i) {
     for(int j = -1; j <= 1; ++j) {
       xOffset[count] = 0.;//1.1*double(i)*lengthScale;
       yOffset[count] = 0.;//1.1*double(j)*lengthScale;
       ++count;
     }
   }

   TLine lineDrawer;
   lineDrawer.SetLineWidth(2);
   TEllipse ellipseDrawer;
   int iColor, iSize;
   double dX1, dY1, dX2, dY2, dJunk;
   /*ifstream ifArborTest("segments_test.txt");
   while( ifArborTest >> dX1 >> dY1 >> dX2 >> dY2 >> iColor >> iSize) {
      lineDrawer.SetLineWidth(max(8-iSize, 1));
      lineDrawer.SetLineColor(3+((iColor+6)%6));
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
   }
   ifArborTest.close();
   */
   ifstream ifArbor(fileName);
   int iOldColor = -1;
   //while( ifArbor >> dX1 >> dY1 >> dJunk >> dX2 >> dY2 >> dJunk >> iColor >> iSize) {
   while( ifArbor >> dX1 >> dY1 >> dX2 >> dY2 >> iColor >> iSize) {
     //iColor -= 10;
     if(0 <= iColor && iColor < 9) {
     lineDrawer.SetLineColor(colorChoice[iSize%8]);
     //lineDrawer.SetLineColor(2+iSize);
     //if( (dX1-dX2)*(dX1-dX2) + (dY1-dY2)*(dY1-dY2) < 0.5 ) {
     //   lineDrawer.SetLineColor(1);
     //} else {
     //   lineDrawer.SetLineColor(2+iSize);
     //}
     lineDrawer.DrawLine(xOffset[iColor] + dX1, yOffset[iColor] + dY1, xOffset[iColor] + dX2, yOffset[iColor] + dY2);
     if( iOldColor != iColor ) {
        iOldColor = iColor;
        ellipseDrawer.SetFillColor(1);//2+((iColor+6)%6));
        ellipseDrawer.DrawEllipse(xOffset[iColor] + dX1, yOffset[iColor] + dY1,
                                  .02*lengthScale, .02*lengthScale, 0., 360., 0.);
     }
     }
   }
   ifArbor.close();

   pCanvas->Update();
   pCanvas->Draw();

}