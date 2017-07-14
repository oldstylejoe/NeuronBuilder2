{
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);

#include <fstream.h>
#include <string>

   TH1D * hArea = new TH1D("hArea", "", 100, 1, 4);
   TH1D * hLength = new TH1D("hLength", "", 100, 1, 4);

   string strJunk;
   double dX, dY, dJunk;
   ifstream ifIn("total_lengths.txt");
   while(ifIn >> strJunk >> dY >> dJunk >> dX) {
      if(dX > 0. && dY > 0.) {
         hArea->Fill(log10(dX));
         hLength->Fill(log10(dY));
      }
   }
   ifIn.close();

   TCanvas * pCanvas = new TCanvas("pCanvas", "", 800, 800);
   pCanvas->Divide(1,2);
   pCanvas->cd(1);
   hArea->Draw();

   pCanvas->cd(2);
   hLength->Draw();

   pCanvas->SetTickx();
   pCanvas->SetTicky();

   pCanvas->Update();
   pCanvas->Draw();
}
