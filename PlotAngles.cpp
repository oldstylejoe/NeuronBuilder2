{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(2);
   gStyle->SetPadColor(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetTitleColor(0);
   gStyle->SetStatColor(0);

   #include <fstream.h>
   #include <vector.h>

   double dX;

   //get the data
   TH1D * h1dData = new TH1D("h1dData", "", 20, 0., 4.);
   ifstream ifExp("angles_all_fish.txt");
   while(ifExp >> dX) {
	   if(dX > 1.e-12) {
		   h1dData->Fill(dX);
	   }
   }
   ifExp.close();
   h1dData->Sumw2();
   h1dData->Scale(1./h1dData->Integral());
   h1dData->SetMarkerStyle(20);
   h1dData->SetMarkerColor(2);
   h1dData->SetLineColor(2);

   //get the simulation data
   TH1D * h1dSim = new TH1D("h1dSim", "", 20, 0., 4.);
   ifstream ifSim("angles_vary0.5.txt");
   while(ifSim >> dX) {
	   if(dX > 1.e-12) {
		   h1dSim->Fill(dX);
	   }
   }
   ifSim.close();
   h1dSim->Sumw2();
   h1dSim->Scale(1./h1dSim->Integral());
   h1dSim->SetMarkerStyle(21);
   h1dSim->SetMarkerColor(4);
   h1dSim->SetLineColor(4);

   THStack * hsAll = new THStack("hAll", "");
   hsAll->Add(h1dData);
   hsAll->Add(h1dSim);

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   hsAll->Draw("nostack");
   hsAll->GetXaxis()->SetTitle("Angle (radians)");
   hsAll->GetXaxis()->SetTitleColor(1);
   hsAll->GetXaxis()->CenterTitle();
   hsAll->GetYaxis()->SetTitle("Probability");
   hsAll->GetYaxis()->CenterTitle();

   pCanvas->Update();
   pCanvas->Draw();

}
