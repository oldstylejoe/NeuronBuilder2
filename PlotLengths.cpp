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
   int iCountData = 0;
   double dpData[3000];  //known amount

   //get the data
   TH1D * h1dData = new TH1D("h1dData", "", 100, 0., 1000.);
   ifstream ifExp("lengths_all_fish.txt");
   while(ifExp >> dX) {
	   if(dX > 1.e-12) {
		   //h1dData->Fill(log(dX)/log(10.));
		   h1dData->Fill(dX);
         if(iCountData < 3000) {
            dpData[iCountData] = dX;
            ++iCountData;
         }
	   }
   }
   ifExp.close();
   h1dData->Sumw2();
   h1dData->Scale(1./h1dData->Integral());
   h1dData->SetMarkerStyle(20);
   h1dData->SetMarkerColor(2);
   h1dData->SetLineColor(2);

   double dShift = 0.;

   //get the simulation data
   int iCountSim = 0;
   TH1D * h1dSim = new TH1D("h1dSim", "", 100, 0., 1000.);
   ifstream ifSim("lengths_vary1.5.txt");
   while(ifSim >> dX) {
	   if(dX > 1.e-12) {
		   h1dSim->Fill(9.*dX);
		   //h1dSim->Fill(log(dX)/log(10.)+dShift);
		   //h1dSim->Fill(dX*pow(10., dShift));
         ++iCountSim;
	   }
   }
   ifSim.close();
   h1dSim->Sumw2();
   h1dSim->Scale(1./h1dSim->Integral());
   h1dSim->SetMarkerStyle(21);
   h1dSim->SetMarkerColor(4);
   h1dSim->SetLineColor(4);

   /*double * dpSim = new double[iCountSim];
   iCountSim = 0;
   ifstream ifSim2("raw_length.txt");
   while(ifSim2 >> dX) {
	   if(dX > 1.e-12) {
         dpSim[iCountSim] = dX;
         ++iCountSim;
	   }
   }
   ifSim2.close();*/

   THStack * hsAll = new THStack("hAll", "");
   hsAll->Add(h1dData);
   hsAll->Add(h1dSim);

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   hsAll->Draw("nostack");
   hsAll->GetXaxis()->SetTitle("Lengths (Log scale #mum)");
   hsAll->GetXaxis()->SetTitleColor(1);
   hsAll->GetXaxis()->CenterTitle();
   hsAll->GetYaxis()->SetTitle("Probability");
   hsAll->GetYaxis()->CenterTitle();

   pCanvas->Update();
   pCanvas->Draw();

   TLegend * leg = new TLegend(0.5, 0.5, 0.6, 0.6);
   leg->AddEntry(h1dData, "Fish RT projections");
   leg->AddEntry(h1dSim, "Model");
   leg->Draw();

   /*TCanvas * pCanvas2 = new TCanvas("pCanvas2", "", 1);

   TGraphQQ * grQQ = new TGraphQQ(iCountData, dpData, iCountSim, dpSim);
   grQQ->Draw("ap");

   pCanvas2->Update();
   pCanvas2->Draw();*/

}
