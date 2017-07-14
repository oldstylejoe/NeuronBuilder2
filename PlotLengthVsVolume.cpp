//Joe Snider

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

   TMultiGraph * mgrData = new TMultiGraph("mgrData", "");

   TGraph * grData = new TGraph();
   grData->SetName("grData");

   ifstream ifData("total_lengths.txt");
   TString strType;
   double dX, dY, dJunk;
   while(strType.ReadToDelim(ifData, ' ')) {
      ifData >> dY >> dX >> dJunk;
	  grData->SetPoint(grData->GetN(), log(dX)/log(10.), log(dY)/log(10.));
   }
   ifData.close();

   grData->SetMarkerColor(2);
   grData->SetMarkerStyle(20);
   grData->SetLineColor(2);

   double dShift = log10(9.);

   TGraph * grSim = new TGraph();
   grSim->SetName("grSim");

   ifstream ifSim("length_volume_vary1.6.txt");
   TString strType;
   while(ifSim >> dY >> dX) {
	  grSim->SetPoint(grSim->GetN(), log(dX)/log(10.)+2.*dShift, log(dY)/log(10.)+dShift);
   }
   ifSim.close();

   grSim->SetMarkerColor(4);
   grSim->SetMarkerStyle(21);
   grSim->SetLineColor(4);

   mgrData->Add(grData, "p");
   mgrData->Add(grSim, "p");

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   mgrData->Draw("a");

   mgrData->GetXaxis()->SetTitle("Area (log scale #mu^{2})");
   mgrData->GetXaxis()->SetTitleColor(1);
   mgrData->GetXaxis()->CenterTitle(true);

   mgrData->GetYaxis()->SetTitle("Length (log scale #mu)");
   mgrData->GetYaxis()->CenterTitle(true);

   TF1 * fTheory = new TF1("fTheory", "0.75*x+.6", 0., 20.);
   fTheory->Draw("same");

   //TLegend * legData = new TLegend(0.5, 0.5, 0.8, 0.8);
   //legData->AddEntry(grAxon, "Axon", "lp");
   //legData->AddEntry(grDendrite, "Dendrite", "lp");
   //legData->AddEntry(grApical, "Apical", "lp");
   //legData->Draw();

   pCanvas->Update();
   pCanvas->Draw();
}
