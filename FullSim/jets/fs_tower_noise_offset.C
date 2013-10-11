//ROOT headers
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TSpectrum.h>
#include <TF1.h>

//STL headers
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

//#define maxHDe 3
#define maxHDe 1
#define maxHDlumi 8
//#define maxHDlumi 2
#define maxHDeta 13
#define maxPrint 2
#define M_PI 3.14159265358979323846

//energy values - global
//Double_t energies[] = {30, 50, 100};
Double_t energies[] = {30};
Double_t lumis[] = {0, 100, 200, 300, 400, 500, 600, 700};
//Double_t lumis[] = {0, 500};
Double_t etas[] = {1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
Int_t year[] = {2017, 2019};

//NB: in this file, "energies" refers to pT

//-------------------------------------------------------
//function to get mean noise from calotowers at given eta
Double_t get_noise(int num, int Lnum, int etabin, bool do_2019, bool do_show, unsigned do_print=0){
	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		return 0;
	}

	//make filenames
	std::stringstream fname, jname, luminame;

	fname << "tree_towernoise_" << year[do_2019] << "_" << energies[num] << "_lumi" << lumis[Lnum] << ".root";
	jname << etas[etabin] << " < #eta < " << etas[etabin+1];
	luminame << "lumi = " << lumis[Lnum] << " fb^{-1}";

	//open file and histo
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TH1F* h_res;
	TTree* totalTree = (TTree*)_file->Get("Total");

	std::stringstream drawname, cutname;
	drawname << "TowerNoisePt>>htemp";
	cutname << "abs(GenJetEta)>" << etas[etabin] << " && abs(GenJetEta)<=" << etas[etabin+1];
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	h_res = (TH1F*)gDirectory->Get("htemp");

	//plotting variables
	TCanvas* can;
	TPaveText* pave;

	//get values from histo
	Double_t m = h_res->GetMean();
	Double_t me = h_res->GetMeanError();
	//Double_t m = h_res->GetBinCenter(h_res->GetMaximumBin()); //peak
	Double_t s = h_res->GetRMS();
	Double_t se = h_res->GetRMSError();
	Int_t N = h_res->GetEntries();

	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};

	//plotting
	if (do_show){
		can = new TCanvas("noise","noise",700,500);
		can->cd();	
		
		//plot histo and fit
		h_res->SetStats(kTRUE);
		gStyle->SetOptStat("emr");
		h_res->SetTitle("");
		h_res->GetXaxis()->SetTitle("CaloTower Noise p_{T} [GeV]");
		h_res->GetYaxis()->SetTitle("");
		h_res->SetLineColor(kBlack);
		h_res->Draw("hist");

		//determine placing of legend and pave
		Double_t xmin;
		if (m/((h_res->GetXaxis()->GetXmax() + h_res->GetXaxis()->GetXmin())/2) < 1) xmin = 0.7;
		else xmin = 0.2;
		
		//pave
		pave = new TPaveText(xmin,0.63,xmin+0.2,0.78,"NDC");
		pave->AddText((jname.str()).c_str());
		pave->AddText(yrnames[do_2019].c_str());
		pave->AddText((luminame.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.04);
		pave->Draw("same");

		std::cout << "Events = " << N << std::endl;
		std::cout << "Mean = " << m << std::endl;
		std::cout << "RMS = " << s << std::endl;

		if(do_print) {
			if(do_print > maxPrint) do_print = 1; //png default
			std::string img[maxPrint] = {"png","eps"};
			std::stringstream oname;
			oname << "tower_" << year[do_2019] << "_noise_" << energies[num] << "gev_lumi" << lumis[Lnum] << "_etabin" << etabin << "." << img[do_print-1];

			can->Print((oname.str()).c_str(),img[do_print-1].c_str());
		}
	}
	else _file->Close();
	
	//return mean of histo
	return m;
}

//---------------------------------------
//function to graph offsets for each lumi
TGraph* get_offsets(int Lnum, bool do_2019, bool do_show, unsigned do_print=0){
	Double_t *xval = new Double_t[maxHDeta];
	Double_t *yval = new Double_t[maxHDeta];
	
	Double_t etastep = (etas[maxHDeta] - etas[0])/maxHDeta;
	
	//get offsets for each eta bin
	for(int i = 0; i < maxHDeta; i++){
		xval[i] = etas[i] + etastep/2.;
		yval[i] = get_noise(0,Lnum,i,do_2019,0);
		yval[i] *= 0.125; //scale from area of ring deta < 0.5 to area of cone dR < 0.5: (pi/4)/(2pi)
	}

	TGraph* graph = new TGraph(maxHDeta,xval,yval);
	
	std::stringstream luminame;
	luminame << "lumi = " << lumis[Lnum] << " fb^{-1}";
	
	//plotting variables
	TCanvas* can;
	TPave* pave;
	
	if(do_show){
	
	
		if(do_print){
		
		}
	}
	
	return graph;
}

//--------------------------------------------------------
//function to plot offsets for all lumis on same pad
void plot_offsets(bool do_2019, unsigned do_print=0){
	TGraph* graphs[maxHDlumi];
	Double_t ymax = 0;
	Double_t ymin = 1e10;
	
	for(int j = 0; j < maxHDlumi; j++){
		//get graph from above function
		graphs[j] = get_offsets(j,do_2019,0);
		Double_t *vals = graphs[j]->GetY();
		
		//check extrema
		for(int i = 0; i < maxHDeta; i++){
			if(ymax < vals[i]) ymax = vals[i];
			if(ymin > vals[i]) ymin = vals[i];
		}
	}
	//manual setting
	ymin = 0;

	//plotting variables
	TCanvas* can = new TCanvas("offsets","offsets",700,500);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.15,0.05,0.15,0.075);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();

	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};
	std::stringstream jname;
	jname << "CaloTowers, dR < 0.5, " << yrnames[do_2019];
	
	double xmin = 0.3;
	TPaveText* pave = new TPaveText(xmin,0.94,xmin+0.5,0.99,"NDC");
	pave->AddText((jname.str()).c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	
	std::stringstream luminames[maxHDlumi];
	TLegend *leg = new TLegend(0.7,0.78,0.9,0.89);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	
	Color_t colors[] = {kBlack, kRed};
	for(int j = 0; j < maxHDlumi; j++){
		//formatting
		graphs[j]->SetTitle("");
		graphs[j]->GetXaxis()->SetTitle("#eta");
		graphs[j]->GetYaxis()->SetTitle("#LTp_{T}#GT Offset [GeV]");
		graphs[j]->SetLineColor(colors[j]);
		graphs[j]->SetMarkerColor(colors[j]);
		graphs[j]->SetMarkerStyle(20);
		graphs[j]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
		
		//more formatting
		graphs[j]->GetXaxis()->SetTitleOffset(0.95);
		graphs[j]->GetYaxis()->SetTitleOffset(1.0);
		graphs[j]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		graphs[j]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		graphs[j]->GetXaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		graphs[j]->GetXaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		graphs[j]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		graphs[j]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		
		//add to legend
		luminames[j] << lumis[j] << " fb^{-1}";
		leg->AddEntry(graphs[j],(luminames[j].str()).c_str(),"pl");
	
		if(j==0) graphs[j]->Draw("APZ");
		else graphs[j]->Draw("PZ same");
	}
	pave->Draw("same");
	leg->Draw("same");
	
	if(do_print){
		if(do_print > maxPrint) do_print = 1; //png default
		std::string img[maxPrint] = {"png","eps"};
		std::stringstream oname;
		oname << "dark_offsets_" << year[do_2019] << "." << img[do_print-1];

		can->Print((oname.str()).c_str(),img[do_print-1].c_str());	
	}
}

//----------------------------------------------------------------
//function to make python file with offsets for each lumi and ieta
void print_offsets(string outname, bool do_2019){
	//open output file
	std::ofstream output(("../python/"+outname).c_str());
	if (!output) {
		std::cerr << "Cannot open the output file " << outname << "\n";
		return;
	}
	
	string s4 = "    "; //python tab
	
	output << "import FWCore.ParameterSet.Config as cms" << std::endl;
	output << std::endl;
	output << "CaloJetPtOffsets = cms.PSet(" << std::endl;
	output << s4 << "etaMin = cms.double(" << etas[0] << ")," << std::endl;
	output << s4 << "etaStep = cms.double(" << (etas[maxHDeta] - etas[0])/maxHDeta << ")," << std::endl;

	output << std::fixed << std::setprecision(6); 
	
	for(int j = 0 ; j < maxHDlumi; j++){
		output << s4 << "offset_" << (int)lumis[j] << " = cms.untracked.vdouble(";
	
		//get offsets for this lumi
		TGraph* gtmp = get_offsets(j,do_2019,0);
		Double_t* ytmp = gtmp->GetY();
		for(int i = 0; i < maxHDeta; i++){
			output << ytmp[i];
			if(i<maxHDeta-1) output << ", ";
		}
		output << ")";
		if(j<maxHDlumi-1) output << ",";
		output << std::endl;
	}
	output << ")" << std::endl;
}
