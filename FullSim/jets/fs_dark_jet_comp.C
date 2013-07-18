//ROOT headers
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
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

//#define maxHDe 3
#define maxHDe 1
#define maxHDlumi 8
//#define maxHDlumi 2
#define maxHDeta 4
#define maxHDqty 3
#define maxPrint 2
#define M_PI 3.14159265358979323846

//energy values - global
//Double_t energies[] = {30, 50, 100};
Double_t energies[] = {30};
Double_t lumis[] = {0, 100, 200, 300, 400, 500, 600, 700};
//Double_t lumis[] = {0, 500};
Double_t etas[] = {1.8, 2.1, 2.4, 2.7, 3.0};
Int_t year[] = {2017, 2019};

//NB: in this file, "energies" refers to pT

//-------------------------------------
//class to store energy resolution data
class energyRes {
	//vars
	private:
		Double_t energy;
		Double_t mu;
		Double_t mu_err;
		Double_t sigma;
		Double_t sigma_err;
		Int_t N_evts;

	public:
		//constructors
		energyRes(Double_t en) : energy(en) {}
		
		//set members
		void setVals(Double_t m, Double_t me, Double_t s, Double_t se, Int_t N){
			mu = m;
			mu_err = me;
			sigma = s;
			sigma_err = se;
			N_evts = N;
		}
		void setEnergy(Double_t en) { energy = en; }
	
		//access members
		Double_t getE(){ return energy; }
		Double_t getMu(){ return mu; }
		Double_t getMuErr(){ return mu_err; }
		Double_t getSigma(){ return sigma; }
		Double_t getSigmaErr(){ return sigma_err; }
		Int_t getN(){ return N_evts; }		

};


//----------------------------------
//function to fit energy resolutions
//qty: 0 = pT, 1 = eta, 2 = phi
std::pair<energyRes*,TH1F*> get_res(int num, int Lnum, int etabin, int qty, bool do_2019, bool do_pf, bool do_fit, bool do_show, unsigned do_print=0){
	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		energyRes* theRes = new energyRes(0);
		theRes->setVals(0,0,0,0,0);
		return std::pair<energyRes*,TH1F*>(theRes,NULL);
	}

	//make filenames
	std::stringstream fname, jname, luminame;

	fname << "tree_jet_" << year[do_2019] << "_" << energies[num] << "_lumi" << lumis[Lnum] << ".root";
	jname << "d jet, p_{T} = " << energies[num] << " GeV, " << etas[etabin] << " < #eta < " << etas[etabin+1];
	luminame << "lumi = " << lumis[Lnum] << " fb^{-1}";

	//open file and histo
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TH1F* h_res;
	TTree* totalTree = (TTree*)_file->Get("Total");

	std::stringstream drawname, cutname;
	string jalgo[] = {"Calo","PF"};
	if(qty==0) drawname << "(" << jalgo[do_pf] << "JetPt/GenJetPt)>>htemp(50,0.0,2.0)";
	else if(qty==1) drawname << "(" << jalgo[do_pf] << "JetEta - GenJetEta)>>htemp(50,-0.5,0.5)";
	else if(qty==2) drawname << "(" << jalgo[do_pf] << "JetPhi - GenJetPhi < -" << M_PI << ")*(" << jalgo[do_pf] << "JetPhi-GenJetPhi+" << 2*M_PI << ")+(" << jalgo[do_pf] << "JetPhi-GenJetPhi>" << M_PI << ")*(" << jalgo[do_pf] << "JetPhi-GenJetPhi-" << 2*M_PI << ")+(" << jalgo[do_pf] << "JetPhi-GenJetPhi<=" << M_PI << " && " << jalgo[do_pf] << "JetPhi-GenJetPhi>=-" << M_PI << ")*(" << jalgo[do_pf] << "JetPhi-GenJetPhi)>>htemp(50,-0.5,0.5)";
	cutname << jalgo[do_pf] << "JetPt>0 && abs(GenJetEta)>" << etas[etabin] << " && abs(GenJetEta)<=" << etas[etabin+1];
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	h_res = (TH1F*)gDirectory->Get("htemp");

	TF1* gfit;

	//plotting variables
	TCanvas* can;
	TLegend* leg;
	TPaveText* pave;

	//create instance of energyRes object
	energyRes* theRes = new energyRes(energies[num]);	
	
	//get values from histo
	Double_t m = h_res->GetMean();
	Double_t me = h_res->GetMeanError();
	//Double_t m = h_res->GetBinCenter(h_res->GetMaximumBin()); //peak
	Double_t s = h_res->GetRMS();
	Double_t se = h_res->GetRMSError();
	Int_t N = h_res->GetEntries();

	//find peak
	TSpectrum *spec = new TSpectrum(5);
	//spec->Search(h_res,6,"nodraw goff");
	spec->Search(h_res,6,"nobackground nodraw goff"); //turn off background removal when nbins too small
	Float_t* xpos = spec->GetPositionX();
	Float_t* ypos = spec->GetPositionY();
	Double_t p = xpos[0];
	Double_t ph = ypos[0];
	if(do_show) std::cout << "peak: " << p << std::endl;

	//names
	std::string ofit;
	if(do_fit) ofit = "fit";
	else ofit = "nofit";

	//setup fitting function
	if (do_fit){
		gfit = new TF1("tot","gaus",h_res->GetXaxis()->GetXmin(),h_res->GetXaxis()->GetXmax());
		gfit->SetParameters((Double_t)ph,p,s);
		if(m > p) gfit->SetRange(p-1.5*s,p+1.0*s); //high tail
		else gfit->SetRange(p-1.0*s,p+1.5*s); //low tail
		gfit->SetLineColor(kRed);
		gfit->SetMarkerColor(kRed);
		gfit->SetLineWidth(2);
	}

	//do fit
	Double_t c_val[3];
	Double_t c_err[3];
	if(do_fit){
		h_res->Fit(gfit,"LNQR");
		for(int j=0;j<3;j++){
			c_val[j] = gfit->GetParameter(j);
			c_err[j] = gfit->GetParError(j);
		}

		//store parameters only for tot (imip = 2)
		theRes->setVals(c_val[1],c_err[1],c_val[2],c_err[2],N);	
	}

	else theRes->setVals(m,me,s,se,N);

	std::stringstream muname, signame, Nname;
	muname.precision(2);
	signame.precision(2);
	if (do_fit) {
		muname << fixed << "#mu = " << c_val[1] << " #pm " << c_err[1];
		signame << fixed << "#sigma = " << c_val[2] << " #pm " << c_err[2];
	}
	else {
		muname << fixed << "Mean = " << m << " #pm " << me;
		signame << fixed << "RMS = " << s << " #pm " << se;	
	}
	Nname << "N = " << N; 

	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};
	std::string qtys[] = {"res","deta","dphi"};
	std::string jalgonames[] = {"calo","pf"};
	std::string qtynames[] = {"p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}","#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}","#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}"};

	//plotting
	if (do_show){
		can = new TCanvas(qtys[qty].c_str(),qtys[qty].c_str(),700,500);
		can->cd();	
		
		//plot histo and fit
		h_res->SetStats(kTRUE);
		gStyle->SetOptStat("mr");
		h_res->SetTitle("");
		h_res->GetXaxis()->SetTitle(qtynames[qty].c_str());
		h_res->GetYaxis()->SetTitle("");
		h_res->SetLineColor(kBlack);
		h_res->Draw("hist");
		if(do_fit) gfit->Draw("same");	

		//determine placing of legend and pave
		Double_t xmin;
		if (m/((h_res->GetXaxis()->GetXmax() + h_res->GetXaxis()->GetXmin())/2) < 1) xmin = 0.7;
		else xmin = 0.2;

		if(do_fit) { //legend
			leg = new TLegend(xmin,0.78,xmin+0.2,0.88);
			leg->AddEntry(h_res,("ak5 " + jalgo[do_pf] + "Jet").c_str());
			leg->AddEntry(gfit,"Fit");
			leg->SetFillColor(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(42);
			leg->SetTextSize(0.04);
			leg->Draw("same");
		}
		
		//pave
		pave = new TPaveText(xmin,0.48,xmin+0.2,0.78,"NDC");
		pave->AddText((jname.str()).c_str());
		pave->AddText(yrnames[do_2019].c_str());
		pave->AddText((luminame.str()).c_str());
		pave->AddText((Nname.str()).c_str());
		pave->AddText((muname.str()).c_str());
		pave->AddText((signame.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.04);
		pave->Draw("same");

		std::cout << "total:" << std::endl;
		if(do_fit){
			std::cout << "mu = " << c_val[1] << " +/- " << c_err[1] << std::endl;
			std::cout << "sigma = " << c_val[2] << " +/- " << c_err[2] << std::endl;
		}
		else{
			std::cout << "Mean = " << m << std::endl;
			std::cout << "RMS = " << s << std::endl;
		}

		if(do_print) {
			if(do_print > maxPrint) do_print = 1; //png default
			std::string img[maxPrint] = {"png","eps"};
			std::stringstream oname;
			oname << "hcal_" << year[do_2019] << "_" << qtys[qty] << "_" << ofit << "_" << energies[num] << "gev_lumi" << lumis[Lnum] << "_etabin" << etabin << "." << img[do_print-1];

			can->Print((oname.str()).c_str(),img[do_print-1].c_str());
		}
	}
	//else _file->Close();
	
	//return data structure with relevant info, and histogram
	return std::pair<energyRes*,TH1F*>(theRes,h_res);
}

//----------------------------------
//function to compare jet pT spectra
void plot_overlay(int num, int etabin, int qty, bool do_2019, bool do_pf, unsigned do_print=0){
	TH1F* the_histos[maxHDlumi];

	Double_t ymax = 0;

	//get histos
	for(int i = 0; i < maxHDlumi; i++){
		the_histos[i] = get_res(num,i,etabin,qty,do_2019,do_pf,0,0).second; //gets histo, discards energyRes object

		//get max
		Double_t max = the_histos[i]->GetBinContent(the_histos[i]->GetMaximumBin());
		
		//check max
		if(ymax < max) ymax = max;
	}
	
	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};
	std::string jalgonames[] = {"calo","pf"};
	std::string qtynames[] = {"p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}","#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}","#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}"};
	Int_t colors[] = {kBlack, kBlue, kMagenta+2, kRed, kCyan+2, kMagenta, kOrange+7, kYellow+3};
	//Int_t styles[] = {1,2,4,3};
	Int_t styles[] = {1,1,1,1,1,1,1,1};
	std::stringstream luminame[maxHDlumi];

	std::string qtys[] = {"res","deta","dphi"};
	TCanvas* can = new TCanvas(qtys[qty].c_str(),qtys[qty].c_str(),700,500);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.15,0.05,0.15,0.075);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();

	double xmin;
	if(qty==0) xmin = 0.65;
	else xmin = 0.2;
	double legy = 0.89;
	TLegend* leg = new TLegend(xmin,legy-0.05*maxHDlumi-0.01,xmin+0.2,legy);

	for(int i = 0; i < maxHDlumi; i++){
		the_histos[i]->SetTitle("");
		the_histos[i]->GetXaxis()->SetTitle(qtynames[qty].c_str());
		the_histos[i]->GetYaxis()->SetTitle("");
		the_histos[i]->SetLineColor(colors[i]);
		the_histos[i]->SetLineStyle(styles[i]);
		the_histos[i]->SetLineWidth(3);
		the_histos[i]->SetMarkerColor(colors[i]);
		the_histos[i]->GetYaxis()->SetRangeUser(0,ymax*1.1);

		//the_histos[i]->GetXaxis()->SetTitleOffset(0.95);
		//the_histos[i]->GetXaxis()->SetLabelOffset(0);
		//the_histos[i]->GetYaxis()->SetTitleOffset(1.2);
		the_histos[i]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_histos[i]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_histos[i]->GetXaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_histos[i]->GetXaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_histos[i]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_histos[i]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		
		luminame[i] << "Lumi: " << lumis[i] << " fb^{-1}";
		leg->AddEntry(the_histos[i],(luminame[i].str()).c_str());
		
		if(i==0) the_histos[i]->Draw("hist");
		else the_histos[i]->Draw("hist same");
	}

	//legend formatting
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->Draw("same");

	std::stringstream jname;
	jname << "d jet, p_{T} = " << energies[num] << " GeV, " << etas[etabin] << " < #eta < " << etas[etabin+1] << ", " << yrnames[do_2019];

	TPaveText* pave = new TPaveText(0.1,0.94,0.9,0.99,"NDC");
	pave->AddText((jname.str()).c_str());
	//pave->AddText(yrnames[do_2019].c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");

	if(do_print){
		if(do_print > maxPrint) do_print = 1; //png default
		std::string img[maxPrint] = {"png","eps"};
		std::stringstream oname;
		oname << "dark_overlay_" << jalgonames[do_pf] << "_" << year[do_2019] << "_" << qtys[qty] << "_" << energies[num] << "gev_etabin" << etabin << "." << img[do_print-1];

		can->Print((oname.str()).c_str(),img[do_print-1].c_str());	
	}


}

//-----------------------------------------------
//function to generate and save all overlay plots
void alltheoverlays(unsigned do_print=1){
	for(int n = 0; n < maxHDe; n++){ //energy
		for(int i = 0; i < maxHDeta; i++){ //eta bins
			int j = 0; //only calo
			//for(int j = 0; j < 2; j++){ //calo or pf
				int y = 1; //only 2019
				//for(int y = 0; y < 2; y++){ //2017 or 2019
					int k = 0; //only pT
					//for(int k = 0; k < maxHDqty; k++){ //qty
						plot_overlay(n,i,k,y,j,do_print);
					//}
				//}
			//}
		}
	}

}

//-------------------------------------------------
//function to compare jet qtys vs. lumi for some pT
TGraphErrors* plot_lumicomp(int num, int etabin, int qty, bool do_2019, bool do_pf, bool do_sigma, bool do_show, unsigned do_print=0){
	bool do_fit = false;
	//if(qty==0) do_fit = true; //fit only for pT qtys
	
	TGraphErrors* the_graphs;
	energyRes* res_temp;
	Double_t vals[maxHDlumi];
	Double_t errs[maxHDlumi];
	Double_t xerrs[maxHDlumi];

	Double_t ymax = 0;
	Double_t ymin = 1e10;
	
	//get means
	for(int i = 0; i < maxHDlumi; i++){
		res_temp = get_res(num,i,etabin,qty,do_2019,do_pf,do_fit,0).first; //gets energyRes, discards histo
		
		Double_t v, ve;
		if(do_sigma && qty == 0){
			v = res_temp->getSigma()/res_temp->getMu();
			ve = res_temp->getSigmaErr()/res_temp->getMu();
		}
		else if(do_sigma){
			v = res_temp->getSigma();
			ve = res_temp->getSigmaErr();
		}
		else {
			v = res_temp->getMu();
			ve = res_temp->getMuErr();
		}

		vals[i] = v;
		errs[i] = ve;
		xerrs[i] = 0;
		
		//check extrema
		if(ymax < vals[i]) ymax = vals[i];
		if(ymin > vals[i]) ymin = vals[i];
	}

	std::string jalgonames[] = {"calo","pf"};
//	std::string qtynames[2][maxHDqty] = { {"#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","<#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}>","<#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}>"},
//								   {"#sigma(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	std::string qtynames[2][maxHDqty] = { {"#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","#LT#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}#GT","#LT#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}#GT"},
								   {"RMS(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	
	std::string qtys[maxHDqty] = {"res","deta","dphi"};
	std::string vtype[2][2] = { {"mean","rms"}, {"mu","sigma"} };
	
	the_graphs = new TGraphErrors(maxHDlumi,lumis,vals,xerrs,errs);
	
	if(do_show){
		TCanvas* can = new TCanvas((qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),(qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),700,500);
		
		//formatting
		the_graphs->SetTitle("");
		the_graphs->GetXaxis()->SetTitle("Luminosity [fb^{-1}]");
		the_graphs->GetYaxis()->SetTitle(qtynames[do_sigma][qty].c_str());
		the_graphs->SetFillColor(0);
		the_graphs->SetMarkerStyle(20);
		the_graphs->SetMarkerColor(kBlack);
		the_graphs->SetMarkerSize(1.5);
		the_graphs->SetLineColor(kBlack);
		the_graphs->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
		
		the_graphs->Draw("APZ");

		if(do_print){
			if(do_print > maxPrint) do_print = 1; //png default
			std::string img[maxPrint] = {"png","eps"};
			std::stringstream oname;
			oname << "dark_lumicomp_" << year[do_2019] << "_" << qtys[qty] << "_" << vtype[do_fit][do_sigma] << "_" << energies[num] << "gev_etabin" << etabin << "." << img[do_print-1];

			can->Print((oname.str()).c_str(),img[do_print-1].c_str());	
		}
	}
	
	return the_graphs;

}

//-------------------------------------------------------------
//function to compare jet qtys vs. lumi for some pT, both years
void comp_lumicomp(int num, int etabin, int qty, bool do_pf, bool do_sigma, unsigned do_print=0){
	bool do_fit = false;
	//if(qty==0) do_fit = true; //fit only for pT qtys
	
	TGraphErrors* the_graphs[2];
	Double_t *vals[2];
	Double_t ymax = 0;
	Double_t ymin = 1e10;
	
	for(int y = 0; y < 2; y++){
		the_graphs[y] = plot_lumicomp(num,etabin,qty,y,do_pf,do_sigma,0);
		vals[y] = the_graphs[y]->GetY();
		
		//check extrema
		for(int i = 0; i < maxHDlumi; i++){
			if(ymax < vals[y][i]) ymax = vals[y][i];
			if(ymin > vals[y][i]) ymin = vals[y][i];
		}
	}

	std::string jalgonames[] = {"calo","pf"};
//	std::string qtynames[2][maxHDqty] = { {"#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","<#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}>","<#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}>"},
//								   {"#sigma(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	std::string qtynames[2][maxHDqty] = { {"#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","#LT#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}#GT","#LT#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}#GT"},
								   {"RMS(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	
								   
	std::string qtys[maxHDqty] = {"res","deta","dphi"};
	std::string vtype[2][2] = { {"mean","rms"}, {"mu","sigma"} };
	
	TCanvas* can = new TCanvas((qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),(qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),700,500);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.15,0.05,0.15,0.075);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();

	std::stringstream jname;
	jname << "d jet, p_{T} = " << energies[num] << " GeV, " << etas[etabin] << " < #eta < " << etas[etabin+1];

	double xmin = 0.3;
	
	TPaveText* pave = new TPaveText(xmin,0.94,xmin+0.5,0.99,"NDC");
	pave->AddText((jname.str()).c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	
	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};
	TLegend *leg = new TLegend(0.2,0.2,0.3,0.3);
	//legend formatting
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	
	Int_t markers[2] = {20, 21};
	Color_t colors[2] = {kBlack, kRed};
	
	for(int y = 0; y < 2; y++){
		//formatting
		the_graphs[y]->SetTitle("");
		the_graphs[y]->GetXaxis()->SetTitle("Luminosity [fb^{-1}]");
		the_graphs[y]->GetYaxis()->SetTitle(qtynames[do_sigma][qty].c_str());
		the_graphs[y]->SetFillColor(0);
		the_graphs[y]->SetMarkerStyle(markers[y]);
		the_graphs[y]->SetMarkerColor(colors[y]);
		the_graphs[y]->SetMarkerSize(1.5);
		the_graphs[y]->SetLineColor(colors[y]);
		the_graphs[y]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
		
		//more formatting
		the_graphs[y]->GetXaxis()->SetTitleOffset(0.95);
		the_graphs[y]->GetXaxis()->SetLabelOffset(0);
		the_graphs[y]->GetYaxis()->SetTitleOffset(1.05);
		the_graphs[y]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[y]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[y]->GetXaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[y]->GetXaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[y]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[y]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		
		//add to legend
		leg->AddEntry(the_graphs[y],yrnames[y].c_str(),"pl");
	}
	the_graphs[0]->Draw("APZ");
	the_graphs[1]->Draw("PZ same");
	pave->Draw("same");
	leg->Draw("same");

	if(do_print){
		if(do_print > maxPrint) do_print = 1; //png default
		std::string img[maxPrint] = {"png","eps"};
		std::stringstream oname;
		oname << "dark_lumicomp_" << jalgonames[do_pf] << "_" << qtys[qty] << "_" << vtype[do_fit][do_sigma] << "_" << energies[num] << "gev_etabin" << etabin << "." << img[do_print-1];

		can->Print((oname.str()).c_str(),img[do_print-1].c_str());	
	}

}

//-----------------------------------------------
//function to generate and save all overlay plots
void allthelumicomps(unsigned do_print=1){
	for(int n = 0; n < maxHDe; n++){ //energy
		for(int i = 0; i < maxHDeta; i++){ //eta bins
			for(int j = 0; j < 2; j++){ //calo or pf
				for(int k = 0; k < maxHDqty; k++){ //qty
					if(k==0) comp_lumicomp(n,i,k,j,0,do_print); //response only for pT
					comp_lumicomp(n,i,k,j,1,do_print); //resolution for all
				}
			}
		}
	}

}

//---------------------------------------------------------------
//function to compare jet qtys vs. lumi for some pT, all eta bins
void comp_etacomp(int num, int qty, bool do_2019, bool do_pf, bool do_sigma, unsigned do_print=0){
	bool do_fit = false;
	//if(qty==0) do_fit = true; //fit only for pT qtys
	
	TGraphErrors* the_graphs[maxHDeta];
	Double_t *vals[maxHDeta];
	Double_t ymax = 0;
	Double_t ymin = 1e10;
	
	for(int i = 0; i < maxHDeta; i++){
		the_graphs[i] = plot_lumicomp(num,i,qty,do_2019,do_pf,do_sigma,0);
		vals[i] = the_graphs[i]->GetY();
		
		//check extrema
		for(int j = 0; j < maxHDlumi; j++){
			if(ymax < vals[i][j]) ymax = vals[i][j];
			if(ymin > vals[i][j]) ymin = vals[i][j];
		}
	}
	
	//manual settings
	ymin = 0;
	if(qty==0 && do_sigma) ymax = 0.58;
	else if (qty==0 && !do_sigma) ymax = 1.0;
	else if (qty==1 && do_sigma) ymax = 0.11;
	else if (qty==2 && do_sigma) ymax = 0.13;

	std::string jalgonames[] = {"calo","pf"};
//	std::string qtynames[2][maxHDqty] = { {"#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","<#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}>","<#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}>"},
//								   {"#sigma(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#mu(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	std::string qtynames[2][maxHDqty] = { {"#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","#LT#eta^{" + jalgonames[do_pf] + "} - #eta^{gen}#GT","#LT#phi^{" + jalgonames[do_pf] + "} - #phi^{gen}#GT"},
								   {"RMS(p_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen})/#LTp_{T}^{" + jalgonames[do_pf] + "}/p_{T}^{gen}#GT","RMS(#eta^{" + jalgonames[do_pf] + "} - #eta^{gen})","RMS(#phi^{" + jalgonames[do_pf] + "} - #phi^{gen})"} };
	
	std::string qtys[maxHDqty] = {"res","deta","dphi"};
	std::string vtype[2][2] = { {"mean","rms"}, {"mu","sigma"} };
	std::string yrnames[] = {"2017 (HPDs)","2019 (SiPMs)"};
	
	TCanvas* can = new TCanvas((qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),(qtys[qty] + "_" + vtype[do_fit][do_sigma]).c_str(),700,500);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.165,0.05,0.15,0.075);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();

	std::stringstream jname;
	jname << "d jet, p_{T} = " << energies[num] << " GeV, " << yrnames[do_2019];

	double xmin = 0.1;
	
	TPaveText* pave = new TPaveText(xmin,0.94,xmin+0.7,0.99,"NDC");
	pave->AddText((jname.str()).c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	
	std::stringstream etanames[maxHDeta];
	TLegend *leg = new TLegend(0.2,0.68,0.4,0.89);
	
	//legend formatting
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	
	Int_t markers[maxHDeta] = {20, 20, 20, 20};
	Color_t colors[maxHDeta] = {kBlack, kBlue, kMagenta, kRed};
	
	for(int i = 0; i < maxHDeta; i++){
		//formatting
		the_graphs[i]->SetTitle("");
		the_graphs[i]->GetXaxis()->SetTitle("Luminosity [fb^{-1}]");
		the_graphs[i]->GetYaxis()->SetTitle(qtynames[do_sigma][qty].c_str());
		the_graphs[i]->SetFillColor(0);
		the_graphs[i]->SetMarkerStyle(markers[i]);
		the_graphs[i]->SetMarkerColor(colors[i]);
		the_graphs[i]->SetMarkerSize(1.5);
		the_graphs[i]->SetLineColor(colors[i]);
		//the_graphs[i]->GetYaxis()->SetRangeUser(ymin*0.8,ymax*1.2);
		the_graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);
		
		//more formatting
		the_graphs[i]->GetXaxis()->SetTitleOffset(0.95);
		the_graphs[i]->GetXaxis()->SetLabelOffset(0);
		the_graphs[i]->GetYaxis()->SetTitleOffset(1.175);
		the_graphs[i]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		
		//add to legend
		etanames[i] << etas[i] << " < #eta < " << etas[i+1];
		leg->AddEntry(the_graphs[i],(etanames[i].str()).c_str(),"pl");
		
		if(i==0) the_graphs[i]->Draw("APZ");
		else the_graphs[i]->Draw("PZ same");
	}
	pave->Draw("same");
	leg->Draw("same");

	if(do_print){
		if(do_print > maxPrint) do_print = 1; //png default
		std::string img[maxPrint] = {"png","eps"};
		std::stringstream oname;
		oname << "dark_etacomp_" << year[do_2019] << "_" << jalgonames[do_pf] << "_" << qtys[qty] << "_" << vtype[do_fit][do_sigma] << "_" << energies[num] << "gev." << img[do_print-1];

		can->Print((oname.str()).c_str(),img[do_print-1].c_str());	
	}

}

//-----------------------------------------------
//function to generate and save all overlay plots
void alltheetacomps(unsigned do_print=1){
	for(int n = 0; n < maxHDe; n++){ //energy
		int y = 1; //only 2019
		//for(int y = 0; y < 2; y++){ //2017, 2019
			int j = 0; //only calo
			//for(int j = 0; j < 2; j++){ //calo or pf
				int k = 0; //only pT
				//for(int k = 0; k < maxHDqty; k++){ //qty
					//if(k==0) comp_etacomp(n,i,k,j,0,do_print); //response only for pT
					comp_etacomp(n,k,y,j,1,do_print); //resolution for all
				//}
			//}
		//}
	}

}


//-------------------------------------
//function to compare jet pT mean & RMS
//TODO: not useful until multiple jet pT samples exist
//void plot_res(int etabin, int qty, bool do_2019, bool do_pf, bool do_sigma, bool do_print){}