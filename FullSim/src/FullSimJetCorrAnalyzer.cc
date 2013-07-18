// Global FWCore clases
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "ForwardCaloUpgrade/FullSim/interface/FullSimJetCorrAnalyzer.h"
#include <sstream>
#include <cmath>
#include <iostream>

//CaloJets
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//PFJets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//GenJets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <Math/VectorUtil.h>

using namespace std;
using namespace edm;
using namespace reco;

double FullSimJetCorrAnalyzer::phi(double x, double y) {
	double phi_ = atan2(y, x);
	return (phi_>=0) ?  phi_ : phi_ + 2*TMath::Pi();
}
double FullSimJetCorrAnalyzer::DeltaPhi(double phi1, double phi2) {
	double phi1_= phi( cos(phi1), sin(phi1) );
	double phi2_= phi( cos(phi2), sin(phi2) );
	double dphi_= phi1_-phi2_;
	if( dphi_> TMath::Pi() ) dphi_-=2*TMath::Pi();
	if( dphi_<-TMath::Pi() ) dphi_+=2*TMath::Pi();

	return dphi_;
}
double FullSimJetCorrAnalyzer::DeltaR(double phi1, double eta1, double phi2, double eta2){
	double dphi = DeltaPhi(phi1,phi2);
	double deta = eta2 - eta1;
	double dR2 = dphi*dphi + deta*deta;
	return sqrt(dR2);
}

FullSimJetCorrAnalyzer::FullSimJetCorrAnalyzer(const edm::ParameterSet& iConfig) { 
	outname = iConfig.getParameter<string>("fileName");
	dRcut = iConfig.getParameter<double>("dRcut");
	
	//information for calojet offset corrections
	lumi = iConfig.getParameter<double>("lumi");
	ParameterSet offsets = iConfig.getParameter<ParameterSet>("offsets");
	etaMin = offsets.getParameter<double>("etaMin");	
	etaStep = offsets.getParameter<double>("etaStep");	
	
	//in case correction is not available, vector will have 0 size
	std::stringstream offname;
	offname << "offset_" << lumi;
	//cout << "vector: " << offname.str() << endl;
	offsetCorr = offsets.getUntrackedParameter<vector<double> >((offname.str()).c_str(), vector<double>(0));
	//cout << "size: " << offsetCorr.size() << endl;
}

FullSimJetCorrAnalyzer::~FullSimJetCorrAnalyzer() { }

// ------------ method called for each event  ------------
void
FullSimJetCorrAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//MATCH JETS
	//----------
	double t_calo_pt_raw, t_calo_pt, t_calo_eta, t_calo_phi, t_calo_area;
	t_calo_pt_raw = t_calo_pt = t_calo_eta = t_calo_phi = t_calo_area = 0;
	double t_pf_pt, t_pf_eta, t_pf_phi, t_pf_area;
	t_pf_pt = t_pf_eta = t_pf_phi = t_pf_area = 0;
	double t_gen_pt, t_gen_eta, t_gen_phi, t_gen_area;
	t_gen_pt = t_gen_eta = t_gen_phi = t_gen_area = 0;

	Handle<CaloJetCollection> h_CaloJets;
	iEvent.getByLabel("ak5CaloJets", h_CaloJets);
	CaloJetCollection::const_iterator caloIt = h_CaloJets->begin();
	CaloJetCollection::const_iterator calo1 = h_CaloJets->end();

	Handle<PFJetCollection> h_PFJets;
	iEvent.getByLabel("ak5PFJets", h_PFJets);
	PFJetCollection::const_iterator pfIt = h_PFJets->begin();
	PFJetCollection::const_iterator pf1 = h_PFJets->end();
	
	Handle<GenJetCollection> h_GenJets;
	iEvent.getByLabel("ak5GenJets", h_GenJets);
	GenJetCollection::const_iterator genIt = h_GenJets->begin();
	GenJetCollection::const_iterator gen1 = h_GenJets->end();
	
	//find leading genjet
	double min_pt = 0;
	
	for(; genIt != h_GenJets->end(); genIt++){
		double curr_pt = genIt->pt();
		if(curr_pt > min_pt) {
			min_pt = curr_pt;
			gen1 = genIt;
		}
	}
	
	//if there is a GenJet, find the leading matched calojet & pfjet
	double calo_min_pt = 0;
	double pf_min_pt = 0;
	if(gen1 != h_GenJets->end()){
		//store the genjet vars
		t_gen_pt = gen1->pt();
		t_gen_eta = gen1->eta();
		t_gen_phi = gen1->phi();
		t_gen_area = gen1->jetArea();
		
		//loop over calo
		for(; caloIt != h_CaloJets->end(); caloIt++){
			//check match within cone
			double dR = ROOT::Math::VectorUtil::DeltaR(gen1->p4(),caloIt->p4());
			if (dR < dRcut) {
				//now check if it could be the leading jet
				double curr_pt = caloIt->pt();
				if(curr_pt > calo_min_pt){
					calo_min_pt = curr_pt;
					calo1 = caloIt;
				}
			}
		}
		
		//loop over pf
		for(; pfIt != h_PFJets->end(); pfIt++){
			//check match within cone
			double dR = ROOT::Math::VectorUtil::DeltaR(gen1->p4(),pfIt->p4());
			if (dR < dRcut) {
				//now check if it could be the leading jet
				double curr_pt = pfIt->pt();
				if(curr_pt > pf_min_pt){
					pf_min_pt = curr_pt;
					pf1 = pfIt;
				}
			}
		}
		
	}

	//if a matching CaloJet was found, store its vars
	//with noise offset correction applied to pT (raw pT also stored)
	if(calo1 != h_CaloJets->end()){
		t_calo_pt_raw = calo1->pt();
		t_calo_eta = calo1->eta();
		t_calo_phi = calo1->phi();
		t_calo_area = calo1->jetArea();
		
		//pT correction
		int ieta = 0;
		double corr = 0;
		//protection for eta too low
		//cout << "calo_eta = " << t_calo_eta << ", etaMin = " << etaMin << endl;
		if(t_calo_eta > etaMin) {
			ieta = abs((int) ((t_calo_eta - etaMin) / etaStep));
			//cout << "ieta = " << ieta << endl;
			//protection for no corrs or eta too high
			int osize = offsetCorr.size();
			if(ieta < osize) {
				//include area scaling: corrs scaled to area Pi/4 (typical cone dR < 0.5)
				corr = offsetCorr[ieta]*t_calo_area/(TMath::Pi()/4);
				
				//cout << "corr = " << offsetCorr[ieta] << "*" << t_calo_area/(TMath::Pi()/4) << " = " << corr << endl;
			}
		}
		
		t_calo_pt = calo1->pt() - corr;
	}
	
	//if a matching PFJet was found, store its vars
	if(pf1 != h_PFJets->end()){
		t_pf_pt = pf1->pt();
		t_pf_eta = pf1->eta();
		t_pf_phi = pf1->phi();
		t_pf_area = pf1->jetArea();
	}
	
	e_gen_pt = t_gen_pt;
	e_gen_eta = t_gen_eta;
	e_gen_phi = t_gen_phi;
	e_gen_area = t_gen_area;
	e_calo_pt_raw = t_calo_pt_raw;
	e_calo_pt = t_calo_pt;
	e_calo_eta = t_calo_eta;
	e_calo_phi = t_calo_phi;
	e_calo_area = t_calo_area;
	e_pf_pt = t_pf_pt;
	e_pf_eta = t_pf_eta;
	e_pf_phi = t_pf_phi;
	e_pf_area = t_pf_area;
	
	tree_tot->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
FullSimJetCorrAnalyzer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void 
FullSimJetCorrAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
FullSimJetCorrAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
	out_file = new TFile(outname.c_str(), "RECREATE");
	tree_tot = new TTree("Total", "Energy Calorimeter info");
	tree_tot->Branch("GenJetPt",&e_gen_pt,"e_gen_pt/D");
	tree_tot->Branch("GenJetEta",&e_gen_eta,"e_gen_eta/D");
	tree_tot->Branch("GenJetPhi",&e_gen_phi,"e_gen_phi/D");
	tree_tot->Branch("GenJetArea",&e_gen_area,"e_gen_area/D");
	tree_tot->Branch("CaloJetPt",&e_calo_pt,"e_calo_pt/D");
	tree_tot->Branch("CaloJetPtRaw",&e_calo_pt_raw,"e_calo_pt_raw/D");
	tree_tot->Branch("CaloJetEta",&e_calo_eta,"e_calo_eta/D");
	tree_tot->Branch("CaloJetPhi",&e_calo_phi,"e_calo_phi/D");
	tree_tot->Branch("CaloJetArea",&e_calo_area,"e_calo_area/D");
	tree_tot->Branch("PFJetPt",&e_pf_pt,"e_pf_pt/D");
	tree_tot->Branch("PFJetEta",&e_pf_eta,"e_pf_eta/D");
	tree_tot->Branch("PFJetPhi",&e_pf_phi,"e_pf_phi/D");
	tree_tot->Branch("PFJetArea",&e_pf_area,"e_pf_area/D");

}

// ------------ method called when ending the processing of a run  ------------
void 
FullSimJetCorrAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) { 
	out_file->cd();
	
	tree_tot->Write();
	
	out_file->Close();
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FullSimJetCorrAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FullSimJetCorrAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FullSimJetCorrAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



