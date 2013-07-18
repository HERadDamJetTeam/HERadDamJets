// Global FWCore clases
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "ForwardCaloUpgrade/FullSim/interface/FullSimTowerNoiseAnalyzer.h"
#include <sstream>
#include <cmath>

//CaloTowers
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

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

double FullSimTowerNoiseAnalyzer::phi(double x, double y) {
	double phi_ = atan2(y, x);
	return (phi_>=0) ?  phi_ : phi_ + 2*TMath::Pi();
}
double FullSimTowerNoiseAnalyzer::DeltaPhi(double phi1, double phi2) {
	double phi1_= phi( cos(phi1), sin(phi1) );
	double phi2_= phi( cos(phi2), sin(phi2) );
	double dphi_= phi1_-phi2_;
	if( dphi_> TMath::Pi() ) dphi_-=2*TMath::Pi();
	if( dphi_<-TMath::Pi() ) dphi_+=2*TMath::Pi();

	return dphi_;
}
double FullSimTowerNoiseAnalyzer::DeltaR(double phi1, double eta1, double phi2, double eta2){
	double dphi = DeltaPhi(phi1,phi2);
	double deta = eta2 - eta1;
	double dR2 = dphi*dphi + deta*deta;
	return sqrt(dR2);
}

FullSimTowerNoiseAnalyzer::FullSimTowerNoiseAnalyzer(const edm::ParameterSet& iConfig) { 
	outname = iConfig.getParameter<string>("fileName");
	dRcut = iConfig.getParameter<double>("dRcut");
}

FullSimTowerNoiseAnalyzer::~FullSimTowerNoiseAnalyzer() { }

// ------------ method called for each event  ------------
void
FullSimTowerNoiseAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//COLLECT CALOTOWERS
	//------------------
	Handle<CaloTowerCollection> CaloTowers;
	bool bCT = iEvent.getByLabel("towerMaker", CaloTowers);
	
	double sum_noise_en, sum_noise_pt;
	sum_noise_en = sum_noise_pt = 0;
	
	//MATCH JETS
	//----------
	double t_gen_eta, t_gen_phi;
	t_gen_eta = t_gen_phi = 0;

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
	
	//if there is a GenJet, find the offset from a cone around negative eta
	if(gen1 != h_GenJets->end()){
		//store the genjet vars
		t_gen_eta = gen1->eta();
		t_gen_phi = gen1->phi();
		
		//loop over calotowers
		if(bCT){
			for (CaloTowerCollection::const_iterator hit = CaloTowers->begin(); hit!=CaloTowers->end(); ++hit) {
				double h_eta = hit->eta();
				//double h_phi = hit->phi();
				double dEtaNoise = fabs(-t_gen_eta - h_eta); //collect noise from negative eta, all phi (ring instead of cone, area = pi)
				double h_en = hit->energy();
				double h_pt = hit->pt();
				if(dEtaNoise < dRcut) { sum_noise_en += h_en; sum_noise_pt += h_pt; }
			}
		}	
		
	}

	e_gen_eta = t_gen_eta;
	e_gen_phi = t_gen_phi;

	e_noise_en = sum_noise_en;
	e_noise_pt = sum_noise_pt;
	
	tree_tot->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
FullSimTowerNoiseAnalyzer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void 
FullSimTowerNoiseAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
FullSimTowerNoiseAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
	out_file = new TFile(outname.c_str(), "RECREATE");
	tree_tot = new TTree("Total", "Energy Calorimeter info");
	tree_tot->Branch("GenJetEta",&e_gen_eta,"e_gen_eta/D");
	tree_tot->Branch("GenJetPhi",&e_gen_phi,"e_gen_phi/D");
	tree_tot->Branch("TowerNoiseEnergy",&e_noise_en,"e_noise_en/D");
	tree_tot->Branch("TowerNoisePt",&e_noise_pt,"e_noise_pt/D");
}

// ------------ method called when ending the processing of a run  ------------
void 
FullSimTowerNoiseAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) { 
	out_file->cd();
	
	tree_tot->Write();
	
	out_file->Close();
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FullSimTowerNoiseAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FullSimTowerNoiseAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FullSimTowerNoiseAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



