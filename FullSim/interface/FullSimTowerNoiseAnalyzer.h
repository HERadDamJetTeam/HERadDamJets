#ifndef FullSimTowerNoiseAnalyzer_h
#define FullSimTowerNoiseAnalyzer_h

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TFile;
class TTree;

class FullSimTowerNoiseAnalyzer : public edm::EDAnalyzer {
	public:
		explicit FullSimTowerNoiseAnalyzer(const edm::ParameterSet&);
		~FullSimTowerNoiseAnalyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
	
		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		double phi(double x, double y);
		double DeltaPhi(double phi1, double phi2);
		double DeltaR(double phi1, double eta1, double phi2, double eta2);				
		
		//member variables
		TFile* out_file;
		TTree* tree_tot;
		double e_gen_eta, e_gen_phi;
		double e_noise_en, e_noise_pt;
		double dRcut;
		std::string outname;

};

//define this as a plug-in
DEFINE_FWK_MODULE(FullSimTowerNoiseAnalyzer);

#endif