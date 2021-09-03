// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"

#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "DataFormats/Math/interface/libminifloat.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

using namespace std;


class ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingNanoAOD(const edm::ParameterSet&);
  ~ScoutingNanoAOD();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
<<<<<<< HEAD
  int getCharge(int pdgId);
  bool jetID(const ScoutingPFJet &pfjet);

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;
=======
  virtual void clearVars();

>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb
  const edm::EDGetTokenT<std::vector<ScoutingMuon> >        muonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingElectron> >  	electronsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPhoton> >  	photonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingParticle> >  	pfcandsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPFJet> >  		pfjetsToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >  	verticesToken;
<<<<<<< HEAD
  const edm::EDGetTokenT<GenEventInfoProduct>               genEvtInfoToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
	
  bool doL1;       
  triggerExpression::Data triggerCache_;
=======
  

  //const edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken;
>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb
      
  // Generator-level information
  // Flags for the different types of triggers used in the analysis
  // For now we are interested in events passing either the single or double lepton triggers


  // Trigger information 
       
  //edm::InputTag                algInputTag_;       
  //edm::EDGetToken              algToken_;
  //l1t::L1TGlobalUtil          *l1GtUtils_;
  //triggerExpression::Data triggerCache_;
  //std::vector<std::string> triggerPathsVector;
  //std::map<std::string, int> triggerPathsMap;
  //const edm::InputTag triggerResultsTag;
  //const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;
  
  bool doL1;       
  
  HLTPrescaleProvider hltPSProv_;
  std::string hltProcess_; //name of HLT process, usually "HLT"

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::Handle<edm::TriggerResults> triggerBits;

  std::vector<std::string>     l1Seeds_;
  std::vector<std::string>     hltSeeds_;
  std::vector<bool>            l1Result_;
<<<<<<< HEAD
       
=======
  std::vector<int>             l1Prescale_;
  std::vector<bool>            hltResult_;


  // From Hardik
  //std::vector<std::string> hltseedsvector;
  //std::vector<pair<string,int>> hltbitmap;
  //std::vector<pair<string,int>> hltprescalemap;
      
  //std::vector<std::string> l1seedsvector;
  //std::vector<pair<string,int>> l1bitmap;
  //std::vector<pair<string,int>> l1prescalemap;
        
>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb
  //Photon
  const static int 	max_pho = 1000;
  UInt_t n_pho;
  vector<Float16_t> 	    Photon_pt;
  vector<Float16_t>        	Photon_eta;
  vector<Float16_t>        	Photon_phi;
  vector<Float16_t>	    	Photon_m;
  vector<Float16_t>	    	Photon_sigmaietaieta;
  vector<Float16_t>	    	Photon_HoE;
  vector<Float16_t>        	Photon_ecaliso;
  vector<Float16_t>	    	Photon_hcaliso;

  //Electron
  const static int 	max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> 	Electron_pt;
  vector<Float16_t>     Electron_eta;
  vector<Float16_t>     Electron_phi;
  vector<Float16_t>	    Electron_m;
  vector<Float16_t>     Electron_d0;
  vector<Float16_t>	    Electron_dz;
  vector<Float16_t>	    Electron_detain;
  vector<Float16_t>	    Electron_dphiin;
  vector<Float16_t>	    Electron_sigmaietaieta;
  vector<Float16_t>	    Electron_HoE;
  vector<Float16_t>	    Electron_ooEMOop;
  vector<Float16_t>	    Electron_mHits;
  vector<Float16_t>     Electron_charge;
  vector<Float16_t>     Electron_ecaliso;
  vector<Float16_t>	    Electron_hcaliso;
  vector<Float16_t>     Electron_tkiso;

  //Muon
  UInt_t n_mu;
  vector<Float16_t> 	Muon_pt;
  vector<Float16_t> 	Muon_eta;
  vector<Float16_t> 	Muon_phi;
  vector<Float16_t> 	Muon_m;
  vector<Float16_t> 	Muon_ecaliso;
  vector<Float16_t> 	Muon_hcaliso;
  vector<Float16_t> 	Muon_trkiso;
  vector<Float16_t> 	Muon_chi2;
  vector<Float16_t> 	Muon_ndof;
  vector<Float16_t> 	Muon_charge;
  vector<Float16_t> 	Muon_dxy;
  vector<Float16_t> 	Muon_dz;
  vector<Float16_t> 	Muon_nvalidmuon_hits;
  vector<Float16_t> 	Muon_nvalidpixelhits;
  vector<Float16_t> 	Muon_nmatchedstations;
  vector<Float16_t>     Muon_type;
  vector<Float16_t>     Muon_nvalidstriphits;
  vector<Float16_t>     Muon_trkqoverp;
  vector<Float16_t>     Muon_trklambda;
  vector<Float16_t>     Muon_trkpt;
  vector<Float16_t>     Muon_trkphi;
  vector<Float16_t>     Muon_trketa;
  vector<Float16_t>     Muon_trkqoverperror;
  vector<Float16_t>     Muon_trklambdaerror;
  vector<Float16_t>     Muon_trkpterror;
  vector<Float16_t>     Muon_trkphierror;
  vector<Float16_t>     Muon_trketaerror;
  vector<Float16_t>     Muon_trkdszerror;
  vector<Float16_t>     Muon_trkdsz;

  //PFJets
  UInt_t n_jet;
  UInt_t n_jetId;
  float ht;
  bool passJetId;
  vector<Float16_t> 	Jet_pt;
  vector<Float16_t>     Jet_eta;
  vector<Float16_t>     Jet_phi;
  vector<Float16_t>	    Jet_m;
  vector<Float16_t>	    Jet_area;
  vector<Float16_t>	    Jet_chargedHadronEnergy;
  vector<Float16_t>     Jet_neutralHadronEnergy;
  vector<Float16_t>	    Jet_photonEnergy;
  vector<Float16_t>	    Jet_electronEnergy;
  vector<Float16_t>	    Jet_muonEnergy;
  vector<Float16_t>	    Jet_HFHadronEnergy;
  vector<Float16_t>	    Jet_HFEMEnergy;
  vector<Float16_t>	    Jet_HOEnergy;
  vector<Float16_t>	    Jet_chargedHadronMultiplicity;
  vector<Float16_t>     Jet_neutralHadronMultiplicity;
  vector<Float16_t>	    Jet_photonMultiplicity;
  vector<Float16_t>	    Jet_electronMultiplicity;
  vector<Float16_t>	    Jet_muonMultiplicity;
  vector<Float16_t>	    Jet_HFHadronMultiplicity;
  vector<Float16_t>	    Jet_HFEMMultiplicity;
  vector<Float16_t> 	Jet_csv;
  vector<Float16_t> 	Jet_mvaDiscriminator;
  vector<Float16_t>  	Jet_nConstituents;
  vector<bool> Jet_passId;

  //PFCand
  UInt_t n_pfcand;
  vector<Float16_t> PFcand_pt;
  vector<Float16_t> PFcand_eta;
  vector<Float16_t> PFcand_phi;
  vector<Float16_t>	PFcand_m;
  vector<Float16_t>	PFcand_pdgid;
  vector<Float16_t>	PFcand_vertex;

  // Fatjets 
  UInt_t n_fatjet;
  vector<Float16_t> FatJet_area;
  vector<Float16_t> FatJet_eta;
  vector<Float16_t> FatJet_n2b1;
  vector<Float16_t> FatJet_n3b1;
  vector<Float16_t> FatJet_phi;
  vector<Float16_t> FatJet_pt;
  vector<Float16_t> FatJet_tau1;
  vector<Float16_t> FatJet_tau2;
  vector<Float16_t> FatJet_tau3;
  vector<Float16_t> FatJet_tau4;
  vector<Float16_t> FatJet_tau21;
  vector<Float16_t> FatJet_tau32;
  vector<Float16_t> FatJet_mass;
  vector<Float16_t> FatJet_msoftdrop;
  vector<Float16_t> FatJet_mtrim;
  vector<Float16_t> FatJet_nconst;

  // Primary vertices
  UInt_t n_pvs;
  vector<Float16_t> Vertex_x;
  vector<Float16_t> Vertex_y;
  vector<Float16_t> Vertex_z;
  vector<Float16_t> Vertex_tracksSize;
  vector<Float16_t> Vertex_chi2;
  vector<Float16_t> Vertex_ndof;
  vector<Float16_t> Vertex_isValidVtx;
        
  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;

};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  muonsToken               (consumes<std::vector<ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken           (consumes<std::vector<ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken             (consumes<std::vector<ScoutingPhoton> >           (iConfig.getParameter<edm::InputTag>("photons"))), 
  pfcandsToken             (consumes<std::vector<ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
  pfjetsToken              (consumes<std::vector<ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))), 
  verticesToken            (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("vertices"))),
//  pileupInfoToken          (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo"))),
//  gensToken                (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("gens"))),
  //genEvtInfoToken          (consumes<GenEventInfoProduct>                    (iConfig.getParameter<edm::InputTag>("geneventinfo"))),    
  doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false),
  hltPSProv_(iConfig,consumesCollector(),*this), //it needs a referernce to the calling module for some reason, hence the *this   
  hltProcess_(iConfig.getParameter<std::string>("hltProcess")),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  l1Seeds_(iConfig.getParameter<std::vector<std::string> >("l1Seeds")),
  hltSeeds_(iConfig.getParameter<std::vector<std::string> >("hltSeeds"))
  //triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
  // triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),

{
 // now do whatever initialization is needed
  usesResource("TFileService");
//  if (doL1) {
//    //algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag"); // might not need
//    //algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_); // might not need
//    //l1GtUtils_ = new l1t::L1TGlobalUtil(iConfig,consumesCollector());	
//  }
//  else {
//    l1Seeds_ = std::vector<std::string>();
//    //l1GtUtils_ = 0;
//  }
//

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree"       , "tree");

  // Event weights
    
  tree->Branch("lumSec"		, &lumSec			 , "lumSec/i" );
  tree->Branch("run"		, &run				 , "run/i" );
    
  // Triggers
<<<<<<< HEAD
  tree->Branch("trig"                 ,&trig       , "trig/b");
  tree->Branch("l1Result"		      ,&l1Result_  );		
=======
  tree->Branch("hltResult"               ,&hltResult_   );              
  tree->Branch("l1Result"		         ,&l1Result_	);		
  tree->Branch("l1Prescale"		         ,&l1Prescale_  );		
>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb
  //Electrons
  tree->Branch("n_ele"            	     ,&n_ele 			, "n_ele/i"		);
  tree->Branch("Electron_pt"             ,&Electron_pt 		 	    );
  tree->Branch("Electron_eta"            ,&Electron_eta 		    );
  tree->Branch("Electron_phi"            ,&Electron_phi 		    );
  tree->Branch("Electron_charge"         ,&Electron_charge          );
  tree->Branch("Electron_m"            	 ,&Electron_m 			    );
  tree->Branch("Electron_tkiso"          ,&Electron_tkiso 		    );
  tree->Branch("Electron_HoE"            ,&Electron_HoE 		    );
  tree->Branch("Electron_sigmaietaieta"  ,&Electron_sigmaietaieta 	);
  tree->Branch("Electron_dphiin"         ,&Electron_dphiin 		    );
  tree->Branch("Electron_detain"         ,&Electron_detain 		    );
  tree->Branch("Electron_mHits"          ,&Electron_mHits 		    );
  tree->Branch("Electron_ooEMOop"        ,&Electron_ooEMOop  	    );

  //Photons
  tree->Branch("n_pho"            	       ,&n_pho 			, "n_pho/i"		);
  tree->Branch("Photon_pt"            	   ,&Photon_pt 			    );
  tree->Branch("Photon_eta"            	   ,&Photon_eta 			);
  tree->Branch("Photon_phi"            	   ,&Photon_phi 			);	
  tree->Branch("Photon_m"            	   ,&Photon_m 			    );
  tree->Branch("Photon_hcaliso"            ,&Photon_hcaliso 		);
  tree->Branch("Photon_ecaliso"            ,&Photon_ecaliso 		);
  tree->Branch("Photon_HoE"            	   ,&Photon_HoE 			);
  tree->Branch("Photon_sigmaietaieta"      ,&Photon_sigmaietaieta	);

  tree->Branch("n_pfcand"            	   ,&n_pfcand 		,"n_pfcand/i"		);	
  tree->Branch("PFcand_pt"        	       ,&PFcand_pt 		 );
  tree->Branch("PFcand_eta"            	   ,&PFcand_eta 	 );
  tree->Branch("PFcand_phi"            	   ,&PFcand_phi		 );
  tree->Branch("PFcand_m"            	   ,&PFcand_m 		 );
  tree->Branch("PFcand_pdgid"              ,&PFcand_pdgid	 );
  tree->Branch("PFcand_vertex"             ,&PFcand_vertex 	 );

  tree->Branch("n_pvs"            	   ,&n_pvs 		,"n_pvs/i"		);	
  tree->Branch("Vertex_x"        	   ,&Vertex_x  		    );
  tree->Branch("Vertex_y"              ,&Vertex_y   	    );
  tree->Branch("Vertex_z"              ,&Vertex_z  		    );
  tree->Branch("Vertex_tracksSize"     ,&Vertex_tracksSize 	);
  tree->Branch("Vertex_chi2"           ,&Vertex_chi2	    );
  tree->Branch("Vertex_ndof"           ,&Vertex_ndof	    );
  tree->Branch("Vertex_isValidVtx"     ,&Vertex_isValidVtx 	);

  tree->Branch("n_mu"            	    , &n_mu 			, "n_mu/i"		);
  tree->Branch("Muon_pt"                , &Muon_pt	);
  tree->Branch("Muon_eta"               , &Muon_eta	);
  tree->Branch("Muon_phi"               , &Muon_phi	);
  tree->Branch("Muon_m"                 , &Muon_m	);
  tree->Branch("Muon_ecaliso"           , &Muon_ecaliso	);
  tree->Branch("Muon_hcaliso"           , &Muon_hcaliso	);
  tree->Branch("Muon_trkiso"            , &Muon_trkiso	);
  tree->Branch("Muon_chi2"              , &Muon_chi2	);
  tree->Branch("Muon_ndof"              , &Muon_ndof	);
  tree->Branch("Muon_charge"            , &Muon_charge	);
  tree->Branch("Muon_dxy"               , &Muon_dxy	);
  tree->Branch("Muon_dz"                , &Muon_dz	);
  tree->Branch("Muon_nvalidmuon_hits"   , &Muon_nvalidmuon_hits	);
  tree->Branch("Muon_validpixelhits"    , &Muon_nvalidpixelhits );
  tree->Branch("Muon_nmatchedstations"  , &Muon_nmatchedstations);
  tree->Branch("Muon_type"              , &Muon_type    );
  tree->Branch("Muon_nvalidstriphits"   , &Muon_nvalidstriphits   );
  tree->Branch("Muon_trkqoverp"         , &Muon_trkqoverp   );
  tree->Branch("Muon_trklambda"         , &Muon_trklambda   );
  tree->Branch("Muon_trkpt"             , &Muon_trkpt    );
  tree->Branch("Muon_trkphi"            , &Muon_trkphi   );
  tree->Branch("Muon_trketa"            , &Muon_trketa   );
  tree->Branch("Muon_trkqoverperror"    , &Muon_trkqoverperror    );
  tree->Branch("Muon_trklambdaerror"    , &Muon_trklambdaerror    );
  tree->Branch("Muon_trkpterror"        , &Muon_trkpterror     );
  tree->Branch("Muon_trkphierror"       , &Muon_trkphierror    );
  tree->Branch("Muon_trketaerror"       , &Muon_trketaerror    );
  tree->Branch("Muon_trkdszerror"       , &Muon_trkdszerror    );
  tree->Branch("Muon_trkdsz"            , &Muon_trkdsz   );


  tree->Branch("ht"                         ,&ht                 );
  tree->Branch("n_jet"            	   	    ,&n_jet 			, "n_jet/i"	  );
  tree->Branch("n_jetId"            	   	,&n_jetId 			, "n_jetId/i" );
  tree->Branch("Jet_pt"            	   	    ,&Jet_pt 				);
  tree->Branch("Jet_eta"            	    ,&Jet_eta 			    );
  tree->Branch("Jet_phi"            	    ,&Jet_phi 			    );
  tree->Branch("Jet_m"            	   	    ,&Jet_m 				);
  tree->Branch("Jet_area"            	    ,&Jet_area			    );
  tree->Branch("Jet_chargedHadronEnergy"    ,&Jet_chargedHadronEnergy 	 );
  tree->Branch("Jet_neutralHadronEnergy"    ,&Jet_neutralHadronEnergy 	 );
  tree->Branch("Jet_photonEnergy"           ,&Jet_photonEnergy 		     );
  tree->Branch("Jet_electronEnergy"         ,&Jet_electronEnergy 		 );
  tree->Branch("Jet_muonEnergy"    		    ,&Jet_muonEnergy 		     );
  tree->Branch("Jet_HFHadronEnergy"         ,&Jet_HFHadronEnergy 		 );
  tree->Branch("Jet_HFEMEnergy"            	,&Jet_HFEMEnergy 		     );
  tree->Branch("Jet_HOEnergy"            	,&Jet_HOEnergy 		         );
  tree->Branch("Jet_chargedHadronMultiplicity"      ,&Jet_chargedHadronMultiplicity 		 );
  tree->Branch("Jet_neutralHadronMultiplicity"      ,&Jet_neutralHadronMultiplicity 		 );
  tree->Branch("Jet_photonMultiplicity"            	,&Jet_photonMultiplicity 		 );
  tree->Branch("Jet_electronMultiplicity"           ,&Jet_electronMultiplicity 		 );
  tree->Branch("Jet_muonMultiplicity"            	,&Jet_muonMultiplicity 		     );
  tree->Branch("Jet_HFHadronMultiplicity"           ,&Jet_HFHadronMultiplicity 		 );
  tree->Branch("Jet_HFEMMultiplicity"            	,&Jet_HFEMMultiplicity 		     );
  tree->Branch("Jet_csv"            	   	,&Jet_csv 		 );
  tree->Branch("Jet_mvaDiscriminator"       ,&Jet_mvaDiscriminator 		 );
  tree->Branch("Jet_nConstituents"           ,&Jet_nConstituents 		 );
  tree->Branch("Jet_passId"           ,&Jet_passId 		 );
  
  tree->Branch("FatJet_area"        ,&FatJet_area   );
  tree->Branch("FatJet_eta"         ,&FatJet_eta    );
  tree->Branch("FatJet_n2b1"        ,&FatJet_n2b1   );
  tree->Branch("FatJet_n3b1"        ,&FatJet_n3b1   );
  tree->Branch("FatJet_phi"         ,&FatJet_phi    );
  tree->Branch("FatJet_pt"          ,&FatJet_pt     );
  tree->Branch("FatJet_tau1"        ,&FatJet_tau1   );
  tree->Branch("FatJet_tau2"        ,&FatJet_tau2   );
  tree->Branch("FatJet_tau3"        ,&FatJet_tau3   );
  tree->Branch("FatJet_tau4"        ,&FatJet_tau4   );
  tree->Branch("FatJet_tau21"       ,&FatJet_tau21  );
  tree->Branch("FatJet_tau32"       ,&FatJet_tau32  );
  tree->Branch("FatJet_mass"        ,&FatJet_mass   );
  tree->Branch("FatJet_msoftdrop"   ,&FatJet_msoftdrop);
  tree->Branch("FatJet_mtrim"       ,&FatJet_mtrim    );
  tree->Branch("FatJet_nconst"      ,&FatJet_nconst   );
  

}


ScoutingNanoAOD::~ScoutingNanoAOD() {
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;
    
  // Handles to the EDM content
  //iEvent.getByToken(triggerBits_, triggerBits);

  //edm::Handle<edm::TriggerResults> triggerResultsH;
  //iEvent.getByToken(triggerResultsToken, triggerResultsH);
    
  Handle<vector<ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);

  Handle<vector<ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  Handle<vector<ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);

  Handle<vector<ScoutingPFJet> > pfjetsH;
  iEvent.getByToken(pfjetsToken, pfjetsH);
    
  Handle<vector<ScoutingParticle> > pfcandsH;
  iEvent.getByToken(pfcandsToken, pfcandsH);

  Handle<vector<ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();


  // Which triggers fired
  hltResult_.clear();

  iEvent.getByToken(triggerBits_, triggerBits);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {                                                          
      const std::string& hltbitName = names.triggerName(i);
      std::string hltpathName = hltbitName;
      bool hltpassFinal = triggerBits->accept(i);

      for(size_t i = 0; i < hltSeeds_.size(); i++){
        TPRegexp pattern(hltSeeds_[i]);
        if( TString(hltpathName).Contains(pattern)){
          hltResult_.push_back(hltpassFinal);
          //std::cout << "HLT Trigger " << hltbitName << " " << hltpassFinal << std::endl;
        }
      }
      
     
  }
  

  // *
  // Electrons here
  // *
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_d0.clear();
  Electron_dz.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_HoE.clear();
  Electron_ooEMOop.clear();
  Electron_mHits.clear();
  Electron_charge.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_tkiso.clear();
  n_ele = 0;
  for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) 
    {
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      Electron_m.push_back(electrons_iter->m());
      Electron_detain.push_back(electrons_iter->dEtaIn());
      Electron_dphiin.push_back(electrons_iter->dPhiIn());
      Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
      Electron_HoE.push_back(electrons_iter->hOverE());	
      Electron_ooEMOop.push_back(electrons_iter->ooEMOop());
      Electron_mHits.push_back(electrons_iter->missingHits());
      Electron_charge.push_back(electrons_iter->charge());
      Electron_tkiso.push_back(electrons_iter->trackIso());
      Electron_ecaliso.push_back(electrons_iter->ecalIso());
      Electron_hcaliso.push_back(electrons_iter->hcalIso());
      n_ele++;
    }

  // *
  // Photons here
  // *
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_m.clear();
  Photon_sigmaietaieta.clear();
  Photon_HoE.clear();
  Photon_ecaliso.clear();
  Photon_hcaliso.clear();
  n_pho = 0;

  for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
    Photon_pt.push_back(photons_iter->pt());
    Photon_eta.push_back(photons_iter->eta());
    Photon_phi.push_back(photons_iter->phi());
    Photon_m.push_back(photons_iter->m());
    Photon_sigmaietaieta.push_back(photons_iter->sigmaIetaIeta());
    Photon_HoE.push_back(photons_iter->hOverE());
    Photon_ecaliso.push_back(photons_iter->ecalIso());
    Photon_hcaliso.push_back(photons_iter->hcalIso());
    
    n_pho++;
  }

  // *
  // Primary vertices
  // * 
  n_pvs = 0;
  Vertex_x.clear();
  Vertex_y.clear();
  Vertex_z.clear();
  Vertex_tracksSize.clear();
  Vertex_chi2.clear();
  Vertex_ndof.clear();
  Vertex_isValidVtx.clear();
    for (auto vertices_iter = verticesH->begin(); vertices_iter != verticesH->end(); ++vertices_iter) {
        Vertex_x.push_back( vertices_iter->x() );
        Vertex_y.push_back( vertices_iter->y() );
        Vertex_z.push_back( vertices_iter->z() );
        Vertex_tracksSize.push_back( vertices_iter->tracksSize() );
        Vertex_chi2.push_back( vertices_iter->chi2() );
        Vertex_ndof.push_back( vertices_iter->ndof() );
        Vertex_isValidVtx.push_back( vertices_iter->isValidVtx() );
        n_pvs++;
    }

  // * 
  // Particle Flow candidates 
  // *
  PFcand_pt.clear();
  PFcand_eta.clear();
  PFcand_phi.clear();
  PFcand_m.clear();
  PFcand_pdgid.clear();
  PFcand_vertex.clear();
  vector<PseudoJet> fj_part;
  n_pfcand = 0;
  for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter) {
    PFcand_pt.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->pt())));
    PFcand_eta.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->eta())));
    PFcand_phi.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->phi())));
    PFcand_m.push_back(pfcands_iter->m());
    PFcand_pdgid.push_back(pfcands_iter->pdgId());
    PFcand_vertex.push_back(pfcands_iter->vertex());

    // Cluster charged PF candidates into fat jets
    if (pfcands_iter->vertex() != 0) continue;
    if (abs(pfcands_iter->eta()) >= 2.4 ) continue;
    if (pfcands_iter->pt() < 1) continue; 
    if (getCharge(pfcands_iter->pdgId()) == 0 ) continue;

    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);
    temp_jet.reset_PtYPhiM(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
    temp_jet.set_user_index(pfcands_iter->pdgId());
    fj_part.push_back(temp_jet);

    n_pfcand++;
  } 

  // 
  // Muons   
  // 
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_chi2.clear();
  Muon_ndof.clear();
  Muon_charge.clear();
  Muon_dxy.clear();
  Muon_dz.clear();
  Muon_nvalidmuon_hits.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nmatchedstations.clear();
  Muon_type.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkpterror.clear();
  Muon_trkphierror.clear();
  Muon_trketaerror.clear();
  Muon_trkdszerror.clear();
  Muon_trkdsz.clear();
  n_mu=0;
  for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
 	Muon_pt.push_back(muons_iter->pt());
   	Muon_eta.push_back(muons_iter->eta());
   	Muon_phi.push_back(muons_iter->phi());
   	Muon_m.push_back(muons_iter->m());
   	Muon_ecaliso.push_back(muons_iter->ecalIso());
   	Muon_hcaliso.push_back(muons_iter->hcalIso());
   	Muon_trkiso.push_back(muons_iter->chi2());
   	Muon_chi2.push_back(muons_iter->ndof());
   	Muon_ndof.push_back(muons_iter->charge());
   	Muon_charge.push_back(muons_iter->dxy());
   	Muon_dxy.push_back(muons_iter->dz());
   	Muon_dz.push_back(muons_iter->nValidMuonHits());
   	Muon_nvalidmuon_hits.push_back(muons_iter->nValidPixelHits());
   	Muon_nvalidpixelhits.push_back(muons_iter->nMatchedStations());
   	Muon_nmatchedstations.push_back(muons_iter->nTrackerLayersWithMeasurement());
    Muon_type.push_back(muons_iter->type());
    Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
    Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
    Muon_trklambda.push_back(muons_iter->trk_lambda());
    Muon_trkpt.push_back(muons_iter->trk_pt());
    Muon_trkphi.push_back(muons_iter->trk_phi());
    Muon_trketa.push_back(muons_iter->trk_eta());
    Muon_trkqoverperror.push_back(muons_iter->dxyError());
    Muon_trklambdaerror.push_back(muons_iter->dzError());
    Muon_trkpterror.push_back(muons_iter->trk_qoverpError());
    Muon_trkphierror.push_back(muons_iter->trk_lambdaError());
    Muon_trketaerror.push_back(muons_iter->trk_phiError());
    Muon_trkdszerror.push_back(muons_iter->trk_dsz());
    Muon_trkdsz.push_back(muons_iter->trk_dszError());
    n_mu++;
 }

  // * 
  // Jets 
  // * 
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_m.clear();
  Jet_area.clear();
  Jet_chargedHadronEnergy.clear();
  Jet_neutralHadronEnergy.clear();
  Jet_photonEnergy.clear();
  Jet_electronEnergy.clear();
  Jet_muonEnergy.clear();
  Jet_HFHadronEnergy.clear();
  Jet_HFEMEnergy.clear();
  Jet_HOEnergy.clear();
  Jet_chargedHadronMultiplicity.clear();
  Jet_neutralHadronMultiplicity.clear();
  Jet_photonMultiplicity.clear();
  Jet_electronMultiplicity.clear();
  Jet_muonMultiplicity.clear();
  Jet_HFHadronMultiplicity.clear();
  Jet_HFEMMultiplicity.clear();
  Jet_csv.clear();
  Jet_mvaDiscriminator.clear();
  Jet_nConstituents.clear();
  Jet_passId.clear();
  n_jet = 0;
  n_jetId = 0;
  ht = 0;
  passJetId = false;
   for (auto pfjet = pfjetsH->begin(); pfjet != pfjetsH->end(); ++pfjet) {

    Jet_pt .push_back( pfjet->pt() );
    Jet_eta.push_back( pfjet->eta());
    Jet_phi.push_back( pfjet->phi());
    Jet_m  .push_back( pfjet->m()  );

    Jet_area.push_back( pfjet->jetArea());

    Jet_chargedHadronEnergy.push_back( pfjet->chargedHadronEnergy());
    Jet_neutralHadronEnergy.push_back( pfjet->neutralHadronEnergy());
    Jet_photonEnergy       .push_back( pfjet->photonEnergy()       );
    Jet_electronEnergy     .push_back( pfjet->electronEnergy()     );
    Jet_muonEnergy         .push_back( pfjet->muonEnergy()     );
    Jet_HFHadronEnergy     .push_back( pfjet->HFHadronEnergy() );
    Jet_HFEMEnergy         .push_back( pfjet->HFEMEnergy()     );
    Jet_HOEnergy           .push_back( pfjet->HOEnergy()       );
    
    Jet_chargedHadronMultiplicity.push_back( pfjet->chargedHadronMultiplicity());
    Jet_neutralHadronMultiplicity.push_back( pfjet->neutralHadronMultiplicity());
    Jet_photonMultiplicity       .push_back( pfjet->photonMultiplicity()       );
    Jet_electronMultiplicity     .push_back( pfjet->electronMultiplicity()     );
    Jet_muonMultiplicity         .push_back( pfjet->muonMultiplicity()         );
    Jet_HFHadronMultiplicity     .push_back( pfjet->HFHadronMultiplicity()     );
    Jet_HFEMMultiplicity         .push_back( pfjet->HFEMMultiplicity()         );

    Jet_csv             .push_back( pfjet->csv() );
    Jet_mvaDiscriminator.push_back( pfjet->mvaDiscriminator()    );
    Jet_nConstituents   .push_back( pfjet->constituents().size() );
    
    n_jet++;

    passJetId = jetID(*pfjet);
    Jet_passId.push_back( passJetId );

    // apply jet ID 
    if ( passJetId == false ) continue; 
    ht += pfjet->pt() ; 
    n_jetId++ ; 
  }

  // * 
  // FatJets 
  // *
  FatJet_area.clear();
  FatJet_eta .clear();
  FatJet_phi .clear();
  FatJet_pt  .clear();
  FatJet_mass.clear();
  FatJet_n2b1.clear();
  FatJet_n3b1.clear();
  FatJet_tau1.clear();
  FatJet_tau2.clear();
  FatJet_tau3.clear();
  FatJet_tau4.clear();
  FatJet_tau21.clear();
  FatJet_tau32.clear();
  FatJet_msoftdrop.clear();
  FatJet_mtrim.clear();
  FatJet_nconst.clear();

  JetDefinition ak15_def = JetDefinition(antikt_algorithm, 1.5);
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
  Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));

  double beta = 1.0;
  Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

  EnergyCorrelatorN2 N2=EnergyCorrelatorN2(1.0);
  EnergyCorrelatorN3 N3=EnergyCorrelatorN3(1.0);

  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  ClusterSequenceArea ak15_cs(fj_part, ak15_def, area_def);
  vector<PseudoJet> ak15_jets = sorted_by_pt(ak15_cs.inclusive_jets(100.0));

  for(auto &j: ak15_jets) {
    FatJet_area.push_back(j.area());
    FatJet_eta .push_back(j.pseudorapidity());
    FatJet_phi .push_back(j.phi_std());
    FatJet_pt  .push_back(j.pt());
    FatJet_mass.push_back(j.m());

    FatJet_nconst.push_back(j.constituents().size());

    PseudoJet sd_ak8 = sd_groomer(j);
    FatJet_msoftdrop.push_back(sd_ak8.m());
    
    PseudoJet trimmed_ak8 = trimmer(j);
    FatJet_mtrim.push_back(trimmed_ak8.m());
    
    // Energy correlation
    FatJet_n2b1.push_back(N2(sd_ak8));
    FatJet_n3b1.push_back(N3(sd_ak8));
    
    // Nsubjettiness, tau 
    FatJet_tau1.push_back(nSub1.result(j));
    FatJet_tau2.push_back(nSub2.result(j));
    FatJet_tau3.push_back(nSub3.result(j));
    FatJet_tau4.push_back(nSub4.result(j));
    FatJet_tau21.push_back(nSub2.result(j)/nSub1.result(j));
    FatJet_tau32.push_back(nSub3.result(j)/nSub2.result(j));
    
    // Jet momentum scaling, rho

    
    // SUEP jet selections 
    //  highest track multiplicty
    //  leading pT
  }
  
 // * 
 // L1 info
 // *
 l1Result_.clear();
 l1Prescale_.clear();

 if (doL1) {

    //I seem to recall this function being slow so perhaps cache for a given lumi 
    //(it only changes on lumi boundaries)  
    //note to the reader, what I'm doing is extremely dangerious (a const cast), never do this!                                     
    //however in this narrow case, it fixes a bug in l1t::L1TGlobalUtil (the method should be const)                                
    //and it is safe for this specific instance                                                                                     
    l1t::L1TGlobalUtil& l1GtUtils = const_cast<l1t::L1TGlobalUtil&> (hltPSProv_.l1tGlobalUtil());


    // For debugging: from https://github.com/Sam-Harper/usercode/blob/09e2252601da473ba02de966930863df57512438/TrigTools/plugins/L1MenuExample.cc
<<<<<<< HEAD
    //std::cout <<"l1 menu: name decisions prescale "<<std::endl;
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    //for(size_t bitNr=0;bitNr<l1GtUtils_->decisionsFinal().size();bitNr++){
    //const std::string& bitName = l1GtUtils_->decisionsFinal()[bitNr].first; // l1GtUtils.decisionsFinal() is of type std::vector<std::pair<std::string,bool> >
=======
    std::cout <<"l1 menu: name decisions prescale "<<std::endl;

    for(size_t bitNr=0;bitNr<l1GtUtils.decisionsFinal().size();bitNr++){
        const std::string& bitName = l1GtUtils.decisionsFinal()[bitNr].first; // l1GtUtils.decisionsFinal() is of type std::vector<std::pair<std::string,bool> >
        bool passInitial = l1GtUtils.decisionsInitial()[bitNr].second; //before masks and prescales, so if we have a 15 GeV electron passing L1_SingleEG10, it will show up as true but will likely not cause a L1 acccept due to the seeds high prescale
        bool passInterm = l1GtUtils.decisionsInterm()[bitNr].second; //after mask (?, unsure what this is)
        bool passFinal = l1GtUtils.decisionsFinal()[bitNr].second; //after masks & prescales, true means it gives a L1 accept to the HLT
        int prescale = l1GtUtils.prescales()[bitNr].second;
        std::cout <<"   "<<bitNr<<" "<<bitName<<" "<<passInitial<<" "<<passInterm<<" "<<passFinal<<" "<<prescale<<std::endl;
        //if (passFinal) std::cout <<"   "<<bitNr<<" "<<bitName<<" "<<passInitial<<" "<<passInterm<<" "<<passFinal<<" "<<prescale<<std::endl;
        for(size_t i = 0; i < l1Seeds_.size(); i++){
          std::string l1Name = l1Seeds_[i];
          std::string pathName = bitName;
          if(bitName.find(l1Name) != std::string::npos ){
             //l1bitmap.push_back(std::make_pair(l1seedsvector[i],passFinal));
             //l1prescalemap.push_back(std::make_pair(l1seedsvector[i],prescale));
             l1Result_  .push_back(passFinal);
             l1Prescale_.push_back(prescale);
          }
        }
    }
>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb


    // Abhijith method
    // Need this! 
    //l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    //
    // not sure what this is for...
    /*	for( int r = 99; r<280; r++){
	string name ("empty");
	bool algoName_ = false;
	algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
	cout << "getAlgNameFromBit = " << algoName_  << endl;
	cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
	}*/

    // Seems like L1 trigger info is messed up...? 
<<<<<<< HEAD
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;	
			
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      l1Result_.push_back( l1htbit );
      }
=======
    //std::cout << "name decision" << std::endl;
    //for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
    //  bool l1htbit = 0;	
	//		
    //  l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
    //  std::cout<<l1Seeds_[iseed]<<"  "<<l1htbit<<std::endl;
    //  l1Result_.push_back( l1htbit );
    //  }
>>>>>>> 19afaa0384fc4f879b8d4f223925d54162886ccb
 }


 tree->Fill();	
	
}


void ScoutingNanoAOD::beginJob() {
  
}

void ScoutingNanoAOD::endJob() {
}

void ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  //triggerPathsVector.push_back("DST_DoubleMu1_noVtx_CaloScouting_v*");
  //triggerPathsVector.push_back("DST_DoubleMu3_noVtx_CaloScouting_v*");
  //triggerPathsVector.push_back("DST_DoubleMu3_noVtx_Mass10_PFScouting_v*");
  //triggerPathsVector.push_back("DST_L1HTT_CaloScouting_PFScouting_v*");
  //triggerPathsVector.push_back("DST_CaloJet40_CaloScouting_PFScouting_v*");
  //triggerPathsVector.push_back("DST_HT250_CaloScouting_v*");
  //triggerPathsVector.push_back("DST_HT410_PFScouting_v*");
  //triggerPathsVector.push_back("DST_HT450_PFScouting_v*");

  //we need to initalise the menu each run (menu can and will change on run boundaries)           
  //HLTConfigProvider hltConfig;
  //bool changedConfig = false;
  //hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

  //for (size_t i = 0; i < triggerPathsVector.size(); i++) {
  //  triggerPathsMap[triggerPathsVector[i]] = -1;
  //}

  //for(size_t i = 0; i < triggerPathsVector.size(); i++){
  //  TPRegexp pattern(triggerPathsVector[i]);
  //  for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
  //    std::string pathName = hltConfig.triggerNames()[j];
  //    if(TString(pathName).Contains(pattern)){
  //  triggerPathsMap[triggerPathsVector[i]] = j;
  //    }
  //  }
  //}


  //we need to initalise the menu each run (menu can and will change on run boundaries)           
  //for L1
  bool changed=false;
  hltPSProv_.init(iRun,iSetup,hltProcess_,changed);
  const l1t::L1TGlobalUtil& l1GtUtils = hltPSProv_.l1tGlobalUtil();
  std::cout <<"l1 menu "<<l1GtUtils.gtTriggerMenuName()<<" version "<<l1GtUtils.gtTriggerMenuVersion()<<" comment "<<std::endl;
  std::cout <<"hlt name "<<hltPSProv_.hltConfigProvider().tableName()<<std::endl;

}

void ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

int ScoutingNanoAOD::getCharge(int pdgId) {
  // following workbook 
  if      (abs(pdgId) == 11 ) return 1; // electron
  else if (abs(pdgId) == 13 ) return 1; // muon
  else if (abs(pdgId) == 211) return 1; // pion
  return 0;
  // 130 = KLong - neutral hadron 
  // 22 = photon 
  // 1 = HF hadron, where HF means forward calo
  // 2 = HF em particle, where HF means forward calo
}
bool ScoutingNanoAOD::jetID(const ScoutingPFJet &pfjet){
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    
    TLorentzVector jet; 
    jet.SetPtEtaPhiM(pfjet.pt(), pfjet.eta(), pfjet.phi(), pfjet.m() );
    
    float NHF  = pfjet.neutralHadronEnergy()/jet.E();
    float NEMF = pfjet.photonEnergy()/jet.E();
    float CHF  = pfjet.chargedHadronEnergy()/jet.E();
    float MUF  = pfjet.muonEnergy()/jet.E();
    float CEMF = pfjet.electronEnergy()/jet.E();
    float NumConst = pfjet.chargedHadronMultiplicity()+pfjet.neutralHadronMultiplicity()+pfjet.photonMultiplicity() + pfjet.electronMultiplicity() + pfjet.muonMultiplicity() + pfjet.HFHadronMultiplicity() + pfjet.HFEMMultiplicity();
    float CHM      = pfjet.chargedHadronMultiplicity() +pfjet.electronMultiplicity() + pfjet.muonMultiplicity(); 
    bool passID = (abs(pfjet.eta())<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );

    return passID;
}

void ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD);
