// -*- C++ -*-
//
// Package:    HiZeeAnalyzer
// Class:      HiZeeAnalyzer
// 
/**\class HiZeeAnalyzer HiZeeAnalyzer.cc based on the code from torsten for Z0 analysis using mu-mu UserCode/tdahms/HiAnalysis/HiZ0/plugins/HiZeeAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Torsten Dahms,40 4-A32,+41227671635,
//         Created:  Mon Nov 29 03:13:35 CET 2010
// $Id: HiZeeAnalyzer.cc,v 1.2 2013/06/18 08:00:21 mgardner Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/EgammaCandidates/interface/Electron.h>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "PatElectron/HiAnalysis/interface/MyCommonHistoManager.h"

// adding Event Plane by dmoon 
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


//
// -------------------------
//

using namespace std;
using namespace reco;
using namespace trigger;
using namespace edm;
using namespace IPTools;
//using namespace math;
//
// class declaration
//

  //bool DEBUG = true;
  bool DEBUG = false;
  //bool TRIG = true;
  bool TRIG = false;
  bool SHORT = false;

class HiZeeAnalyzer : public edm::EDAnalyzer {
public:
  explicit HiZeeAnalyzer(const edm::ParameterSet&);
  ~HiZeeAnalyzer();

  typedef math::XYZTLorentzVector LorentzVector;
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void InitEvent();
  void InitTree();

  void makeCuts(int sign) ;
  bool checkCuts(const pat::CompositeCandidate* cand, const pat::Electron* electron1,  const pat::Electron* electron2, bool(HiZeeAnalyzer::* callFunc1)(const pat::Electron*), bool(HiZeeAnalyzer::* callFunc2)(const pat::Electron*)); 

  void fillGenInfo();

  void fillRecoElectrons( );
  bool isElectronInAccept(const pat::Electron* aElectron);
  bool isElectronInAccept(const reco::GsfElectron* aElectron);

  bool selGlobalElectron(const pat::Electron* aElectron);

  void fillRecoHistos(int lastSign);

  void fillTreeElectron(const pat::Electron* electron, int trigBits);
  void fillTreeZ(int iSign, int count);

  void checkTriggers(const pat::CompositeCandidate* aZCand);

  TLorentzVector lorentzMomentum(const reco::Candidate::LorentzVector& p);
  // ----------member data ---------------------------
  enum StatBins {
    BIN_nEvents = 0,
    BIN_HLT_HIPhoton15_Photon20 = 1,
    BIN_HLT_HISinglePhoton15 = 2
  };

  // TFile
  TFile* fOut;

  // TTree
  TTree* myTree;

  TClonesArray* Reco_ele_4mom;
  TClonesArray* Reco_ele_3vec;

  TClonesArray* Reco_QQ_4mom;
  TClonesArray* Reco_QQ_SC_4mom;
  TClonesArray* Reco_QQ_Tk_4mom;
  TClonesArray* Reco_QQ_ClTk_4mom;

  TClonesArray* Reco_QQ_elepl_4mom;
  TClonesArray* Reco_QQ_elepl_SC_4mom;
  TClonesArray* Reco_QQ_elepl_Tk_4mom;
  TClonesArray* Reco_QQ_elepl_ClTk_4mom;

  TClonesArray* Reco_QQ_elemi_4mom;
  TClonesArray* Reco_QQ_elemi_SC_4mom;
  TClonesArray* Reco_QQ_elemi_Tk_4mom;
  TClonesArray* Reco_QQ_elemi_ClTk_4mom;

  TClonesArray* Gen_ele_4mom;
  TClonesArray* Gen_ele_3vec;
  TClonesArray* Gen_QQ_4mom;
  TClonesArray* Gen_QQ_elepl_4mom;
  TClonesArray* Gen_QQ_elemi_4mom;

  static const int Max_QQ_size = 10000;
  static const int Max_ele_size = 100;

  int Gen_QQ_size; // number of generated Z0
  
  int Gen_ele_size; // number of generated electrons
  int Gen_ele_charge[10000]; // electron charge

  int Reco_QQ_size;       // Number of reconstructed Z0 
  int Reco_QQ_sign[10000];   /* Ele Ele combinations sign:
           0 = +/- (signal)
           1 = +/+
           2 = -/- 
        */
  int Reco_QQ_trig[10000];      // Vector of trigger bits matched to the Z0
  float Reco_QQ_VtxProb[10000]; // chi2 probability of vertex fitting 

  int Reco_QQ_elepl_charge[10000];
  float Reco_QQ_elepl_he[10000]; //  h/e: keep < 0.15 (pp); keep < 0.2 (PbPb)
  float Reco_QQ_elepl_sigmaietaieta[10000]; // shower shape variable: keep < 0.01 EB; keep < 0.035 EE (pPb)
  float Reco_QQ_elepl_deltaetain[10000]; // difference in eta between supercluster and inner tracker: keep < 0.03
  float Reco_QQ_elepl_deltaphiin[10000]; // difference in phi between supercluster and inner tracker: keep < 0.15
  float Reco_QQ_elepl_eseedpout[10000];     // Vector of eseedpout
  float Reco_QQ_elepl_ep[10000];    // Vector of ep
  float Reco_QQ_elepl_eseedp[10000];     // Vector of eseedp
  float Reco_QQ_elepl_eelepout[10000];    // Vector of eelepout
  float Reco_QQ_elepl_ecalE[10000];    // Vector of ecalE
  float Reco_QQ_elepl_trackP[10000];    // Vector of trackP
  float Reco_QQ_elepl_iso03Tk[10000];    // Vector of iso03Tk
  float Reco_QQ_elepl_iso03Ecal[10000];    // Vector of iso03Ecal
  float Reco_QQ_elepl_iso03Hcal[10000];    // Vector of iso03Hcal
  float Reco_QQ_elepl_dxy[10000];    // Vector of dxy
  float Reco_QQ_elepl_dz[10000];    // Vector of dz
  
  int Reco_QQ_elemi_charge[10000];
  float Reco_QQ_elemi_he[10000]; //  h/e: keep < 0.15 (pp); keep < 0.2 (PbPb)
  float Reco_QQ_elemi_sigmaietaieta[10000]; // shower shape variable: keep < 0.01 EB; keep < 0.035 EE (pPb)
  float Reco_QQ_elemi_deltaetain[10000]; // difference in eta between supercluster and inner tracker: keep < 0.03
  float Reco_QQ_elemi_deltaphiin[10000]; // difference in phi between supercluster and inner tracker: keep < 0.15
  float Reco_QQ_elemi_eseedpout[10000];     // Vector of eseedpout
  float Reco_QQ_elemi_ep[10000];    // Vector of ep
  float Reco_QQ_elemi_eseedp[10000];     // Vector of eseedp
  float Reco_QQ_elemi_eelepout[10000];    // Vector of eelepout
  float Reco_QQ_elemi_ecalE[10000];    // Vector of ecalE
  float Reco_QQ_elemi_trackP[10000];    // Vector of trackP
  float Reco_QQ_elemi_iso03Tk[10000];    // Vector of iso03Tk
  float Reco_QQ_elemi_iso03Ecal[10000];    // Vector of iso03Ecal
  float Reco_QQ_elemi_iso03Hcal[10000];    // Vector of iso03Hcal
  float Reco_QQ_elemi_dxy[10000];    // Vector of dxy
  float Reco_QQ_elemi_dz[10000];    // Vector of dz

  int Reco_ele_size;           // Number of reconstructed electrons
  int Reco_ele_trig[100];      // Vector of trigger bits matched to the electrons
  float Reco_ele_ptErr[100];   // Vector of err on pt of electrons
  float Reco_ele_phiErr[100];  // Vector of err on phi of electrons
  float Reco_ele_etaErr[100];  // Vector of err on eta of electrons
  float Reco_ele_d0[100];      // Vector of d0 of electrons
  float Reco_ele_d0err[100];   // Vector of d0err of electrons
  float Reco_ele_dz[100];      // Vector of dz of electrons
  float Reco_ele_dzerr[100];   // Vector of dzerr of electrons
  int Reco_ele_charge[100];  // Vector of charge of electrons
  float Reco_ele_normChi2[100];   // Vector of chi2/ndof of electrons
  int Reco_ele_nhitsCSC[100];    // Vector of number of valid hits of electrons
  int Reco_ele_nhitsDT[100];    // Vector of number of valid hits of electrons
  int Reco_ele_nhitsTrack[100];    // Vector of number of valid hits of electrons
  float Reco_ele_caloComp[100];    // Vector of calorimeter compatibilities
  float Reco_ele_segmComp[100];    // Vector of electron segment compatibilities 
  float Reco_ele_iso[100];    // Vector of isolations (NOW ONLY SUMPt OF TRACKS) 
  int Reco_ele_nhitsStrip[100];  // Vectors of strip/pixel hits
  int Reco_ele_nhitsPixB[100];
  int Reco_ele_nhitsPixE[100];
  int Reco_ele_nhitsPix1Hit[100];
  int Reco_ele_nhitsPix1HitBE[100];

  float Reco_ele_he[100];    // Vector of he
  float Reco_ele_sigmaietaieta[100];    // Vector of sigmaietaieta
  float Reco_ele_eseedpout[100];     // Vector of eseedpout
  float Reco_ele_ep[100];    // Vector of ep
  float Reco_ele_eseedp[100];     // Vector of eseedp
  float Reco_ele_eelepout[100];    // Vector of eelepout
  float Reco_ele_deltaetain[100];   // Vector of deltaetain
  float Reco_ele_deltaphiin[100];   // Vector of deltaphiin
  float Reco_ele_sigmaetaeta[100];   // Vector of sigmaetaeta
  float Reco_ele_e15[100];    // Vector of e15
  float Reco_ele_e25max[100];    // Vector of e25max
  float Reco_ele_e55[100];    // Vector of e55
  float Reco_ele_fbrem[100];    // Vector of fbrem
  float Reco_ele_mva[100];   // Vector of mva

  int Reco_ele_isbarrel[100];   // Vector of isbarrel
  int Reco_ele_isendcap[100];   // Vector of isendcap
  
  // event counters
  TH1F* hStats;

  // centrality
  TH1F *hCent;

  // number of primary vertices
  TH1F* hPileUp;

  // z vertex distribution
  TH1F* hZVtx;

  // centrality
  CentralityProvider* centrality_;
  int centBin;

  // handles
  edm::Handle<pat::CompositeCandidateCollection> collZ;
  edm::Handle<pat::ElectronCollection> collElectron;
  edm::Handle<reco::GsfElectronCollection> recoElectron;

  edm::Handle<reco::GenParticleCollection> collGenParticles;

  // data members
  edm::InputTag       _patElectron;
  edm::InputTag       _patZ;
  edm::InputTag       _genParticle;
  edm::InputTag       _recoElectron;
  edm::InputTag       _thePVs;
  std::string         _histfilename;

  bool           _applycuts;
  bool           _storeSs;
  bool           _fillTree;
  bool           _theMinimumFlag;
  bool           _fillSingleElectrons;
  bool           _isHI;
  bool           _isMC;

  int _z0PDG;

  double _maxEta;
  double _maxPt;

  std::vector<std::string> theTriggerNames;
  std::vector<std::string> HLTLastFilters;
  int _doubleTrigNum;
  
  std::vector<const pat::CompositeCandidate*>   _thePassedCands[3];

  // number of events
  unsigned int nEvents;

  unsigned int runNb;
  unsigned int eventNb;
  unsigned int lumiSection;

  math::XYZPoint RefVtx;
  float zVtx;
  float nPV;

 // Triger stuff
  // PUT HERE THE *LAST FILTERS* OF THE BITS YOU LIKE
  static const unsigned int sNTRIGGERS = 8;
  unsigned int NTRIGGERS;
  // MC 8E29
  bool isTriggerMatched[sNTRIGGERS];
  bool alreadyFilled[sNTRIGGERS];
  int HLTriggers;

  const edm::ParameterSet _iConfig;
};

//
// constructors and destructor
//
HiZeeAnalyzer::HiZeeAnalyzer(const edm::ParameterSet& iConfig):
  _patElectron(iConfig.getParameter<edm::InputTag>("srcElectron")), // usually filled with patElectronsWithTrigger
  _patZ(iConfig.getParameter<edm::InputTag>("src")),
  _genParticle(iConfig.getParameter<edm::InputTag>("genParticles")),
  _recoElectron(iConfig.getParameter<edm::InputTag>("recoElectron")),
  _thePVs(iConfig.getParameter<edm::InputTag>("primaryVertexTag")),
  _histfilename(iConfig.getParameter<std::string>("histFileName")),
  _applycuts(iConfig.getParameter<bool>("applyCuts")),
  _storeSs(iConfig.getUntrackedParameter<bool>("storeSameSign",false)),
  _fillTree(iConfig.getParameter<bool>("fillTree")),  
  _theMinimumFlag(iConfig.getParameter<bool>("minimumFlag")),  
  _fillSingleElectrons(iConfig.getParameter<bool>("fillSingleElectrons")),
  _isHI(iConfig.getUntrackedParameter<bool>("isHI",true) ),
  _isMC(iConfig.getUntrackedParameter<bool>("isMC",false) ),
  _z0PDG(iConfig.getParameter<int>("z0PDG")),
  _maxEta(iConfig.getParameter<double>("maxEta")),
  _maxPt(iConfig.getParameter<double>("maxPt")),
  theTriggerNames(iConfig.getParameter<vector <std::string> >("trigPath")),
  HLTLastFilters(iConfig.getParameter<vector <std::string> >("trigFilter")),
  _doubleTrigNum(iConfig.getParameter<int>("doubleTrigNum")),
  NTRIGGERS(iConfig.getParameter<uint32_t>("NumberOfTriggers")),
  _iConfig(iConfig)
{
   //now do what ever initialization is needed
  nEvents = 0;
  centrality_ = 0;

  isTriggerMatched[0]=true; // first entry 'hardcoded' true to accept "all" events
  for (unsigned int iTr = 1; iTr<NTRIGGERS; ++iTr) {
    isTriggerMatched[iTr] = false;
  }

}


HiZeeAnalyzer::~HiZeeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HiZeeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   using namespace edm;
  InitEvent();

  nEvents++;
  hStats->Fill(BIN_nEvents);
   
  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lumiSection = iEvent.luminosityBlock();
  
  edm::Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(_thePVs, privtxs);
  reco::VertexCollection::const_iterator privtx;

  nPV = privtxs->size();
   
  if ( privtxs->begin() != privtxs->end() ) {
    privtx=privtxs->begin();
    RefVtx = privtx->position();
  } else {
    RefVtx.SetXYZ(0.,0.,0.);
  }

  zVtx = RefVtx.Z();

  hZVtx->Fill(zVtx);
  if (fabs(zVtx) > _iConfig.getParameter< double > ("maxAbsZ")) return;
  hPileUp->Fill(nPV);

  if (_isHI){
    if(!centrality_) centrality_ = new CentralityProvider(iSetup);
    centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
    centBin = centrality_->getBin();
  } else centBin = 0;

  hCent->Fill(centBin);

  iEvent.getByLabel(_patZ,collZ); 
  iEvent.getByLabel(_patElectron,collElectron);
  iEvent.getByLabel(_recoElectron,recoElectron);

  if (_isMC) {
    iEvent.getByLabel(_genParticle,collGenParticles);
    fillGenInfo();
  }
  
  // APPLY CUTS
  int lastSign = 0;
  makeCuts(0);
  if (_storeSs) {
    makeCuts(1);
    makeCuts(2);
    lastSign = 2;
  }

  if (_fillSingleElectrons)
    fillRecoElectrons( );

  fillRecoHistos(lastSign);

  if (_fillTree)
    myTree->Fill();

  return;
}

void
HiZeeAnalyzer::fillRecoHistos(int lastSign) {
  for (int iSign = 0; iSign <= lastSign; ++iSign) {
    if (DEBUG) cout << "there were " << _thePassedCands[iSign].size() << " passed Cands" << endl;
    for( unsigned int count = 0; count < _thePassedCands[iSign].size(); count++) { 
      const pat::CompositeCandidate* aZCand = _thePassedCands[iSign].at(count); 

      checkTriggers(aZCand);
      if (_fillTree) fillTreeZ(iSign, count);
    }
  }

  return;
}

void
HiZeeAnalyzer::fillTreeElectron(const pat::Electron* electron, int trigBits) {
  if (Reco_ele_size >= Max_ele_size) {
    std::cout << "Too many electrons: " << Reco_ele_size << std::endl;
    std::cout << "Maximum allowed: " << Max_ele_size << std::endl;
    return;
  }

  Reco_ele_charge[Reco_ele_size] = electron->charge();

  // if PAT doesn't include these, you can do:
  // most methods not included in reco::Electron, so need GsfElectron
  const reco::GsfElectron *r_electron = dynamic_cast<const reco::GsfElectron *>(electron->originalObject());

  // E_T for the supercluster > 4, is already implemented at the reconstruction stage.
  Reco_ele_he[Reco_ele_size] = r_electron->hadronicOverEm(); //  h/e: keep < 0.15 (pp); keep < 0.2 (PbPb)
  Reco_ele_sigmaietaieta[Reco_ele_size] = r_electron->sigmaIetaIeta(); // shower shape variable: keep < 0.01 EB; keep < 0.035 EE (pPb)
  Reco_ele_eseedpout[Reco_ele_size] = r_electron->eSeedClusterOverPout();
  Reco_ele_ep[Reco_ele_size] = r_electron->eSuperClusterOverP();
  Reco_ele_eseedp[Reco_ele_size] = r_electron->eSeedClusterOverP();
  Reco_ele_eelepout[Reco_ele_size] = r_electron->eEleClusterOverPout();
  Reco_ele_deltaetain[Reco_ele_size] = r_electron->deltaEtaSuperClusterTrackAtVtx(); // difference in eta between supercluster and inner tracker: keep < 0.03
  Reco_ele_deltaphiin[Reco_ele_size] = r_electron->deltaPhiSuperClusterTrackAtVtx(); // difference in phi between supercluster and inner tracker: keep < 0.15
  Reco_ele_sigmaetaeta[Reco_ele_size] = r_electron->sigmaEtaEta();
  Reco_ele_e15[Reco_ele_size] = r_electron->e1x5();
  Reco_ele_e25max[Reco_ele_size] = r_electron->e2x5Max();
  Reco_ele_e55[Reco_ele_size] = r_electron->e5x5();
  Reco_ele_fbrem[Reco_ele_size] = r_electron->fbrem();
  Reco_ele_mva[Reco_ele_size] = r_electron->mva();

  Reco_ele_isbarrel[Reco_ele_size] = r_electron->isEB();
  Reco_ele_isendcap[Reco_ele_size] = r_electron->isEE();
  
  TLorentzVector vElectron = lorentzMomentum(electron->p4());
  new((*Reco_ele_4mom)[Reco_ele_size])TLorentzVector(vElectron);

  Reco_ele_trig[Reco_ele_size] = trigBits;

  Reco_ele_size++;
  return;
}

void
HiZeeAnalyzer::fillTreeZ(int iSign, int count) {
  if (Reco_QQ_size >= Max_QQ_size) {
    std::cout << "Too many dielectrons: " << Reco_QQ_size << std::endl;
    std::cout << "Maximum allowed: " << Max_QQ_size << std::endl;
    return;
  }

  const pat::CompositeCandidate* aZCand = _thePassedCands[iSign].at(count);

  const pat::Electron* electron1 = dynamic_cast<const pat::Electron*>(aZCand->daughter("electron1"));
  const pat::Electron* electron2 = dynamic_cast<const pat::Electron*>(aZCand->daughter("electron2"));

  int trigBits=0;
  for (unsigned int iTr=1; iTr<NTRIGGERS; ++iTr) {
    if (isTriggerMatched[iTr]) {
      trigBits += pow(2,iTr-1);
    }
  }

  Reco_QQ_sign[Reco_QQ_size] = iSign;

  Reco_QQ_trig[Reco_QQ_size] = trigBits;

  const reco::GsfTrackRef track1 = electron1->gsfTrack();
  const reco::GsfTrackRef track2 = electron2->gsfTrack();

  const reco::TrackRef cl_track1 = electron1->closestCtfTrackRef();
  const reco::TrackRef cl_track2 = electron2->closestCtfTrackRef();

  double cl_px1, cl_py1, cl_pz1, cl_p1;
  double cl_px2, cl_py2, cl_pz2, cl_p2;

  if (cl_track1.isNonnull()) {
    cl_px1 = cl_track1->px();
    cl_py1 = cl_track1->py();
    cl_pz1 = cl_track1->pz();
    cl_p1 = cl_track1->p();
  } else {
    cl_px1 = track1->px();
    cl_py1 = track1->py();
    cl_pz1 = track1->pz();
    cl_p1 = track1->p();
  }

  if (cl_track2.isNonnull()) {
    cl_px2 = cl_track2->px();
    cl_py2 = cl_track2->py();
    cl_pz2 = cl_track2->pz();
    cl_p2 = cl_track2->p();
  } else {
    cl_px2 = track2->px();
    cl_py2 = track2->py();
    cl_pz2 = track2->pz();
    cl_p2 = track2->p();
  }

  math::XYZPoint v(0, 0, 0);

  math::XYZVector sc1 = electron1->superCluster()->energy() * (electron1->superCluster()->position() - v).unit();
  double t = sqrt(pow(0.0005,2) + sc1.mag2());

  math::XYZVector sc2 = electron2->superCluster()->energy() * (electron2->superCluster()->position() - v).unit();
  double u = sqrt(pow(0.0005,2) + sc2.mag2());

  TLorentzVector vElectron1 = lorentzMomentum(electron1->p4());
  TLorentzVector vElectron1_SC(sc1.x(),sc1.y(),sc1.z(),t);
  TLorentzVector vElectron1_Tk(track1->px(),track1->py(),track1->pz(),track1->p());
  TLorentzVector vElectron1_ClTk(cl_px1,cl_py1,cl_pz1,cl_p1);

  TLorentzVector vElectron2 = lorentzMomentum(electron2->p4());
  TLorentzVector vElectron2_SC(sc2.x(),sc2.y(),sc2.z(),u);
  TLorentzVector vElectron2_Tk(track2->px(),track2->py(),track2->pz(),track2->p());
  TLorentzVector vElectron2_ClTk(cl_px2,cl_py2,cl_pz2,cl_p2);

  const reco::GsfElectron *electron_pl;
  const reco::GsfElectron *electron_mi;
  
  if (electron1->charge() > electron2->charge()) {

    new((*Reco_QQ_elepl_4mom)[Reco_QQ_size])TLorentzVector(vElectron1);
    new((*Reco_QQ_elepl_SC_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_SC);
    new((*Reco_QQ_elepl_Tk_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_Tk);
    new((*Reco_QQ_elepl_ClTk_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_ClTk);

    new((*Reco_QQ_elemi_4mom)[Reco_QQ_size])TLorentzVector(vElectron2);
    new((*Reco_QQ_elemi_SC_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_SC);
    new((*Reco_QQ_elemi_Tk_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_Tk);
    new((*Reco_QQ_elemi_ClTk_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_ClTk);

    Reco_QQ_elepl_charge[Reco_QQ_size] = electron1->charge();
    Reco_QQ_elemi_charge[Reco_QQ_size] = electron2->charge();
    
    electron_pl = dynamic_cast<const reco::GsfElectron *>(electron1->originalObject());
    electron_mi = dynamic_cast<const reco::GsfElectron *>(electron2->originalObject());
  }
  else {
    new((*Reco_QQ_elepl_4mom)[Reco_QQ_size])TLorentzVector(vElectron2);
    new((*Reco_QQ_elepl_SC_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_SC);
    new((*Reco_QQ_elepl_Tk_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_Tk);
    new((*Reco_QQ_elepl_ClTk_4mom)[Reco_QQ_size])TLorentzVector(vElectron2_ClTk);

    new((*Reco_QQ_elemi_4mom)[Reco_QQ_size])TLorentzVector(vElectron1);
    new((*Reco_QQ_elemi_SC_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_SC);
    new((*Reco_QQ_elemi_Tk_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_Tk);
    new((*Reco_QQ_elemi_ClTk_4mom)[Reco_QQ_size])TLorentzVector(vElectron1_ClTk);

    Reco_QQ_elepl_charge[Reco_QQ_size] = electron2->charge();
    Reco_QQ_elemi_charge[Reco_QQ_size] = electron1->charge();
    
    electron_pl = dynamic_cast<const reco::GsfElectron *>(electron2->originalObject());
    electron_mi = dynamic_cast<const reco::GsfElectron *>(electron1->originalObject());
  }

  Reco_QQ_elepl_he[Reco_QQ_size] = electron_pl->hadronicOverEm(); //  h/e: keep < 0.15 (pp); keep < 0.2 (PbPb)
  Reco_QQ_elepl_sigmaietaieta[Reco_QQ_size] = electron_pl->sigmaIetaIeta(); // shower shape variable: keep < 0.01 EB; keep < 0.035 EE (pPb)
  Reco_QQ_elepl_deltaetain[Reco_QQ_size] = electron_pl->deltaEtaSuperClusterTrackAtVtx(); // difference in eta between supercluster and inner tracker: keep < 0.03
  Reco_QQ_elepl_deltaphiin[Reco_QQ_size] = electron_pl->deltaPhiSuperClusterTrackAtVtx(); // difference in phi between supercluster and inner tracker: keep < 0.15
  Reco_QQ_elepl_eseedpout[Reco_QQ_size] = electron_pl->eSeedClusterOverPout();
  Reco_QQ_elepl_ep[Reco_QQ_size] = electron_pl->eSuperClusterOverP();
  Reco_QQ_elepl_eseedp[Reco_QQ_size] = electron_pl->eSeedClusterOverP();
  Reco_QQ_elepl_eelepout[Reco_QQ_size] = electron_pl->eEleClusterOverPout();  
  Reco_QQ_elepl_ecalE[Reco_QQ_size] = electron_pl->ecalEnergy();
  Reco_QQ_elepl_trackP[Reco_QQ_size] = electron_pl->ecalEnergy() / electron_pl->eSuperClusterOverP();
  Reco_QQ_elepl_iso03Tk[Reco_QQ_size] = electron_pl->dr03TkSumPt();
  Reco_QQ_elepl_iso03Ecal[Reco_QQ_size] = electron_pl->dr03EcalRecHitSumEt();
  Reco_QQ_elepl_iso03Hcal[Reco_QQ_size] = electron_pl->dr03HcalTowerSumEt();
  Reco_QQ_elepl_dxy[Reco_QQ_size] = electron_pl->gsfTrack()->dxy(RefVtx);
  Reco_QQ_elepl_dz[Reco_QQ_size] = electron_pl->gsfTrack()->dz(RefVtx);

  Reco_QQ_elemi_he[Reco_QQ_size] = electron_mi->hadronicOverEm(); //  h/e: keep < 0.15 (pp); keep < 0.2 (PbPb)
  Reco_QQ_elemi_sigmaietaieta[Reco_QQ_size] = electron_mi->sigmaIetaIeta(); // shower shape variable: keep < 0.01 EB; keep < 0.035 EE (pPb)
  Reco_QQ_elemi_deltaetain[Reco_QQ_size] = electron_mi->deltaEtaSuperClusterTrackAtVtx(); // difference in eta between supercluster and inner tracker: keep < 0.03
  Reco_QQ_elemi_deltaphiin[Reco_QQ_size] = electron_mi->deltaPhiSuperClusterTrackAtVtx(); // difference in phi between supercluster and inner tracker: keep < 0.15
  Reco_QQ_elemi_eseedpout[Reco_QQ_size] = electron_mi->eSeedClusterOverPout();
  Reco_QQ_elemi_ep[Reco_QQ_size] = electron_mi->eSuperClusterOverP();
  Reco_QQ_elemi_eseedp[Reco_QQ_size] = electron_mi->eSeedClusterOverP();
  Reco_QQ_elemi_eelepout[Reco_QQ_size] = electron_mi->eEleClusterOverPout();
  Reco_QQ_elemi_ecalE[Reco_QQ_size] = electron_mi->ecalEnergy();
  Reco_QQ_elemi_trackP[Reco_QQ_size] = electron_mi->ecalEnergy() / electron_pl->eSuperClusterOverP();
  Reco_QQ_elemi_iso03Tk[Reco_QQ_size] = electron_mi->dr03TkSumPt();
  Reco_QQ_elemi_iso03Ecal[Reco_QQ_size] = electron_mi->dr03EcalRecHitSumEt();
  Reco_QQ_elemi_iso03Hcal[Reco_QQ_size] = electron_mi->dr03HcalTowerSumEt();
  Reco_QQ_elemi_dxy[Reco_QQ_size] = electron_mi->gsfTrack()->dxy(RefVtx);
  Reco_QQ_elemi_dz[Reco_QQ_size] = electron_mi->gsfTrack()->dz(RefVtx);
  
  TLorentzVector vZ = lorentzMomentum(aZCand->p4());
  TLorentzVector vZ_SC = TLorentzVector(vElectron1_SC) + TLorentzVector(vElectron2_SC);
  TLorentzVector vZ_Tk = TLorentzVector(vElectron1_Tk) + TLorentzVector(vElectron2_Tk);
  TLorentzVector vZ_ClTk = TLorentzVector(vElectron1_ClTk) + TLorentzVector(vElectron2_ClTk);

  new((*Reco_QQ_4mom)[Reco_QQ_size])TLorentzVector(vZ);
  new((*Reco_QQ_SC_4mom)[Reco_QQ_size])TLorentzVector(vZ_SC);
  new((*Reco_QQ_Tk_4mom)[Reco_QQ_size])TLorentzVector(vZ_Tk);
  new((*Reco_QQ_ClTk_4mom)[Reco_QQ_size])TLorentzVector(vZ_ClTk);

  Reco_QQ_VtxProb[Reco_QQ_size] = aZCand->userFloat("vProb");

  Reco_QQ_size++;
  return;
}

void
HiZeeAnalyzer::checkTriggers(const pat::CompositeCandidate* aZCand) {
  const pat::Electron* electron1 = dynamic_cast<const pat::Electron*>(aZCand->daughter("electron1"));
  const pat::Electron* electron2 = dynamic_cast<const pat::Electron*>(aZCand->daughter("electron2"));

  // Trigger passed
  for (unsigned int iTr = 1; iTr<NTRIGGERS; ++iTr) {
    const pat::TriggerObjectStandAloneCollection ele1HLTMatchesFilter = electron1->triggerObjectMatchesByFilter( HLTLastFilters.at(iTr) );
    const pat::TriggerObjectStandAloneCollection ele2HLTMatchesFilter = electron2->triggerObjectMatchesByFilter( HLTLastFilters.at(iTr) );
    
    const pat::TriggerObjectStandAloneCollection ele1HLTMatchesPath = electron1->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );
    const pat::TriggerObjectStandAloneCollection ele2HLTMatchesPath = electron2->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );
    
    bool pass1 = false;
    bool pass2 = false;

    pass1 = (ele1HLTMatchesPath.size() > 0);
    pass2 = (ele2HLTMatchesPath.size() > 0);
    //pass1 = (ele1HLTMatchesPath.size()||ele1HLTMatchesFilter.size()) > 0;
    //pass2 = (ele2HLTMatchesPath.size()||ele2HLTMatchesFilter.size()) > 0;

    if ((signed)iTr <= _doubleTrigNum) {  // double triggers here
      isTriggerMatched[iTr] = pass1 && pass2;
    } else {        // single triggers here
      isTriggerMatched[iTr] = pass1 || pass2;
    }
  }

  for (unsigned int iTr=1;iTr<NTRIGGERS;++iTr) {
    if (isTriggerMatched[iTr]) {
      // fill event counting histogram only once per event, also if several electrons fired trigger
      if (alreadyFilled[iTr]) continue;
      hStats->Fill(iTr);
      HLTriggers += pow(2,iTr-1);
      alreadyFilled[iTr]=true;
    }
  }
  
  return;
}

void
HiZeeAnalyzer::makeCuts(int sign) { //need to check why there are so few dielectrons
  if (collZ.isValid()) {
    for(std::vector<pat::CompositeCandidate>::const_iterator it=collZ->begin(); it!=collZ->end(); ++it) {
      
      const pat::CompositeCandidate* cand = &(*it);  
      if (fabs(cand->rapidity()) >= _maxEta) continue;

      const pat::Electron* electron1 = dynamic_cast<const pat::Electron*>(cand->daughter("electron1"));
      const pat::Electron* electron2 = dynamic_cast<const pat::Electron*>(cand->daughter("electron2"));

      if (fabs(electron1->rapidity()) >= _maxEta || fabs(electron2->rapidity()) >= _maxEta) continue;
      
      bool thisSign = ( (sign == 0 && electron1->charge() + electron2->charge() == 0) || (sign == 1 && electron1->charge() + electron2->charge() == 2) || (sign == 2 && electron1->charge() + electron2->charge() == -2) );

     if (DEBUG) cout << "thisSign was " << thisSign << endl;

      if (thisSign) {
        // global + global?
        if (checkCuts(cand,electron1,electron2,&HiZeeAnalyzer::selGlobalElectron,&HiZeeAnalyzer::selGlobalElectron)){
          _thePassedCands[sign].push_back(cand);
          continue;
        }
      }
    }
  }
  
  return;
}


bool
HiZeeAnalyzer::checkCuts(const pat::CompositeCandidate* cand, const pat::Electron* electron1,  const pat::Electron* electron2, bool(HiZeeAnalyzer::* callFunc1)(const pat::Electron*), bool(HiZeeAnalyzer::* callFunc2)(const pat::Electron*)) {
  if (DEBUG) cout << "Check dielectron cuts: (this->*callFunc1)(electron1) " << (this->*callFunc1)(electron1) << " (this->*callFunc2)(electron2) " << (this->*callFunc2)(electron2) << " _applycuts" << _applycuts << " cand->userFloat('vProb') " << cand->userFloat("vProb") << endl;

  if ( (  (this->*callFunc1)(electron1) &&  (this->*callFunc2)(electron2) ) && (!_applycuts ||cand->userFloat("vProb") > -0.1) ) //need to figure out what is wrong with vProb, which is being set to 0! will have to look back at creation of z composite candidate
    return true;
  else
    return false;
}


bool
HiZeeAnalyzer::isElectronInAccept(const pat::Electron* aElectron) {
  return (fabs(aElectron->eta()) < _maxEta && aElectron->pt() > _maxPt);
}

bool
HiZeeAnalyzer::isElectronInAccept(const reco::GsfElectron* aElectron) {
  return (fabs(aElectron->eta()) < _maxEta && aElectron->pt() > _maxPt);
}

// can't find Global for electrons
bool
HiZeeAnalyzer::selGlobalElectron(const pat::Electron* aElectron) {
  
  if(!_applycuts) return true;

  const reco::GsfElectron *r_electron = dynamic_cast<const reco::GsfElectron *>(aElectron->originalObject());

  if (DEBUG) cout << "Electron_HoverE" << r_electron->hadronicOverEm() << endl;

  // Z tuned as of 2011-03-18
  return isElectronInAccept(aElectron); // need to check whether hadronicOverEm is a good cut
  //return (isElectronInAccept(aElectron) && r_electron->hadronicOverEm() < 0.2);
}

void
HiZeeAnalyzer::InitEvent()
{
  for (unsigned int iTr=1;iTr<NTRIGGERS;++iTr) {
    alreadyFilled[iTr]=false;
  }
  HLTriggers = 0;

  _thePassedCands[0].clear();
  _thePassedCands[1].clear();
  _thePassedCands[2].clear();

  Reco_QQ_size = 0;
  Reco_ele_size = 0;

  Gen_QQ_size = 0;
  Gen_ele_size = 0;

  Reco_QQ_4mom->Clear();
  Reco_QQ_SC_4mom->Clear();
  Reco_QQ_Tk_4mom->Clear();
  Reco_QQ_ClTk_4mom->Clear();

  Reco_QQ_elepl_4mom->Clear();
  Reco_QQ_elepl_SC_4mom->Clear();
  Reco_QQ_elepl_Tk_4mom->Clear();
  Reco_QQ_elepl_ClTk_4mom->Clear();

  Reco_QQ_elemi_4mom->Clear();
  Reco_QQ_elemi_SC_4mom->Clear();
  Reco_QQ_elemi_Tk_4mom->Clear();
  Reco_QQ_elemi_ClTk_4mom->Clear();

  Reco_ele_4mom->Clear();
  Reco_ele_3vec->Clear();

  if (_isMC) {
    Gen_QQ_4mom->Clear();
    Gen_QQ_elepl_4mom->Clear();
    Gen_QQ_elemi_4mom->Clear();
    Gen_ele_4mom->Clear();
    Gen_ele_3vec->Clear();
  }

  return;
}

void
HiZeeAnalyzer::fillGenInfo()
{
  if (collGenParticles.isValid()) {
    for(std::vector<reco::GenParticle>::const_iterator it=collGenParticles->begin(); it!=collGenParticles->end();++it) {

      if (Gen_QQ_size >= Max_QQ_size) {
        std::cout << "Too many dielectrons: " << Gen_QQ_size << std::endl;
        std::cout << "Maximum allowed: " << Max_QQ_size << std::endl;
        return;
      }

      if (Gen_ele_size >= Max_ele_size) {
        std::cout << "Too many electrons: " << Gen_ele_size << std::endl;
        std::cout << "Maximum allowed: " << Max_ele_size << std::endl;
        return;
      }

      const reco::GenParticle* gen = &(*it);
      if (DEBUG) {
        if (abs(gen->pdgId()) == 11) {
          cout << "Electron " << gen->pdgId() << " status " << gen->status() << " pT " << gen->pt() << " eta " << gen->eta() << " phi "  << gen->phi() << " mass " << gen->mass() << endl;
          cout << "Electron with " << gen->numberOfDaughters() << " children." << endl;
          if (gen->numberOfDaughters() > 0) {
            cout << "Daught1 pdgId " << gen->daughter(0)->pdgId() << " status " << gen->daughter(0)->status() << " pT " << gen->daughter(0)->pt() << " eta " << gen->daughter(0)->eta() << " phi "  << gen->daughter(0)->phi() << " mass " << gen->daughter(0)->mass() << endl;
          }
          if (gen->numberOfDaughters() > 1) {
            cout << "Daught2 pdgId " << gen->daughter(1)->pdgId() << " status " << gen->daughter(1)->status() << " mass " << gen->daughter(0)->mass() << " pT " << gen->daughter(1)->pt() << " eta " << gen->daughter(1)->eta() << " phi "  << gen->daughter(1)->phi() << " mass " << gen->daughter(1)->mass() << endl;
          }
          if (gen->numberOfDaughters() > 2) {
            cout << "Daught3 pdgId " << gen->daughter(2)->pdgId() << " status " << gen->daughter(2)->status() << " pT " << gen->daughter(2)->pt() << " eta " << gen->daughter(2)->eta() << " phi "  << gen->daughter(2)->phi() << " mass " << gen->daughter(2)->mass() << endl;
          }
          if (gen->numberOfDaughters() > 3) {
            cout << "Daught4 pdgId " << gen->daughter(3)->pdgId() << " status " << gen->daughter(3)->status() << " pT " << gen->daughter(3)->pt() << " eta " << gen->daughter(3)->eta() << " phi "  << gen->daughter(3)->phi() << " mass " << gen->daughter(3)->mass() << endl;
          }
        } else if (gen->pdgId() == 23) {
          cout << "Z with " << gen->numberOfDaughters() << " children." << endl;
          cout << "Z pT " << gen->pt() << " eta " << gen->eta() << " phi "  << gen->phi() << " mass " << gen->mass() << endl;
          if (gen->numberOfDaughters() > 0) {
            cout << "Daught1 pdgId " << gen->daughter(0)->pdgId() << " status " << gen->daughter(0)->status() << " pT " << gen->daughter(0)->pt() << " eta " << gen->daughter(0)->eta() << " phi "  << gen->daughter(0)->phi() << " mass " << gen->daughter(0)->mass() << endl;
          }
          if (gen->numberOfDaughters() > 1) {
            cout << "Daught2 pdgId " << gen->daughter(1)->pdgId() << " status " << gen->daughter(1)->status() << " mass " << gen->daughter(0)->mass() << " pT " << gen->daughter(1)->pt() << " eta " << gen->daughter(1)->eta() << " phi "  << gen->daughter(1)->phi() << " mass " << gen->daughter(1)->mass() << endl;
            TLorentzVector vElectron1 = lorentzMomentum(gen->daughter(0)->p4());
            TLorentzVector vElectron2 = lorentzMomentum(gen->daughter(1)->p4());
            TLorentzVector vZ0 = vElectron1 + vElectron2;
               
            cout << "Z2ee pT " << vZ0.Pt() << " Z0 eta " << vZ0.Eta() << " Z0 phi "  << vZ0.Phi() << " Z0 mass " << vZ0.M() << endl;
          }
          if (gen->numberOfDaughters() > 2) {
            cout << "Daught3 pdgId " << gen->daughter(2)->pdgId() << " status " << gen->daughter(2)->status() << " pT " << gen->daughter(2)->pt() << " eta " << gen->daughter(2)->eta() << " phi "  << gen->daughter(2)->phi() << " mass " << gen->daughter(2)->mass() << endl;
          }
          if (gen->numberOfDaughters() > 3) {
            cout << "Daught4 pdgId " << gen->daughter(3)->pdgId() << " status " << gen->daughter(3)->status() << " pT " << gen->daughter(3)->pt() << " eta " << gen->daughter(3)->eta() << " phi "  << gen->daughter(3)->phi() << " mass " << gen->daughter(3)->mass() << endl;
          }
        }
      }

      // if (abs(gen->pdgId()) == _z0PDG  && gen->status() == 2) {
      if (abs(gen->pdgId()) == _z0PDG  && gen->status() == 3) { // so, with the generation of Zs we have right now, the generated status of the Zs that decay into electrons is 3
        TLorentzVector vZ = lorentzMomentum(gen->p4());
        new((*Gen_QQ_4mom)[Gen_QQ_size])TLorentzVector(vZ);
        if (gen->numberOfDaughters() >= 2) {
          const reco::Candidate* genElectron1 = gen->daughter(0);
          const reco::Candidate* genElectron2 = gen->daughter(1);
          if ( abs(genElectron1->pdgId()) == 11 && abs(genElectron2->pdgId()) == 11 && (genElectron1->status() == 1 || genElectron1->status() == 3) && (genElectron2->status() == 1 || genElectron2->status() == 3) ) {
            TLorentzVector vElectron1 = lorentzMomentum(genElectron1->p4());
            TLorentzVector vElectron2 = lorentzMomentum(genElectron2->p4());
      
            if (genElectron1->charge() > genElectron2->charge()) {
              new((*Gen_QQ_elepl_4mom)[Gen_QQ_size])TLorentzVector(vElectron1);
              new((*Gen_QQ_elemi_4mom)[Gen_QQ_size])TLorentzVector(vElectron2);
            }
            else {
              new((*Gen_QQ_elepl_4mom)[Gen_QQ_size])TLorentzVector(vElectron2);
              new((*Gen_QQ_elemi_4mom)[Gen_QQ_size])TLorentzVector(vElectron1);
            }
          }
        }
        Gen_QQ_size++;
      }

      if (abs(gen->pdgId()) == 11  && gen->status() == 3) {
        Gen_ele_charge[Gen_ele_size] = gen->charge();

        TLorentzVector vElectron = lorentzMomentum(gen->p4());
        new((*Gen_ele_4mom)[Gen_ele_size])TLorentzVector(vElectron);

        Gen_ele_size++;
      }
    }
  }

  return;
}

void
HiZeeAnalyzer::fillRecoElectrons()
{
  if (collElectron.isValid()) {
    if (DEBUG||SHORT) {
      if (collElectron->size()!=0) cout << collElectron->size() << " electrons in collElectron" << endl;
      int num = 0;
      for(vector<pat::Electron>::const_iterator it=collElectron->begin();it!=collElectron->end();++it) {
       const pat::Electron* electron = &(*it);
       if (isElectronInAccept(electron)) cout << "Electron " << num++ << " has pt " << electron->pt() << " eta " << electron->eta() << " and phi " << electron->phi() << endl;
      }
      
      if (recoElectron.isValid()) {
        // now look at the gsf electron collection:
        if (recoElectron->size()!=0) cout << recoElectron->size() << " electrons in recoElectron" << endl;
        num = 0;
        for(vector<reco::GsfElectron>::const_iterator it=recoElectron->begin();it!=recoElectron->end();++it) {
         const reco::GsfElectron* electron = &(*it);
         if (isElectronInAccept(electron)) cout << "Electron " << num++ << " has pt " << electron->pt() << " eta " << electron->eta() << " and phi " << electron->phi() << endl;
        }
      } else cout << "UNABLE TO FIND RECO COLLECTION!!!" << endl;
    }
    
    for(vector<pat::Electron>::const_iterator it=collElectron->begin();it!=collElectron->end();++it) {
      const pat::Electron* electron = &(*it);

      if (selGlobalElectron(electron)) {
        int trigBits=0;
        for (unsigned int iTr=1; iTr<NTRIGGERS; ++iTr) {
          if (TRIG) cout << "Filter to check match: " << HLTLastFilters.at(iTr) << endl;
          const pat::TriggerObjectStandAloneCollection eleHLTMatchesFilter = electron->triggerObjectMatchesByFilter(  HLTLastFilters.at(iTr) );
          if (TRIG) cout << "Path to check match: " << theTriggerNames.at(iTr)<< endl;
          const pat::TriggerObjectStandAloneCollection eleHLTMatchesPath = electron->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );
          if (TRIG) cout << "Number of matches to the filter " << eleHLTMatchesFilter.size() << endl;
          if (SHORT) {
            const pat::TriggerObjectStandAloneCollection test = electron->triggerObjectMatchesByFilter(  "*" );
            for (unsigned int i = 0; i < test.size(); i++) {
              std::vector< std::string > filter = test[i].filterLabels();
              std::vector< std::string > path = test[i].pathNames();
              for (unsigned int j = 0; j < filter.size(); j++) {
                cout << "Filter " << filter[j] << endl;
              }
              for (unsigned int j = 0; j < path.size(); j++) {
                cout << "Path " << path[j] << endl;
              }
            }
          }
    // apparently matching by path gives false positives so we use matching by filter for all triggers for which we know the filter name
          if ( eleHLTMatchesPath.size() > 0 ) {
            trigBits += pow(2,iTr-1);
          }
        }
        if (_fillTree) fillTreeElectron(electron, trigBits);
      }
    }
  }
  
  return;
}

void
HiZeeAnalyzer::InitTree()
{
  Reco_ele_4mom = new TClonesArray("TLorentzVector", 100);
  Reco_ele_3vec = new TClonesArray("TVector3", 100);

  Reco_QQ_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_SC_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_Tk_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_ClTk_4mom = new TClonesArray("TLorentzVector",10000);

  Reco_QQ_elepl_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elepl_SC_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elepl_Tk_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elepl_ClTk_4mom = new TClonesArray("TLorentzVector",10000);

  Reco_QQ_elemi_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elemi_SC_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elemi_Tk_4mom = new TClonesArray("TLorentzVector",10000);
  Reco_QQ_elemi_ClTk_4mom = new TClonesArray("TLorentzVector",10000);

  if (_isMC) {
    Gen_ele_4mom = new TClonesArray("TLorentzVector", 2);
    Gen_ele_3vec = new TClonesArray("TVector3", 2);
    Gen_QQ_4mom = new TClonesArray("TLorentzVector", 2);
    Gen_QQ_elepl_4mom = new TClonesArray("TLorentzVector", 2);
    Gen_QQ_elemi_4mom = new TClonesArray("TLorentzVector", 2);
  }

  myTree = new TTree("myTree","My TTree of dielectrons");
  
  myTree->Branch("eventNb", &eventNb,   "eventNb/i");
  myTree->Branch("runNb",   &runNb,     "runNb/i");
  myTree->Branch("LS",      &lumiSection, "LS/i"); 
  myTree->Branch("zVtx",    &zVtx,        "zVtx/F"); 
  myTree->Branch("HLTriggers", &HLTriggers, "HLTriggers/I");
  myTree->Branch("Centrality", &centBin, "Centrality/I");

  myTree->Branch("Reco_QQ_size", &Reco_QQ_size,  "Reco_QQ_size/I");
  myTree->Branch("Reco_QQ_sign", Reco_QQ_sign,   "Reco_QQ_sign[Reco_QQ_size]/I");

  myTree->Branch("Reco_QQ_4mom", "TClonesArray", &Reco_QQ_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_SC_4mom", "TClonesArray", &Reco_QQ_SC_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_Tk_4mom", "TClonesArray", &Reco_QQ_Tk_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_ClTk_4mom", "TClonesArray", &Reco_QQ_ClTk_4mom, 32000, 0);

  myTree->Branch("Reco_QQ_elepl_4mom", "TClonesArray", &Reco_QQ_elepl_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elepl_SC_4mom", "TClonesArray", &Reco_QQ_elepl_SC_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elepl_Tk_4mom", "TClonesArray", &Reco_QQ_elepl_Tk_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elepl_ClTk_4mom", "TClonesArray", &Reco_QQ_elepl_ClTk_4mom, 32000, 0);

  myTree->Branch("Reco_QQ_elemi_4mom", "TClonesArray", &Reco_QQ_elemi_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elemi_SC_4mom", "TClonesArray", &Reco_QQ_elemi_SC_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elemi_Tk_4mom", "TClonesArray", &Reco_QQ_elemi_Tk_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_elemi_ClTk_4mom", "TClonesArray", &Reco_QQ_elemi_ClTk_4mom, 32000, 0);

  myTree->Branch("Reco_QQ_elepl_charge", Reco_QQ_elepl_charge, "Reco_QQ_elepl_charge[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_elepl_he", Reco_QQ_elepl_he, "Reco_QQ_elepl_he[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_sigmaietaieta", Reco_QQ_elepl_sigmaietaieta, "Reco_QQ_elepl_sigmaietaieta[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_deltaetain", Reco_QQ_elepl_deltaetain, "Reco_QQ_elepl_deltaetain[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_deltaphiin", Reco_QQ_elepl_deltaphiin, "Reco_QQ_elepl_deltaphiin[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_eseedpout",Reco_QQ_elepl_eseedpout,"Reco_QQ_elepl_eseedpout[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_ep",Reco_QQ_elepl_ep,"Reco_QQ_elepl_ep[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_eseedp",Reco_QQ_elepl_eseedp,"Reco_QQ_elepl_eseedp[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_eelepout",Reco_QQ_elepl_eelepout,"Reco_QQ_elepl_eelepout[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_ecalE",Reco_QQ_elepl_ecalE,"Reco_QQ_elepl_ecalE[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_trackP",Reco_QQ_elepl_trackP,"Reco_QQ_elepl_trackP[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_iso03Tk",Reco_QQ_elepl_iso03Tk,"Reco_QQ_elepl_iso03Tk[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_iso03Ecal",Reco_QQ_elepl_iso03Ecal,"Reco_QQ_elepl_iso03Ecal[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_iso03Hcal",Reco_QQ_elepl_iso03Hcal,"Reco_QQ_elepl_iso03Hcal[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_dxy",Reco_QQ_elepl_dxy,"Reco_QQ_elepl_dxy[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elepl_dz",Reco_QQ_elepl_dz,"Reco_QQ_elepl_dz[Reco_QQ_size]/F");
  
  myTree->Branch("Reco_QQ_elemi_charge", Reco_QQ_elemi_charge, "Reco_QQ_elemi_charge[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_elemi_he", Reco_QQ_elemi_he, "Reco_QQ_elemi_he[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_sigmaietaieta", Reco_QQ_elemi_sigmaietaieta, "Reco_QQ_elemi_sigmaietaieta[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_deltaetain", Reco_QQ_elemi_deltaetain, "Reco_QQ_elemi_deltaetain[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_deltaphiin", Reco_QQ_elemi_deltaphiin, "Reco_QQ_elemi_deltaphiin[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_eseedpout",Reco_QQ_elemi_eseedpout,"Reco_QQ_elemi_eseedpout[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_ep",Reco_QQ_elemi_ep,"Reco_QQ_elemi_ep[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_eseedp",Reco_QQ_elemi_eseedp,"Reco_QQ_elemi_eseedp[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_eelepout",Reco_QQ_elemi_eelepout,"Reco_QQ_elemi_eelepout[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_ecalE",Reco_QQ_elemi_ecalE,"Reco_QQ_elemi_ecalE[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_trackP",Reco_QQ_elemi_trackP,"Reco_QQ_elemi_trackP[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_iso03Tk",Reco_QQ_elemi_iso03Tk,"Reco_QQ_elemi_iso03Tk[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_iso03Ecal",Reco_QQ_elemi_iso03Ecal,"Reco_QQ_elemi_iso03Ecal[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_iso03Hcal",Reco_QQ_elemi_iso03Hcal,"Reco_QQ_elemi_iso03Hcal[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_dxy",Reco_QQ_elemi_dxy,"Reco_QQ_elemi_dxy[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_elemi_dz",Reco_QQ_elemi_dz,"Reco_QQ_elemi_dz[Reco_QQ_size]/F");
  
  myTree->Branch("Reco_QQ_trig", Reco_QQ_trig,   "Reco_QQ_trig[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_VtxProb", Reco_QQ_VtxProb,   "Reco_QQ_VtxProb[Reco_QQ_size]/F");

  myTree->Branch("Reco_ele_size", &Reco_ele_size,  "Reco_ele_size/I");
  myTree->Branch("Reco_ele_charge", Reco_ele_charge,   "Reco_ele_charge[Reco_ele_size]/I");
  myTree->Branch("Reco_ele_4mom", "TClonesArray", &Reco_ele_4mom, 32000, 0);
  myTree->Branch("Reco_ele_trig", Reco_ele_trig,   "Reco_ele_trig[Reco_ele_size]/I");

  //electron 1
  myTree->Branch("Reco_ele_he",Reco_ele_he,"Reco_ele_he[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_sigmaietaieta",Reco_ele_sigmaietaieta,"Reco_ele_sigmaietaieta[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_eseedpout",Reco_ele_eseedpout,"Reco_ele_eseedpout[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_ep",Reco_ele_ep,"Reco_ele_ep[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_eseedp",Reco_ele_eseedp,"Reco_ele_eseedp[Reco_ele_size]/D");  
  myTree->Branch("Reco_ele_eelepout",Reco_ele_eelepout,"Reco_ele_eelepout[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_deltaetain",Reco_ele_deltaetain,"Reco_ele_deltaetain[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_deltaphiin",Reco_ele_deltaphiin,"Reco_ele_deltaphiin[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_sigmaetaeta",Reco_ele_sigmaetaeta,"Reco_ele_sigmaetaeta[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_e15",Reco_ele_e15,"Reco_ele_e15[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_e25max",Reco_ele_e25max,"Reco_ele_e25max[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_e55",Reco_ele_e55,"Reco_ele_e55[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_fbrem",Reco_ele_fbrem,"Reco_ele_fbrem[Reco_ele_size]/D");
  myTree->Branch("Reco_ele_mva",Reco_ele_mva,"Reco_ele_mva[Reco_ele_size]/D");

  myTree->Branch("Reco_ele_isbarrel",Reco_ele_isbarrel,"Reco_ele_isbarrel[Reco_ele_size]/I");
  myTree->Branch("Reco_ele_isendcap",Reco_ele_isendcap,"Reco_ele_isendcap[Reco_ele_size]/I");    
    
    if (_isMC) {
      myTree->Branch("Gen_QQ_size",      &Gen_QQ_size,    "Gen_QQ_size/I");
    myTree->Branch("Gen_QQ_4mom",      "TClonesArray", &Gen_QQ_4mom, 32000, 0);
    myTree->Branch("Gen_QQ_elepl_4mom", "TClonesArray", &Gen_QQ_elepl_4mom, 32000, 0);
    myTree->Branch("Gen_QQ_elemi_4mom", "TClonesArray", &Gen_QQ_elemi_4mom, 32000, 0);

    myTree->Branch("Gen_ele_size",   &Gen_ele_size,  "Gen_ele_size/I");
    myTree->Branch("Gen_ele_charge", Gen_ele_charge, "Gen_ele_charge[Gen_ele_size]/I");
    myTree->Branch("Gen_ele_4mom",   "TClonesArray", &Gen_ele_4mom, 32000, 0);
  }

  if (!_theMinimumFlag) {
    myTree->Branch("Reco_ele_phiErr",   Reco_ele_phiErr,  "Reco_ele_phiErr[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_etaErr",   Reco_ele_etaErr,  "Reco_ele_etaErr[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_ptErr",    Reco_ele_ptErr,   "Reco_ele_ptErr[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_d0",       Reco_ele_d0,      "Reco_ele_d0[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_d0err",    Reco_ele_d0err,   "Reco_ele_d0err[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_dz",       Reco_ele_dz,      "Reco_ele_dz[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_dzerr",    Reco_ele_dzerr,   "Reco_ele_dzerr[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_normChi2",     Reco_ele_normChi2,    "Reco_ele_normChi2[Reco_ele_size]/F");
    myTree->Branch("Reco_ele_nhitsTrack",    Reco_ele_nhitsTrack,   "Reco_ele_nhitsTrack[Reco_ele_size]/I");      
    myTree->Branch("Reco_ele_nhitsStrip",    Reco_ele_nhitsStrip,   "Reco_ele_nhitsStrip[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsPixB",    Reco_ele_nhitsPixB,   "Reco_ele_nhitsPixB[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsPixE",    Reco_ele_nhitsPixE,   "Reco_ele_nhitsPixE[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsPix1Hit",    Reco_ele_nhitsPix1Hit,   "Reco_ele_nhitsPix1Hit[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsPix1HitBE",    Reco_ele_nhitsPix1HitBE,   "Reco_ele_nhitsPix1HitBE[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsDT",    Reco_ele_nhitsDT,   "Reco_ele_nhitsDT[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_nhitsCSC",    Reco_ele_nhitsCSC,   "Reco_ele_nhitsCSC[Reco_ele_size]/I");
    myTree->Branch("Reco_ele_caloComp",   Reco_ele_caloComp,  "Reco_ele_caloComp[Reco_ele_size]/F"); 
    myTree->Branch("Reco_ele_segmComp",   Reco_ele_segmComp,  "Reco_ele_segmComp[Reco_ele_size]/F"); 
    myTree->Branch("Reco_ele_iso",   Reco_ele_iso,  "Reco_ele_iso[Reco_ele_size]/F");  
  }


}

// ------------ method called once each job just before starting event loop  ------------
void 
HiZeeAnalyzer::beginJob()
{
  fOut = new TFile(_histfilename.c_str(), "RECREATE");
  InitTree();

  // book histos
  hStats = new TH1F("hStats","hStats;;Number of Events",20,0,20);
  hStats->GetXaxis()->SetBinLabel(1,"All");
  for (int i=2; i< (int) theTriggerNames.size()+1; ++i) {
    hStats->GetXaxis()->SetBinLabel(i,theTriggerNames.at(i-1).c_str());
  }
  hStats->Sumw2();

  hCent = new TH1F("hCent","hCent;centrality bin;Number of Events",40,0,40);
  hCent->Sumw2();

  hPileUp = new TH1F("hPileUp","Number of Primary Vertices;n_{PV};counts", 50, 0, 50);
  hPileUp->Sumw2();

  hZVtx = new TH1F("hZVtx","Primary z-vertex distribution;z_{vtx} [cm];counts", 120, -30, 30);
  hZVtx->Sumw2();

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiZeeAnalyzer::endJob() {
  std::cout << "Total number of events = " << nEvents << std::endl;

  fOut->cd();
  hStats->Write();
  hCent->Write();
  hPileUp->Write();
  hZVtx->Write();

  if (_fillTree)
    myTree->Write();

  return;
}

TLorentzVector
HiZeeAnalyzer::lorentzMomentum(const reco::Candidate::LorentzVector& p) {
  TLorentzVector res;
  res.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());

  return res;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiZeeAnalyzer);
