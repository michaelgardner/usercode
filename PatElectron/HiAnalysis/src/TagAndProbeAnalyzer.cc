// -*- C++ -*-
//
// Package:    TagAndProbeAnalyzer
// Class:      TagAndProbeAnalyzer
// 
/**\class TagAndProbeAnalyzer TagAndProbeAnalyzer.cc based on the code from torsten for Z0 analysis using mu-mu UserCode/tdahms/HiAnalysis/HiZ0/plugins/TagAndProbeAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Torsten Dahms,40 4-A32,+41227671635,
//         Created:  Mon Nov 29 03:13:35 CET 2010
// $Id: TagAndProbeAnalyzer.cc,v 1.2 2013/06/18 08:00:21 mgardner Exp $
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

// For H/E - Iso on SC
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"

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

  bool INFO = false;

class TagAndProbeAnalyzer : public edm::EDAnalyzer {
public:
  explicit TagAndProbeAnalyzer(const edm::ParameterSet&);
  ~TagAndProbeAnalyzer();

  typedef math::XYZTLorentzVector LorentzVector ;
  typedef edm::View<reco::Track> trackCollection ;
  
private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void InitTree();
  void InitEvent();

  double FindCentWeight(int bin);

  bool inAcceptance(const pat::Electron* ele);
  bool inAcceptanceSC(const TLorentzVector &ele);

  bool matchSingle(const pat::Electron* ele);
  bool matchDouble(const pat::Electron* ele);

  bool isElectron(const pat::Electron* ele);
  bool isGoodSC(const reco::SuperCluster &probe);
  
  bool ZinRange(const TLorentzVector &mom_v);
  bool ZinRangeSC(const TLorentzVector &mom_v);

  bool reconstructedSC(const reco::SuperClusterRef probe, const edm::Handle<pat::ElectronCollection> electronColl);
  
  bool matchSC(const pat::Electron* tag, const reco::SuperClusterRef probe);

  void fillTree(const TLorentzVector &mom_v, const TLorentzVector &tag_v, const TLorentzVector &probe_v);
  
  // TFile
  TFile* fOut;

  // TTree
  TTree* eleId;
  TTree* eleRec;
  TTree* eleTrg;
  TTree* fitter_tree;

  CentralityProvider* centrality_;
  int centBin;
  int theCentralityBin;

  // handles
  edm::Handle<CaloTowerCollection> * towersH_;
  edm::Handle<pat::ElectronCollection> tagElectronsDblTrg;
  edm::Handle<pat::ElectronCollection> tagElectronsSglTrg;
  edm::Handle<pat::ElectronCollection> probeElectronsID;
  edm::Handle<pat::ElectronCollection> probeElectronsTrig;
  edm::Handle<reco::SuperClusterCollection> probeElectronsReco;

  std::string  _outputName; 
  int _centMin;
  int _centMax;
  double _ptWeight;
  bool  _isMC;
  bool  _isHI;
  double   max_eta;
  double  min_pt;
  double  min_mass;
  double  max_mass;
  double  min_mass_SC;
  double  max_mass_SC;
  std::string  _dblEleTrg;
  std::string  _sglEleTrg;
  
  edm::InputTag _tagSglTrg;
  edm::InputTag _tagDblTrg;
  edm::InputTag _probeTrg;
  edm::InputTag _probeEleId;
  edm::InputTag _probeRec;
  edm::InputTag _towerMaker;
  
  // data members

  // number of events
  float abseta;
  float eta;
  float pt;

  float weight;
  float mass;
  float tag_eta;
  float tag_pt;
  float pair_abseta;
  float pair_absy;
  float pair_eta;
  float pair_pt;
  float pair_y;

  int pair_PassID;
  int tag_PassID;
  int PassID;

  int pair_SCMatch;
  int tag_SCMatch;
  int SCMatch;

  int pair_HLT;
  int tag_HLT;
  int HLT;
};

//
// constructors and destructor
//
TagAndProbeAnalyzer::TagAndProbeAnalyzer(const edm::ParameterSet& iConfig):
    _outputName(iConfig.getParameter<std::string>("outputName")),
    _centMin(iConfig.getParameter<int>("centMin")),
    _centMax(iConfig.getParameter<int>("centMax")),
    _ptWeight(iConfig.getParameter<double>("ptWeight")),
    _isMC(iConfig.getParameter<bool>("isMC")),
    _isHI(iConfig.getParameter<bool>("isHI")),
    max_eta(iConfig.getParameter<double>("maxEta")),
    min_pt(iConfig.getParameter<double>("minPt")),
    min_mass(iConfig.getParameter<double>("minMass")),
    max_mass(iConfig.getParameter<double>("maxMass")),
    min_mass_SC(iConfig.getParameter<double>("minMassSC")),
    max_mass_SC(iConfig.getParameter<double>("maxMassSC")),
                _dblEleTrg(iConfig.getParameter<std::string>("dblEleTrg")),
                _sglEleTrg(iConfig.getParameter<std::string>("sglEleTrg")),
    _tagSglTrg(iConfig.getParameter<edm::InputTag>("tagSglTrg")),
    _tagDblTrg(iConfig.getParameter<edm::InputTag>("tagDblTrg")),
    _probeTrg(iConfig.getParameter<edm::InputTag>("probeTrg")),
    _probeEleId(iConfig.getParameter<edm::InputTag>("probeEleId")),
    _probeRec(iConfig.getParameter<edm::InputTag>("probeRec")),
    _towerMaker(iConfig.getParameter<edm::InputTag>("towerMaker"))
{
  centrality_ = 0;

   abseta = 0;
  eta= 0;
  pt= 0;
  weight= 0;
  mass= 0;
  tag_eta= 0;
  tag_pt= 0;
  pair_abseta= 0;
  pair_absy= 0;
  pair_eta= 0;
  pair_pt= 0;
  pair_y= 0;

  pair_PassID= 0;
  tag_PassID= 0;
  PassID= 0;

  pair_SCMatch= 0;
  tag_SCMatch= 0;
  SCMatch= 0;

  pair_HLT= 0;
  tag_HLT= 0;
  HLT= 0;
}


TagAndProbeAnalyzer::~TagAndProbeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TagAndProbeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (INFO) cout << "Begin Analyzing" << endl;
  //   using namespace edm;
  InitEvent();

  if (INFO) cout << "A1" << endl;

  cout << "Event with centrality of ";

  if (_isHI) {
    if(!centrality_) centrality_ = new CentralityProvider(iSetup);
    if (INFO) cout << "A2" << endl;
    centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
    if (INFO) cout << "A3" << endl;
    if (_centMin<0) centBin = 0;
    else centBin = centrality_->getBin();

    cout << centBin << endl;

    if ((centBin < _centMin) || (centBin >= _centMax)) {
      if (INFO) cout << centBin << "Not in cent range " << _centMin << " to " << _centMax << endl;
      return;
    }
  } else centBin = 0;

  if (INFO) cout << "A5" << endl;
    
  towersH_ = new edm::Handle<CaloTowerCollection>();
  
  iEvent.getByLabel(_tagSglTrg,tagElectronsSglTrg);
  iEvent.getByLabel(_tagDblTrg,tagElectronsDblTrg);
  iEvent.getByLabel(_probeEleId,probeElectronsID);
  iEvent.getByLabel(_probeTrg,probeElectronsTrig);
  iEvent.getByLabel(_probeRec,probeElectronsReco);
  iEvent.getByLabel(_towerMaker,*towersH_);

  if (INFO) cout << "A" << endl;
  
  if (tagElectronsSglTrg.isValid()) { 
    for(vector<pat::Electron>::const_iterator it=tagElectronsSglTrg->begin();it!=tagElectronsSglTrg->end();++it) {
      const pat::Electron* tag = &(*it);
      if (!(inAcceptance(tag))) continue;
      if (!(matchSingle(tag))) continue;//get rid of those tags not matched to trigger
      if (!(isElectron(tag))) continue;//get rid of those tags that do not pass ElectronID cuts
      
      // Tag and Probe - Trigger
      if (probeElectronsTrig.isValid() && matchSingle(tag)) {
        for(vector<pat::Electron>::const_iterator ip=probeElectronsTrig->begin();ip!=probeElectronsTrig->end();++ip) {

          const pat::Electron* probe = &(*ip);
          if (!(inAcceptance(probe))) continue;
          if (!(isElectron(probe))) continue;//get rid of those probes that do not pass ElectronID cuts

          TLorentzVector tag_v, probe_v, mom_v;
          tag_v.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(),0.0005);
          probe_v.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(),0.0005);
          mom_v = tag_v+probe_v;
          
          if (!(ZinRange(mom_v))) continue;
          fillTree(mom_v, tag_v, probe_v);

          if (matchDouble(probe)) HLT= 1;

          eleTrg->Fill();
          InitEvent();
        }
      }
    }
  }
        
  if (INFO) cout << "B" << endl;
  if (tagElectronsDblTrg.isValid()) { 
    for(vector<pat::Electron>::const_iterator it=tagElectronsDblTrg->begin();it!=tagElectronsDblTrg->end();++it) {
      const pat::Electron* tag = &(*it);
      if (!(inAcceptance(tag))) continue;
      if (!(matchDouble(tag))) continue;//get rid of those tags not matched to trigger
      if (!(isElectron(tag))) continue;//get rid of those tags that do not pass ElectronID cuts
      
    // Tag and Probe - ElectronID
      if (probeElectronsID.isValid()) {
        for(vector<pat::Electron>::const_iterator ip=probeElectronsID->begin();ip!=probeElectronsID->end();++ip) {

          const pat::Electron* probe = &(*ip);
          if (!(inAcceptance(probe))) continue;
          
          TLorentzVector tag_v, probe_v, mom_v;
          tag_v.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(),0.0005);
          probe_v.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(),0.0005);
          mom_v = tag_v+probe_v;
          
          if (!(ZinRange(mom_v))) continue;
          fillTree(mom_v,tag_v, probe_v);

          if (isElectron(probe)) PassID= 1;
          
          eleId->Fill();
          InitEvent();
        }
      }
      
  if (INFO) cout << "C" << endl;
    // Tag and Probe - Reco
      if (probeElectronsReco.isValid()) {
        for(unsigned int ip = 0; ip < probeElectronsReco->size(); ip++) {
          const reco::SuperCluster & probe_v1 = (*probeElectronsReco)[ip];
          reco::SuperClusterRef probe(probeElectronsReco,ip);
          if (matchSC(tag,probe)) continue;
  if (INFO) cout << "C2" << endl;
          //if (!(isGoodSC(probe))) continue;
          if (!(isGoodSC(probe_v1))) continue;
  if (INFO) cout << "C3" << endl;

          math::XYZPoint v(0, 0, 0); // this should be taken from something else...
          math::XYZVector P = probe->energy() * (probe->position() - v).unit();
          double t = sqrt(pow(0.0005,2) + P.mag2());
          TLorentzVector probe_v(P.x(), P.y(), P.z(), t);
          
          TLorentzVector tag_v,mom_v;
          tag_v.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(),0.0005);
          mom_v=tag_v+probe_v;
          
          if (!(inAcceptanceSC(probe_v))) continue;
          if (!(ZinRangeSC(mom_v))) continue;
          fillTree(mom_v, tag_v, probe_v);

          if (reconstructedSC(probe,probeElectronsID)) SCMatch= 1;

          eleRec->Fill();
          InitEvent();
        }
      }
    }
  } else if (INFO) cout << "A1" << endl;
  if (INFO) cout << "D" << endl;

  return;
}

void TagAndProbeAnalyzer::fillTree(const TLorentzVector &mom_v, const TLorentzVector &tag_v, const TLorentzVector &probe_v)
{
  abseta = fabs(probe_v.Eta());
  eta= probe_v.Eta();
  pt= probe_v.Pt();
  if (_isMC && _isHI) weight = _ptWeight*FindCentWeight(centBin);
  else if (_isMC) weight = _ptWeight;
  else weight = 1;

  mass= mom_v.M();

  tag_eta= tag_v.Eta();
  tag_pt= tag_v.Pt();

    if (INFO) cout << "E" << endl;

  pair_abseta= fabs(mom_v.Eta());
  pair_absy= fabs(mom_v.Rapidity());
  pair_eta= mom_v.Eta();
  pair_pt= mom_v.Pt();
  pair_y= mom_v.Rapidity();

  if (mom_v.M()) {
    tag_PassID= 1;
    tag_SCMatch= 1;
    tag_HLT= 1;

    pair_PassID= 1;
    pair_SCMatch= 1;
    pair_HLT= 1;
  }
}

bool TagAndProbeAnalyzer::ZinRange(const TLorentzVector &mom_v)
{
  return (mom_v.M()>min_mass && mom_v.M()<max_mass && fabs(mom_v.Rapidity())<max_eta);
}

bool TagAndProbeAnalyzer::ZinRangeSC(const TLorentzVector &mom_v)
{
  return (mom_v.M()>min_mass_SC && mom_v.M()<max_mass_SC && fabs(mom_v.Rapidity())<max_eta);
}

bool TagAndProbeAnalyzer::inAcceptance(const pat::Electron* ele)
{
  return (ele->pt()>min_pt && fabs(ele->eta())<max_eta);
}

bool TagAndProbeAnalyzer::inAcceptanceSC(const TLorentzVector &ele)
{
  return (ele.Et()>min_pt && fabs(ele.Eta())<max_eta);
}

//bool TagAndProbeAnalyzer::isGoodSC(const reco::SuperClusterRef probe)
bool TagAndProbeAnalyzer::isGoodSC(const reco::SuperCluster &probe)
{
   if (INFO) cout << "C22" << endl;
  EgammaTowerIsolation * towerIso1_ = new EgammaTowerIsolation(0.15,0.,0.,1,towersH_->product());
  if (INFO) cout << "C23" << endl;
  EgammaTowerIsolation * towerIso2_ = new EgammaTowerIsolation(0.15,0.,0.,2,towersH_->product());
  if (INFO) cout << "C24" << endl;

  //double HoE = towerIso1_->getTowerESum(&(*probe))) + towerIso2_->getTowerESum(&(*probe));
  double HoE = towerIso1_->getTowerESum(&probe) + towerIso2_->getTowerESum(&probe);
  if (INFO) cout << "C25" << endl;
  HoE /= probe.energy();

  if (INFO) cout << "C26" << endl;
  return (HoE < 0.2);
 
 // return true;  
}

bool TagAndProbeAnalyzer::isElectron(const pat::Electron* ele)
{
  return (fabs(ele->hadronicOverEm())<0.2 && fabs(ele->sigmaIetaIeta())<0.011 && fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.03 && fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.15 && fabs(ele->dB()) < 0.02 && fabs(1./ele->ecalEnergy()-ele->eSuperClusterOverP()/ele->ecalEnergy()) < 0.1);
}

bool TagAndProbeAnalyzer::matchSingle(const pat::Electron* ele)
{
//  return !(ele->triggerObjectMatchesByPath("HLT_HISinglePhoton15_v*").empty());
//  return !(ele->triggerObjectMatchesByPath("HLT_HISinglePhoton20_v*").empty());
  return !(ele->triggerObjectMatchesByPath(_sglEleTrg).empty());
}

bool TagAndProbeAnalyzer::matchDouble(const pat::Electron* ele)
{
//  return !(ele->triggerObjectMatchesByPath("HLT_HIPhoton15_Photon20_v*").empty());
  return !(ele->triggerObjectMatchesByPath(_dblEleTrg).empty());
}

bool TagAndProbeAnalyzer::reconstructedSC(const reco::SuperClusterRef probe, const edm::Handle<pat::ElectronCollection> electronColl)
{
  bool match = false;

  if (INFO) cout << "F" << endl;

  for(vector<pat::Electron>::const_iterator ip=electronColl->begin();ip!=electronColl->end();++ip) {
    const pat::Electron* electron = &(*ip);
    if (matchSC(electron,probe)) {
      match = true;
      break;
    }
  }

  if (INFO) cout << "G" << endl;

  return match;
}

bool TagAndProbeAnalyzer::matchSC(const pat::Electron* tag, const reco::SuperClusterRef probe)
{
  math::XYZPoint v(0, 0, 0); // this should be taken from something else...

  math::XYZVector T = tag->superCluster()->energy() * (tag->superCluster()->position() - v).unit();
  math::XYZVector P = probe->energy() * (probe->position() - v).unit();

  double t = sqrt(pow(0.0005,2) + T.mag2());
  double p = sqrt(pow(0.0005,2) + P.mag2());

  TLorentzVector tag_v(T.x(), T.y(), T.z(), t);
  TLorentzVector probe_v(P.x(), P.y(), P.z(), p);

  return (tag_v.DeltaR(probe_v) < 0.0001 && fabs(tag_v.Pt() - probe_v.Pt()) < 0.0001);
}

double TagAndProbeAnalyzer::FindCentWeight(int bin)
{
 double NCollArray[40]={1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,521.9120,456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695};
 return(NCollArray[bin]);
}


void
TagAndProbeAnalyzer::InitEvent()
{
  if (INFO) cout << "Initializing Event" << endl;

  abseta = 0;
  eta= 0;
  pt= 0;
  weight= 0;
  mass= 0;
  tag_eta= 0;
  tag_pt= 0;
  pair_abseta= 0;
  pair_absy= 0;
  pair_eta= 0;
  pair_pt= 0;
  pair_y= 0;

  pair_PassID= 0;
  tag_PassID= 0;
  PassID= 0;

  pair_SCMatch= 0;
  tag_SCMatch= 0;
  SCMatch= 0;

  pair_HLT= 0;
  tag_HLT= 0;
  HLT= 0;

  if (INFO) cout << "Done Initializing Event" << endl;
  
  return;
}

void
TagAndProbeAnalyzer::InitTree()
{
  cout << "Initialize Tree" << endl;
  fitter_tree = new TTree("fitter_tree","fitter_tree");

  eleId = new TTree("eleId","My TTree of dielectrons");
  
  eleId->Branch("abseta", &abseta, "abseta/F");
  eleId->Branch("eta", &eta, "eta/F");
  eleId->Branch("pt", &pt, "pt/F");
  eleId->Branch("PassID", &PassID, "PassID/I");
  eleId->Branch("weight", &weight, "weight/F");
  eleId->Branch("mass", &mass, "mass/F");
  eleId->Branch("tag_eta", &tag_eta, "tag_eta/F");
  eleId->Branch("tag_pt", &tag_pt, "tag_pt/F");
  eleId->Branch("tag_PassID", &tag_PassID, "tag_PassID/I");
  eleId->Branch("pair_abseta", &pair_abseta, "pair_abseta/F");
  eleId->Branch("pair_absy", &pair_absy, "pair_absy/F");
  eleId->Branch("pair_eta", &pair_eta, "pair_eta/F");
  eleId->Branch("pair_pt", &pair_pt, "pair_pt/F");
  eleId->Branch("pair_y", &pair_y, "pair_y/F");
  eleId->Branch("pair_PassID", &pair_PassID, "pair_PassID/I");

  eleRec = new TTree("eleRec","My TTree of dielectrons");
  
  eleRec->Branch("abseta", &abseta, "abseta/F");
  eleRec->Branch("eta", &eta, "eta/F");
  eleRec->Branch("pt", &pt, "pt/F");
  eleRec->Branch("SCMatch", &SCMatch, "SCMatch/I");
  eleRec->Branch("weight", &weight, "weight/F");
  eleRec->Branch("mass", &mass, "mass/F");
  eleRec->Branch("tag_eta", &tag_eta, "tag_eta/F");
  eleRec->Branch("tag_pt", &tag_pt, "tag_pt/F");
  eleRec->Branch("tag_SCMatch", &tag_SCMatch, "tag_SCMatch/I");
  eleRec->Branch("pair_abseta", &pair_abseta, "pair_abseta/F");
  eleRec->Branch("pair_absy", &pair_absy, "pair_absy/F");
  eleRec->Branch("pair_eta", &pair_eta, "pair_eta/F");
  eleRec->Branch("pair_pt", &pair_pt, "pair_pt/F");
  eleRec->Branch("pair_y", &pair_y, "pair_y/F");
  eleRec->Branch("pair_SCMatch", &pair_SCMatch, "pair_SCMatch/I");

  eleTrg = new TTree("eleTrg","My TTree of dielectrons");
  
  eleTrg->Branch("abseta", &abseta, "abseta/F");
  eleTrg->Branch("eta", &eta, "eta/F");
  eleTrg->Branch("pt", &pt, "pt/F");
  eleTrg->Branch("HLT", &HLT, "HLT/I");
  eleTrg->Branch("weight", &weight, "weight/F");
  eleTrg->Branch("mass", &mass, "mass/F");
  eleTrg->Branch("tag_eta", &tag_eta, "tag_eta/F");
  eleTrg->Branch("tag_pt", &tag_pt, "tag_pt/F");
  eleTrg->Branch("tag_HLT", &tag_HLT, "tag_HLT/I");
  eleTrg->Branch("pair_abseta", &pair_abseta, "pair_abseta/F");
  eleTrg->Branch("pair_absy", &pair_absy, "pair_absy/F");
  eleTrg->Branch("pair_eta", &pair_eta, "pair_eta/F");
  eleTrg->Branch("pair_pt", &pair_pt, "pair_pt/F");
  eleTrg->Branch("pair_y", &pair_y, "pair_y/F");
  eleTrg->Branch("pair_HLT", &pair_HLT, "pair_HLT/I");
  
  cout << "Done Initializing Tree" << endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TagAndProbeAnalyzer::beginJob()
{
  if (INFO) cout << "Beginning" << endl;
  fOut = new TFile(_outputName.c_str(), "RECREATE");
  InitTree();

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TagAndProbeAnalyzer::endJob() {

  if (INFO) cout << "ZZ1" << endl;

  fOut->cd();
  fOut->mkdir("EleID");
  fOut->cd("EleID");

  eleId->SetName("fitter_tree");
  eleId->Write("fitter_tree");
  
  if (INFO) cout << "ZZ2" << endl;

  fOut->cd();
  fOut->mkdir("EleRec");
  fOut->cd("EleRec");
  
  eleRec->SetName("fitter_tree");
  eleRec->Write("fitter_tree");

  if (INFO) cout << "ZZ3" << endl;

  fOut->cd();
  fOut->mkdir("EleTrg");
  fOut->cd("EleTrg");

  eleTrg->SetName("fitter_tree");
  eleTrg->Write("fitter_tree");

  if (INFO) cout << "ZZ4" << endl;

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagAndProbeAnalyzer);
