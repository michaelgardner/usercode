#ifndef PatElectron_HiAnalysis_HiZ02EleElePAT_h
#define PatElectron_HiAnalysis_HiZ02EleElePAT_h

// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Electron.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
//
// class decleration
//

class HiZ02EleElePAT : public edm::EDProducer {
 public:
  explicit HiZ02EleElePAT(const edm::ParameterSet&);
  ~HiZ02EleElePAT();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
 private:
  edm::InputTag electrons_;
  edm::InputTag thebeamspot_;
  edm::InputTag thePVs_;
  StringCutObjectSelector<pat::Electron> higherPuritySelection_;
  StringCutObjectSelector<pat::Electron> lowerPuritySelection_; 
  StringCutObjectSelector<reco::Candidate, true> dielectronSelection_;
  bool addMCTruth_;
  GreaterByPt<pat::CompositeCandidate> pTComparator_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

#endif
