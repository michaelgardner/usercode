// ///////////////////////////////////////////////////
// This was created by Michael Gardner, to parallel 
// what is done for muons. Creates dielectron object.
// Last updated May 9th, 2013
// ///////////////////////////////////////////////////

#include "PatElectron/HiAnalysis/interface/HiZ02EleElePAT.h"

//Headers for the data items
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/EgammaCandidates/interface/ElectronFwd.h>
#include <DataFormats/EgammaCandidates/interface/Electron.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//Headers for services and tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "PatElectron/HiAnalysis/interface/VertexReProducer.h"

HiZ02EleElePAT::HiZ02EleElePAT(const edm::ParameterSet& iConfig):
  electrons_(iConfig.getParameter<edm::InputTag>("electrons")),
  thebeamspot_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
  thePVs_(iConfig.getParameter<edm::InputTag>("primaryVertexTag")),
  higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
  lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
  dielectronSelection_(iConfig.existsAs<std::string>("dielectronSelection") ? iConfig.getParameter<std::string>("dielectronSelection") : ""),
  addMCTruth_(iConfig.getParameter<bool>("addMCTruth"))
{  
    produces<pat::CompositeCandidateCollection>();  
}


HiZ02EleElePAT::~HiZ02EleElePAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiZ02EleElePAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::auto_ptr<pat::CompositeCandidateCollection> z0Output(new pat::CompositeCandidateCollection);
  
  Vertex thePrimaryV;
  Vertex theBeamSpotV; 

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  Handle< View<pat::Electron> > electrons;
  iEvent.getByLabel(electrons_,electrons);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter;
  
  // Z0 candidates only from electrons
  for(View<pat::Electron>::const_iterator it = electrons->begin(), itend = electrons->end(); it != itend; ++it){
    // both must pass low quality
    if(!lowerPuritySelection_(*it)) continue; 
    for(View<pat::Electron>::const_iterator it2 = it+1; it2 != itend;++it2){
      // both must pass low quality
      if(!lowerPuritySelection_(*it2)) continue; 
      // one must pass tight quality
      if (!(higherPuritySelection_(*it) || higherPuritySelection_(*it2))) continue;

      pat::CompositeCandidate myCand;

      // ---- no explicit order defined ----
      myCand.addDaughter(*it, "electron1");
      myCand.addDaughter(*it2,"electron2");	

      // ---- define and set candidate's 4momentum  ----  
      LorentzVector z0 = it->p4() + it2->p4();
      myCand.setP4(z0);
      myCand.setCharge(it->charge()+it2->charge());

      // ---- apply the dielectron cut ----
      if(!dielectronSelection_(myCand)) continue;

      if (it->gsfTrack().isNonnull() && it2->gsfTrack().isNonnull()) {
        vector<TransientTrack> t_tks;
 	      t_tks.push_back(theTTBuilder->build(*it->gsfTrack()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(*it2->gsfTrack())); // 
        TransientVertex myVertex = vtxFitter.vertex(t_tks);
        if (myVertex.isValid()) {
          float vChi2 = myVertex.totalChiSquared();
          float vNDF  = myVertex.degreesOfFreedom();
          float vProb(TMath::Prob(vChi2,(int)vNDF));

          myCand.addUserFloat("vChi2",vChi2);
          myCand.addUserFloat("vNDF",vNDF);
          myCand.addUserFloat("vNChi2",vChi2/vNDF);
          myCand.addUserFloat("vProb",vProb);
        }
      }

			// ---- MC Truth, if enabled ----
			if (addMCTruth_) {
				reco::GenParticleRef genEle1 = it->genParticleRef();
				reco::GenParticleRef genEle2 = it2->genParticleRef();
				if (genEle1.isNonnull() && genEle2.isNonnull()) {
					if (genEle1->numberOfMothers()>0 && genEle2->numberOfMothers()>0) {
						reco::GenParticleRef mom1 = genEle1->motherRef();
						reco::GenParticleRef mom2 = genEle2->motherRef();
						if (mom1.isNonnull() && (mom1 == mom2)) {
							myCand.setGenParticleRef(mom1); // set
							myCand.embedGenParticle();      // and embed
						}
					} else { // if the electrons aren't matched to the pat, we will just fill the Z's properties.
						Handle<GenParticleCollection>theGenParticles;
						iEvent.getByLabel("hiGenParticles", theGenParticles);
						if (theGenParticles.isValid()){
							for(size_t iGenParticle=0; iGenParticle<theGenParticles->size();++iGenParticle) {
								const Candidate & genCand = (*theGenParticles)[iGenParticle];
								if (genCand.pdgId()==23) {
									reco::GenParticleRef mom1(theGenParticles,iGenParticle);
									myCand.setGenParticleRef(mom1);
									myCand.embedGenParticle();
								}
							}
						}
					}
				}
      }

      // ---- Push back output ----  
      z0Output->push_back(myCand);
    }
  }

  std::sort(z0Output->begin(),z0Output->end(),pTComparator_);

  iEvent.put(z0Output);
	
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiZ02EleElePAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiZ02EleElePAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiZ02EleElePAT);
