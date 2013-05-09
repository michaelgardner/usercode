import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *

def z02EleElePAT(process, GlobalTag='STARTHI44_V12::All', MC=False, HLT='HLT', Filter=True):
    # Setup the process
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        # fileMode = cms.untracked.string('MERGE'),
    )
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 100
    process.load('Configuration.StandardSequences.GeometryExtended_cff')
    process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    ## Standard PAT Configuration File
    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    ## load tau sequences up to selectedPatElectrons
    process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cff")
    process.load("PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi")
    # produce missing l1extraParticles
    process.load('Configuration.StandardSequences.L1Reco_cff')

    process.GlobalTag.globaltag = GlobalTag
		
    # Drop the DQM stuff on input
    process.source = cms.Source("PoolSource",
        inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop *_hiEvtPlane_*_*"),
        fileNames = cms.untracked.vstring()
    )
    
    # Merge CaloElectrons and General Tracks into the collection of reco::Electrons
    # Prune generated particles to electrons and their parents
    process.genElectrons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("hiGenParticles"),
        select = cms.vstring(
            "drop  *  ",                     # this is the default
            "++keep abs(pdgId) = 11",        # keep electrons and their parents
            "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
        )
    )

	# BSC or HF coincidence (masked unprescaled L1 bits)
    process.load('L1Trigger.Skimmer.l1Filter_cfi')
    process.bscOrHfCoinc = process.l1Filter.clone(
		    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
		    )

    # Common offline event selection
    process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

    # Pat elements
    from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

    # Make PAT Electrons
    import PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi
    process.patTrigger = PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi.patTrigger.clone()
    #process.patTrigger.onlyStandAlone = True

    pathTriggerEle ='(path("HLT_HIPhoton15_Photon20_v*",0,0) && filter("hltL1sL1DoubleEG5BptxAND")) || (path("HLT_HISinglePhoton20_v*",0,0) && filter("hltL1sL1SingleEG5BptxAND"))'
    
    process.eleTriggerMatchHLT = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
        src     = cms.InputTag( "selectedPatElectrons" ),
        matched = cms.InputTag( "patTrigger" ),
        matchedCuts = cms.string(pathTriggerEle),
        #matchedCuts = cms.string(""),
        maxDPtRel   = cms.double( 10.0 ),
        maxDeltaR   = cms.double( 0.5 ),
        resolveAmbiguities    = cms.bool( False ), # changed from True
        resolveByMatchQuality = cms.bool( False ) # changed from True
    )
    
    process.patElectronsWithTrigger = cms.EDProducer("PATTriggerMatchElectronEmbedder",
       src     = cms.InputTag( "selectedPatElectrons" ),
       matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'))
    )

    # change the values used in the pat creation.
    process.patElectrons.electronSource = cms.InputTag("ecalDrivenGsfElectrons")
    process.patElectrons.addElectronID = cms.bool(False)
    process.patElectrons.pvSrc = cms.InputTag("hiSelectedVertex")

    process.electronMatch.src = cms.InputTag("ecalDrivenGsfElectrons")
    process.electronMatch.matched = cms.InputTag("hiGenParticles")
    process.electronMatch.mcStatus = cms.vint32(3)

    # Make a sequence
    process.pat_step = cms.Sequence(
        process.makePatElectrons * process.selectedPatElectrons
        )

    process.pat_trigger_step = cms.Sequence(
        process.patTrigger * process.eleTriggerMatchHLT * process.patElectronsWithTrigger
        )
      
    # Make dielectron candidates
    process.z02EleElePat = cms.EDProducer('HiZ02EleElePAT',
        electrons = cms.InputTag("patElectronsWithTrigger"),
        beamSpotTag = cms.InputTag("offlineBeamSpot"),
        primaryVertexTag = cms.InputTag("hiSelectedVertex"),
        # At least one electron must pass this selection
        higherPuritySelection = cms.string("abs(eta)<2.4 && pt>10"),
        # BOTH electrons must pass this selection
        lowerPuritySelection  = cms.string("abs(eta)<2.4 && pt>10"),
        dielectronSelection  = cms.string("20 < mass"), ## The dielectron must pass this selection before vertexing
        addMCTruth = cms.bool(MC), ## Add the common MC mother of the two electrons, if any
        resolvePileUpAmbiguity = cms.bool(False)   ## Order PVs by their vicinity to the Z vertex, not by sumPt
    )

    # check if there is at least one (inclusive) tracker+tracker di-electron
    process.z02EleElePatFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('z02EleElePat'),
        minNumber = cms.uint32(1),
    )

    process.patElectronSequence = cms.Sequence(
        process.pat_step *
        process.pat_trigger_step
    )
	
    # the z02EleEle path
    process.Z02EleElePAT = cms.Path(
        process.pat_step *
        process.pat_trigger_step *
        process.z02EleElePat *
        process.z02EleElePatFilter
    )
