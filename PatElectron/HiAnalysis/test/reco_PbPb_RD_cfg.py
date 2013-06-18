# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO2')

# import of standard configurations
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#Reconstruction			

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')
process.load("RecoHI.HiEgammaAlgos.HiElectronSequence_cff") # gsf electrons

process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('hiReco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# reprocess pixel seeds for electron reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)

process.hiPixelSeedReco = cms.Sequence(process.rechits+process.hiPixelClusterVertex+process.hiPixel3ProtoTracks+process.hiPixelMedianVertex+process.hiSelectedProtoTracks+process.hiPixelAdaptiveVertex+process.hiBestAdaptiveVertex+process.hiSelectedVertex+process.hiPixel3PrimTracks+process.hiPixelTrackSeeds)

#   Before selection is made, merge the Barrel and EndCap SC's.
process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag('correctedIslandBarrelSuperClusters'), cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
	 )

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("superClusterMerger"),
    particleType = cms.int32(11),
 )

#   Get the above SC's Candidates and place a cut on their Et.
process.goodSuperClusters = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("superClusterCands"),
    cut = cms.string('et > 15.0 && abs(eta)<1.5'),
    filter = cms.bool(True)
    )

process.electronSuperClusterCounter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodSuperClusters"),
    minNumber = cms.uint32(2),
    filter = cms.bool(True)
    )

# Other statements
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V27A::All'

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

process.source = cms.Source("PoolSource",
    inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop *_hiEvtPlane_*_*"),
    fileNames = cms.untracked.vstring()
    )

process.source.fileNames = cms.untracked.vstring(#'file:/afs/cern.ch/work/m/mgardner/private/Z_ee/Data_RAW/pickevents.root')
'file:/afs/cern.ch/work/m/mgardner/private/Z_ee/Data_RAW/4E347442-1010-E111-9DD9-003048F11CF0.root')

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root'),
#    outputCommands = cms.untracked.vstring('keep *'),
    outputCommands = cms.untracked.vstring('drop *','keep *_towerMaker_*_*','keep *_ecalRecHit_*_*','keep *_TriggerResults_*_*','keep *_hltTriggerSummaryAOD_*_*','keep *_offlineBeamSpot_*_*','keep *_hiCentrality_*_*','keep *_hiSelectedVertex_*_*','keep *_l1extraParticles_*_*','keep *_ecalDrivenGsfElectrons_*_*','keep *_ecalDrivenGsfElectronCores_*_*','keep *_electronGsfTracks_*_*','keep recoSuperClusters_*_*_*','keep *_hltL1GtObjectMap_*_*','keep *_gtDigis_*_*'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('reconstruction_step'))
    )

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hlt_DiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt_DiEle.HLTPaths = ["HLT_HISinglePhoton15_v1","HLT_HIPhoton15_Photon20_v1"]
process.hlt_DiEle.throw = False

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1000)
    input = cms.untracked.int32(-1)
)

# Path and EndPath definitions
process.filter = cms.Sequence(process.hlt_DiEle+
                              process.superClusterMerger+
                              process.superClusterCands+
                              process.goodSuperClusters+
                              process.electronSuperClusterCounter+
                              process.collisionEventSelection)

process.reconstruction_step = cms.Path(process.filter+
                                       process.hiPixelSeedReco+
                                       process.hiElectronSequence)
#process.reconstruction_step = cms.Path(process.hiElectronSequence)

process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.output_step)
