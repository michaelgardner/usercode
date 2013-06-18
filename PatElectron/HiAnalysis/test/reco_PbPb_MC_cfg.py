# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
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

process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')
process.load("RecoHI.HiEgammaAlgos.HiElectronSequence_cff")                 # gsf electrons

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('hiReco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

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

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    splitLevel = cms.untracked.int32(0),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECODEBUG'),
        filterName = cms.untracked.string('')),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    #fileName = cms.untracked.string('file:/afs/cern.ch/work/m/mgardner/private/Z_ee/Data_RAW/SD_PhotonHI_MC.root')
    fileName = cms.untracked.string('SD_PhotonHI_MC.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filter_step'))
    )

process.output.outputCommands.extend(['keep *_towerMaker_*_*','keep *_ecalRecHit_*_*','keep *_TriggerResults_*_*','keep *_hltTriggerSummaryAOD_*_*','keep *_offlineBeamSpot_*_*','keep *_hiCentrality_*_*','keep *_hiSelectedVertex_*_*','keep *_l1extraParticles_*_*','keep *_ecalDrivenGsfElectrons_*_*','keep *_ecalDrivenGsfElectronCores_*_*','keep *_electronGsfTracks_*_*','keep recoSuperClusters_*_*_*','keep *_hltL1GtObjectMap_*_*','keep *_gtDigis_*_*'])

# Other statements
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI44_V12::All'

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

process.source = cms.Source("PoolSource",
    inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop *_hiEvtPlane_*_*"),
    fileNames = cms.untracked.vstring()
    )

process.source.fileNames = cms.untracked.vstring(
# 'root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/mgardner/NoDYGenSimZee2040/v3/GenSimMixedZee_9_1_YNh.root'
'/store/user/lamia/HidjetQuenchedMinBias/GenSimZeeNoDY_Hydjet1p8_40_60/f66175446ef3029f23597f878c9d9ae4/GenSimMixedZee_100_2_q0A.root'
#     'file:/afs/cern.ch/user/m/mgardner/Zee_production/CMSSW_4_4_5_patch1/src/Configuration/GenProduction/test/GenSimMixedZee.root'
#'/store/group/phys_heavyions/dileptons/mgardner/Z_ee/Gen_4060/mgardner/HidjetQuenchedMinBias/NoDYGenSimZee_Hydjet1p8_40_60/f66175446ef3029f23597f878c9d9ae4/GenSimMixedZee_9_1_0n2.root'
)

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hlt_DiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt_DiEle.HLTPaths = ["HLT_HIPhoton15_Photon20_v*","HLT_HISinglePhoton15_v*"] # not positive about version
process.hlt_DiEle.throw = False

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(10)
)

process.filter = cms.Sequence(process.hlt_DiEle+
                              process.superClusterMerger+
                              process.superClusterCands+
                              process.goodSuperClusters+
                              process.electronSuperClusterCounter+
                              process.collisionEventSelection)

# Path and EndPath definitions
process.filter_step = cms.Path(process.filter)
process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withPF)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withRegitMu)

process.output_step = cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.L1Reco_step,process.ProdEvtPlane,process.Z02EleElePAT,process.p)
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.filter_step,process.output_step)
