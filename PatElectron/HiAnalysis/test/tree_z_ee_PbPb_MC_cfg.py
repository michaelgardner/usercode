# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('ZEE')

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

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltZ0HI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltZ0HI.HLTPaths = ["HLT_HIPhoton15_Photon20_v*","HLT_HISinglePhoton15_v*"] # not positive about version
#process.hltZ0HI.HLTPaths = ["HLT_HISinglePhoton15_v1"] # not positive about version

from PatElectron.HiAnalysis.z02EleElePAT_cff import *

z02EleElePAT(process, GlobalTag=process.GlobalTag.globaltag, MC=True, HLT="HLT", Filter=False)

process.source.fileNames = cms.untracked.vstring(
'/store/group/phys_heavyions/dileptons/mgardner/Z_ee/Reco/4060/SD_PhotonHI_MC_9_1_B2M.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.makePatElectrons = cms.Sequence(process.patElectrons)

process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)

recoCollection = "ecalDrivenGsfElectrons"
vertexCollection = "hiSelectedVertex"
genCollection = "hiGenParticles"

# change the values used in the pat creation.
process.patElectrons.electronSource = cms.InputTag(recoCollection)
process.patElectrons.addElectronID = cms.bool(False)
process.patElectrons.pvSrc = cms.InputTag(vertexCollection)

process.electronMatch.src = cms.InputTag(recoCollection)
process.electronMatch.matched = cms.InputTag(genCollection)
process.electronMatch.mcStatus = cms.vint32(3)

process.z02EleElePat.primaryVertexTag = cms.InputTag(vertexCollection)

process.eleTriggerMatchHLT.matchedCuts = cms.string('(path("HLT_HIPhoton15_Photon20_v*",0,0) && filter("hltL1sL1DoubleEG5BptxAND")) || (path("HLT_HISinglePhoton20_v*",0,0) && filter("hltL1sL1SingleEG5BptxAND"))')

process.hiz0 = cms.EDAnalyzer('HiZeeAnalyzer',
    srcElectron = cms.InputTag("patElectronsWithTrigger"),
    src = cms.InputTag("z02EleElePat"),
    genParticles = cms.InputTag(genCollection),
    primaryVertexTag = cms.InputTag(vertexCollection),

    #-- Reco Details
    recoElectron = cms.InputTag(recoCollection),
    histFileName = cms.string("Z_ee_Hist.root"),
    applyCuts = cms.bool(True),
    storeSameSign = cms.untracked.bool(True),
    isHI = cms.untracked.bool(True),

    #-- Histogram configuration
    fillTree = cms.bool(True),
    minimumFlag = cms.bool(True),
    fillSingleElectrons = cms.bool(True),

    #-- Gen Details
    z0PDG = cms.int32(23),
    isMC = cms.untracked.bool(True),

    #--
    maxEta = cms.double(1.44),
    maxPt = cms.double(20.),
				
    NumberOfTriggers = cms.uint32(3),

    trigFilter = cms.vstring("","hltL1sL1DoubleEG*","hltL1sL1SingleEG*"),

    trigPath = cms.vstring("NoTrigger","HLT_HIPhoton15_Photon20_v*","HLT_HISinglePhoton15_v*"),

    # first set should be double, and second should be single. doubleTrigNum is number of double triggers you have, and controls where to go from checking double triggers (both elements must be matched) to single triggers
    doubleTrigNum = cms.int32(1),

    maxAbsZ = cms.double(24.0),
    )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withPF)
process.Z02EleElePAT = cms.Path(
    process.collisionEventSelection *    
    process.pat_step *
    process.pat_trigger_step *
    process.z02EleElePat *
    process.z02EleElePatFilter
)
process.p = cms.Path(process.hiz0)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.L1Reco_step,process.ProdEvtPlane,process.Z02EleElePAT,process.p)
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.L1Reco_step,process.Z02EleElePAT,process.p)
