# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

# import of standard configurations
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load("Configuration.StandardSequences.RawToDigi_cff")		    # RawToDigi
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

process.options = cms.untracked.PSet()

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

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string(""),
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
process.hlt_DiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt_DiEle.HLTPaths = ["HLT_HIPhoton15_Photon20_v*"]

from PatElectron.HiAnalysis.z02EleElePAT_cff import *

z02EleElePAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)

process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/mgardner/private/Z_ee/Data_RAW/4E347442-1010-E111-9DD9-003048F11CF0.root')

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
    isMC = cms.untracked.bool(False),

    #--
    maxEta = cms.double(1.44),
    maxPt = cms.double(20.),

    NumberOfTriggers = cms.uint32(3),

    trigFilter = cms.vstring("","hltL1sL1DoubleEG*","hltL1sL1SingleEG*"),

    trigPath = cms.vstring("NoTrigger","HLT_HIPhoton15_Photon20_v*","HLT_HISinglePhoton15_v*"),

    # the first set should be double, second should be single. doubleTrigNum is number of double triggers you have, and controls where to go from checking double triggers (both elements must be matched) to single triggers
    doubleTrigNum = cms.int32(1),
				
    maxAbsZ = cms.double(24.0),
    )

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withPF)
#process.reconstruction_step = cms.Path(process.hiPixelSeedReco+process.reconstructionHeavyIons_withPF)

process.filter = cms.Sequence(process.hlt_DiEle+
                              process.superClusterMerger+
                              process.superClusterCands+
                              process.goodSuperClusters+
                              process.electronSuperClusterCounter+
                              process.collisionEventSelection)

process.reconstruction_step = cms.Path(process.hiPixelSeedReco + process.hiElectronSequence)

process.full_reco = cms.Path(process.reconstructionHeavyIons_HcalNZS_withPF)

process.p = cms.Path(process.hiz0)
#process.output_step = cms.EndPath(process.outTnP)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.L1Reco_step,process.Z02EleElePAT,process.p)
#process.schedule = cms.Schedule(process.reconstruction_step,process.L1Reco_step,process.ProdEvtPlane,process.Z02EleElePAT,process.TagAndProbeSta,process.TagAndProbeEleID,process.TagAndProbeTrig,process.output_step,process.p)
#process.schedule = cms.Schedule(process.reconstruction_step,process.Z02EleElePAT,process.output_step)

for path in process.paths:
    getattr(process,path)._seq = process.filter + getattr(process,path)._seq
