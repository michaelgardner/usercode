# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

# import of standard configurations
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load("Configuration.StandardSequences.RawToDigi_cff")		    # RawToDigi
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
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

# Other statements
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V43D::All'

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowersPlusTrunc"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("pACentrality"),
    pPbRunFlip = cms.untracked.uint32(211313)
    )

# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

from PatElectron.HiAnalysis.z02EleElePAT_cff import *

z02EleElePAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)

process.source.fileNames = cms.untracked.vstring('/store/data/Run2013A/PPPhoton/RECO/PromptReco-v1/000/211/822/00000/681D2556-0C78-E211-8618-E0CB4E4408E3.root')

#process.load("PatElectron.HiAnalysis.dataset_list.RD_AAAA_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.makePatElectrons = cms.Sequence(process.patElectrons)

process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)

recoCollection = "gsfElectrons"
vertexCollection = "offlinePrimaryVertices"
genCollection = "hiGenParticles"

# change the values used in the pat creation.
process.patElectrons.electronSource = cms.InputTag(recoCollection)
process.patElectrons.addElectronID = cms.bool(False)
process.patElectrons.pvSrc = cms.InputTag(vertexCollection)

process.electronMatch.src = cms.InputTag(recoCollection)
process.electronMatch.matched = cms.InputTag(genCollection)
process.electronMatch.mcStatus = cms.vint32(3)

process.z02EleElePat.primaryVertexTag = cms.InputTag(vertexCollection)

process.eleTriggerMatchHLT.matchedCuts = cms.string('(path("HLT_PAPhoton*_Photon*_*_v*",0,0) && filter("hltL1sL1DoubleEG*")) || (path("HLT_PAPhoton*_*_v*",0,0) && filter("hltL1sL1SingleEG*"))')

process.mergedSC =  cms.EDProducer("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'), cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
)

process.mergedSCCands = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("mergedSC"),
    particleType = cms.int32(11),
)

#   Get the above SC's Candidates and place a cut on their Et.
process.goodSC = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("mergedSCCands"),
    cut = cms.string('et > 15.0 && abs(eta)<1.5'),
    filter = cms.bool(True)
    )

process.electronSCCounter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodSC"),
    minNumber = cms.uint32(2),
    filter = cms.bool(True)
    )

ELE_ID = "(abs(hadronicOverEm)<0.12 && abs(sigmaIetaIeta)<0.011 && abs(deltaEtaSuperClusterTrackAtVtx)<0.007 && abs(deltaPhiSuperClusterTrackAtVtx)<0.15 && abs(dB) < 0.02 && abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05)"
TRG1 = "!triggerObjectMatchesByPath('HLT_PAPhoton20_NoCaloIdVL_v*').empty()"
TRG2 = "!triggerObjectMatchesByPath('HLT_PAPhoton20_Photon15_NoCaloIdVL_v*').empty()"
IN_ACCEPTANCE = "(abs(eta)<1.44 && pt>20)"
SC_CUT = "energy > 20"

process.tagElectronsSglTrg = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(ELE_ID + " && " + IN_ACCEPTANCE + " && " + TRG1)
)

process.tagElectronsDblTrg = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(ELE_ID + " && " + IN_ACCEPTANCE + " && " + TRG2)
)

# must be STA, so we can measure inner tracking efficiency
process.probeElectronsID = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(IN_ACCEPTANCE)
)

# must be a GLB to measure trigger efficiency
process.probeElectronsTrig = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(ELE_ID + " && " + IN_ACCEPTANCE)
)

process.probeElectronsReco = cms.EDFilter("SuperClusterSelector",
    src = cms.InputTag("mergedSC"),
    cut = cms.string(SC_CUT)
)

process.probeElectronsSC = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("mergedSCCands"),
    cut = cms.string("et > 15.0 && abs(eta)<2.5")
)

process.TagAndProbe = cms.Path(
    process.mergedSC *
    process.mergedSCCands *
    process.tagElectronsSglTrg *
    process.tagElectronsDblTrg *
    process.probeElectronsID *
    process.probeElectronsTrig *
    process.probeElectronsReco *
    process.probeElectronsSC #this is not used right now, but we'll see.
)

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hlt_DiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt_DiEle.HLTPaths = ["HLT_PAPhoton20_Photon15_NoCaloIdVL_v1","HLT_PAPhoton20_NoCaloIdVL_v*"]

process.filter = cms.Sequence(process.hlt_DiEle+
                              process.mergedSC+
                              process.mergedSCCands+
                              process.goodSC+
                              process.electronSCCounter)
						
process.tnp = cms.EDAnalyzer("TagAndProbeAnalyzer",
    tagSglTrg = cms.InputTag("tagElectronsSglTrg"),
    tagDblTrg = cms.InputTag("tagElectronsDblTrg"),
    probeTrg = cms.InputTag("probeElectronsTrig"),
    probeEleId = cms.InputTag("probeElectronsID"),
    probeRec = cms.InputTag("probeElectronsReco"),
    towerMaker = cms.InputTag("towerMaker"), # used for H/E for SC

    isMC = cms.bool(False),
    isHI = cms.bool(False),

    #--
    maxEta = cms.double(1.44),
    minPt = cms.double(20.),
    minMass = cms.double(60.),
    maxMass = cms.double(120.),
    minMassSC = cms.double(40.),
    maxMassSC = cms.double(140.),

    dblEleTrg = cms.string("HLT_PAPhoton20_Photon15_NoCaloIdVL_v*"),

    sglEleTrg = cms.string("HLT_PAPhoton20_NoCaloIdVL_v*"),

    outputName = cms.string('pat_RD.root'),
    centMin = cms.int32(-1000),
    centMax = cms.int32(1000),
    ptWeight = cms.double(1)
    )

process.p = cms.Path(process.tnp)

# Schedule definition
#process.schedule = cms.Schedule(process.Z02EleElePAT,process.TagAndProbe,process.output_step,process.p)
process.schedule = cms.Schedule(process.Z02EleElePAT,process.TagAndProbe,process.p)

for path in process.paths:
    getattr(process,path)._seq = process.filter + getattr(process,path)._seq
