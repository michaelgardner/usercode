# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions STARTHI44_V7::All --datatier GEN-SIM-RECODEBUG --eventcontent FEVTDEBUGHLT --geometry Extended --filein=input.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

process.options = cms.untracked.PSet()

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
# pile up rejection
process.load("Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff")
# Centrality for pPb
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

process.GlobalTag.globaltag = 'GR_P_V43D::All'

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality"),
  pPbRunFlip = cms.untracked.uint32(211313)
  )

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hlt_DiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt_DiEle.HLTPaths = ["HLT_PAPhoton20_Photon15_NoCaloIdVL_v*","HLT_PAPhoton20_Photon15_TightCaloIdVL_v*"]

from PatElectron.HiAnalysis.z02EleElePAT_cff import *

z02EleElePAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)

process.source.fileNames = cms.untracked.vstring('/store/data/Run2013A/PPPhoton/RECO/PromptReco-v1/000/211/822/00000/681D2556-0C78-E211-8618-E0CB4E4408E3.root')

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(100)
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

				trigPath = cms.vstring("NoTrigger","HLT_PAPhoton20_Photon15_NoCaloIdVL_v*","HLT_PAPhoton20_NoCaloIdVL_v*"),

				# the first set of triggers should be double, and the second should be single. doubleTrigNum is the number of double triggers you have, and controls where to go from checking double triggers (both elements must be matched) to single triggers
        doubleTrigNum = cms.int32(1),
				
        maxAbsZ = cms.double(24.0),
			 )

# Path and EndPath definitions
process.filter = cms.Sequence(process.hlt_DiEle)

process.pACentrality_step = cms.Sequence(process.siPixelRecHits+process.pACentrality)

process.p = cms.Path(process.PAcollisionEventSelection * process.pACentrality_step * process.hiz0)
#process.p = cms.Path(process.PAcollisionEventSelection * process.pileupVertexFilterCutGplus * process.pACentrality_step * process.hiz0)
#process.p = cms.Path(process.PAcollisionEventSelection * process.pileupVertexFilterCutGplus * process.pACentrality_step * process.hiz0)
#process.output_step = cms.EndPath(process.outTnP)

# Schedule definition
process.schedule = cms.Schedule(process.L1Reco_step,process.Z02EleElePAT,process.p)
#process.schedule = cms.Schedule(process.L1Reco_step,process.ProdEvtPlane,process.Z02EleElePAT,process.p)

for path in process.paths:
    getattr(process,path)._seq = process.filter + getattr(process,path)._seq
