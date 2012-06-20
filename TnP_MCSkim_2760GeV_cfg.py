# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("Onia2MuMuPAT")

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI44_V7::All'

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

from CmsHi.Analysis2010.CommonFunctions_cff import * 
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
    centralitySrc = cms.InputTag("hiCentrality"),
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass")
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
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = ["HLT_HIL1DoubleMu0_HighQ_v1",
                              "HLT_HIL2Mu3_NHitQ_v1",
                              "HLT_HIL2Mu7_v1","HLT_HIL2Mu15_v1",
                              "HLT_HIL2DoubleMu0_NHitQ_v1",
                              "HLT_HIL2DoubleMu3_v1",
                              "HLT_HIL3Mu3_v1",
                              "HLT_HIL3DoubleMuOpen_v1","HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1"
                              ]
#process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HISIGNAL")
process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True

from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *

#onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=True, HLT="HLT", Filter=False)
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)

process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = False
process.onia2MuMuPatGlbGlb.resolvePileUpAmbiguity = False
#process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.source.fileNames = cms.untracked.vstring(
    "dummy.root")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.outOnia2MuMu.fileName = cms.untracked.string( 'onia2MuMuPAT_Data.root' )
process.outTnP.fileName = cms.untracked.string( 'tnp_MC_JPsi.root' )

process.outOnia2MuMu.outputCommands.extend(cms.untracked.vstring('keep *_generator_*_*'))

process.outTnP.outputCommands.extend(cms.untracked.vstring('keep *_generator_*_*'))

process.e = cms.EndPath(process.outTnP)

process.patMuonSequence = cms.Sequence(
    process.genMuons *
    process.patMuonsWithTriggerSequence
    )

process.schedule = cms.Schedule(process.L1Reco_step,
                                process.TagAndProbeSta, process.TagAndProbeMuID, process.TagAndProbeTrig,
                                process.e)


