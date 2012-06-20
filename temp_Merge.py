import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('STARTHI44_V7::All')

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.load('CmsHi.Analysis2010.CommonFunctions_cff')
process.HeavyIonGlobalParameters = cms.PSet(
        centralityVariable = cms.string("HFtowers"),
        nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
        centralitySrc = cms.InputTag("hiCentrality")
)

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
#process.centralityFilter.selectedBins = [CCCC] # 0 - 10 %
#process.centralityFilter.selectedBins = [0,1,2,3] # 0 - 10 %
#process.centralityFilter.selectedBins = [4,5,6,7] # 10 - 20 %
#process.centralityFilter.selectedBins = [8,9,10,11,12,13,14,15,16,17,18,19] # 20 - 50 %
#process.centralityFilter.selectedBins = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # 50 - 100 %
process.centralityFilter.selectedBins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # MinBias

process.load("MuonAnalysis.TagAndProbe.MC_JPsi_AAAA_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
)

TRACK_CUTS = "track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0"
GLB_CUTS = "isGlobalMuon && globalTrack.normalizedChi2 < 20  && abs(dB) < 3 && abs(track.dz) < 15"#move id cut to tracking efficiency
QUALITY_CUTS =  GLB_CUTS + ' && ' + TRACK_CUTS
IN_ACCEPTANCE = '( abs(eta) < 2.4 && pt >= 4 )'

# CANNOT GET MC WORKING WITH THE PT > 4 condition!!! Have no idea why not!
# Maybe the PATMuonSelector does not keep the correct information?
#####################
process.probeMuonsHighPt = cms.EDFilter("PATMuonSelector",
	src = cms.InputTag("probeMuons"),
	cut = cms.string("pt>4")
)
process.probeMuonsTrkHighPt = cms.EDFilter("PATMuonSelector",
	src = cms.InputTag("probeMuonsTrk"),
	cut = cms.string("pt>4")
)
process.probeMuonsStaHighPt = cms.EDFilter("PATMuonSelector",
	src = cms.InputTag("probeMuonsSta"),
	cut = cms.string("pt>4")
)
#################
# in temp_Trg.py
process.tpPairsTrigNew = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string('2.6 < mass < 3.5 && pt > 6.5'),
	decay = cms.string('tagMuonsSglTrg@+ probeMuons@-')
)

process.tpPairsTrigNewHighPt = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string('2.6 < mass < 3.5 && pt > 6.5'),
	decay = cms.string('tagMuonsSglTrg@+ probeMuonsHighPt@-')
)
# in temp_MuId.py
process.tpPairsMuIDNew = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string('2.6 < mass < 3.5 && pt > 6.5'),
	decay = cms.string('tagMuonsDblTrg@+ probeMuonsTrk@-')
)

process.tpPairsMuIDNewHighPt = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string('2.6 < mass < 3.5 && pt > 6.5'),
	decay = cms.string('tagMuonsDblTrg@+ probeMuonsTrkHighPt@-')
)
# in temp_Trk.py
process.tpPairsStaNew = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string(' 2 < mass < 5 && pt > 6.5'),
	decay = cms.string('tagMuonsDblTrg@+ probeMuonsSta@-')
)

process.tpPairsStaNewHighPt = cms.EDProducer("CandViewShallowCloneCombiner",
	cut = cms.string(' 2 < mass < 5 && pt > 6.5'),
	decay = cms.string('tagMuonsDblTrg@+ probeMuonsStaHighPt@-')
)
#####################

# in temp_Merge.py
# Make the fit tree and save it in the "MuonID" directory
process.MuonTrg = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNew"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
#        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()||!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
        HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
#        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()||!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
        HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      GLB = cms.string("isGlobalMuon"),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsSglTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(True),
#    checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuons"),
		eventWeight = cms.double(DDDD*BBBB)
)

process.MuonTrgHighPt = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNewHighPt"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
#        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()||!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
        HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
#        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()||!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()||!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
        HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
        HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      GLB = cms.string("isGlobalMuon"),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsSglTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(False),
#    checkMotherInUnbiasEff = cms.bool(False),
		allProbes     = cms.InputTag("probeMuonsHighPt"),
		eventWeight = cms.double(DDDD*BBBB)
)

#in temp_Trk.py
process.MuonIDsta = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNew"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsDblTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsStaMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(True),
#    checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
		eventWeight = cms.double(DDDD*BBBB)
)

process.MuonIDstaHighPt = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNewHighPt"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
        isGlb = cms.string("isGlobalMuon &&" + TRACK_CUTS),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsDblTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsStaMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(False),
#    checkMotherInUnbiasEff = cms.bool(False),
    allProbes     = cms.InputTag("probeMuonsStaHighPt"),
		eventWeight = cms.double(DDDD*BBBB)
)

# in temp_MuId.py
# Make the fit tree and save it in the "MuonID" directory
process.MuonID = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIDNew"),
    #tagProbePairs = cms.InputTag("tpPairsTracks"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
        #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
        #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsDblTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsTrkMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(True),
#    checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
		eventWeight = cms.double(DDDD*BBBB)
)

# in temp_MuId.py
# Make the fit tree and save it in the "MuonID" directory
process.MuonIDHighPt = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIDNewHighPt"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
        #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
        #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
        PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("tagMuonsDblTrgMCMatch"),
#    probeMatches  = cms.InputTag("probeMuonsTrkMCMatch"),
#    motherPdgId = cms.int32(443),
#    makeMCUnbiasTree = cms.bool(False),
#    checkMotherInUnbiasEff = cms.bool(False),
    allProbes     = cms.InputTag("probeMuonsTrkHighPt"),
		eventWeight = cms.double(DDDD*BBBB)
)

process.probe_producer = cms.Sequence(
				process.probeMuonsTrkHighPt +
				process.probeMuonsStaHighPt +
				process.probeMuonsHighPt
)

process.pairs_producer = cms.Sequence(
				process.tpPairsMuIDNew +
				process.tpPairsMuIDNewHighPt +
				process.tpPairsStaNew +
				process.tpPairsStaNewHighPt +
				process.tpPairsTrigNew +
				process.tpPairsTrigNewHighPt
)

process.tree_producer = cms.Sequence(
				process.MuonID +
				process.MuonIDHighPt +
				process.MuonIDsta +
				process.MuonIDstaHighPt +
				process.MuonTrg +
				process.MuonTrgHighPt
)

process.tnpSimpleSequence = cms.Sequence(
				process.probe_producer *
				process.pairs_producer *
				process.tree_producer
)

process.tagAndProbe = cms.Path(
        process.centralityFilter *
        process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_PMerge_MC_JPsi_AAAA_CentCCCC.root"))
