import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

PDFName = "cbPlusExpo"

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring("tnp_JPsi_MC_Merge_ALL.root"),
    InputDirectoryName = cms.string("MuonTrgHighPt"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("tnp_JPsi_Ana_MC_WT_Trg_ALL_total.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    binsForMassPlots = cms.uint32(52),
    WeightVariable = cms.string("weight"),
    
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        #mass = cms.vstring("Tag-Probe Mass", "2.8", "3.3", "GeV/c^{2}"),
        mass = cms.vstring("Tag-Probe Mass", "2.6", "3.5", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        weight = cms.vstring("weight","0.0","10000.0",""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        HLTL1 = cms.vstring("HLTL1", "dummy[true=1,false=0]"),
    ),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
		PDFs = cms.PSet(
			cbPlusExpo = cms.vstring(
				"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.01,0.1], alpha[1.0, 0.2, 3.0], n[2, 0.5, 100.])",
				"Exponential::backgroundPass(mass, lp[0,-5,5])",
				"Exponential::backgroundFail(mass, lf[0,-5,5])",
				"efficiency[0.9,0,1]",
				"signalFractionInPassing[0.9]"
			),
		),
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        HLTL1_pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("HLTL1","true"),
            UnbinnedVariables = cms.vstring("mass","weight"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(0, 3,	6,	9,	12,	20),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        HLTL1_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("HLTL1","true"),
            UnbinnedVariables = cms.vstring("mass","weight"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.4, -1.6, -0.8, 0.8, 1.6, 2.4),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        HLTL1_1bin_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("HLTL1","true"),
            UnbinnedVariables = cms.vstring("mass","weight"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.4, 2.4),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        HLTL1_1bin_pt_eta = cms.PSet(
          EfficiencyCategoryAndState = cms.vstring("HLTL1","true"),
          UnbinnedVariables = cms.vstring("mass","weight"),
          BinnedVariables = cms.PSet(
            pt = cms.vdouble(0, 20),
            eta = cms.vdouble(-2.4, 2.4),
          ),  
          BinToPDFmap = cms.vstring(PDFName)
        ),
    )
)

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

