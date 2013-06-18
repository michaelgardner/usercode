# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/Pyquen_GammaJet_pt30_cfi.py -n 10 -s GEN:hiSignal,SIM,DIGI,L1,DIGI2RAW,HLT:HIon --conditions auto:starthi --scenario HeavyIons --geometry Extended --himix --datatier GEN-SIM-RAWDEBUG --eventcontent=RAWDEBUG --filein=input.root --processName HISIGNAL --no_exec

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.register ('randomNumber',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.register ('skipEvents',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.skipEvents = 0
ivars.randomNumber = 1
#ivars.inputFiles = "/store/user/mnguyen//Hydjet1p8_Winter2012/Hydjet1p8_Winter2012/36b8e428366a7e9103136b276c0cd872/hydjet18_GEN_SIM_9_2_rS7.root"
ivars.inputFiles = "file:/afs/cern.ch/work/m/mgardner/private/Z_ee/HYDJET_Background/hydjet18_GEN_SIM_9_2_rS7.root"
ivars.outputFile = 'GenSimMixedZee.root'

#ivars.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.HiEventMixing_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealisticHI2011Collision_cfi')
process.load('SimGeneral.MixingModule.himixGEN_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.Sim_cff')
process.load('SimGeneral.MixingModule.himixSIMExtended_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('SimGeneral.MixingModule.himixDIGI_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_HIon_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.RandomNumberGeneratorService.hiSignal.initialSeed = ivars.randomNumber
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Input source
#process.source = cms.Source("EmptySource")

process.source = cms.Source("PoolSource",
                            #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),                            
                            skipEvents=cms.untracked.uint32(ivars.skipEvents),
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(
        ivars.inputFiles
        ),
                            inputCommands = cms.untracked.vstring('drop *', 
                                                                  'keep *_generator_*_*', 
                                                                  'keep *_g4SimHits_*_*'),
                            dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
                            )


process.options = cms.untracked.PSet( )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.5 $'),
    annotation = cms.untracked.string('Zmumu embedded in Hydjet 1.8 at sqrt(s) = 2.76TeV'),
    name = cms.untracked.string('$Source:  $')
    )

# Output definition

process.RAWSIMHLToutput = cms.OutputModule("PoolOutputModule",
#process.RAWDEBUGoutput = cms.OutputModule("PoolOutputModule",
#process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMHLTEventContent.outputCommands,
    fileName = cms.untracked.string(ivars.outputFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAWDEBUG')
        ),
                                          SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
        )
                                          )

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'STARTHI44_V7::All'
process.GlobalTag.globaltag = 'STARTHI44_V12::All'

process.load('Configuration.GenProduction.Pyquen_zeeNoDY_ptHat020_cfi')

process.ProductionFilterSequence = cms.Sequence(process.hiSignal)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_himix)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
#process.RAWDEBUGoutput_step = cms.EndPath(process.RAWDEBUGoutput)
process.RAWSIMHLToutput_step = cms.EndPath(process.RAWSIMHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
process.schedule.extend([process.endjob_step,process.RAWSIMHLToutput_step])


# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq
    
    from Configuration.GenProduction.customiseCaloDigisNZS import *
    customiseHcalNZS(process)
    





