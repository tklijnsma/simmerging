import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing("analysis")
options.parseArguments()

def add_debug_module(process, module_name='DoFineCalo'):
    process.MessageLogger.cerr.threshold = "DEBUG"
    process.MessageLogger.cerr.FwkReport.limit = 0
    process.MessageLogger.cerr.FwkSummary.limit = 0
    process.MessageLogger.cerr.default.limit = 0
    process.MessageLogger.debugModules.append(module_name)
    setattr(
        process.MessageLogger.cerr, module_name,
        cms.untracked.PSet(limit = cms.untracked.int32(10000000))
        )

from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
process = cms.Process('simmerger', Phase2C11)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.Geometry.GeometryExtended2026D71_cff")
process.load('Configuration.Geometry.GeometryExtended2026D71Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles))
output_file = options.inputFiles[0].replace('SIM', 'SIMMERGED')

add_debug_module(process, 'SimMerging')

# process.simulation_step = cms.Path(process.psim)

process.simmerger = cms.EDProducer("simmerger")
process.simmerger_step = cms.Path(process.simmerger)
process.end_step = cms.EndPath(process.endOfProcess)

process.load('Configuration.EventContent.EventContent_cff')
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    # SelectEvents = cms.untracked.PSet(
    #     SelectEvents = cms.vstring('simulation_step')
    #     ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
        ),
    fileName = cms.untracked.string('file:{}'.format(output_file)),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
    )
process.FEVTDEBUGoutput.outputCommands.extend([
    'drop *',
    'keep *_simmerger_*_*',
    'keep SimTracks_*_*_*',
    'keep SimVertexs_*_*_*',
    'keep *_genParticles_*_*',
    ])
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(
    # process.simulation_step,
    process.simmerger_step,
    process.end_step,
    process.FEVTDEBUGoutput_step,
    )
