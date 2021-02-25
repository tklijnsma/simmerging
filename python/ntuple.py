import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing("analysis")
options.parseArguments()
from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
process = cms.Process('ntupler', Phase2C11)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.Geometry.GeometryExtended2026D71_cff")
process.load('Configuration.Geometry.GeometryExtended2026D71Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles))
output_file = options.inputFiles[0].replace('SIM', 'NTUPLE')
if output_file == options.inputFiles[0]:
    raise Exception('About to overwrite input!')
process.TFileService = cms.Service("TFileService", fileName=cms.string(output_file))
process.ntupler = cms.EDAnalyzer("hgcalfinecalontupler")
process.step = cms.Path(process.ntupler)
process.end_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.step, process.end_step)