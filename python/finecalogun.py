import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

def parse_options():
    from FWCore.ParameterSet.VarParsing import VarParsing
    options = VarParsing("analysis")
    # Shortcuts for options
    register_float = lambda name, default, help="": options.register(
        name, default, VarParsing.multiplicity.singleton, VarParsing.varType.float, help
        )
    register_bool = lambda name, default, help="": options.register(
        name, default, VarParsing.multiplicity.singleton, VarParsing.varType.bool, help
        )
    register_int = lambda name, default, help="": options.register(
        name, default, VarParsing.multiplicity.singleton, VarParsing.varType.int, help
        )
    # Actual options
    register_float('pt', 100., 'Pt in GeV (only for single particle guns)')
    register_float('minpt', -1., 'Pt in GeV (only for single particle guns)')
    register_float('maxpt', -1., 'Pt in GeV (only for single particle guns)')
    register_int('pdgid', -13, 'PdgID of the single particle gun (default -13 muon)')
    register_int('seed', 1001, 'RNG seed')
    register_bool('debug', False, 'Enable debug logging')
    register_bool('dofinecalo', True, 'Enables fine calo for HGCAL')
    register_bool('ntuple', False, 'Runs ntupler instead of saving SIM level')
    register_bool('profiling', False, 'Enables memory and timing profiling')
    options.maxEvents = 1
    options.outputFile = 'default.root'
    options.parseArguments()
    return options


def init_process(options):
    from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
    if options.dofinecalo:
        from Configuration.ProcessModifiers.fineCalo_cff import fineCalo
        process = cms.Process('SIM', Phase2C11, fineCalo)
    else:
        process = cms.Process('SIM', Phase2C11)
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.load('FWCore.MessageService.MessageLogger_cfi')
    # process.load("Configuration.Geometry.GeometryExtended2026D71_cff")
    # process.load('Configuration.Geometry.GeometryExtended2026D71Reco_cff')
    process.load("Configuration.Geometry.GeometryExtended2026D84_cff")
    process.load('Configuration.Geometry.GeometryExtended2026D84Reco_cff')
    process.load("SimGeneral.MixingModule.mixNoPU_cfi")
    process.load('Configuration.StandardSequences.MagneticField_cff')
    process.load('Configuration.StandardSequences.Generator_cff')
    process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
    process.load('GeneratorInterface.Core.genFilterSummary_cff')
    process.load('Configuration.StandardSequences.SimIdeal_cff')
    process.load('Configuration.StandardSequences.Digi_cff')
    process.load('Configuration.StandardSequences.SimL1Emulator_cff')
    process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
    process.load('Configuration.StandardSequences.DigiToRaw_cff')
    process.load('HLTrigger.Configuration.HLT_Fake2_cff')
    process.load('Configuration.StandardSequences.EndOfProcess_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

    # reset all random numbers to ensure statistically distinct but reproducible jobs
    from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
    randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
    randHelper.resetSeeds(options.seed)

    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(options.maxEvents),
        output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
        )

    process.source = cms.Source("EmptySource")

    process.options = cms.untracked.PSet(
        FailPath = cms.untracked.vstring(),
        IgnoreCompletely = cms.untracked.vstring(),
        Rethrow = cms.untracked.vstring(),
        SkipEvent = cms.untracked.vstring(),
        allowUnscheduled = cms.obsolete.untracked.bool,
        canDeleteEarly = cms.untracked.vstring(),
        emptyRunLumiMode = cms.obsolete.untracked.string,
        eventSetup = cms.untracked.PSet(
            forceNumberOfConcurrentIOVs = cms.untracked.PSet(),
            numberOfConcurrentIOVs = cms.untracked.uint32(1)
            ),
        fileMode = cms.untracked.string('FULLMERGE'),
        forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
        makeTriggerResults = cms.obsolete.untracked.bool,
        numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
        numberOfConcurrentRuns = cms.untracked.uint32(1),
        numberOfStreams = cms.untracked.uint32(0),
        numberOfThreads = cms.untracked.uint32(1),
        printDependencies = cms.untracked.bool(False),
        sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
        throwIfIllegalParameter = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False)
        )

    # Production Info
    process.configurationMetadata = cms.untracked.PSet(
        annotation = cms.untracked.string('SingleMuPt1000_pythia8_cfi nevts:10'),
        name = cms.untracked.string('Applications'),
        version = cms.untracked.string('$Revision: 1.19 $')
        )

    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = autoCond['phase2_realistic']


    # _____________________________
    # Main customisations

    if options.minpt == -1. and options.maxpt == -1.:
        minpt = options.pt - 0.01
        maxpt = options.pt + 0.01
    else:
        minpt = options.minpt
        maxpt = options.maxpt
    print 'Using pdgid={}, minpt={}, maxpt={}'.format(options.pdgid, minpt, maxpt)
    add_single_particle_gun(process, options.pdgid, minpt=minpt, maxpt=maxpt)
    if options.debug: add_debug_module(process, 'DoFineCalo')

    # _____________________________
    # Finalizing

    # The standard steps
    process.generation_step = cms.Path(process.pgen)
    process.simulation_step = cms.Path(process.psim)
    process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
    process.endjob_step = cms.EndPath(process.endOfProcess)

    # Make nicely formatted output root file
    if options.outputFile.startswith('default'):
        from time import strftime
        outputFile = '{outputlvl}_event{seed}_pdgid{pdgid}_{pt}GeV_{date}_{finecalo}_n{nevents}.root'.format(
            outputlvl = 'NTUPLE' if options.ntuple else 'SIM',
            seed = options.seed,
            pdgid = abs(options.pdgid),
            pt = int(options.pt) if abs(options.pdgid) != 6 else 14000,
            date = strftime('%b%d'),
            finecalo = 'finecalo' if options.dofinecalo else 'nofine',
            nevents = options.maxEvents
            )
    else:
        outputFile = options.outputFile

    # Figure out the output format
    if options.ntuple:
        process.TFileService = cms.Service("TFileService",
            fileName = cms.string(outputFile)
            )
        process.ntupler = cms.EDAnalyzer('HistoryNTupler')
        process.ntupler_step = cms.Path(process.HistoryNTupler)
        process.schedule = cms.Schedule(
            process.generation_step,
            process.genfiltersummary_step,
            process.simulation_step,
            process.ntupler_step,
            process.endjob_step,
            )
    else:
        process.load('Configuration.EventContent.EventContent_cff')
        process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
            SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('generation_step')
                ),
            dataset = cms.untracked.PSet(
                dataTier = cms.untracked.string('GEN-SIM'),
                filterName = cms.untracked.string('')
                ),
            fileName = cms.untracked.string('file:{}'.format(outputFile)),
            outputCommands = process.FEVTDEBUGEventContent.outputCommands,
            splitLevel = cms.untracked.int32(0)
            )
        process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)
        process.schedule = cms.Schedule(
            process.generation_step,
            process.genfiltersummary_step,
            process.simulation_step,
            process.endjob_step,
            process.FEVTDEBUGoutput_step,
            )

    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)
    # filter all path with the production filter sequence
    for path in process.paths:
        getattr(process,path).insert(0, process.generator)

    if options.profiling:
        # customisation function from Validation.Performance.TimeMemoryInfo
        from Validation.Performance.TimeMemoryInfo import customise 
        process = customise(process)

    if options.dofinecalo:
        for _pset in ["CaloSD", "CaloTrkProcessing", "TrackingAction"]:
            getattr(process.g4SimHits,_pset).UseFineCalo = [2]

    # Add early deletion of temporary data products to reduce peak memory need
    from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
    process = customiseEarlyDelete(process)
    return process


def add_single_particle_gun(
    process, pdgid=-13, minpt=None, maxpt=None
    ):
    pythia_parameters = cms.PSet(parameterSets = cms.vstring())
    if 1 <= abs(pdgid) <= 5:
        pythia_parameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CUEP8M1SettingsBlock,
            parameterSets = cms.vstring('pythia8CommonSettings','pythia8CUEP8M1Settings')
            )
    if abs(pdgid) == 6:
        # process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
        process.generator = cms.EDFilter("Pythia8GeneratorFilter",
            PythiaParameters = cms.PSet(
                parameterSets = cms.vstring(
                    'pythia8CommonSettings',
                    'pythia8CP5Settings',
                    'processParameters'
                ),
                processParameters = cms.vstring(
                    'Top:gg2ttbar = on ',
                    'Top:qqbar2ttbar = on ',
                    '6:m0 = 175 '
                ),
                pythia8CP5Settings = cms.vstring(
                    'Tune:pp 14',
                    'Tune:ee 7',
                    'MultipartonInteractions:ecmPow=0.03344',
                    'MultipartonInteractions:bProfile=2',
                    'MultipartonInteractions:pT0Ref=1.41',
                    'MultipartonInteractions:coreRadius=0.7634',
                    'MultipartonInteractions:coreFraction=0.63',
                    'ColourReconnection:range=5.176',
                    'SigmaTotal:zeroAXB=off',
                    'SpaceShower:alphaSorder=2',
                    'SpaceShower:alphaSvalue=0.118',
                    'SigmaProcess:alphaSvalue=0.118',
                    'SigmaProcess:alphaSorder=2',
                    'MultipartonInteractions:alphaSvalue=0.118',
                    'MultipartonInteractions:alphaSorder=2',
                    'TimeShower:alphaSorder=2',
                    'TimeShower:alphaSvalue=0.118',
                    'SigmaTotal:mode = 0',
                    'SigmaTotal:sigmaEl = 21.89',
                    'SigmaTotal:sigmaTot = 100.309',
                    'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118'
                ),
                pythia8CommonSettings = cms.vstring(
                    'Tune:preferLHAPDF = 2',
                    'Main:timesAllowErrors = 10000',
                    'Check:epTolErr = 0.01',
                    'Beams:setProductionScalesFromLHEF = off',
                    'SLHA:minMassSM = 1000.',
                    'ParticleDecays:limitTau0 = on',
                    'ParticleDecays:tau0Max = 10',
                    'ParticleDecays:allowPhotonRadiation = on'
                )
            ),
            comEnergy = cms.double(14000.0),
            filterEfficiency = cms.untracked.double(1.0),
            maxEventsToPrint = cms.untracked.int32(0),
            pythiaHepMCVerbosity = cms.untracked.bool(False),
            pythiaPylistVerbosity = cms.untracked.int32(0)
        )
    else:
        # Single particle gun (+anti)
        process.generator = cms.EDFilter("Pythia8PtGun",
            PGunParameters = cms.PSet(
                AddAntiParticle = cms.bool(True),
                MinEta = cms.double(1.6),
                MaxEta = cms.double(2.8),
                MinPhi = cms.double(-3.14159265359),
                MaxPhi = cms.double(3.14159265359),
                MinPt = cms.double(minpt),
                MaxPt = cms.double(maxpt),
                ParticleID = cms.vint32(pdgid),
                ),
            PythiaParameters = pythia_parameters,
            Verbosity = cms.untracked.int32(0),
            firstRun = cms.untracked.uint32(1),
            psethack = cms.string('single {} pt {}to{}'.format(pdgid, minpt, maxpt))
            )


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

process = init_process(parse_options())
print 'Final process.schedule:'
print process.schedule
