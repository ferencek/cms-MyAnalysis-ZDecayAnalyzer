import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runHadronic', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run hadronic Z decays"
)

options.parseArguments()

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/f/ferencek/MC_production/CMSSW_5_3_7_patch4/src/Z1JetToEE_Hadronizer_TuneZ2star_8TeV_generic_LHE_pythia_cff_py_GEN.root'
    )
)

if options.runHadronic:
    process.source.fileNames = [
        'file:/afs/cern.ch/work/f/ferencek/MC_production/CMSSW_5_3_7_patch4/src/Z1JetToQQ_Hadronizer_TuneZ2star_8TeV_generic_LHE_pythia_cff_py_GEN.root'
    ]

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('histos_' + ('hadronic' if options.runHadronic else 'leptonic') + '.root')
)

process.leptonicZDecays = cms.EDAnalyzer('ZDecayAnalyzer',
    GenParticleTag            = cms.InputTag('genParticles'),
    BosonPdgId                = cms.int32(23),
    BosonDecayProdPdgIds      = cms.vint32(11)
)

process.hadronicZDecays = cms.EDAnalyzer('ZDecayAnalyzer',
    GenParticleTag            = cms.InputTag('genParticles'),
    BosonPdgId                = cms.int32(23),
    BosonDecayProdPdgIds      = cms.vint32(1,2,3,4,5)
)

process.hadronicZDecaysNoB = cms.EDAnalyzer('ZDecayAnalyzer',
    GenParticleTag            = cms.InputTag('genParticles'),
    BosonPdgId                = cms.int32(23),
    BosonDecayProdPdgIds      = cms.vint32(1,2,3,4)
)

process.p = cms.Path(process.leptonicZDecays)

if options.runHadronic:
    process.p = cms.Path(process.hadronicZDecays + process.hadronicZDecaysNoB)
