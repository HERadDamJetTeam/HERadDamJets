import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGBLTG::All', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# this inputs the input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part1.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part2.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part3.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part4.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part5.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part6.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part7.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part8.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part9.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part10.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part11.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part12.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part13.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part14.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part15.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part16.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part17.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part18.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part19.root',
    'file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_part20.root'
    ),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.demo = cms.EDAnalyzer('FullSimTowerNoiseAnalyzer',
    fileName = cms.string("tree_towernoise_20YEAR_ENERGYIN_lumiLUMIDRK.root"),
	
	dRcut = cms.double(0.5)
	
)

process.p = cms.Path(process.demo)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#SOURCE process.load('#SRC')
