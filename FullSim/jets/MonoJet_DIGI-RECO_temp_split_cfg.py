# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: FastSimulation/Configuration/python/DiJets_cfi.py --python_filename MonoJet_2017_DIGI-RECO_test_cfg.py -s DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --conditions STAR17_61_V1A::All --geometry Extended2017 --datatier GEN-SIM-RECO -n 10 --eventcontent FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2017,SLHCUpgradeSimulations/Configuration/combinedCustoms.noCrossing --filein file:monojet_gensim_2017.root --fileout file:monojet_digireco_2017.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended20YEARReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

#random number seed
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(1PNUMBER)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:jet_gensim_20YEAR_ENERGYIN_partPNUMBER.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('FastSimulation/Configuration/python/DiJets_cfi.py nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = cms.untracked.vstring('drop *',
									   'keep PCaloHits_g4SimHits_EcalHitsEB_*',
									   'keep PCaloHits_g4SimHits_EcalHitsEE_*',
									   'keep PCaloHits_g4SimHits_HcalHits_*',
									   'keep EcalRecHitsSorted_ecalRecHit_*_*',
									   'keep HBHERecHitsSorted_*_*_*',
									   'keep HFRecHitsSorted_*_*_*',
									   'keep HORecHitsSorted_horeco_*_*',
									   'keep CaloTowersSorted_towerMaker_*_*',
									   'keep recoCaloJets_ak5CaloJets_*_*',
									   'keep recoGenJets_ak5GenJets_*_*',
									   'keep recoPFJets_ak5PFJets_*_*'
									   ),	
    fileName = cms.untracked.string('file:jet_digireco_20YEAR_ENERGYIN_lumiLUMIDRK_partPNUMBER.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'STAR17_61_V1A::All', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGBLTG::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_20YEAR,noCrossing,ageHcal,agePixel

#call to customisation function cust_20YEAR imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_20YEAR(process)

#call to customisation function noCrossing imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = noCrossing(process)

#aging customizations
process = ageHcal(process,LUMIDRK)
process = agePixel(process,LUMIDRK)

#ECAL customizations
if not hasattr(process.GlobalTag,'toGet'):
    process.GlobalTag.toGet=cms.VPSet()

#globaltag conditions
if LUMIDRK==0:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),
                 tag = cms.string("EcalTPGLinearizationConst_beamv5_startup_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_2011_V3_Bon_start_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_v02_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_p1p2p3_v2_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_beam2012_option1_v00_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_startup_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
    )
else:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),
                 tag = cms.string("EcalTPGLinearizationConst_TLLUMIDRK_ILINSTLUMI_v1_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_TLLUMIDRK_ILINSTLUMI_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_TLLUMIDRK_ILINSTLUMI_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_TLLUMIDRK_ILINSTLUMI_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_42X_ECAL_LAS")
                 ),
    	cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TLLUMIDRK_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 ),
    	cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
    )

# End of customisation functions
