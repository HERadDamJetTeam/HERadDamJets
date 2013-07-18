import FWCore.ParameterSet.Config as cms

CaloJetPtOffsets = cms.PSet(
    etaMin = cms.double(1.7),
    etaStep = cms.double(0.1),
    #offset_0 = cms.untracked.vdouble(0.228516, 0.235295, 0.237238, 0.229460, 0.201133, 0.169657, 0.150337, 0.122920, 0.239256, 0.377470, 0.447478, 0.468063, 0.570233),
    #offset_500 = cms.untracked.vdouble(3.666318, 2.857824, 1.652837, 1.188886, 1.266323, 1.193647, 1.266500, 1.277904, 1.343364, 1.397737, 1.404945, 1.344286, 1.277148)
)
