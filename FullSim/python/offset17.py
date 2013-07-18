import FWCore.ParameterSet.Config as cms

CaloJetPtOffsets = cms.PSet(
    etaMin = cms.double(1.7),
    etaStep = cms.double(0.1),
    #offset_0 = cms.untracked.vdouble(0.406489, 0.391783, 0.414526, 0.399961, 0.333607, 0.279294, 0.233064, 0.199558, 0.160209, 0.123127, 0.103338, 0.081801, 0.064306),
    #offset_500 = cms.untracked.vdouble(4.086449, 3.177738, 1.990882, 1.530543, 1.626077, 1.516009, 1.542726, 1.590620, 1.452601, 1.335337, 1.156529, 1.093183, 0.949487)
)
