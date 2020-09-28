import FWCore.ParameterSet.Config as cms

myOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                               src = cms.InputTag('offlinePrimaryVertices'),
                                               cut = cms.string("isValid & ndof > 4 & tracksSize > 0 & abs(z) <= 24 & abs(position.Rho) <= 2."),
                                               filter = cms.bool(False)
)

myOfflinePrimaryVerticesWithBS = myOfflinePrimaryVertices.clone()
myOfflinePrimaryVerticesWithBS.src = cms.InputTag('offlinePrimaryVerticesWithBS')

myVtxAnalysis = cms.EDAnalyzer('MyVtx',
                                use_only_charged_tracks = cms.untracked.bool(True),
                                verbose = cms.untracked.bool(False),
                                trackingParticleCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                                trackingVertexCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                                trackAssociatorMap = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
                                vertexAssociator = cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
                                vertexRecoCollection = cms.untracked.InputTag("myOfflinePrimaryVertices"),
				tracks = cms.untracked.InputTag("generalTracks")
)
from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(myVtxAnalysis,
    trackingParticleCollection = "mixData:MergedTrackTruth",
    trackingVertexCollection = "mixData:MergedTrackTruth",
)

myVtxAnalysisTrackingOnly = myVtxAnalysis.clone()

##########

myVtxAnalysisSelection = cms.Sequence(
    cms.ignore(myOfflinePrimaryVertices)
    + cms.ignore(myOfflinePrimaryVerticesWithBS)
)

##########

myVtxAnalysisSequence = cms.Sequence(
    myVtxAnalysisSelection
    + myVtxAnalysis
)

myVtxAnalysisSequenceTrackingOnly = cms.Sequence(
    myVtxAnalysisSelection
    + myVtxAnalysisTrackingOnly
)
