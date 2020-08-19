import FWCore.ParameterSet.Config as cms
myTrkNtpAnalysis = cms.EDAnalyzer('myTrkNtp')
myTrkNtpAnalysisSequence = cms.Sequence(myTrkNtpAnalysis)
