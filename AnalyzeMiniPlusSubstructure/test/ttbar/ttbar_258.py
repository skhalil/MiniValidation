import FWCore.ParameterSet.Config as cms

###############################################
useMiniAOD = True

# AOD
pfcandidates          = 'particleFlow'
chsstring             = 'pfNoPileUpJME'
genjetparticles       = 'genParticles'
importantgenparticles = 'genParticles'
tracks                = 'generalTracks'
vertices              = 'offlinePrimaryVertices'
mergedvertices        = 'inclusiveMergedVertices' 
mergedvertices2       = '' 
primaryvertices       = 'offlinePrimaryVertices'

#miniAOD
if useMiniAOD:
  pfcandidates          = 'packedPFCandidates'
  genjetparticles       = 'packedGenParticles'
  importantgenparticles = 'prunedGenParticles'
  tracks                = 'unpackedTracksAndVertices'
  vertices              = 'unpackedTracksAndVertices'
  mergedvertices        = 'unpackedTracksAndVertices'
  mergedvertices2       = 'secondary'
  primaryvertices       = 'offlineSlimmedPrimaryVertices'


print 'useMiniAOD = '+str(useMiniAOD)
print ' pfcandidates          = '+pfcandidates         
print ' genjetparticles       = '+genjetparticles      
print ' importantgenparticles = '+importantgenparticles
print ' tracks                = '+tracks               
print ' vertices              = '+vertices             
print ' mergedvertices        = '+mergedvertices       
print ' mergedvertices2       = '+mergedvertices2      
print ' primaryvertices       = '+primaryvertices 


###############################################
# SETUP
process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) , allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


###############################################
# SOURCE
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F4E3856E-6E75-E411-9CFA-0025907DC9B4.root'


  )
)

###############################################
# ANA
process.demo = cms.EDAnalyzer("AnalyzeMiniPlusSubstructure",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger")
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ttbar258.root"),
      closeFileFast = cms.untracked.bool(True)
  )

###############################################
# RECO AND GEN SETUP
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag ='PHYS14_25_V2'
#'START70_V6::All'
#'START70_V6::All'

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
#process.fixedGridRhoFastjetAll.pfCandidatesTag = pfcandidates
process.fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'
process.fixedGridRhoAll.pfCandidatesTag        = 'packedPFCandidates'
# process.fixedGridRhoAll.pfCandidatesTag = .InputTag("packedPFCandidates")
# process.fixedGridRhoFastjetAll = fixedGridRhoFastjetAll.clone( pfCandidatesTag = cms.InputTag("packedPFCandidates"))
# process.fixedGridRhoAll = fixedGridRhoAll.clone( pfCandidatesTag = cms.InputTag("packedPFCandidates"))


from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.caTopTaggers_cff import *




###############################################

process.content = cms.EDAnalyzer("EventContentAnalyzer")


process.p = cms.Path(
  #process.fixedGridRhoFastjetAll
  process.demo
  )
