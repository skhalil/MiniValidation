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
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.FwkJob.limit=1
process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )

###############################################
# SOURCE
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D3EAF1-3174-E411-A5B2-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02EF3EFC-0475-E411-A9DB-002590DB9166.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0434E222-1C75-E411-B4D4-0025907FD34C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04829E8D-6174-E411-97F0-0025901D4940.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04C50B58-7B76-E411-9101-003048D437D2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/060D9464-7174-E411-8EC9-0025907DC9BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06180A42-DA75-E411-B64E-00266CF2AE10.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/063CF1C6-C176-E411-A90A-003048F0E83A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0661643D-4F76-E411-862B-002590AC4E28.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0A2155A4-D376-E411-A1AD-00266CFFA764.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0AE61EFC-3074-E411-854B-002481E101DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0AF8FAF3-3174-E411-A5DE-0025907DCA7E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0C9118BE-0776-E411-A71D-00215AD4D6E2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CA23B27-0275-E411-9638-003048F0E3AC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CB9FD3A-0675-E411-835E-002590AC4FEC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CF2A9A8-3174-E411-A8C5-0025901D4D76.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CF8FFF5-3274-E411-8936-00266CFFA764.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E62BB2E-4174-E411-9C28-002590494C94.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EA5C14E-BC76-E411-8BBA-0025907DC9D0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/102DAE83-0776-E411-B487-002590DB9296.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/104CB0EF-3074-E411-880B-0025901D4B22.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/10A792C6-1076-E411-979E-0025901D4C92.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/10C070BF-0776-E411-8AC8-003048D436EA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1231C32D-3574-E411-BEA5-00266CFFA1AC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1245FBA8-3174-E411-87A0-002590AC4FEC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/14176E42-0575-E411-9E29-0025907DC9C4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/14B0F7DC-B976-E411-A96D-002590DB9296.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/161E8356-0A75-E411-B87B-002590AC5012.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/16FD4024-0475-E411-8A46-002481E0DC4C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1A5B3CE4-3774-E411-855A-0025901D4AF0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1AC8E581-2D75-E411-AE89-003048F0E5A4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1C6B1F46-9F74-E411-8D2B-003048D437EC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1C7C9202-3774-E411-99AB-002481E0DDE8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1EB0B07C-3674-E411-B066-002590AC4BF6.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1EB94895-C576-E411-BD66-0025907FD40C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1EB99435-D575-E411-8D19-003048D436B2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/20010105-3374-E411-B2D6-0025901D481E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2006F5CA-C176-E411-87C7-002590AC4C56.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2067F954-3474-E411-B8C3-002590DB9286.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2243C8A1-3574-E411-BEDD-002590DB041E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2283E236-0876-E411-AA99-002590AC5066.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/22C2A84A-1276-E411-A9A5-002590AC4CB8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/22FFAEE9-3174-E411-AC80-002590AC4C6C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24A8DC0C-3174-E411-821C-002590AC4CEC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24CBA008-3574-E411-8959-003048F0E522.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28A6B754-6F75-E411-AB30-0025901D4138.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A4E17CC-E375-E411-B36E-002590AC50C2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A6C7E9F-1375-E411-8EE8-002590DB927A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A9A30D7-3374-E411-8740-002590AC500E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2C84A2C5-C176-E411-87D8-0025904B578E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2C968AB7-B776-E411-99F6-0025901D42C0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/307294F6-3674-E411-BCF4-003048D4DEAE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/30FD7404-3174-E411-9D77-0025907FD2D8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3283B3DE-1076-E411-BEB7-002590A2CCE6.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/32BF9687-0375-E411-A710-002590AC4C08.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/36716EE5-3174-E411-864E-0025901D4090.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/38DF80DF-A876-E411-B16D-002590AC4CC4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/38F153C4-3174-E411-9F9E-00266CFFA3FC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3ABB400D-0275-E411-91F7-002590AC5068.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C0205B1-6574-E411-984C-0025907DC9C4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C426275-DA75-E411-B1F9-002590DB92A8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C5717C2-0976-E411-BDED-0025907DCA4A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E0F0447-7174-E411-ABF7-002481E101DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E50A413-0D76-E411-BD22-0025901D4C92.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E929E21-0976-E411-9130-002590AC4CEC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3EBD12F2-3174-E411-BF1E-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/405C27A4-0575-E411-B0E8-0025907DCA72.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/408A2472-7575-E411-A071-0025907DC9F8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/40EECAAC-D575-E411-9E48-003048F0E826.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/40F0574A-B776-E411-B423-0025907DCA9C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4243C89C-0575-E411-B658-002481E0D9BC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/424463AC-7C75-E411-B3D6-002590A2CD68.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/42994324-7174-E411-9A8C-00266CFFA044.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/42CC517E-7575-E411-A19B-002590AC539E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/442A5AB8-2D75-E411-8F46-00266CFFA2EC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/44530924-D575-E411-B1E5-002590DB925E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4458A7FD-0C76-E411-80B7-0025901D4C98.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/463B61F0-3174-E411-9DCC-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/482239B8-3574-E411-B7DE-002590DB91F4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/48305281-0C76-E411-B950-003048F30422.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4881873E-BE76-E411-BDE9-0025901D4764.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4A5A8838-0375-E411-AA2E-003048F0E812.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4ACAA13D-BE76-E411-A23D-0025901D4764.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C4171AC-3174-E411-A403-003048D4DF08.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C59912E-1C75-E411-97CD-003048D4DFAA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C8DDF0D-DD75-E411-A5AA-0025901D4B20.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C94DAFB-3474-E411-82DF-002590DB91CE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4CD83CF4-3174-E411-9FB0-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4CF7FF46-6174-E411-9944-002481E94BCA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4EA7942C-BF76-E411-87B2-003048D43982.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4ECB3BF1-3174-E411-91E4-002590DB91F4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4ED11EF4-3174-E411-BD03-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5044277F-2375-E411-9661-002590AC5082.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/52FCA641-B776-E411-A2FB-00266CFFA704.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/541838A6-3274-E411-B7B4-0025901D4C3C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5465E663-0A75-E411-AF2F-002590AC4FC8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/54E93D3F-EE75-E411-83FE-002590DB9358.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/560C5F02-C176-E411-821E-003048D439BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/56612700-1175-E411-A50B-0025904B1428.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/569FFEE7-3274-E411-B872-0025904B578E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/56C3AA36-0876-E411-B7FB-00266CF2AACC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/58A0F116-3174-E411-BAFD-0025907DCA4A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/58FB1F2A-ED75-E411-B019-003048F0E832.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A12CA60-0C76-E411-8A70-002590AC4C08.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A7DFA05-3774-E411-8EBC-00266CF32A00.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A8B1611-DD75-E411-9AE5-00266CFFA0B0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5ACE92DB-0C76-E411-9E62-003048F0E82C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5C62E0C4-1375-E411-9F4A-00266CFFA25C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5E297F32-BF76-E411-9625-002590DB91BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5ED65C61-0B76-E411-9D3D-00259081A2C8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6009D345-7B76-E411-9950-002590DB91FA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62517940-0A75-E411-AB83-002590AC5482.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62B34D39-0375-E411-87A5-0025907FD24A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64196AE5-3774-E411-B868-00266CF9B420.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/641A7B60-BC76-E411-8221-0025904B578E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64E57633-6F75-E411-AE68-00266CFFA604.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66203DDB-3775-E411-9E9E-002590AC4C9E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66689298-3874-E411-BBFB-0025904B1340.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66A6025E-D275-E411-B485-002590DB9188.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68548FD6-7C75-E411-8EC4-002590AC4B3E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68B439F0-3174-E411-A50C-0025904B144E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6A99352C-9F74-E411-8244-002590A2CD44.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6C725C73-0A75-E411-91E3-003048D4364C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6E47D20A-0276-E411-9773-00215AD4D6E2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/70221D45-C276-E411-B9CF-002590DB9258.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7284FA91-3174-E411-AB57-00266CF33118.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/72C0782B-D275-E411-ABBE-0025901D4090.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/72C087A0-C076-E411-B02A-002590DB92C4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/747FCFF4-C276-E411-A4BF-0025904B5FAE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/76E53E7C-6F75-E411-86D7-00266CFFA704.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/788520A5-3274-E411-B300-002590DB91C8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7A0034B0-3274-E411-A5A3-002590A52B4A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7ACD4BB1-7D75-E411-A921-002590AC5074.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7AE2C6CA-2F74-E411-8237-00266CF422D8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7C44A5C7-3174-E411-9923-00266CF25320.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7C8801E6-3174-E411-BF9E-0025901D4AF0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7CCE517B-0575-E411-877E-002590DB9214.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7CF49B55-6F75-E411-811E-003048D439A0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7E803F15-1075-E411-A96C-00266CF32EAC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/80003E1C-0675-E411-BD36-003048D4399E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8097C8AB-A876-E411-90C4-003048F0E424.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/82209B13-EE75-E411-8C38-003048D439C0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/822245B8-BF76-E411-8B53-002590AC5074.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8431D528-1475-E411-9177-0025901D4D76.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8454E20A-0275-E411-A450-002590DB917C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/84B4F842-DA75-E411-8148-002481E0D3C0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/86A67A1B-0976-E411-9673-003048F0E82C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/86DD8F22-7174-E411-B0B3-003048F0E820.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/88A4F21E-0275-E411-A264-003048F0E832.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8A731319-ED75-E411-8D97-002590AC4C7C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8A790641-3274-E411-8F0D-0025904B5FAE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C399FB5-1076-E411-BC85-002590A2CDA8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C65D0B1-3174-E411-971E-002590A52B4A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C6B9438-0375-E411-A391-003048F0E1AE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8CE103E3-0975-E411-B9B2-003048D4399E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8E10ACC4-3174-E411-9B4C-002590AC4B54.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8E27CA78-D875-E411-B20C-0025907DC9D2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/906D84DE-1076-E411-BF61-002590A2CCE6.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/909D0488-0375-E411-8A6B-002590AC4CC4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/90D3A80B-4F76-E411-9889-002590AC4B5A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/946DCAC3-4C74-E411-8753-002590DB91CA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/94EC2B61-D275-E411-B825-0025901D4D76.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/987048F1-BF76-E411-8C8B-00259074B2BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/98807FEB-3574-E411-9D05-002590AB38DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/988386C1-1D75-E411-AEC3-0025904B1452.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/98D86D1D-DA75-E411-95A4-002590DB041A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9ADAD08B-3275-E411-A7CB-003048D437DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9AF7EE5E-3474-E411-B3C6-0025901D4AF0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C459463-3574-E411-ADF4-002590DB9296.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C771AFB-3574-E411-B1FE-00266CFFA148.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C84DAD8-A374-E411-A16D-0025907DC9B0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C9B85EE-3074-E411-9D16-002590A7A2E0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E0E1BA7-DA75-E411-B6C4-003048F0E84E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E2F6B01-3875-E411-8186-002590AC5066.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A20D12CF-C176-E411-B194-002590A2CD44.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A22845E0-0975-E411-858E-0025901D4B4C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A24A3F93-3174-E411-80A7-003048F0E3BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A2B871BA-1D75-E411-B0BE-0025901D4C46.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A40F4C6F-DF75-E411-8FE6-00266CF32BC4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A4623817-0276-E411-A9EA-00266CFFA7A8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A47DEDA3-6574-E411-BE4B-002481E0DC66.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A4F62A81-0776-E411-955B-002590DB9188.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A888DEF0-C576-E411-9FEC-002590494CB2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AA2325BC-B776-E411-870A-003048D4399C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AC685927-D275-E411-BBB1-0025907DC9C8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/ACB4EF86-C176-E411-A346-00266CFF9FFC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE393589-0776-E411-977B-002590DB9258.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE3E7AAE-1076-E411-BECD-002590AC4BF6.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE416221-DB75-E411-94D6-003048D43996.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE839598-C076-E411-9931-002590DB0640.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AECF34F2-3274-E411-AE88-002481E94C7A.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AEF96180-0A76-E411-ABED-002590DB9286.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B2297A5D-6E75-E411-8F08-00266CF330D8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B26D69D8-C376-E411-876B-0025907FD2DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B4070E36-0876-E411-8BCD-002590DB92BE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B49DF412-1475-E411-91CF-002590DB9196.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B67338C9-4C74-E411-AFED-002590AC4FEA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B6C3101B-0D76-E411-B21C-00266CF1074C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B81D2D7D-0A76-E411-B4E3-0025901D4ADE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B8C3665D-0A75-E411-B0A2-002590AC4C6C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BA1AF2B2-BF76-E411-B27C-003048D437EC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BA53B3A8-3174-E411-A862-0025907DCA0C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BABEFAFE-0475-E411-84F6-003048F02CB4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BAF6AC29-4174-E411-9AC5-0025904B1372.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BE2D6273-C376-E411-9F78-003048F0E194.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BE494C81-0A76-E411-B7E6-002590A2CD44.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C0F681DC-E375-E411-8AB8-0025901D4B04.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C0FC2A5D-6E75-E411-AE87-003048F0E82C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C20DB61D-0675-E411-8B59-0025901D4B4C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C2DB8CE4-3174-E411-8F5C-00266CF422D8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C2FD308C-C176-E411-B4CB-002590DB9296.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4A4BECC-0976-E411-8D3D-002590AC4CAA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4EC08D5-2F74-E411-A823-00266CF32F90.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C8A17219-0475-E411-8FAD-002590DB91B4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C8FC1057-3274-E411-9F12-00266CF33054.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CA49FC20-C276-E411-85F5-002590DB9302.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CAE5FA96-C076-E411-A7AC-002590DB91C6.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CC99B1C0-0776-E411-AA16-002590AB3A70.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CEBD11F4-3174-E411-8D59-00266CFFA768.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CEBF9165-3474-E411-9A79-0025907FD34C.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CEC3829C-BF76-E411-A6B8-003048F0E194.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CEFB3A83-DF75-E411-B8CF-00266CFFA6A8.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D0126D62-D875-E411-A5CF-002590DB91FA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D099D827-0B77-E411-A87C-002590494CC2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D20124EA-3174-E411-9887-003048F0E5B4.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D2281AE0-7A77-E411-8FAB-0025901D42C0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D23C3675-0776-E411-A78C-002590AC4C48.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D27B4D96-C176-E411-BA4B-00266CFFA2EC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D2A395C1-C176-E411-BD3F-002590DB05DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D2EE68AE-D575-E411-929A-002481E0D3C0.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D4205ED7-BE76-E411-863C-0025907FD442.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D4C7FD87-2375-E411-8982-00266CFFA090.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DA4DDAFF-3474-E411-8428-002590AC4C9E.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DA885D3B-0375-E411-A203-0025901D4124.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DAA3270D-1075-E411-8535-002590A2CD44.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC1035C0-3174-E411-83E8-003048D43788.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC44095D-3274-E411-831E-003048D3C886.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DE49746D-0475-E411-8B81-003048F0E524.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DE8CA3D0-B976-E411-8389-002590DB9296.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E055A2FF-3174-E411-9E30-00266CFFA3FC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E072FA5C-0475-E411-A213-002481E76372.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E07437F5-3774-E411-A089-002590AC4B76.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2B1B3FF-3074-E411-91A7-002590DB05DA.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2CBC38C-3275-E411-A3C3-0025904B5FAE.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E4E1455D-0475-E411-AF60-00266CFFA750.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E66FB37E-3674-E411-AD05-0025904B1366.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E6EDB1E4-A374-E411-BB3C-003048F0E812.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E8244374-DA75-E411-87C2-00266CFFA7BC.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E8585D1C-DA75-E411-BAFB-002481E0D5E2.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E85C60F8-3474-E411-9610-0025907DCA72.root',
'root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E868CCCC-3275-E411-8317-0025901D4D20.root'
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
      fileName = cms.string("test.root"),
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
