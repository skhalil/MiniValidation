// system include files
#include <memory>
#include <vector>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"


// TFile
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// root

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TUnixSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include "TH1.h"
#include "TH2.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TFormula.h>
#include <TMath.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <TLatex.h>
#include <fstream>
#include <TText.h>
#include "TGraph.h"
#include <TChain.h>
#include <TChainElement.h>
#include <algorithm>




//
// class declaration
//


// global stuff

const int NUM_CUTFLOW = 12;
std::string  cut_flow_case[NUM_CUTFLOW] = {"_LepPtCase", "_LepEtaCase", "_METCase", "_MTCase", "_bJetPtCase", "_LepbJetDPhiCase","_LepJetMinDRCase", "_ForwardJetCase" , "_ST500Cut", "_ST600Cut", "_ST700Cut","_ST800Cut"  };


class AnalyzeMiniPlusSubstructure : public edm::EDAnalyzer {
public:
    explicit AnalyzeMiniPlusSubstructure(const edm::ParameterSet&);
    ~AnalyzeMiniPlusSubstructure();
    bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
    static bool compareByPt(const TLorentzVector a, const TLorentzVector b ){
        return a.Pt() > b.Pt();
    }
   // void make1DHist(const char* hname, const char* htitle, int nbins, double xlow, double xhigh, edm::Service<TFileService> &fs  );
   // void create1DHist(const char* hname);
private:

    std::map< std::string, TH1D* > hists;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    
    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    
    // ----------root ---------------------------
    
    TH1D * NVTX_GOOD            ;
    TH1D * NVTX_ALL             ;
    TH1D * RHO                  ;
    
    
    TH1D * MU_PT                ;
    TH1D * MU_ETA               ;
    TH1D * MU_PHI               ;
    TH1D * MU_dzPV              ;
    TH1D * MU_dxyPV             ;
    TH1D * MU_ISO               ;
    TH1D * MU_nStation          ;
    
    TH1D * EL_PT                ;
    TH1D * EL_ETA               ;
    TH1D * EL_PHI               ;
    TH1D * EL_ISO               ;
    TH1D * EL_SIGMAETAETA       ;
    TH1D * EL_HoverE            ;
    TH1D * EL_DO                ;
    TH1D * EL_DZ                ;
    
    
    
    TH1D * AK4PF_PT             ;
    TH1D * AK4PF_PHI            ;
    TH1D * AK4PF_ETA            ;
    TH1D * AK4PF_RAP            ;
    TH1D * AK4PF_MASS           ;
    TH1D * AK4PF_NCONST         ;
    TH1D * AK4PF_AREA           ;
    TH1D * AK4PF_CH_MULT        ;
    TH1D * AK4PF_NE_MULT        ;
    TH1D * AK4PF_CHEF           ;
    TH1D * AK4PF_NHEF           ;
    TH1D * AK4PF_CEEF           ;
    TH1D * AK4PF_NEEF           ;
    TH1D * AK4PF_CMEF           ;
    TH1D * AK4PF_CSV             ;
    TH1D * AK4PF_PUID            ;

    
    TH1D * AK8PF_PT             ;
    TH1D * AK8PF_PHI            ;
    TH1D * AK8PF_ETA            ;
    TH1D * AK8PF_RAP            ;
    TH1D * AK8PF_MASS           ;
    TH1D * AK8PF_MASSPRUNED    ;
    TH1D * AK8PF_MASSTRIMMED   ;
    TH1D * AK8PF_MASSFILTERED   ;
    TH1D * AK8PF_MASSTOP       ;
    TH1D * AK8PF_NCONST       ;
    TH1D * AK8PF_AREA         ;
    TH1D * AK8PF_CH_MULT      ;
    TH1D * AK8PF_NE_MULT      ;
    TH1D * AK8PF_CHEF         ;
    TH1D * AK8PF_NHEF         ;
    TH1D * AK8PF_CEEF         ;
    TH1D * AK8PF_NEEF         ;
    TH1D * AK8PF_CMEF         ;
    
    TH1D * MET_PT             ;
    TH1D * MET_PHI            ;
    TH1D * MET_SUMET          ;
    TH1D * MET_GENMET         ;
    

    TH1D * AK4PF_MULTIPLICITY       ;
    TH1D * AK4PF_LEADINGPT          ;
    TH1D * AK4PF_LEADINGETA         ;
    TH1D * AK4PF_SECLEADINGPT       ;
    TH1D * AK4PF_SECLEADINGETA      ;
    TH1D * AK4PF_THIRDLEADINGPT     ;
    TH1D * AK4PF_THIRDLEADINGETA    ;
    TH1D * AK4PF_FOURTHLEADINGPT    ;
    TH1D * AK4PF_FOURTHLEADINGETA   ;
    
    TH1D * AK4PF_CSV_LMULTIPLICITY  ;
    TH1D * AK4PF_CSV_MMULTIPLICITY  ;
    
    
    TH1D * AK4PF_MCSV_LEADINGPT     ;
    TH1D * AK4PF_MCSV_LEADINGETA    ;
    TH1D * AK4PF_MCSV_SECLEADINGPT  ;
    TH1D * AK4PF_MCSV_SECLEADINGETA ;
    

    TH1D * AK8PF_tau2Overtau1;
    TH1D * MT;
    
    
    TH1D * AK4PF_MULTIPLICITY_SLF       ;
    TH1D * AK4PF_LEADINGPT_SLF          ;
    TH1D * AK4PF_LEADINGETA_SLF         ;
    TH1D * AK4PF_SECLEADINGPT_SLF       ;
    TH1D * AK4PF_SECLEADINGETA_SLF      ;
    TH1D * AK4PF_THIRDLEADINGPT_SLF     ;
    TH1D * AK4PF_THIRDLEADINGETA_SLF    ;
    TH1D * AK4PF_FOURTHLEADINGPT_SLF    ;
    TH1D * AK4PF_FOURTHLEADINGETA_SLF   ;

    
    
    TH1D * AK4PF_CSV_LMULTIPLICITY_SLF  ;
    TH1D * AK4PF_CSV_MMULTIPLICITY_SLF  ;
    
    
    TH1D * AK4PF_MCSV_LEADINGPT_SLF     ;
    TH1D * AK4PF_MCSV_LEADINGETA_SLF    ;
    TH1D * AK4PF_MCSV_SECLEADINGPT_SLF  ;
    TH1D * AK4PF_MCSV_SECLEADINGETA_SLF ;
    TH1D * MT_SLF;
    
    TH1D * JetClesestToMuPt;
    TH1D * ClosestJetToMuMinDR ;
    
    TH1D * JetClesestToElePt ;
    TH1D * ClosestJetToEleMinDR;
    
    TH1D * ST_SLF ;
    TH1D * HT_SLF;
    
    TH1D * LeadingBJetEleDPhi;
    TH1D * LeadingBJetMuDPhi;
    
    TH1D * pfMET_SLF;
    TH1D * TwoLeadingJetsDR;
    TH1D * JetsDR;
    
    TH1D * EL_PT_SLF                ;
    TH1D * MU_PT_SLF                ;
    TH1D * TriMass_SLF              ;
    TH1D * ForwardJetPt             ;
    TH1D * ForwardJetEta            ;
    
    TH1D * ForwardGenJetEta         ;
    TH1D * ForwardGenJetPt          ;
    
    TH1D * Forward_Gen_Reco_DR      ;
    TH1D * AK8MatchedToB_MinDR;
    
    TH1D * MET_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * MT_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LepPt_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LepEta_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * BJetPt_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LepJetMinDR_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LepBJetDPhi_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * DiJetMass_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * TriMass_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * ST_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * ForwardJetEta_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * ForwardJetPt_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK4PF_LEADINGPT_SLF_CUTFLOW[NUM_CUTFLOW]          ;
    TH1D * AK4PF_SECLEADINGPT_SLF_CUTFLOW[NUM_CUTFLOW]       ;
    TH1D * AK8LeadingBDR_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK8PrunedMass_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK8PF_WTagged_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK4HT_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK4Multiplicity_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK4HTEta2p5_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LooseBJetPt_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * LooseBJetsDR_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * DRLeadingbClosestjet_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * PtClosestjetToB_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK8MatchedToBTau1_SLF_CUTFLOW[NUM_CUTFLOW];
    TH1D * AK8MatchedToBTau2_SLF_CUTFLOW[NUM_CUTFLOW];
    
};

AnalyzeMiniPlusSubstructure::AnalyzeMiniPlusSubstructure(const edm::ParameterSet& iConfig):
vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
{
    

   
    
    edm::Service<TFileService> fs;
    
    
    
    MT = fs->make<TH1D>("MT"      ,"MT", 100,0,250  );
    
    NVTX_GOOD             = fs->make<TH1D>("NVTX_GOOD"      ,"", 80,0,80     );
    NVTX_ALL              = fs->make<TH1D>("NVTX_ALL"       ,"", 80,0,80     );
    RHO                   = fs->make<TH1D>("RHO"            ,"", 200,0,100   );
    
    MU_PT                 = fs->make<TH1D>("MU_PT"        ,"MU_PT",         100,0,500      );
    MU_ETA                = fs->make<TH1D>("MU_ETA"       ,"MU_ETA",        100,-3,3       );
    MU_PHI                = fs->make<TH1D>("MU_PHI"       ,"MU_PHI",        100,-3.2,3.2   );
    MU_dzPV               = fs->make<TH1D>("MU_dzPV"      ,"MU_dzPV",       100,0,1        );
    MU_dxyPV              = fs->make<TH1D>("MU_dxyPV"     ,"MU_dxyPV",      100,0,1        );
    MU_ISO                = fs->make<TH1D>("MU_ISO"       ,"MU_ISO",        100,0,1       );
    MU_nStation           = fs->make<TH1D>("MU_nStation"  ,"MU_nStation",   5, -0.5,4.5    );
    
    EL_PT                 = fs->make<TH1D>("EL_PT"                ,"EL_PT",         100,0,500       );
    EL_ETA                = fs->make<TH1D>("EL_ETA"               ,"EL_ETA",        100,-3,3        );
    EL_PHI                = fs->make<TH1D>("EL_PHI"               ,"EL_PHI",        100,-3.2,3.2    );
    EL_ISO                = fs->make<TH1D>("EL_ISO"               ,"EL_ISO",        100,0,1         );
    EL_HoverE             = fs->make<TH1D>("EL_HoverE"            ,"EL_HoverE",     100,0,0.25      );
    EL_SIGMAETAETA        = fs->make<TH1D>("EL_SIGMAETAETA"       ,"EL_SIGMAETAETA",100,0,0.04      );
    EL_DO                 = fs->make<TH1D>("EL_DO"                ,"EL_DO"         ,200,-1,1      );
    EL_DZ                 = fs->make<TH1D>("EL_DZ"                ,"EL_DZ"         ,200,-1,1      );
    
    
    AK4PF_MULTIPLICITY      = fs->make<TH1D>("AK4PF_MULTIPLICITY"              ,"AK4PF_MULTIPLICITY",   10,-0.5,9.5     );
    AK4PF_LEADINGPT         = fs->make<TH1D>("AK4PF_LEADINGPT"                 ,"AK4PF_LEADINGPT",      100,0,1000      );
    AK4PF_LEADINGETA        = fs->make<TH1D>("AK4PF_LEADINGETA"                ,"AK4PF_LEADINGETA",     100,-5,5        );
    AK4PF_SECLEADINGPT      = fs->make<TH1D>("AK4PF_SECLEADINGPT"              ,"AK4PF_SECLEADINGPT",   100,0,1000      );
    AK4PF_SECLEADINGETA     = fs->make<TH1D>("AK4PF_SECLEADINGETA"             ,"AK4PF_SECLEADINGETA",  100,-5,5        );
    AK4PF_THIRDLEADINGPT    = fs->make<TH1D>("AK4PF_THIRDLEADINGPT"            ,"AK4PF_THIRDLEADINGPT", 100,0,1000      );
    AK4PF_THIRDLEADINGETA   = fs->make<TH1D>("AK4PF_THIRDLEADINGETA"           ,"AK4PF_THIRDLEADINGETA",100,-5,5        );
    AK4PF_FOURTHLEADINGPT   = fs->make<TH1D>("AK4PF_FOURTHLEADINGPT"           ,"AK4PF_FOURTHLEADINGPT",100,0,1000      );
    AK4PF_FOURTHLEADINGETA   = fs->make<TH1D>("AK4PF_FOURTHLEADINGETA"           ,"AK4PF_FOURTHLEADINGETA",100,-5,5     );
    
    AK4PF_CSV_LMULTIPLICITY = fs->make<TH1D>("AK4PF_CSV_LMULTIPLICITY"               ,"AK4PF_CSV_LMULTIPLICITY", 5,-0.5,4.5   );
    AK4PF_CSV_MMULTIPLICITY = fs->make<TH1D>("AK4PF_CSV_MMULTIPLICITY"               ,"AK4PF_CSV_MMULTIPLICITY", 5,-0.5,4.5   );
    
    
    AK4PF_MCSV_LEADINGPT = fs->make<TH1D>("AK4PF_MCSV_LEADINGPT"               ,"AK4PF_MCSV_LEADINGPT", 100,0,1000   );
    AK4PF_MCSV_LEADINGETA = fs->make<TH1D>("AK4PF_MCSV_LEADINGETA"               ,"AK4PF_MCSV_LEADINGETA", 100,-3,3  );
    AK4PF_MCSV_SECLEADINGPT = fs->make<TH1D>("AK4PF_MCSV_SECLEADINGPT"               ,"AK4PF_MCSV_SECLEADINGPT", 100,0,1000   );
    AK4PF_MCSV_SECLEADINGETA = fs->make<TH1D>("AK4PF_MCSV_SECLEADINGETA"               ,"AK4PF_MCSV_SECLEADINGETA",100,-3,3    );
    
    
    AK4PF_PT              = fs->make<TH1D>("AK4PF_PT"             ,"", 100,0,250   );
    AK4PF_PHI             = fs->make<TH1D>("AK4PF_PHI"            ,"",100,-3.2,3.2 );
    AK4PF_ETA             = fs->make<TH1D>("AK4PF_ETA"            ,"", 200,-5,5  );
    AK4PF_RAP             = fs->make<TH1D>("AK4PF_RAP"            ,"",100,-4,4     );
    AK4PF_MASS            = fs->make<TH1D>("AK4PF_MASS"           ,"",100,0,500   );
    AK4PF_NCONST          = fs->make<TH1D>("AK4PF_NCONST"         ,"",1000,0,1000  );
    AK4PF_AREA            = fs->make<TH1D>("AK4PF_AREA"           ,"",100,0,1     );
    AK4PF_CH_MULT         = fs->make<TH1D>("AK4PF_CH_MULT"        ,"",500,0,500   );
    AK4PF_NE_MULT         = fs->make<TH1D>("AK4PF_NE_MULT"        ,"",500,0,500   );
    AK4PF_CHEF            = fs->make<TH1D>("AK4PF_CHEF"           ,"",200,0,1     );
    AK4PF_NHEF            = fs->make<TH1D>("AK4PF_NHEF"           ,"",200,0,1     );
    AK4PF_CEEF            = fs->make<TH1D>("AK4PF_CEEF"           ,"",200,0,1     );
    AK4PF_NEEF            = fs->make<TH1D>("AK4PF_NEEF"           ,"",200,0,1     );
    AK4PF_CMEF            = fs->make<TH1D>("AK4PF_CMEF"           ,"",200,0,1     );
    AK4PF_CSV             = fs->make<TH1D>("AK4PF_CSV"            ,"",100,0,1     );
    AK4PF_PUID            = fs->make<TH1D>("AK4PF_PUID"           ,"",100,-1,1     );

    AK8PF_PT              = fs->make<TH1D>("AK8PF_PT"             ,"",100,0,1000   );
    AK8PF_PHI             = fs->make<TH1D>("AK8PF_PHI"            ,"",100,-3.2,3.2 );
    AK8PF_ETA             = fs->make<TH1D>("AK8PF_ETA"            ,"",200,-5,5  );
    AK8PF_RAP             = fs->make<TH1D>("AK8PF_RAP"            ,"",100,-4,4     );
    AK8PF_MASS            = fs->make<TH1D>("AK8PF_MASS"           ,"",100,0,500   );
    AK8PF_MASSPRUNED      = fs->make<TH1D>("AK8PF_MASSPRUNED"     ,"",100,0,500   );
    AK8PF_MASSTRIMMED     = fs->make<TH1D>("AK8PF_MASSTRIMMED"    ,"",100,0,500   );
    AK8PF_MASSFILTERED    = fs->make<TH1D>("AK8PF_MASSFILTERED"   ,"",100,0,500   );
    AK8PF_MASSTOP         = fs->make<TH1D>("AK8PF_MASSTOP"        ,"",100,0,500   );
    AK8PF_NCONST          = fs->make<TH1D>("AK8PF_NCONST"         ,"",1000,0,1000  );
    AK8PF_AREA            = fs->make<TH1D>("AK8PF_AREA"           ,"",100,0,4     );
    AK8PF_CH_MULT         = fs->make<TH1D>("AK8PF_CH_MULT"        ,"",500,0,500   );
    AK8PF_NE_MULT         = fs->make<TH1D>("AK8PF_NE_MULT"        ,"",500,0,500   );
    AK8PF_CHEF            = fs->make<TH1D>("AK8PF_CHEF"           ,"",100,0,1     );
    AK8PF_NHEF            = fs->make<TH1D>("AK8PF_NHEF"           ,"",100,0,1     );
    AK8PF_CEEF            = fs->make<TH1D>("AK8PF_CEEF"           ,"",100,0,1     );
    AK8PF_NEEF            = fs->make<TH1D>("AK8PF_NEEF"           ,"",100,0,1     );
    AK8PF_CMEF            = fs->make<TH1D>("AK8PF_CMEF"           ,"",100,0,1     );
    
    MET_PT                = fs->make<TH1D>("MET_PT"              ,"",100,0,1000  );
    MET_PHI               = fs->make<TH1D>("MET_PHI"             ,"",100,-3.2,3.2);
    MET_SUMET             = fs->make<TH1D>("MET_SUMET"           ,"",200,1000,9000  );
    MET_GENMET            = fs->make<TH1D>("MET_GENMET"          ,"",100,0,1000  );
    


    
    
    AK8PF_tau2Overtau1 = fs->make<TH1D>("AK8PF_tau2Overtau1",   "AK8PF_tau2Overtau1",100, 0, 1);
    
    
    
    // Single lepton final state
    
    AK4PF_MULTIPLICITY_SLF      = fs->make<TH1D>("AK4PF_MULTIPLICITY_SLF"              ,"AK4PF_MULTIPLICITY_SLF",   10,-0.5,9.5     );
    AK4PF_LEADINGPT_SLF         = fs->make<TH1D>("AK4PF_LEADINGPT_SLF"                 ,"AK4PF_LEADINGPT_SLF",      100,0,1000      );
    AK4PF_LEADINGETA_SLF        = fs->make<TH1D>("AK4PF_LEADINGETA_SLF"                ,"AK4PF_LEADINGETA_SLF",     100,-5,5        );
    AK4PF_SECLEADINGPT_SLF      = fs->make<TH1D>("AK4PF_SECLEADINGPT_SLF"              ,"AK4PF_SECLEADINGPT_SLF",   100,0,1000      );
    AK4PF_SECLEADINGETA_SLF     = fs->make<TH1D>("AK4PF_SECLEADINGETA_SLF"             ,"AK4PF_SECLEADINGETA_SLF",  100,-5,5        );
    AK4PF_THIRDLEADINGPT_SLF    = fs->make<TH1D>("AK4PF_THIRDLEADINGPT_SLF"            ,"AK4PF_THIRDLEADINGPT_SLF", 100,0,1000      );
    AK4PF_THIRDLEADINGETA_SLF   = fs->make<TH1D>("AK4PF_THIRDLEADINGETA_SLF"           ,"AK4PF_THIRDLEADINGETA_SLF",100,-5,5        );
    AK4PF_FOURTHLEADINGPT_SLF   = fs->make<TH1D>("AK4PF_FOURTHLEADINGPT_SLF"           ,"AK4PF_FOURTHLEADINGPT_SLF",100,0,1000      );
    AK4PF_FOURTHLEADINGETA_SLF   = fs->make<TH1D>("AK4PF_FOURTHLEADINGETA_SLF"           ,"AK4PF_FOURTHLEADINGETA_SLF",100,-5,5     );
    
    AK4PF_CSV_LMULTIPLICITY_SLF = fs->make<TH1D>("AK4PF_CSV_LMULTIPLICITY_SLF"               ,"AK4PF_CSV_LMULTIPLICITY_SLF", 5,-0.5,4.5   );
    AK4PF_CSV_MMULTIPLICITY_SLF = fs->make<TH1D>("AK4PF_CSV_MMULTIPLICITY_SLF"               ,"AK4PF_CSV_MMULTIPLICITY_SLF", 5,-0.5,4.5   );
    
    
    AK4PF_MCSV_LEADINGPT_SLF = fs->make<TH1D>("AK4PF_MCSV_LEADINGPT_SLF"               ,"AK4PF_MCSV_LEADINGPT_SLF", 100,0,1000   );
    AK4PF_MCSV_LEADINGETA_SLF = fs->make<TH1D>("AK4PF_MCSV_LEADINGETA_SLF"               ,"AK4PF_MCSV_LEADINGETA_SLF", 100,-3,3  );
    AK4PF_MCSV_SECLEADINGPT_SLF = fs->make<TH1D>("AK4PF_MCSV_SECLEADINGPT_SLF"               ,"AK4PF_MCSV_SECLEADINGPT_SLF", 100,0,1000   );
    AK4PF_MCSV_SECLEADINGETA_SLF = fs->make<TH1D>("AK4PF_MCSV_SECLEADINGETA_SLF"               ,"AK4PF_MCSV_SECLEADINGETA_SLF",100,-3,3    );
    
    MT_SLF  = fs->make<TH1D>("MT_SLF","MT_SLF", 100,0,250  );
    pfMET_SLF = fs->make<TH1D>("pfMET_SLF"      ,"pfMET_SLF", 100,0,500  );
    
    JetClesestToMuPt = fs->make<TH1D>("JetClesestToMuPt","JetClesestToMuPt", 100,0,500  );
    ClosestJetToMuMinDR = fs->make<TH1D>("ClosestJetToMuMinDR","ClosestJetToMuMinDR", 100,0,5  );
    
    JetClesestToElePt = fs->make<TH1D>("JetClesestToElePt","JetClesestToElePt", 100,0,500  );
    ClosestJetToEleMinDR = fs->make<TH1D>("ClosestJetToEleMinDR","ClosestJetToEleMinDR", 100,0,5  );
    
    ST_SLF = fs->make<TH1D>("ST_SLF","ST_SLF", 150,0,1500  );
    HT_SLF = fs->make<TH1D>("HT_SLF","HT_SLF", 100,0,1000  );
    
    LeadingBJetEleDPhi = fs->make<TH1D>("LeadingBJetEleDPhi","LeadingBJetEleDPhi", 100,-5,5  );
    LeadingBJetMuDPhi = fs->make<TH1D>("LeadingBJetMuDPhi","LeadingBJetMuDPhi", 100,-5,5  );
    
    TwoLeadingJetsDR = fs->make<TH1D>("TwoLeadingJetsDR","TwoLeadingJetsDR", 100,0,10  );
    JetsDR = fs->make<TH1D>("JetsDR","JetsDR", 100,0,10  );
    
    
    EL_PT_SLF                 = fs->make<TH1D>("EL_PT_SLF"                ,"EL_PT_SLF",         100,0,500       );
    MU_PT_SLF                 = fs->make<TH1D>("MU_PT_SLF"                ,"MU_PT_SLF",         100,0,500       );
    
    TriMass_SLF = fs->make<TH1D>("TriMass_SLF"                ,"TriMass_SLF",         150,0,1500       );
    
    ForwardJetEta = fs->make<TH1D>("ForwardJetEta"              ,"ForwardJetEta",          100,-5,5     );
    ForwardJetPt = fs->make<TH1D>("ForwardJetPt"                ,"ForwardJetPt",        100, 0 , 500    );
    
    ForwardGenJetEta = fs->make<TH1D>("ForwardGenJetEta"              ,"ForwardGenJetEta",       100,-5,5        );
    ForwardGenJetPt  = fs->make<TH1D>("ForwardGenJetPt"               ,"ForwardGenJetPt",        100, 0 , 500    );
    
    Forward_Gen_Reco_DR = fs->make<TH1D>("Forward_Gen_Reco_DR"        ,"Forward_Gen_Reco_DR",        100, 0, 10    );
    
    
    AK8MatchedToB_MinDR = fs->make<TH1D>("AK8MatchedToB_MinDR"        ,"AK8MatchedToB_MinDR",        100, 0, 1    );
    
    // cutflow histos
    for(int iCutFlow =0; iCutFlow < NUM_CUTFLOW; iCutFlow++){
        
        std::string iCut = cut_flow_case[iCutFlow];
        
        std::string name = ("MET"+iCut).c_str();
        MET_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("MT"+iCut).c_str();
        MT_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 500);
        
        name = ("LepPt"+iCut).c_str();
        LepPt_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("LepEta"+iCut).c_str();
        LepEta_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(),  60, -3, 3);
        
        name = ("BJetPt"+iCut).c_str();
        BJetPt_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("LepJetMinDR"+iCut).c_str();
        LepJetMinDR_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 5);
        
        name = ("LepBJetDPhi"+iCut).c_str();
        LepBJetDPhi_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 35, 0, 3.5);
        
        name = ("DiJetMass"+iCut).c_str();
        DiJetMass_SLF_CUTFLOW[iCutFlow]= fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 250);
        
        name = ("TriMass"+iCut).c_str();
        TriMass_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(),  150, 0, 1500);
        
        name = ("ST"+iCut).c_str();
        ST_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 150, 0, 1500);
        
        name = ("ForwardJetEta"+iCut).c_str();
        ForwardJetEta_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, -5.0, 5.0);
        
        name = ("ForwardJetPt"+iCut).c_str();
        ForwardJetPt_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 500);
        
        name = ("LeadingJetPt"+iCut).c_str();
        AK4PF_LEADINGPT_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("SecLeadingJetPt"+iCut).c_str();
        AK4PF_SECLEADINGPT_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 500);
        
        name = ("AK8MatchedToLeadingBDR"+iCut).c_str();
        AK8LeadingBDR_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 0.5);
        
        name = ("AK8PrunedMass"+iCut).c_str();
        AK8PrunedMass_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 250);
        
        name = ("AK8PFWTagged"+iCut).c_str();
        AK8PF_WTagged_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 5, -0.5, 4.5);
        
        name = ("AK4HT"+iCut).c_str();
        AK4HT_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("AK4HTEta2p5"+iCut).c_str();
        AK4HTEta2p5_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 1000);
        
        name = ("AK4Multiplicity"+iCut).c_str();
        AK4Multiplicity_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 10, -0.5, 9.5);
        
        name = ("LooseBTagPt"+iCut).c_str();
        LooseBJetPt_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 500);
        
        name = ("LooseBJetsDR"+iCut).c_str();
        LooseBJetsDR_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 10);
        
        name = ("LeadingbClosestjetDR"+iCut).c_str();
        DRLeadingbClosestjet_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 10);
        
        name = ("PtClosestjetToB"+iCut).c_str();
        PtClosestjetToB_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 100, 0, 500);
        
        name = ("AK8MatchedToBTau1"+iCut).c_str();
        AK8MatchedToBTau1_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 10, 0, 1);
        
        name = ("AK8MatchedToBTau2"+iCut).c_str();
        AK8MatchedToBTau2_SLF_CUTFLOW[iCutFlow] = fs->make<TH1D>(name.c_str(), name.c_str(), 10, 0, 1);

    }
}

AnalyzeMiniPlusSubstructure::~AnalyzeMiniPlusSubstructure()
{
}

bool AnalyzeMiniPlusSubstructure::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
    //particle is already the ancestor
    if(ancestor == particle ) return true;
    
    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++)
    {
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}


void AnalyzeMiniPlusSubstructure::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //  bool verbose = false;
    using namespace edm;
    using namespace std;
    using namespace reco;
    //using namespace pat;
    //using namespace fastjet;
    
    
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();
    //int size = vertices->size();
    int count_vertex = 0;
    int count_good_vertex = 0;
    for(reco::VertexCollection::const_iterator v=vertices->begin();v!=vertices->end(); ++v)
    {
        count_vertex++;
        if ( v->ndof() > 4 && fabs(v->z()) < 24 && fabs(v->position().rho()) < 2 ) count_good_vertex++;
    }
    NVTX_GOOD->Fill(count_good_vertex);
    NVTX_ALL->Fill(count_vertex);
    
    edm::Handle<double> h_rho;
    // iEvent.getByLabel( "fixedGridRhoAll", h_rho );  //ex. kt6PFJets
    iEvent.getByLabel( "fixedGridRhoFastjetAll", h_rho );  //ex. kt6PFJets
    double rho = *h_rho;
    RHO->Fill(rho);
    
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    MET_PT        ->Fill( met.pt()             );
    MET_PHI       ->Fill( met.phi()            );
    MET_SUMET     ->Fill( met.sumEt()          );
    MET_GENMET    ->Fill( met.genMET()->pt()   );
    
    
    // Reco muons
    std::vector<TLorentzVector> tightMuons;
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    for (const pat::Muon &mu : *muons) {
        if ( mu.pt()  < 10 ) continue;
        if ( TMath::Abs(mu.eta()) > 2.4 )continue;  // 2.5
        if ( !mu.isTightMuon(PV)  ) continue;
        
        double chargedHadronIso = mu.pfIsolationR03().sumChargedHadronPt;
        double neutralHadronIso = mu.pfIsolationR03().sumNeutralHadronEt;
        double photonIso = mu.pfIsolationR03().sumPhotonEt;
        double beta = mu.pfIsolationR03().sumPUPt;
        double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu.pt() ;
        
        MU_ISO ->Fill( pfRelIso );
        if(pfRelIso > 0.15) continue;
        
        MU_PT               ->Fill( mu.pt() );
        MU_ETA              ->Fill( mu.eta() );
        MU_PHI              ->Fill( mu.phi() );
        MU_nStation         ->Fill( mu.numberOfMatchedStations() );
        MU_dzPV             ->Fill( fabs(mu.muonBestTrack()->dz(PV.position())) );
        MU_dxyPV            ->Fill( fabs(mu.muonBestTrack()->dxy(PV.position())) );
        
        
        TLorentzVector passedMu;
        passedMu.SetPtEtaPhiE(mu.pt(), mu.eta(), mu.phi(), mu.energy());
        tightMuons.push_back(passedMu);
    }
    
    
    // Reco electrons
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    std::vector<TLorentzVector> tightElectrons;
    
    
    for(const pat::Electron &el : *electrons){
        
        if(fabs(el.eta()) > 2.4)continue;
        if(el.pt() <= 10) continue;
        if(fabs(el.superCluster()->eta())>1.479 && fabs(el.superCluster()->eta())<1.566) continue;  // exclude the gap region

        GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
        // Compute isolation with delta beta correction for PU
        float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
        float relIsoWithDBeta_ = absiso/el.pt();
        
        float d0 = (-1) * el.gsfTrack()->dxy(PV.position());
        float dz = el.gsfTrack()->dz(PV.position());
        
        EL_SIGMAETAETA      ->Fill( el.full5x5_sigmaIetaIeta() );
        EL_HoverE           ->Fill( el.hadronicOverEm() );
        EL_DO               ->Fill( d0 );
        EL_DZ               ->Fill( dz );
        EL_ISO              ->Fill( relIsoWithDBeta_ );
    
        //barrel
        if( fabs(el.superCluster()->eta()) <  1.479){
            if( TMath::Abs(el.deltaPhiSuperClusterTrackAtVtx()) > 0.031 ) continue;  // 0.15
            if( TMath::Abs(el.deltaEtaSuperClusterTrackAtVtx()) > 0.0091 ) continue;  //0.007
            if( TMath::Abs(el.full5x5_sigmaIetaIeta()) > 0.0106 ) continue;  //0.01
            if( TMath::Abs(el.hadronicOverEm()) > 0.0532 ) continue;  //recommended is 0.12 but HLT applies 0.1
            if(TMath::Abs(d0) >  0.0126                  ) continue;
            if(TMath::Abs(dz) >  0.0116                  ) continue;
            if(TMath::Abs(relIsoWithDBeta_) >  0.1649    ) continue;
        }
        else{
            if( TMath::Abs(el.deltaPhiSuperClusterTrackAtVtx()) > 0.0359 ) continue;  // 0.1
            if( TMath::Abs(el.deltaEtaSuperClusterTrackAtVtx()) > 0.0106 ) continue; // 0.009
            if( TMath::Abs(el.full5x5_sigmaIetaIeta()) > 0.0305 ) continue;
            if( TMath::Abs(el.hadronicOverEm()) > 0.0835 ) continue;   /// at the HLT 0.075  recommended is 0.1
            if(TMath::Abs(d0) >   0.0163                 ) continue;
            if(TMath::Abs(dz) >   0.5999                 ) continue;
            if(TMath::Abs(relIsoWithDBeta_) >   0.2075   ) continue;
        }
        
        
        EL_PT               ->Fill( el.pt() );
        EL_ETA              ->Fill( el.eta() );
        EL_PHI              ->Fill( el.phi() );

        TLorentzVector passedEle;
        passedEle.SetPtEtaPhiE(el.pt(), el.eta(), el.phi(), el.energy());
        tightElectrons.push_back(passedEle);
    }
    
    TLorentzVector GENMET;
    GENMET.SetPtEtaPhiM(met.genMET()->pt(),0., met.genMET()->phi(), 0);
    
    TLorentzVector MET;
    MET.SetPtEtaPhiM(met.pt(),0., met.phi(), 0);
    
    double reco_mt = 0;
    if((tightElectrons.size() + tightMuons.size()) ==1 ){
        if(tightElectrons.size()==1){
            reco_mt = TMath::Sqrt( 2*tightElectrons[0].Pt() * met.pt() * ( 1 - TMath::Cos(tightElectrons[0].DeltaPhi(MET) ) ) );
        }
        else if(tightMuons.size()==1){
            reco_mt = TMath::Sqrt( 2*tightMuons[0].Pt() * met.pt() * ( 1 - TMath::Cos(tightMuons[0].DeltaPhi(MET) ) ) );
        }
    }
    
    if(reco_mt != 0)MT       ->Fill(reco_mt);
    
    
    // Pruned particles are the one containing "important" stuff
    Handle<edm::View<reco::GenParticle> > pruned;
    iEvent.getByToken(prunedGenToken_,pruned);
    
    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    Handle<edm::View<pat::PackedGenParticle> > packed;
    iEvent.getByToken(packedGenToken_,packed);
    
    
    
    std::vector<TLorentzVector> jetsAtMatrixElement;
    
    for(size_t i=0; i<pruned->size();i++){
        
        if((abs((*pruned)[i].status()) == 23 || abs((*pruned)[i].status()) == 33 ) &&   ( abs((*pruned)[i].pdgId()) == 1   ||
                                                                                          abs((*pruned)[i].pdgId()) == 2   ||
                                                                                          abs((*pruned)[i].pdgId()) == 3   ||
                                                                                          abs((*pruned)[i].pdgId()) == 4   ||
                                                                                          abs((*pruned)[i].pdgId()) == 5   ||
                                                                                          abs((*pruned)[i].pdgId()) == 21 )){
           
            const Candidate * motherInPrunedCollection = (*pruned)[i].mother(0) ;
            if(abs(motherInPrunedCollection->pdgId()) == 24 || abs(motherInPrunedCollection->pdgId()) == 6000006 || abs(motherInPrunedCollection->pdgId()) == 23 )continue;
            
            const Candidate * GenJet = &(*pruned)[i];
            TLorentzVector tmp;
            tmp.SetPtEtaPhiE(GenJet->pt(), GenJet->eta(), GenJet->phi(), GenJet->energy());
            jetsAtMatrixElement.push_back(tmp);
            
        }
    }
    
    
    
    
    
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    std::vector<TLorentzVector> Jets;
    std::vector<TLorentzVector> LooseBJets;
    std::vector<TLorentzVector> MediumBJets;

    bool isLepJetOverlap1p0 = false;
    for (const pat::Jet &j : *jets) {
        
        if(j.pt()< 30) continue;
        TLorentzVector passedJet;
        passedJet.SetPtEtaPhiE(j.pt(), j.eta(), j.phi(), j.energy());
        
        bool isOverlapped = false;
        for(unsigned int j = 0; j < tightMuons.size(); j++){
            
            double dr = tightMuons[j].DeltaR(passedJet);
            if( dr < 0.4) isOverlapped = true;
            
            if( dr < 1 && dr > 0.4 ) isLepJetOverlap1p0 = true;
        }
        
        for(unsigned int j = 0; j < tightElectrons.size(); j++){
            
            double dr = tightElectrons[j].DeltaR(passedJet);
            if( dr < 0.3) isOverlapped = true;
            if( dr < 1 && dr > 0.4 ) isLepJetOverlap1p0 = true;
        }
        
        
        if(isOverlapped)continue;
        

        Jets.push_back(passedJet);
        
        
        AK4PF_PT       ->Fill( j.pt()                          );
        AK4PF_PHI      ->Fill( j.phi()                         );
        AK4PF_ETA      ->Fill( j.eta()                         );
        AK4PF_RAP      ->Fill( j.rapidity()                    );
        AK4PF_MASS     ->Fill( j.mass()                        );
        AK4PF_NCONST   ->Fill( j.numberOfDaughters()           );
        AK4PF_AREA     ->Fill( j.jetArea()                     );
        AK4PF_CH_MULT  ->Fill( j.chargedMultiplicity()         );
        AK4PF_NE_MULT  ->Fill( j.neutralMultiplicity()         );
        AK4PF_CHEF     ->Fill( j.chargedHadronEnergyFraction() );
        AK4PF_NHEF     ->Fill( j.neutralHadronEnergyFraction() );
        AK4PF_CEEF     ->Fill( j.chargedEmEnergyFraction()     );
        AK4PF_NEEF     ->Fill( j.neutralEmEnergyFraction()     );
        AK4PF_CMEF     ->Fill( j.chargedMuEnergyFraction()     );
        AK4PF_CSV      ->Fill( std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"))              );
        AK4PF_PUID     ->Fill( j.userFloat("pileupJetId:fullDiscriminant")                                    );
        
        
        if(TMath::Abs(j.eta())  < 2.5){
            if(std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")) > 0.423) LooseBJets.push_back(passedJet);
            if(std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")) > 0.814) MediumBJets.push_back(passedJet);  // Medium WP: 0.814
        }
    }

    
    AK4PF_MULTIPLICITY->Fill(Jets.size());
    if(Jets.size() > 0)AK4PF_LEADINGPT->Fill(Jets[0].Pt());
    if(Jets.size() > 0)AK4PF_LEADINGETA->Fill(Jets[0].Eta());
    if(Jets.size() > 1)AK4PF_SECLEADINGPT->Fill(Jets[1].Pt());
    if(Jets.size() > 1)AK4PF_SECLEADINGETA->Fill(Jets[1].Eta());
    if(Jets.size() > 2)AK4PF_THIRDLEADINGPT->Fill(Jets[2].Pt());
    if(Jets.size() > 2)AK4PF_THIRDLEADINGETA->Fill(Jets[2].Eta());
    if(Jets.size() > 3)AK4PF_FOURTHLEADINGPT->Fill(Jets[3].Pt());
    if(Jets.size() > 3)AK4PF_FOURTHLEADINGETA->Fill(Jets[3].Eta());
    
    AK4PF_CSV_LMULTIPLICITY->Fill( LooseBJets.size() );
    AK4PF_CSV_MMULTIPLICITY->Fill( MediumBJets.size() );
    
    
    if(MediumBJets.size() > 0)AK4PF_MCSV_LEADINGPT->Fill(MediumBJets[0].Pt());
    if(MediumBJets.size() > 0)AK4PF_MCSV_LEADINGETA->Fill(MediumBJets[0].Eta());
    if(MediumBJets.size() > 1)AK4PF_MCSV_SECLEADINGPT->Fill(MediumBJets[1].Pt());
    if(MediumBJets.size() > 1)AK4PF_MCSV_SECLEADINGETA->Fill(MediumBJets[1].Eta());
    
    
    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken(fatjetToken_, fatjets);

    std::pair<TLorentzVector, double> AK8MatchedToB;
    AK8MatchedToB.second = 999.;
    double dr_min = 0.3;
    int numWtag   = 0;
    std::pair<double, double> numSubJets;
    for (const pat::Jet &j : *fatjets) {
        
        if( j.pt() < 200 )       continue;
        if( fabs(j.eta()) > 2.5) continue;
        
        double tau1 = j.userFloat("NjettinessAK8:tau1");    //
        double tau2 = j.userFloat("NjettinessAK8:tau2");    //  Access the n-subjettiness variables
    
        AK8PF_tau2Overtau1       ->Fill( tau2/tau1                                     );
        AK8PF_PT       ->Fill( j.pt()                                                  );
        AK8PF_PHI      ->Fill( j.phi()                                                 );
        AK8PF_ETA      ->Fill( j.eta()                                                 );
        AK8PF_RAP      ->Fill( j.rapidity()                                            );
        AK8PF_MASS     ->Fill( j.mass()                                                );
        AK8PF_MASSPRUNED     ->Fill( j.userFloat("ak8PFJetsCHSPrunedLinks")            );
        AK8PF_MASSTRIMMED    ->Fill( j.userFloat("ak8PFJetsCHSTrimmedLinks")           );
        AK8PF_MASSFILTERED   ->Fill( j.userFloat("ak8PFJetsCHSFilteredLinks")          );
        AK8PF_MASSTOP        ->Fill( j.userFloat("cmsTopTagPFJetsCHSLinksAK8")         );
        AK8PF_NCONST   ->Fill( j.numberOfDaughters()                                   );
        AK8PF_AREA     ->Fill( j.jetArea()                                             );
        AK8PF_CH_MULT  ->Fill( j.chargedMultiplicity()                                 );
        AK8PF_NE_MULT  ->Fill( j.neutralMultiplicity()                                 );
        AK8PF_CHEF     ->Fill( j.chargedHadronEnergyFraction()                         );
        AK8PF_NHEF     ->Fill( j.neutralHadronEnergyFraction()                         );
        AK8PF_CEEF     ->Fill( j.chargedEmEnergyFraction()                             );
        AK8PF_NEEF     ->Fill( j.neutralEmEnergyFraction()                             );
        AK8PF_CMEF     ->Fill( j.chargedMuEnergyFraction()                             );
        
        double prunedMass  = j.userFloat("ak8PFJetsCHSPrunedLinks");
        if(prunedMass > 60 && prunedMass < 100)numWtag++;
        
        if(MediumBJets.size() > 0){
            
            TLorentzVector tmp;
            tmp.SetPtEtaPhiE(j.pt() , j.eta(), j.phi(), j.energy());
            double dr = MediumBJets[0].DeltaR(tmp);
            if(dr < dr_min){
                dr_min = dr;
                AK8MatchedToB.first  = tmp;
                AK8MatchedToB.second = j.userFloat("ak8PFJetsCHSPrunedLinks");
                numSubJets.first = tau1;
                numSubJets.second = tau2;
            }
        }
    }
    
    
    

    
    
    
    // Single lepton final state  (*_SLF)
    
    if((tightElectrons.size() + tightMuons.size()) == 1 ){
        
        if(MediumBJets.size() > 0 && MediumBJets[0].Pt() > 200){
            AK8MatchedToB_MinDR->Fill(dr_min);
        }
        
        MT_SLF       ->Fill(reco_mt);
        AK4PF_MULTIPLICITY_SLF->Fill(Jets.size());
        if(Jets.size() > 0)AK4PF_LEADINGPT_SLF->Fill(Jets[0].Pt());
        if(Jets.size() > 0)AK4PF_LEADINGETA_SLF->Fill(Jets[0].Eta());
        if(Jets.size() > 1)AK4PF_SECLEADINGPT_SLF->Fill(Jets[1].Pt());
        if(Jets.size() > 1)AK4PF_SECLEADINGETA_SLF->Fill(Jets[1].Eta());
        if(Jets.size() > 2)AK4PF_THIRDLEADINGPT_SLF->Fill(Jets[2].Pt());
        if(Jets.size() > 2)AK4PF_THIRDLEADINGETA_SLF->Fill(Jets[2].Eta());
        if(Jets.size() > 3)AK4PF_FOURTHLEADINGPT_SLF->Fill(Jets[3].Pt());
        if(Jets.size() > 3)AK4PF_FOURTHLEADINGETA_SLF->Fill(Jets[3].Eta());
        
        AK4PF_CSV_LMULTIPLICITY_SLF->Fill( LooseBJets.size() );
        AK4PF_CSV_MMULTIPLICITY_SLF->Fill( MediumBJets.size() );
        
        
        if(MediumBJets.size() > 0)AK4PF_MCSV_LEADINGPT_SLF->Fill(MediumBJets[0].Pt());
        if(MediumBJets.size() > 0)AK4PF_MCSV_LEADINGETA_SLF->Fill(MediumBJets[0].Eta());
        if(MediumBJets.size() > 1)AK4PF_MCSV_SECLEADINGPT_SLF->Fill(MediumBJets[1].Pt());
        if(MediumBJets.size() > 1)AK4PF_MCSV_SECLEADINGETA_SLF->Fill(MediumBJets[1].Eta());
        
        
        if(Jets.size() > 1){
            
            TwoLeadingJetsDR->Fill(Jets[0].DeltaR(Jets[1]));
        }
        
        for(unsigned int i = 0; i < Jets.size(); i++){
            for(unsigned int j = i + 1; j < Jets.size(); j++){
                double dr = Jets[i].DeltaR(Jets[j]);
                JetsDR->Fill(dr);
            }
        }
        
        pfMET_SLF->Fill(met.pt());
        
        
        double ST = 0;
        double TriMass = 0;

        
        // leptons
        
        double deltaphi_bjet_lep = 999;
        double lepton_pt = 0;
        double lepton_eta = 999.;
        double min_dr_jet_lep = 999;
       if(tightMuons.size() == 1){
            lepton_pt = tightMuons[0].Pt();
            lepton_eta = tightMuons[0].Eta();
            MU_PT_SLF->Fill(tightMuons[0].Pt());
            if(MediumBJets.size() > 0){
                TriMass = (MediumBJets[0] + tightMuons[0] + MET ).M();
                LeadingBJetMuDPhi->Fill(tightMuons[0].DeltaPhi(MediumBJets[0]));
                ST = MediumBJets[0].Pt() + tightMuons[0].Pt() + met.pt();
                
                deltaphi_bjet_lep = TMath::Abs( tightMuons[0].DeltaPhi(MediumBJets[0]) );
            }
            double min_dr_jetmu  = 999;
            TLorentzVector JetClesestToMu;
            for(unsigned int i = 0; i < Jets.size(); i++){
                
                double dr = Jets[i].DeltaR(tightMuons[0]);
                if(dr < min_dr_jetmu){
                    min_dr_jetmu = dr;
                    JetClesestToMu = Jets[i];
                    min_dr_jet_lep = dr;
                }
            }
            if(min_dr_jetmu != 999){
                JetClesestToMuPt->Fill(JetClesestToMu.Pt());
                ClosestJetToMuMinDR->Fill(min_dr_jetmu);
            }
        }
        else if(tightElectrons.size() == 1){
            lepton_pt = tightElectrons[0].Pt();
            lepton_eta = tightElectrons[0].Eta();
            EL_PT_SLF->Fill(tightElectrons[0].Pt());
            if(MediumBJets.size() > 0){
                TriMass = (MediumBJets[0] + tightElectrons[0] + MET ).M();
                LeadingBJetEleDPhi->Fill(tightElectrons[0].DeltaPhi(MediumBJets[0]));
                ST = MediumBJets[0].Pt() + tightElectrons[0].Pt() + met.pt();
                deltaphi_bjet_lep = TMath::Abs( tightElectrons[0].DeltaPhi(MediumBJets[0]) );
            }
            double min_dr_jetele = 999;
            TLorentzVector JetClesestToEle;
            for(unsigned int i = 0; i < Jets.size(); i++){
                
                double dr = Jets[i].DeltaR(tightElectrons[0]);
                if(dr < min_dr_jetele){
                    min_dr_jetele = dr;
                    JetClesestToEle = Jets[i];
                    min_dr_jet_lep = dr;
                }
            }
            if(min_dr_jetele != 999){
                JetClesestToElePt->Fill(JetClesestToEle.Pt());
                ClosestJetToEleMinDR->Fill(min_dr_jetele);
            }
        }
        if(ST)ST_SLF ->Fill(ST);
        if(TriMass) TriMass_SLF->Fill(TriMass);
        
        double forward_jet_eta = 999.;
        double forward_jet_pt  = 0.;
        double find_forward_eta = 0.;
        TLorentzVector forward_jet_tl;
        for(unsigned int i = 0; i < Jets.size(); i++){
            
            if(MediumBJets.size() > 0 && Jets[i] == MediumBJets[0])continue;
            
            if(TMath::Abs(Jets[i].Eta()) > find_forward_eta ){
                find_forward_eta = TMath::Abs(Jets[i].Eta());
                forward_jet_eta = Jets[i].Eta();
                forward_jet_pt  = Jets[i].Pt();
                forward_jet_tl = Jets[i];
            }
        }
 
        if(TMath::Abs(forward_jet_eta) < 5){
            ForwardJetEta->Fill(forward_jet_eta);
            ForwardJetPt->Fill(forward_jet_pt);
        }
        
        double forward_genjet_eta = 999.;
        double forward_genjet_pt  = 0.;
        find_forward_eta = 0;
        TLorentzVector forward_genjet_tl;
        for(unsigned int i = 0; i < jetsAtMatrixElement.size(); i++){
            
            if(jetsAtMatrixElement[i].Pt() < 30)continue;
            if(TMath::Abs(jetsAtMatrixElement[i].Eta()) > find_forward_eta ){
                find_forward_eta = TMath::Abs(jetsAtMatrixElement[i].Eta());
                forward_genjet_eta = jetsAtMatrixElement[i].Eta();
                forward_genjet_pt  = jetsAtMatrixElement[i].Pt();
                forward_genjet_tl = jetsAtMatrixElement[i];
            }
        }
        
        if(TMath::Abs(forward_genjet_eta) < 5){
            ForwardGenJetEta->Fill(forward_genjet_eta);
            ForwardGenJetPt->Fill(forward_genjet_pt);
        }
        double forward_gen_reco_dr = -1.;
        if(forward_genjet_tl.Pt() > 0 && forward_jet_tl.Pt() > 0) forward_gen_reco_dr= forward_genjet_tl.DeltaR(forward_jet_tl);
        if(forward_gen_reco_dr)Forward_Gen_Reco_DR->Fill(forward_gen_reco_dr);
        
        
        // cut flow
        
        int cut_flow = -1;
        if(lepton_pt > 50){
            cut_flow++;
            if(TMath::Abs(lepton_eta) < 2.5){
                cut_flow++;
                if(MET.Pt() > 80){
                    cut_flow++;
                    if(reco_mt < 130){
                        cut_flow++;
                        if(MediumBJets.size() > 0 && MediumBJets[0].Pt() > 200){
                            cut_flow++;
                            if(deltaphi_bjet_lep > 2 && deltaphi_bjet_lep != 999){
                                cut_flow++;
                                if(!isLepJetOverlap1p0){
                                    cut_flow++;
                                    if(TMath::Abs(forward_jet_eta) != 999.  && TMath::Abs(forward_jet_eta) > 2){
                                        cut_flow++;
                                        if(ST > 500){
                                            cut_flow++;
                                            if(ST > 600){
                                                cut_flow++;
                                                if(ST > 700){
                                                    cut_flow++;
                                                    if(ST > 800){
                                                        cut_flow++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        double dijet_mass = 0.;
        double find_w = 999;
        for(unsigned int i = 0; i < Jets.size(); i++){
            for(unsigned int j = 0; j < Jets.size(); j++){
            
                if(i == j) continue;
                if(MediumBJets.size() > 0){
                    if(Jets[i] == MediumBJets[0] || Jets[j] == MediumBJets[0])continue;
                }
                double mass = (Jets[i]+Jets[j]).M();
                
                if(TMath::Abs(mass - 80.) < find_w){
                    find_w = TMath::Abs(mass - 80.);
                    dijet_mass = mass;
                }
            }
        }
        
        
        double dr_leadingb_closestjet = 999.;
        TLorentzVector tl_leadingb_closestjet;
        if(MediumBJets.size() > 0 && Jets.size() > 1 ){

            for(unsigned int i = 0; i < Jets.size(); i++){
                
                if(Jets[i].Pt() < 50)continue;
                if(Jets[i] == MediumBJets[0])continue;
                double dr = MediumBJets[0].DeltaR(Jets[i]);
                if(dr < dr_leadingb_closestjet){
                    dr_leadingb_closestjet = dr;
                    tl_leadingb_closestjet = Jets[i];
                }
            }
        }
        
        
        // HT calculation
        double HT = 0;
        double HT_eta_2p5 = 0;
        for(unsigned int i = 0; i < Jets.size(); i++){
            
            if(MediumBJets.size() > 0){
                if(Jets[i] == MediumBJets[0])continue;
            }
            HT = HT + Jets[i].Pt();
            
            if(TMath::Abs(Jets[i].Eta()) < 2.5)HT_eta_2p5 = HT_eta_2p5 + Jets[i].Pt();
        }
        
        HT_SLF ->Fill(HT);

        
        // Filling cutflow histos
        for( int iCutFlow = 0; iCutFlow <= cut_flow; iCutFlow++){
            
            MET_SLF_CUTFLOW[iCutFlow]->Fill(MET.Pt());
            MT_SLF_CUTFLOW[iCutFlow]->Fill(reco_mt);
            LepPt_SLF_CUTFLOW[iCutFlow]->Fill(lepton_pt);
            LepEta_SLF_CUTFLOW[iCutFlow]->Fill(lepton_eta);
            if(MediumBJets.size() > 0)BJetPt_SLF_CUTFLOW[iCutFlow]->Fill(MediumBJets[0].Pt());
            if(min_dr_jet_lep != 999)LepJetMinDR_SLF_CUTFLOW[iCutFlow]->Fill(min_dr_jet_lep);
            if(deltaphi_bjet_lep != 999)LepBJetDPhi_SLF_CUTFLOW[iCutFlow]->Fill(deltaphi_bjet_lep);
            if(dijet_mass != 0.)DiJetMass_SLF_CUTFLOW[iCutFlow]->Fill(dijet_mass);
            if(TriMass != 0)TriMass_SLF_CUTFLOW[iCutFlow]->Fill(TriMass);
            if(ST != 0)ST_SLF_CUTFLOW[iCutFlow]->Fill(ST);
            if(TMath::Abs(forward_jet_eta) != 999)ForwardJetEta_SLF_CUTFLOW[iCutFlow]->Fill(forward_jet_eta);
            if(forward_jet_pt > 0)ForwardJetPt_SLF_CUTFLOW[iCutFlow]->Fill(forward_jet_pt);
            if(MediumBJets.size() > 0)AK8LeadingBDR_SLF_CUTFLOW[iCutFlow]->Fill(MediumBJets[0].DeltaR(AK8MatchedToB.first));
            if(AK8MatchedToB.second != 999.){
                AK8PrunedMass_SLF_CUTFLOW[iCutFlow]->Fill(AK8MatchedToB.second);
                AK8MatchedToBTau1_SLF_CUTFLOW[iCutFlow]->Fill(numSubJets.first);
                AK8MatchedToBTau2_SLF_CUTFLOW[iCutFlow]->Fill(numSubJets.second);
            }
            if(Jets.size() > 0)AK4PF_LEADINGPT_SLF_CUTFLOW[iCutFlow]->Fill(Jets[0].Pt());
            if(Jets.size() > 1)AK4PF_SECLEADINGPT_SLF_CUTFLOW[iCutFlow]->Fill(Jets[1].Pt());
            AK8PF_WTagged_SLF_CUTFLOW[iCutFlow]->Fill(numWtag);
            AK4HT_SLF_CUTFLOW[iCutFlow]->Fill(HT);
            AK4HTEta2p5_SLF_CUTFLOW[iCutFlow]->Fill(HT_eta_2p5);
            AK4Multiplicity_SLF_CUTFLOW[iCutFlow]->Fill(Jets.size());
            if(LooseBJets.size() > 1)LooseBJetPt_SLF_CUTFLOW[iCutFlow]->Fill(LooseBJets[1].Pt());
            if(LooseBJets.size() > 1)LooseBJetsDR_SLF_CUTFLOW[iCutFlow]->Fill(LooseBJets[0].DeltaR(LooseBJets[1]));
            if(dr_leadingb_closestjet != 999.){
                DRLeadingbClosestjet_SLF_CUTFLOW[iCutFlow]->Fill(dr_leadingb_closestjet);
                PtClosestjetToB_SLF_CUTFLOW[iCutFlow]->Fill(tl_leadingb_closestjet.Pt());
            }
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeMiniPlusSubstructure);







