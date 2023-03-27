// ----------------------------------------------------------------------------
// 'SCorrelatorJetTree.system.h'
// Derek Anderson
// 01.18.2023
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Methods relevant to various
// internal operations are
// collected here.
//
// Derived from code by Antonio
// Silva (thanks!!)
// ----------------------------------------------------------------------------

#pragma once

using namespace std;
using namespace findNode;



// system methods -------------------------------------------------------------

void SCorrelatorJetTree::InitVariables() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::InitVariables() Initializing class members..." << endl;
  }

  // initialize class members as needed
  m_hm                   = nullptr;
  m_particleflow_mineta  = -1.1;
  m_particleflow_maxeta  = 1.1;
  m_track_minpt          = 0.;
  m_track_maxpt          = 9999.;
  m_track_mineta         = -1.1;
  m_track_maxeta         = 1.1;
  m_EMCal_cluster_minpt  = 0.;
  m_EMCal_cluster_maxpt  = 9999.;
  m_EMCal_cluster_mineta = -1.1;
  m_EMCal_cluster_maxeta = 1.1;
  m_HCal_cluster_minpt   = 0.;
  m_HCal_cluster_maxpt   = 9999.;
  m_HCal_cluster_mineta  = -1.1;
  m_HCal_cluster_maxeta  = 1.1;
  m_MC_particle_minpt    = 0.;
  m_MC_particle_maxpt    = 9999.;
  m_MC_particle_mineta   = -1.1;
  m_MC_particle_maxeta   = 1.1;
  m_add_particleflow     = true;
  m_add_tracks           = false;
  m_add_EMCal_clusters   = false;
  m_add_HCal_clusters    = false;
  m_jetr                 = 0.4;
  m_jetType              = 0;
  m_jetalgo              = antikt_algorithm;
  m_recomb_scheme        = pt_scheme;
  m_doQualityPlots       = true;
  m_save_dst             = false;
  m_recoNumJets           = 0;
  m_recoPartonID[0]       = -9999;
  m_recoPartonID[1]       = -9999;
  m_recoPartonMomX[0]     = -9999.;
  m_recoPartonMomX[1]     = -9999.;
  m_recoPartonMomY[0]     = -9999.;
  m_recoPartonMomY[1]     = -9999.;
  m_recoPartonMomZ[0]     = -9999.;
  m_recoPartonMomZ[1]     = -9999.;
  m_trueNumJets           = 0;
  m_truePartonID[0]       = -9999;
  m_truePartonID[1]       = -9999;
  m_truePartonMomX[0]     = -9999.;
  m_truePartonMomX[1]     = -9999.;
  m_truePartonMomY[0]     = -9999.;
  m_truePartonMomY[1]     = -9999.;
  m_truePartonMomZ[0]     = -9999.;
  m_truePartonMomZ[1]     = -9999.;
  m_recoJetNCst.clear();
  m_recoJetId.clear();
  m_recoJetTruId.clear();
  m_recoJetE.clear();
  m_recoJetPt.clear();
  m_recoJetEta.clear();
  m_recoJetPhi.clear();
  m_recoJetArea.clear();
  m_trueJetNCst.clear();
  m_trueJetId.clear();
  m_trueJetTruId.clear();
  m_trueJetE.clear();
  m_trueJetPt.clear();
  m_trueJetEta.clear();
  m_trueJetPhi.clear();
  m_trueJetArea.clear();
  m_recoCstZ.clear();
  m_recoCstDr.clear();
  m_recoCstE.clear();
  m_recoCstJt.clear();
  m_recoCstEta.clear();
  m_recoCstPhi.clear();
  m_trueCstZ.clear();
  m_trueCstDr.clear();
  m_trueCstE.clear();
  m_trueCstJt.clear();
  m_trueCstEta.clear();
  m_trueCstPhi.clear();
  m_outFile = new TFile();
  return;

}  // end 'InitVariables()'



void SCorrelatorJetTree::InitHists() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::InitHists() Initializing QA histograms..." << endl;
  }

  // binning
  const unsigned long nNumBins(500);
  const unsigned long nEneBins(200);
  const unsigned long nPhiBins(63);
  const unsigned long nEtaBins(40);
  const unsigned long nPtBins(200);
  const unsigned long nAreaBins(1000);
  const double        rNumBins[NRange]  = {0.,    500.};
  const double        rPhiBins[NRange]  = {-3.15, 3.15};
  const double        rEtaBins[NRange]  = {-2.,   2.};
  const double        rEneBins[NRange]  = {0.,    100.};
  const double        rPtBins[NRange]   = {0.,    100.};
  const double        rAreaBins[NRange] = {0.,    10.};

  m_outFile -> cd();
  // no. of objects in acceptance
  m_hNumObject[OBJECT::TRACK]             = new TH1D("hNumTrks",       "N_{accept} (tracks)",       nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::ECLUST]            = new TH1D("hNumEClust",     "N_{accept} (EMCal clust.)", nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::HCLUST]            = new TH1D("hNumHClust",     "N_{accept} (HCal clust.)",  nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::FLOW]              = new TH1D("hNumFlow",       "N_{accept} (flow)",         nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::PART]              = new TH1D("hNumPart",       "N_{accept} (par.s)",        nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::TJET]              = new TH1D("hNumTruthJet",   "N_{jet} (truth)",           nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::RJET]              = new TH1D("hNumRecoJets",   "N_{jet} (reco.)",           nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::TCST]              = new TH1D("hNumTruthCst",   "N_{cst} (truth)",           nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumObject[OBJECT::RCST]              = new TH1D("hNumRecoCst",    "N_{cst} (reco.)",           nNumBins,  rNumBins[0],  rNumBins[1]);
  // sum of cst. energies
  m_hSumCstEne[CST_TYPE::TRACK_CST]       = new TH1D("hSumTrackEne",   "#SigmaE (cst. track)",      nEneBins,  rEneBins[0],  rEneBins[1]);
  m_hSumCstEne[CST_TYPE::CALO_CST]        = new TH1D("hSumCaloEne",    "#SigmaE (cst. calo.)",      nEneBins,  rEneBins[0],  rEneBins[1]);
  m_hSumCstEne[CST_TYPE::FLOW_CST]        = new TH1D("hSumFlowEne",    "#SigmaE (cst. flow)",       nEneBins,  rEneBins[0],  rEneBins[1]);
  m_hSumCstEne[CST_TYPE::PART_CST]        = new TH1D("hSumPartEne",    "#SigmaE (cst. par.s)",      nEneBins,  rEneBins[0],  rEneBins[1]);
  // track QA
  m_hObjectQA[OBJECT::TRACK][INFO::PT]    = new TH1D("hTrackPt",       "p_{T} (tracks)",            nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::TRACK][INFO::ETA]   = new TH1D("hTrackEta",      "#eta (tracks)",             nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::TRACK][INFO::PHI]   = new TH1D("hTrackPhi",      "#phi (tracks)",             nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::TRACK][INFO::ENE]   = new TH1D("hTrackEne",      "E (tracks)",                nEneBins,  rEneBins[0],  rEneBins[1]);
  // emcal cluster QA
  m_hObjectQA[OBJECT::ECLUST][INFO::PT]   = new TH1D("hEClustPt",      "p_{T} (EMCal clust.)",      nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::ECLUST][INFO::ETA]  = new TH1D("hEClustEta",     "#eta (EMCal clust.)",       nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::ECLUST][INFO::PHI]  = new TH1D("hEClustPhi",     "#phi (EMCal clust.)",       nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::ECLUST][INFO::ENE]  = new TH1D("hEClustEne",     "E (EMCal clust.)",          nEneBins,  rEneBins[0],  rEneBins[1]);
  // hcal cluster QA
  m_hObjectQA[OBJECT::HCLUST][INFO::PT]   = new TH1D("hHClustPt",      "p_{T} (HCal clust.)",       nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::HCLUST][INFO::ETA]  = new TH1D("hHClustEta",     "#eta (HCal clust.)",        nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::HCLUST][INFO::PHI]  = new TH1D("hHClustPhi",     "#phi (HCal clust.)",        nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::HCLUST][INFO::ENE]  = new TH1D("hHClustEne",     "E (HCal clust.)",           nEneBins,  rEneBins[0],  rEneBins[1]);
  // particle flow QA
  m_hObjectQA[OBJECT::FLOW][INFO::PT]     = new TH1D("hFlowPt",        "p_{T} (flow)",              nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::FLOW][INFO::ETA]    = new TH1D("hFlowEta",       "#eta (flow)",               nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::FLOW][INFO::PHI]    = new TH1D("hFlowPhi",       "#phi (flow)",               nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::FLOW][INFO::ENE]    = new TH1D("hFlowEne",       "E (flow)",                  nEneBins,  rEneBins[0],  rEneBins[1]);
  // particle QA
  m_hObjectQA[OBJECT::PART][INFO::PT]     = new TH1D("hPartPt",        "p_{T} (par.s)",             nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::PART][INFO::ETA]    = new TH1D("hPartEta",       "#eta (par.s)",              nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::PART][INFO::PHI]    = new TH1D("hPartPhi",       "#phi (par.s)",              nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::PART][INFO::ENE]    = new TH1D("hPartEne",       "E (par.s)",                 nEneBins,  rEneBins[0],  rEneBins[1]);
  // truth jet QA
  m_hObjectQA[OBJECT::TJET][INFO::PT]     = new TH1D("hTruthJetPt",    "p_{T} (truth jet)",         nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::TJET][INFO::ETA]    = new TH1D("hTruthJetEta",   "#eta (truth jet)",          nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::TJET][INFO::PHI]    = new TH1D("hTruthJetPhi",   "#phi (truth jet)",          nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::TJET][INFO::ENE]    = new TH1D("hTruthJetEne",   "E (truth jet)",             nEneBins,  rEneBins[0],  rEneBins[1]);
  m_hJetArea[0]                           = new TH1D("hTruthJetArea",  "Area (truth jet)",          nAreaBins, rAreaBins[0], rAreaBins[1]);
  m_hJetNumCst[0]                         = new TH1D("hTruthJetNCst",  "N_{cst} (truth jet)",       nNumBins,  rNumBins[0],  rNumBins[1]);
  // reco jet QA
  m_hObjectQA[OBJECT::RJET][INFO::PT]     = new TH1D("hRecoJetPt",     "p_{T} (reco. jet)",         nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::RJET][INFO::ETA]    = new TH1D("hRecoJetEta",    "#eta (reco. jet)",          nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::RJET][INFO::PHI]    = new TH1D("hRecoJetPhi",    "#phi (reco. jet)",          nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::RJET][INFO::ENE]    = new TH1D("hRecoJetEne",    "E (reco. jet)",             nEneBins,  rEneBins[0],  rEneBins[1]);
  m_hJetArea[1]                           = new TH1D("hRecoJetArea",   "Area (reco. jet)",          nAreaBins, rAreaBins[0], rAreaBins[1]);
  m_hJetNumCst[1]                         = new TH1D("hRecoJetNCst",   "N_{cst} (reco. jet)",       nNumBins,  rNumBins[0],  rNumBins[1]);
  // truth cst. QA
  m_hObjectQA[OBJECT::TCST][INFO::PT]     = new TH1D("hTruthCstPt",    "p_{T} (truth cst.)",        nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::TCST][INFO::ETA]    = new TH1D("hTruthCstEta",   "#eta (truth cst.)",         nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::TCST][INFO::PHI]    = new TH1D("hTruthCstPhi",   "#phi (truth cst.)",         nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::TCST][INFO::ENE]    = new TH1D("hTruthCstEne",   "E (truth cst.)",            nEneBins,  rEneBins[0],  rEneBins[1]);
  // reco cst. QA
  m_hObjectQA[OBJECT::RCST][INFO::PT]     = new TH1D("hRecoCstPt",     "p_{T} (reco. cst.)",        nPtBins,   rPtBins[0],   rPtBins[1]);
  m_hObjectQA[OBJECT::RCST][INFO::ETA]    = new TH1D("hRecoCstEta",    "#eta (reco. cst.)",         nEtaBins,  rEtaBins[0],  rEtaBins[1]);
  m_hObjectQA[OBJECT::RCST][INFO::PHI]    = new TH1D("hRecoCstPhi",    "#phi (reco. cst.)",         nPhiBins,  rPhiBins[0],  rPhiBins[1]);
  m_hObjectQA[OBJECT::RCST][INFO::ENE]    = new TH1D("hRecoCstEne",    "E (reco. cst.)",            nEneBins,  rEneBins[0],  rEneBins[1]);
  // no. of cst.s
  m_hNumCstAccept[CST_TYPE::TRACK_CST][0] = new TH1D("hNumTrkCstTot",  "N_{cst}^{trk} total",       nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::TRACK_CST][1] = new TH1D("hNumTrkCstAcc",  "N_{cst}^{trk} accepted",    nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::CALO_CST][0]  = new TH1D("hNumCaloCstTot", "N_{cst}^{clust} total",     nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::CALO_CST][1]  = new TH1D("hNumCaloCstAcc", "N_{cst}^{clust} accepted",  nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::FLOW_CST][0]  = new TH1D("hNumFlowCstTot", "N_{cst}^{flow} total",      nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::FLOW_CST][1]  = new TH1D("hNumFlowCstAcc", "N_{cst}^{flow} accepted",   nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::PART_CST][0]  = new TH1D("hNumPartCstTot", "N_{cst}^{par} total",       nNumBins,  rNumBins[0],  rNumBins[1]);
  m_hNumCstAccept[CST_TYPE::PART_CST][1]  = new TH1D("hNumPartCstAcc", "N_{cst}^{par} accepted",    nNumBins,  rNumBins[0],  rNumBins[1]);

}  // end 'InitHists()'



void SCorrelatorJetTree::InitTrees() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::InitTrees() Initializing output trees..." << endl;
  }

  // initialize output trees
  m_recoTree = new TTree("RecoJetTree", "A tree of reconstructed jets");
  m_recoTree -> Branch("EvtNumJets",   &m_recoNumJets,       "NumJets/I");
  m_recoTree -> Branch("Parton3_ID",   &m_recoPartonID[0],   "Parton3_ID/I");
  m_recoTree -> Branch("Parton4_ID",   &m_recoPartonID[1],   "Parton4_ID/I");
  m_recoTree -> Branch("Parton3_MomX", &m_recoPartonMomX[0], "Parton3_MomX/D");
  m_recoTree -> Branch("Parton3_MomY", &m_recoPartonMomY[0], "Parton3_MomY/D");
  m_recoTree -> Branch("Parton3_MomZ", &m_recoPartonMomZ[0], "Parton3_MomZ/D");
  m_recoTree -> Branch("Parton4_MomX", &m_recoPartonMomX[1], "Parton4_MomX/D");
  m_recoTree -> Branch("Parton4_MomY", &m_recoPartonMomY[1], "Parton4_MomY/D");
  m_recoTree -> Branch("Parton4_MomZ", &m_recoPartonMomZ[1], "Parton4_MomZ/D");
  m_recoTree -> Branch("JetNumCst",    &m_recoJetNCst);
  m_recoTree -> Branch("JetID",        &m_recoJetId);
  m_recoTree -> Branch("JetTruthID",   &m_recoJetTruId);
  m_recoTree -> Branch("JetEnergy",    &m_recoJetE);
  m_recoTree -> Branch("JetPt",        &m_recoJetPt);
  m_recoTree -> Branch("JetEta",       &m_recoJetEta);
  m_recoTree -> Branch("JetPhi",       &m_recoJetPhi);
  m_recoTree -> Branch("JetArea",      &m_recoJetArea);
  m_recoTree -> Branch("CstZ",         &m_recoCstZ);
  m_recoTree -> Branch("CstDr",        &m_recoCstDr);
  m_recoTree -> Branch("CstEnergy",    &m_recoCstE);
  m_recoTree -> Branch("CstJt",        &m_recoCstJt);
  m_recoTree -> Branch("CstEta",       &m_recoCstEta);
  m_recoTree -> Branch("CstPhi",       &m_recoCstPhi);

  m_trueTree = new TTree("TruthJetTree", "A tree of truth jets");
  m_trueTree -> Branch("EvtNumJets",   &m_trueNumJets,       "NumJets/I");
  m_trueTree -> Branch("Parton3_ID",   &m_truePartonID[0],   "Parton3_ID/I");
  m_trueTree -> Branch("Parton4_ID",   &m_truePartonID[1],   "Parton4_ID/I");
  m_trueTree -> Branch("Parton3_MomX", &m_truePartonMomX[0], "Parton3_MomX/D");
  m_trueTree -> Branch("Parton3_MomY", &m_truePartonMomY[0], "Parton3_MomY/D");
  m_trueTree -> Branch("Parton3_MomZ", &m_truePartonMomZ[0], "Parton3_MomZ/D");
  m_trueTree -> Branch("Parton4_MomX", &m_truePartonMomX[1], "Parton4_MomX/D");
  m_trueTree -> Branch("Parton4_MomY", &m_truePartonMomY[1], "Parton4_MomY/D");
  m_trueTree -> Branch("Parton4_MomZ", &m_truePartonMomZ[1], "Parton4_MomZ/D");
  m_trueTree -> Branch("JetNumCst",    &m_trueJetNCst);
  m_trueTree -> Branch("JetID",        &m_trueJetId);
  m_trueTree -> Branch("JetTruthID",   &m_trueJetTruId);
  m_trueTree -> Branch("JetEnergy",    &m_trueJetE);
  m_trueTree -> Branch("JetPt",        &m_trueJetPt);
  m_trueTree -> Branch("JetEta",       &m_trueJetEta);
  m_trueTree -> Branch("JetPhi",       &m_trueJetPhi);
  m_trueTree -> Branch("JetArea",      &m_trueJetArea);
  m_trueTree -> Branch("CstZ",         &m_trueCstZ);
  m_trueTree -> Branch("CstDr",        &m_trueCstDr);
  m_trueTree -> Branch("CstEnergy",    &m_trueCstE);
  m_trueTree -> Branch("CstJt",        &m_trueCstJt);
  m_trueTree -> Branch("CstEta",       &m_trueCstEta);
  m_trueTree -> Branch("CstPhi",       &m_trueCstPhi);
  return;

}  // end 'InitTrees()'



void SCorrelatorJetTree::FillTrees() {

  /* filling trees will go here */
  return;

}  // end 'FillTrees()'



void SCorrelatorJetTree::SaveOutput() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::SaveOutput() Saving output trees and histograms..." << endl;
  }

  // save QA histograms if need be
  const string  sQuality[NDirectory + 1] = {"Tracks", "CaloClusters", "ParticleFlow", "Particles", "TruthJets", "RecoJets", "QA"};
  TDirectory   *dQuality[NDirectory + 1];
  if (m_doQualityPlots) {

    // create QA directories
    dQuality[NDirectory] = (TDirectory*) m_outFile -> mkdir(sQuality[NDirectory].data());
    for (size_t iDirect = 0; iDirect < NDirectory; iDirect++) {
      dQuality[iDirect] = (TDirectory*) dQuality[NDirectory] -> mkdir(sQuality[iDirect].data());
    }

    // save object-specific QA hists
    for (size_t iObj = OBJECT::TRACK; iObj < NObjType; iObj++) {
      switch (iObj) {
        case OBJECT::TRACK:
          dQuality[0] -> cd();
          break;
        case OBJECT::HCLUST:
          dQuality[1] -> cd();
          break;
        case OBJECT::ECLUST:
          dQuality[1] -> cd();
          break;
        case OBJECT::FLOW:
          dQuality[2] -> cd();
          break;
        case OBJECT::PART:
          dQuality[3] -> cd();
          break;
        case OBJECT::TJET:
          dQuality[4] -> cd();
          break;
        case OBJECT::RJET:
          dQuality[5] -> cd();
          break;
        case OBJECT::TCST:
          dQuality[4] -> cd();
          break;
        case OBJECT::RCST:
          dQuality[5] -> cd();
          break;
      }
      m_hNumObject[iObj] -> Write();
      for (size_t iInfo = INFO::PT; iInfo < NInfoQA; iInfo++) {
        m_hObjectQA[iObj][iInfo] -> Write();
      }
    }  // end object loop

    // save cst-specific histograms
    for (size_t iCst = CST_TYPE::TRACK_CST; iCst < NCstType; iCst++) {
      switch (iCst) {
        case CST_TYPE::TRACK_CST:
          dQuality[0] -> cd();
          break;
        case CST_TYPE::CALO_CST:
          dQuality[1] -> cd();
          break;
        case CST_TYPE::FLOW_CST:
          dQuality[2] -> cd();
          break;
        case CST_TYPE::PART_CST:
          dQuality[3] -> cd();
          break;
      }
      m_hNumCstAccept[iCst][0] -> Write();
      m_hNumCstAccept[iCst][1] -> Write();
      m_hSumCstEne[iCst]       -> Write();
    }  // end cst loop

    // save jet-specific histograms
    dQuality[4]     -> cd();
    m_hJetArea[0]   -> Write();
    m_hJetNumCst[0] -> Write();
    dQuality[5]     -> cd();
    m_hJetArea[1]   -> Write();
    m_hJetNumCst[1] -> Write();
  }

  // save output trees
  m_outFile -> cd();
  m_recoTree -> Write();
  m_trueTree -> Write();
  return;

}  // end 'SaveOutput()'



void SCorrelatorJetTree::ResetTreeVariables() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::ResetTreeVariables() Resetting tree variables..." << endl;
  }
  m_recoNumJets       = 0;
  m_recoPartonID[0]   = -9999;
  m_recoPartonID[1]   = -9999;
  m_recoPartonMomX[0] = -9999.;
  m_recoPartonMomX[1] = -9999.;
  m_recoPartonMomY[0] = -9999.;
  m_recoPartonMomY[1] = -9999.;
  m_recoPartonMomZ[0] = -9999.;
  m_recoPartonMomZ[1] = -9999.;
  m_trueNumJets       = 0;
  m_truePartonID[0]   = -9999;
  m_truePartonID[1]   = -9999;
  m_truePartonMomX[0] = -9999.;
  m_truePartonMomX[1] = -9999.;
  m_truePartonMomY[0] = -9999.;
  m_truePartonMomY[1] = -9999.;
  m_truePartonMomZ[0] = -9999.;
  m_truePartonMomZ[1] = -9999.;
  m_recoJetNCst.clear();
  m_recoJetId.clear();
  m_recoJetTruId.clear();
  m_recoJetE.clear();
  m_recoJetPt.clear();
  m_recoJetEta.clear();
  m_recoJetPhi.clear();
  m_trueJetNCst.clear();
  m_trueJetId.clear();
  m_trueJetTruId.clear();
  m_trueJetE.clear();
  m_trueJetPt.clear();
  m_trueJetEta.clear();
  m_trueJetPhi.clear();
  m_recoCstZ.clear();
  m_recoCstDr.clear();
  m_recoCstE.clear();
  m_recoCstJt.clear();
  m_recoCstEta.clear();
  m_recoCstPhi.clear();
  m_trueCstZ.clear();
  m_trueCstDr.clear();
  m_trueCstE.clear();
  m_trueCstJt.clear();
  m_trueCstEta.clear();
  m_trueCstPhi.clear();
  return;

}  // end 'ResetTreeVariables()



int SCorrelatorJetTree::CreateJetNode(PHCompositeNode* topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::CreateJetNode(PHCompositeNode*) Creating jet node..." << endl;
  }

  // create iterator & DST node
  PHNodeIterator   iter(topNode);
  PHCompositeNode *lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!lowerNode) {
    lowerNode = new PHCompositeNode("DST");
    topNode   -> addNode(lowerNode);
    cout << "DST node added" << endl;
  }

  // construct jet tree name
  string baseName;
  string jetNodeName;
  string jetNodeNameMC;
  if (m_jetcontainer_name.empty()) {
    baseName = "JetTree";
  } else {
    baseName = m_jetcontainer_name;
  }

  // can't have forward slashes in DST or else you make a subdirectory on save!!!
  string undrscr = "_";
  string nothing = "";

  // define good strings to replace bad ones
  map<string, string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings) {
    size_t pos;
    while ((pos = baseName.find(badString)) != string::npos) {
      baseName.replace(pos, 1, goodString);
    }
  }

  // construct jet node name
  jetNodeName   = baseName + "_Jet_Container";
  jetNodeNameMC = baseName + "_MC_Jet_Container";

  // construct jet maps
  m_jetMap = new JetMapv1();
  if (m_ismc && m_save_truth_dst) {
    m_truth_jetMap = new JetMapv1();
  }

  // add jet node
  PHIODataNode<PHObject>* jetNode = new PHIODataNode<PHObject>(m_jetMap, jetNodeName.c_str(), "PHObject");
  lowerNode -> addNode(jetNode);
  cout << jetNodeName << " node added" << endl;

  // save truth DST if needed
  if(m_ismc && m_save_truth_dst) {
    PHIODataNode<PHObject> *jetNodeMC = new PHIODataNode<PHObject>(m_truth_jetMap, jetNodeNameMC.c_str(), "PHObject");
    lowerNode -> addNode(jetNodeMC);
    cout << jetNodeNameMC << " node added" << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'CreateJetNode(PHCompositeNode*)'

// end ------------------------------------------------------------------------
