// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.sys.h'
// Derek Anderson
// 01.18.2023
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Initially derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#pragma once

using namespace std;
using namespace findNode;



namespace SColdQcdCorrelatorAnalysis {

  // system methods -----------------------------------------------------------

  void SCorrelatorJetTreeMaker::InitVariables() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::InitVariables() Initializing class members..." << endl;
    }
    // initialize parton and other variables
    m_trueVtx      = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_recoVtx      = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_partonID[0]  = -9999;
    m_partonID[1]  = -9999;
    m_partonMom[0] = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_partonMom[1] = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_recoJets.clear();
    m_trueJets.clear();
    m_vecEvtsToGrab.clear();
    m_mapCstToEmbedID.clear();
    return;

  }  // end 'InitVariables()'



  void SCorrelatorJetTreeMaker::InitFuncs() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::InitFuncs() Initializing functions for internal calculations..." << endl;
    }

    // pt range of functions
    const pair<float, float> ptRange = {0., 100.};

    // initialize functions
    m_fSigDcaXY = new TF1("fSigmaDcaXY", "[0]+[1]/x+[2]/(x*x)", ptRange.first, ptRange.second);
    m_fSigDcaZ  = new TF1("fSigmaDcaZ",  "[0]+[1]/x+[2]/(x*x)", ptRange.first, ptRange.second);
    for (uint8_t iParam = 0; iParam < CONST::NParam; iParam++) {
      m_fSigDcaXY -> SetParameter(iParam, m_parSigDcaXY[iParam]);
      m_fSigDcaZ  -> SetParameter(iParam, m_parSigDcaZ[iParam]);
    }
    return;

  }  // end 'InitFuncs()'



  void SCorrelatorJetTreeMaker::InitTrees() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::InitTrees() Initializing output trees..." << endl;
    }

    // initialize trees
    m_trueTree = new TTree("TruthJetTree", "A tree of truth jets");
    m_recoTree = new TTree("RecoJetTree",  "A tree of reconstructed jets");
    return;

  }  // end 'InitTrees()'



  void SCorrelatorJetTreeMaker::InitEvals(PHCompositeNode* topNode) {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::InitEvals(PHCompositeNode*) Initializing evaluators..." << endl;
    }

    m_evalStack = new SvtxEvalStack(topNode);
    if (!m_evalStack) {
      cerr << "SCorrelatorJetTreeMaker::InitEvals(PHCompositeNode*) PANIC: couldn't grab SvtxEvalStack! Aborting!" << endl;
      assert(m_evalStack);
    } else {
      m_evalStack -> next_event(topNode);
    }

    m_trackEval = m_evalStack -> get_track_eval();
    if (!m_trackEval) {
      cerr << "SCorrelatorJetTreeMaker::InitEvals(PHCompositeNode*) PANIC: couldn't grab track evaluator! Aborting!" << endl;
      assert(m_trackEval);
    }
    return;

  }  // end 'InitEvals(PHCompositeNode*)'



  void SCorrelatorJetTreeMaker::FillTrueTree() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::FillTrueTree() Filling truth jet tree..." << endl;
    }

    // prepare vectors for filling
    m_trueJetNCst.clear();
    m_trueJetID.clear();
    m_trueJetE.clear();
    m_trueJetPt.clear();
    m_trueJetEta.clear();
    m_trueJetPhi.clear();
    m_trueJetArea.clear();
    m_trueCstID.clear();
    m_trueCstEmbedID.clear();
    m_trueCstZ.clear();
    m_trueCstDr.clear();
    m_trueCstE.clear();
    m_trueCstPt.clear();
    m_trueCstEta.clear();
    m_trueCstPhi.clear();

    // declare vectors to storing constituents
    vector<int>    vecTruCstID;
    vector<int>    vecTruCstEmbedID;
    vector<double> vecTruCstZ;
    vector<double> vecTruCstDr;
    vector<double> vecTruCstE;
    vector<double> vecTruCstPt;
    vector<double> vecTruCstEta;
    vector<double> vecTruCstPhi;
    vecTruCstID.clear();
    vecTruCstZ.clear();
    vecTruCstDr.clear();
    vecTruCstE.clear();
    vecTruCstPt.clear();
    vecTruCstEta.clear();
    vecTruCstPhi.clear();

    // fill jets/constituent variables
    unsigned int nTruJet(0);
    unsigned int nTruCst(0);
    for (unsigned int iTruJet = 0; iTruJet < m_trueJets.size(); ++iTruJet) {

      // get jet info
      const unsigned int jetNCst  = m_trueJets[iTruJet].constituents().size();
      const unsigned int jetTruID = iTruJet;
      const double       jetPhi   = m_trueJets[iTruJet].phi_std();
      const double       jetEta   = m_trueJets[iTruJet].pseudorapidity();
      const double       jetArea  = 0.;  // FIXME: jet area needs to be defined
      const double       jetE     = m_trueJets[iTruJet].E();
      const double       jetPt    = m_trueJets[iTruJet].perp();
      const double       jetPx    = m_trueJets[iTruJet].px();
      const double       jetPy    = m_trueJets[iTruJet].py();
      const double       jetPz    = m_trueJets[iTruJet].pz();
      const double       jetP     = sqrt((jetPx * jetPx) + (jetPy * jetPy) + (jetPz * jetPz));

      // clear constituent vectors
      vecTruCstID.clear();
      vecTruCstEmbedID.clear();
      vecTruCstZ.clear();
      vecTruCstDr.clear();
      vecTruCstE.clear();
      vecTruCstPt.clear();
      vecTruCstEta.clear();
      vecTruCstPhi.clear();

      // loop over constituents
      vector<fastjet::PseudoJet> trueCsts = m_trueJets[iTruJet].constituents();
      for (unsigned int iTruCst = 0; iTruCst < trueCsts.size(); ++iTruCst) {

        // get constituent info
        const double cstPhi = trueCsts[iTruCst].phi_std();
        const double cstEta = trueCsts[iTruCst].pseudorapidity();
        const double cstE   = trueCsts[iTruCst].E();
        const double cstPt  = trueCsts[iTruCst].perp();
        const double cstPx  = trueCsts[iTruCst].px();
        const double cstPy  = trueCsts[iTruCst].py();
        const double cstPz  = trueCsts[iTruCst].pz();
        const double cstP   = ((cstPx * cstPx) + (cstPy * cstPy) + (cstPz * cstPz));
        const double cstZ   = cstP / jetP;
        const double cstDf  = cstPhi - jetPhi;
        const double cstDh  = cstEta - jetEta;
        const double cstDr  = sqrt((cstDf * cstDf) + (cstDh * cstDh));

        // get barcode and embedding ID
        const int cstID   = trueCsts[iTruCst].user_index();
        const int embedID = m_mapCstToEmbedID[cstID];

        // add csts to vectors
        vecTruCstID.push_back(abs(cstID));
        vecTruCstEmbedID.push_back(embedID);
        vecTruCstZ.push_back(cstZ);
        vecTruCstDr.push_back(cstDr);
        vecTruCstE.push_back(cstE);
        vecTruCstPt.push_back(cstPt);
        vecTruCstEta.push_back(cstEta);
        vecTruCstPhi.push_back(cstPhi);

        // fill QA histograms and increment counters
        m_hObjectQA[OBJECT::TCST][INFO::PT]  -> Fill(cstPt);
        m_hObjectQA[OBJECT::TCST][INFO::ETA] -> Fill(cstEta);
        m_hObjectQA[OBJECT::TCST][INFO::PHI] -> Fill(cstPhi);
        m_hObjectQA[OBJECT::TCST][INFO::ENE] -> Fill(cstE);
        ++nTruCst;
      }  // end constituent loop

      // store jet/cst output
      m_trueJetNCst.push_back(jetNCst);
      m_trueJetID.push_back(jetTruID);
      m_trueJetE.push_back(jetE);
      m_trueJetPt.push_back(jetPt);
      m_trueJetEta.push_back(jetEta);
      m_trueJetPhi.push_back(jetPhi);
      m_trueJetArea.push_back(jetArea);
      m_trueCstID.push_back(vecTruCstID);
      m_trueCstEmbedID.push_back(vecTruCstEmbedID);
      m_trueCstZ.push_back(vecTruCstZ);
      m_trueCstDr.push_back(vecTruCstDr);
      m_trueCstE.push_back(vecTruCstE);
      m_trueCstPt.push_back(vecTruCstPt);
      m_trueCstEta.push_back(vecTruCstEta);
      m_trueCstPhi.push_back(vecTruCstPhi);

      // fill QA histograms and increment counters
      m_hJetArea[0]                        -> Fill(jetArea);
      m_hJetNumCst[0]                      -> Fill(jetNCst);
      m_hObjectQA[OBJECT::TJET][INFO::PT]  -> Fill(jetPt);
      m_hObjectQA[OBJECT::TJET][INFO::ETA] -> Fill(jetEta);
      m_hObjectQA[OBJECT::TJET][INFO::PHI] -> Fill(jetPhi);
      m_hObjectQA[OBJECT::TJET][INFO::ENE] -> Fill(jetE);
      ++nTruJet;
    }  // end jet loop

    // fill QA histograms
    m_hNumObject[OBJECT::TJET] -> Fill(nTruJet);
    m_hNumObject[OBJECT::TCST] -> Fill(nTruCst);

    // store evt info
    m_trueNumJets       = nTruJet;
    m_truePartonID[0]   = m_partonID[0];
    m_truePartonID[1]   = m_partonID[1];
    m_truePartonMomX[0] = m_partonMom[0].x();
    m_truePartonMomX[1] = m_partonMom[1].x();
    m_truePartonMomY[0] = m_partonMom[0].y();
    m_truePartonMomY[1] = m_partonMom[1].y();
    m_truePartonMomZ[0] = m_partonMom[0].z();
    m_truePartonMomZ[1] = m_partonMom[1].z(); 
    m_trueVtxX          = m_trueVtx.x();
    m_trueVtxY          = m_trueVtx.y();
    m_trueVtxZ          = m_trueVtx.z();

    // fill output tree
    m_trueTree -> Fill();
    return;

  }  // end 'FillTrueTree()'



  void SCorrelatorJetTreeMaker::FillRecoTree() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::FillRecoTree() Filling reco jet tree..." << endl;
    }

    // prepare vectors for filling
    m_recoJetNCst.clear();
    m_recoJetID.clear();
    m_recoJetE.clear();
    m_recoJetPt.clear();
    m_recoJetEta.clear();
    m_recoJetPhi.clear();
    m_recoJetArea.clear();
    m_recoCstMatchID.clear();
    m_recoCstZ.clear();
    m_recoCstDr.clear();
    m_recoCstE.clear();
    m_recoCstPt.clear();
    m_recoCstEta.clear();
    m_recoCstPhi.clear();

    // declare vectors for storing constituents
    vector<int>    vecRecCstMatchID;
    vector<double> vecRecCstZ;
    vector<double> vecRecCstDr;
    vector<double> vecRecCstE;
    vector<double> vecRecCstPt;
    vector<double> vecRecCstEta;
    vector<double> vecRecCstPhi;
    vecRecCstMatchID.clear();
    vecRecCstZ.clear();
    vecRecCstDr.clear();
    vecRecCstE.clear();
    vecRecCstPt.clear();
    vecRecCstEta.clear();
    vecRecCstPhi.clear();

    // fill jet/constituent variables
    unsigned long nRecJet(0);
    unsigned long nRecCst(0);
    for (unsigned int iJet = 0; iJet < m_recoJets.size(); ++iJet) {

      // get jet info
      const unsigned int jetNCst  = m_recoJets[iJet].constituents().size();
      const unsigned int jetRecID = iJet;
      const double       jetPhi   = m_recoJets[iJet].phi_std();
      const double       jetEta   = m_recoJets[iJet].pseudorapidity();
      const double       jetArea  = 0.;  // FIXME: jet area needs to be defined
      const double       jetE     = m_recoJets[iJet].E();
      const double       jetPt    = m_recoJets[iJet].perp();
      const double       jetPx    = m_recoJets[iJet].px();
      const double       jetPy    = m_recoJets[iJet].py();
      const double       jetPz    = m_recoJets[iJet].pz();
      const double       jetP     = sqrt((jetPx * jetPx) + (jetPy * jetPy) + (jetPz * jetPz));

      // clear constituent vectors
      vecRecCstMatchID.clear();
      vecRecCstZ.clear();
      vecRecCstDr.clear();
      vecRecCstE.clear();
      vecRecCstPt.clear();
      vecRecCstEta.clear();
      vecRecCstPhi.clear();

      // loop over constituents
      vector<fastjet::PseudoJet> recoCsts = m_recoJets[iJet].constituents();
      for (unsigned int iCst = 0; iCst < recoCsts.size(); ++iCst) {

        // get constituent info
        const double cstMatchID = recoCsts[iCst].user_index();
        const double cstPhi     = recoCsts[iCst].phi_std();
        const double cstEta     = recoCsts[iCst].pseudorapidity();
        const double cstE       = recoCsts[iCst].E();
        const double cstPt      = recoCsts[iCst].perp();
        const double cstPx      = recoCsts[iCst].px();
        const double cstPy      = recoCsts[iCst].py();
        const double cstPz      = recoCsts[iCst].pz();
        const double cstP       = ((cstPx * cstPx) + (cstPy * cstPy) + (cstPz * cstPz));
        const double cstZ       = cstP / jetP;
        const double cstDf      = cstPhi - jetPhi;
        const double cstDh      = cstEta - jetEta;
        const double cstDr      = sqrt((cstDf * cstDf) + (cstDh * cstDh));

        // add csts to vectors
        vecRecCstMatchID.push_back(cstMatchID);
        vecRecCstZ.push_back(cstZ);
        vecRecCstDr.push_back(cstDr);
        vecRecCstE.push_back(cstE);
        vecRecCstPt.push_back(cstPt);
        vecRecCstEta.push_back(cstEta);
        vecRecCstPhi.push_back(cstPhi);

        // fill QA histograms and increment counters
        m_hObjectQA[OBJECT::RCST][INFO::PT]  -> Fill(cstPt);
        m_hObjectQA[OBJECT::RCST][INFO::ETA] -> Fill(cstEta);
        m_hObjectQA[OBJECT::RCST][INFO::PHI] -> Fill(cstPhi);
        m_hObjectQA[OBJECT::RCST][INFO::ENE] -> Fill(cstE);
        ++nRecCst;
      }  // end constituent loop

      // store jet/cst output
      m_recoJetNCst.push_back(jetNCst);
      m_recoJetID.push_back(jetRecID);
      m_recoJetE.push_back(jetE);
      m_recoJetPt.push_back(jetPt);
      m_recoJetEta.push_back(jetEta);
      m_recoJetPhi.push_back(jetPhi);
      m_recoJetArea.push_back(jetArea);
      m_recoCstMatchID.push_back(vecRecCstMatchID);
      m_recoCstZ.push_back(vecRecCstZ);
      m_recoCstDr.push_back(vecRecCstDr);
      m_recoCstE.push_back(vecRecCstE);
      m_recoCstPt.push_back(vecRecCstPt);
      m_recoCstEta.push_back(vecRecCstEta);
      m_recoCstPhi.push_back(vecRecCstPhi);

      // fill QA histograms and increment counters
      m_hJetArea[1]                        -> Fill(jetArea);
      m_hJetNumCst[1]                      -> Fill(jetNCst);
      m_hObjectQA[OBJECT::RJET][INFO::PT]  -> Fill(jetPt);
      m_hObjectQA[OBJECT::RJET][INFO::ETA] -> Fill(jetEta);
      m_hObjectQA[OBJECT::RJET][INFO::PHI] -> Fill(jetPhi);
      m_hObjectQA[OBJECT::RJET][INFO::ENE] -> Fill(jetE);
      ++nRecJet;
    }  // end jet loop

    // fill QA histograms
    m_hNumObject[OBJECT::RJET] -> Fill(nRecJet);
    m_hNumObject[OBJECT::RCST] -> Fill(nRecCst);

    // store event info
    m_recoNumJets = nRecJet;
    m_recoVtxX    = m_recoVtx.x();
    m_recoVtxY    = m_recoVtx.y();
    m_recoVtxZ    = m_recoVtx.z();

    // fill object tree
    m_recoTree -> Fill();
    return;

  }  // end 'FillRecoTree()'



  void SCorrelatorJetTreeMaker::SaveOutput() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::SaveOutput() Saving output trees and histograms..." << endl;
    }

    // save QA histograms if need be
    const string sQuality[CONST::NDirectory + 1] = {"Tracks", "CaloClusters", "ParticleFlow", "Particles", "TruthJets", "RecoJets", "QA"};
    TDirectory*  dQuality[CONST::NDirectory + 1];
    if (m_doQualityPlots) {

      // create QA directories
      dQuality[CONST::NDirectory] = (TDirectory*) m_outFile -> mkdir(sQuality[CONST::NDirectory].data());
      for (size_t iDirect = 0; iDirect < CONST::NDirectory; iDirect++) {
        dQuality[iDirect] = (TDirectory*) dQuality[CONST::NDirectory] -> mkdir(sQuality[iDirect].data());
      }

      // save object-specific QA hists
      for (size_t iObj = OBJECT::TRACK; iObj < CONST::NObjType; iObj++) {
        switch (iObj) {
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
          default:
            /* do nothing */
            break;
        }
        m_hNumObject[iObj] -> Write();
        for (size_t iInfo = INFO::PT; iInfo < CONST::NInfoQA; iInfo++) {
          m_hObjectQA[iObj][iInfo] -> Write();
        }
      }  // end object loop

      // save cst-specific histograms
      for (size_t iCst = CST_TYPE::TRACK_CST; iCst < CONST::NCstType; iCst++) {
        switch (iCst) {
          case CST_TYPE::ECAL_CST:
            dQuality[1] -> cd();
            break;
          case CST_TYPE::HCAL_CST:
            dQuality[1] -> cd();
            break;
          case CST_TYPE::FLOW_CST:
            dQuality[2] -> cd();
            break;
          case CST_TYPE::PART_CST:
            dQuality[3] -> cd();
            break;
          default:
            /* do nothing */
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

    // save QA tuples
    dQuality[0] -> cd();
    m_ntTrkQA   -> Write();

    // save output trees
    m_outFile  -> cd();
    m_recoTree -> Write();
    if (m_isMC) {
      m_trueTree -> Write();
    }
    return;

  }  // end 'SaveOutput()'



  void SCorrelatorJetTreeMaker::ResetVariables() {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::ResetTreeVariables() Resetting tree variables..." << endl;
    }

    // reset fastjet members
    m_trueJetDef = NULL;
    m_recoJetDef = NULL;
    m_trueClust  = NULL;
    m_recoClust  = NULL;

    // reset parton and other variables
    m_partonID[0]  = -9999;
    m_partonID[1]  = -9999;
    m_partonMom[0] = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_partonMom[1] = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_vecEvtsToGrab.clear();
    m_mapCstToEmbedID.clear();

    // reset truth (inclusive) tree variables
    m_trueVtx           = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_trueNumJets       = 0;
    m_trueNumChrgPars   = -9999;
    m_trueSumPar        = -9999.;
    m_truePartonID[0]   = -9999;
    m_truePartonID[1]   = -9999;
    m_truePartonMomX[0] = -9999.;
    m_truePartonMomX[1] = -9999.;
    m_truePartonMomY[0] = -9999.;
    m_truePartonMomY[1] = -9999.;
    m_truePartonMomZ[0] = -9999.;
    m_truePartonMomZ[1] = -9999.;
    m_trueJets.clear();
    m_trueJetNCst.clear();
    m_trueJetID.clear();
    m_trueJetE.clear();
    m_trueJetPt.clear();
    m_trueJetEta.clear();    
    m_trueJetPhi.clear();
    m_trueJetArea.clear();
    m_trueCstID.clear();
    m_trueCstEmbedID.clear();
    m_trueCstZ.clear();
    m_trueCstDr.clear();
    m_trueCstE.clear();
    m_trueCstPt.clear();
    m_trueCstEta.clear();
    m_trueCstPhi.clear();

    // reset reco tree variables
    m_recoVtx     = CLHEP::Hep3Vector(-9999., -9999., -9999.);
    m_recoNumJets = 0;
    m_recoNumTrks = -9999;
    m_recoSumECal = -9999.;
    m_recoSumHCal = -9999.;
    m_recoJets.clear();
    m_recoJetNCst.clear();
    m_recoJetID.clear();
    m_recoJetE.clear();
    m_recoJetPt.clear();
    m_recoJetEta.clear();
    m_recoJetPhi.clear();
    m_recoJetArea.clear();
    m_recoCstZ.clear();
    m_recoCstDr.clear();
    m_recoCstE.clear();
    m_recoCstPt.clear();
    m_recoCstEta.clear();
    m_recoCstPhi.clear();
    return;

  }  // end 'ResetTreeVariables()



  void SCorrelatorJetTreeMaker::DetermineEvtsToGrab(PHCompositeNode* topNode) {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::DetermineEvtsToGrab() Determining which subevents to grab..." << endl;
    }

    // make sure vector is clear
    m_vecEvtsToGrab.clear();

    // if not embedding, grab signal event
    // otherwise add all subevents
    if (!m_isEmbed) {
      m_vecEvtsToGrab.push_back(1);
    } else {
      PHHepMCGenEventMap* mapMcEvts = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      for (PHHepMCGenEventMap::ConstIter itEvt = mapMcEvts -> begin(); itEvt != mapMcEvts -> end(); ++itEvt) {
        m_vecEvtsToGrab.push_back(itEvt -> second -> get_embedding_id());
      }
    }
    return;

  }  // end 'DetermineEvtsToGrab()'



  int SCorrelatorJetTreeMaker::CreateJetNode(PHCompositeNode* topNode) {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::CreateJetNode(PHCompositeNode*) Creating jet node..." << endl;
    }

    // create iterator & DST node
    PHNodeIterator   iter(topNode);
    PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!lowerNode) {
      lowerNode = new PHCompositeNode("DST");
      topNode   -> addNode(lowerNode);
      cout << "DST node added" << endl;
    }

    // construct jet tree name
    string baseName;
    string recoNodeName;
    string trueNodeName;
    if (m_jetTreeName.empty()) {
      baseName = "JetTree";
    } else {
      baseName = m_jetTreeName;
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
    recoNodeName = baseName + "_RecoJets";
    trueNodeName = baseName + "_TruthJets";

    // construct jet maps
    m_recoJetMap = new JetMapv1();
    if  (m_isMC && m_saveDST) {
     m_trueJetMap = new JetMapv1();
    }

    // add jet node
    if (m_saveDST) {
      PHIODataNode<PHObject>* recoJetNode = new PHIODataNode<PHObject>(m_recoJetMap, recoNodeName.c_str(), "PHObject");
      lowerNode -> addNode(recoJetNode);
      cout << recoNodeName << " node added" << endl;
    }

    // save truth DSTs if needed
    if(m_isMC && m_saveDST) {
      PHIODataNode<PHObject> *trueJetNode = new PHIODataNode<PHObject>(m_trueJetMap, trueNodeName.c_str(), "PHObject");
      lowerNode -> addNode(trueJetNode);
      cout << trueNodeName << " node added" << endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'CreateJetNode(PHCompositeNode*)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
