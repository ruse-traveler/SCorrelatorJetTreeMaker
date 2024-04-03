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

  void SCorrelatorJetTreeMaker::InitTrees() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::InitTrees() Initializing output trees..." << endl;
    }

    // initialize reco tree
    m_recoTree = new TTree("RecoJetTree",  "A tree of reconstructed jets");
    if (m_config.isLegacy) {
      m_recoLegacy.SetTreeAddresses(m_recoTree);
    } else {
      m_recoOutput.SetTreeAddresses(m_recoTree);
    }

    // initialize truth tree
    if (m_config.isSimulation) {
      m_trueTree = new TTree("TruthJetTree", "A tree of truth jets");
      if (m_config.isLegacy) {
        m_trueLegacy.SetTreeAddresses(m_trueTree);
      } else {
        m_trueOutput.SetTreeAddresses(m_trueTree);
      }
    }
    return;

  }  // end 'InitTrees()'



  void SCorrelatorJetTreeMaker::InitEvals(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
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



  void SCorrelatorJetTreeMaker::FillTrees() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::FillTrueTree() Filling jet trees..." << endl;
    }

    m_recoTree -> Fill();
    if (m_config.isSimulation) {
      m_trueTree -> Fill();
    }
    return;

  }  // end'FillTrees()'



  void SCorrelatorJetTreeMaker::SaveOutput() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::SaveOutput() Saving output trees and histograms..." << endl;
    }

    // save output trees
    m_outFile  -> cd();
    m_recoTree -> Write();
    if (m_config.isSimulation) {
      m_trueTree -> Write();
    }
    return;

  }  // end 'SaveOutput()'



  void SCorrelatorJetTreeMaker::ResetSysVariables() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 2)) {
      cout << "SCorrelatorJetTreeMaker::ResetSysVariables() Resetting system variables..." << endl;
    }

    m_vecEvtsToGrab.clear();
    m_mapCstToEmbedID.clear();
    return;

  }  // end 'ResetSysVariables()'



  void SCorrelatorJetTreeMaker::ResetJetVariables() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 2)) {
      cout << "SCorrelatorJetTreeMaker::ResetJetVariables() Resetting jet variables..." << endl;
    }

    m_recoJets.clear();
    m_trueJets.clear();
    m_recoJetInput.clear();
    m_trueJetInput.clear();
    m_recoSourceMap.clear();
    m_trueSourceMap.clear();
    return;

  }  // end 'ResetJetVariables()'



  int SCorrelatorJetTreeMaker::CreateJetNode(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::CreateJetNode(PHCompositeNode*, string) Creating jet node..." << endl;
    }

    // construct jet tree name
    string baseName;
    if (m_config.jetTreeName.empty()) {
      baseName = "JetTree";
    } else {
      baseName = m_config.jetTreeName;
    }

    // construct jet node name
    recoNodeName = baseName + "_RecoJets";
    trueNodeName = baseName + "_TruthJets";

    // construct jet maps
    //   - FIXME I don't think the jet maps are actually being
    //     filled? And I want my jet trees in there anyways...
    m_recoJetMap = new JetMapv1();
    if (m_config.isSimulation) {
      m_trueJetMap = new JetMapv1();
    }

    // create nodes
    Interfaces::CreateNode(topNode, recoNodeName, m_recoJetMap);
    if (m_config.isSimulation) {
      Interfaces::CreateNode(topNode, trueNodeName, m_trueJetMap);
    }
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'CreateJetNode(PHCompositeNode*)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
