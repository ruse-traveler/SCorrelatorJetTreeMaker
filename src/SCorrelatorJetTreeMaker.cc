// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.cc'
// Derek Anderson
// 12.04.2022
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Initially derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#define SCORRELATORJETTREE_CC

// user includes
#include "SCorrelatorJetTreeMaker.h"
#include "SCorrelatorJetTreeMaker.io.h"
#include "SCorrelatorJetTreeMaker.ana.h"
#include "SCorrelatorJetTreeMaker.evt.h"
#include "SCorrelatorJetTreeMaker.sys.h"

using namespace std;
using namespace fastjet;
using namespace findNode;



namespace SColdQcdCorrelatorAnalysis {

  // ctor/dtor ----------------------------------------------------------------

  SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(const string& name, const bool debug) : SubsysReco(name) {

    if (debug) {
      cout << "SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(string&, bool) Calling ctor" << endl;
    }
    InitVariables();  // TODO remove

  }  // end ctor(string&, bool)



  SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(SCorrelatorJetTreeMaker& config) : SubsysReco(config.moduleName) {

    m_config = config;
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(SCorrelatorJetTreeMakerConfig&) Calling ctor" << endl;
    }
    InitVariables();  // TODO remove

  }  // end ctor(SCorrelatorJetTreeMakerConfig&)


  SCorrelatorJetTreeMaker::~SCorrelatorJetTreeMaker() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::~SCorrelatorJetTreeMaker() Calling dtor" << endl;
    }

    // clean up dangling pointers
    //   - FIXME use smart pointers instead
    if (m_histMan) {
      delete m_histMan;
      m_histMan = NULL;
    }
    if (m_evalStack) {
      delete m_evalStack;
      m_evalStack = NULL;
      m_trackEval = NULL;
    }
    if (m_trueJetDef) {
      delete m_trueJetDef;
      m_trueJetDef = NULL;
    }
    if (m_recoJetDef) {
      delete m_recoJetDef;
      m_recoJetDef = NULL;
    }
    if (m_trueClust) {
      delete m_trueClust;
      m_trueClust = NULL;
    }
    if (m_recoClust) {
      delete m_recoClust;
      m_recoClust = NULL;
    }

  }  // end dtor



  // F4A methods --------------------------------------------------------------

  int SCorrelatorJetTreeMaker::Init(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::Init(PHCompositeNode*) Initializing..." << endl;
    }

    // intitialize output file
    m_outFile = new TFile(m_config.outFileName.c_str(), "RECREATE");
    if (!m_outFile) {
      cerr << "PANIC: couldn't open SCorrelatorJetTreeMaker output file!" << endl;
    }

    // create node for jet-tree
    if (m_saveDST) {
      CreateJetNode(topNode);
    }

    // initialize QA histograms/tuples, output trees, and functions
    InitHists();
    InitTrees();
    InitFuncs();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'Init(PHcompositeNode*)'



  int SCorrelatorJetTreeMaker::process_event(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::process_event(PHCompositeNode*) Processing Event..." << endl;
    }

    // reset event-wise variables & members
    ResetVariables();

    // initialize evaluator & determine subevts to grab for event
    if (m_config.isMC) {
      InitEvals(topNode);
      DetermineEvtsToGrab(topNode);
    }

    // get event-wise variables
    GetEventVariables(topNode);
    if (m_config.isMC) {
      GetPartonInfo(topNode);
    }

    // check if reconstructed vertex is in in acceptance
    bool isGoodEvt = true;
    if (m_config.doVtxCut) {
      isGoodEvt = IsGoodVertex(m_recoVtx);
    }

    // set event status
    int eventStatus = Fun4AllReturnCodes::EVENT_OK;
    if (m_config.doVtxCut && !isGoodEvt) {
      eventStatus = Fun4AllReturnCodes::DISCARDEVENT;
    } else {
      eventStatus = Fun4AllReturnCodes::EVENT_OK;
    }

    // if event is good, continue processing
    if (isGoodEvt) {

      // find jets
      MakeRecoJets(topNode);
      if (m_config.isMC) {
        MakeTrueJets(topNode);
      }

      // fill output trees
      FillRecoTree();
      if (m_config.isMC) {
        FillTrueTree();
      }
    }
    return eventStatus;

  }  // end 'process_event(PHCompositeNode*)'



  int SCorrelatorJetTreeMaker::End(PHCompositeNode* topNode) {

    // print debug statements
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::End(PHCompositeNode*) This is the End..." << endl;
    }

    // save output and close
    SaveOutput();
    m_outFile -> cd();
    m_outFile -> Close();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'End(PHcompositeNode*)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
