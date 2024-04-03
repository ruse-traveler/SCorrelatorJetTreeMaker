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
#include "SCorrelatorJetTreeMaker.sys.h"

using namespace std;
using namespace fastjet;
using namespace findNode;



namespace SColdQcdCorrelatorAnalysis {

  // ctor/dtor ----------------------------------------------------------------

  SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(const string& name, const bool debug) : SubsysReco(name) {

    if (debug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(string&, bool) Calling ctor" << endl;
    }

    // clear vector/map members
    ResetSysVariables();
    ResetJetVariables();

  }  // end ctor(string&, bool)



  SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(SCorrelatorJetTreeMaker& config) : SubsysReco(config.moduleName) {

    m_config = config;
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::SCorrelatorJetTreeMaker(SCorrelatorJetTreeMakerConfig&) Calling ctor" << endl;
    }

    // clear vector/map members
    ResetSysVariables();
    ResetJetVariables();

  }  // end ctor(SCorrelatorJetTreeMakerConfig&)



  SCorrelatorJetTreeMaker::~SCorrelatorJetTreeMaker() {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::~SCorrelatorJetTreeMaker() Calling dtor" << endl;
    }

    // clean up dangling pointers
    //   - FIXME use smart pointers instead
    if (m_evalStack) {
      delete m_evalStack;
      m_evalStack = NULL;
      m_trackEval = NULL;
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
    if (m_config.saveDST) {
      CreateJetNode(topNode);
    }

    // initialize output trees
    InitTrees();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'Init(PHcompositeNode*)'



  int SCorrelatorJetTreeMaker::process_event(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::process_event(PHCompositeNode*) Processing Event..." << endl;
    }

    // initialize evaluator & determine subevts to grab for event
    if (m_config.isSimulation) {
      InitEvals(topNode);
      DetermineEvtsToGrab(topNode);
    }

    // get event-wise variables
    GetEventVariables(topNode);

    // if needed, check reconstructed vtx
/* FIXME finish!
    if (m_config.doVtxCut && !IsGoodVertex()) {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
*/

    // find jets
    MakeRecoJets(topNode);
    if (m_config.isMC) {
      MakeTrueJets(topNode);
    }

    /* TODO add jet collector function */

    // fill output trees
    FillTrees();
    return eventStatus;

  }  // end 'process_event(PHCompositeNode*)'



  int SCorrelatorJetTreeMaker::End(PHCompositeNode* topNode) {

    // print debug statements
    if (m_config.isDebugOn) {
      cout << "SCorrelatorJetTreeMaker::End(PHCompositeNode*) This is the End..." << endl;
    }

    SaveOutput();

    // close file and exit
    m_outFile -> cd();
    m_outFile -> Close();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'End(PHcompositeNode*)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
