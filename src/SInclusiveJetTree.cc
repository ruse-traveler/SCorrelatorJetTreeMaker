// 'SInclusiveJetTree.cc'
// Derek Anderson
// 12.04.202
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)

#define SINCLUSIVEJETTREE_CC

// user include
#include "SInclusiveJetTree.h"
#include "SInclusiveJetTree.io.h"
// f4a includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
// phool includes
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
// tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
// calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>
// g4 includes
#include <g4jets/Jet.h>
#include <g4jets/Jetv1.h>
#include <g4jets/JetMap.h>
#include <g4jets/FastJetAlgo.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
// fastjet includes
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
// particle flow includes
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>
// hepmc includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
// root includes
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
// standard c includes
#include <string>
#include <cassert>
#include <sstream>

using namespace std;
using namespace fastjet;

// global constants
static const bool Debug(true);



// ctor/dtor ------------------------------------------------------------------

SInclusiveJetTree::SInclusiveJetTree(const string &name, const string &outfile) : SubsysReco(name) {
  m_outfilename = outfile;
  initializeVariables();
  initializeTrees();
  if (Debug) {
    std::cout << "SInclusiveJetTree::SInclusiveJetTree(const std::string &name) Calling ctor" << std::endl;
  }
}  // end ctor(string, string)



SInclusiveJetTree::~SInclusiveJetTree() {
  if (Debug) {
    std::cout << "SInclusiveJetTree::~SInclusiveJetTree() Calling dtor" << std::endl;
  }
}  // end dtor



// F4A methods ----------------------------------------------------------------

int SInclusiveJetTree::Init(PHCompositeNode *topNode) {
  std::cout << "SInclusiveJetTree::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}



int SInclusiveJetTree::process_event(PHCompositeNode *topNode) {
  std::cout << "SInclusiveJetTree::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}



int SInclusiveJetTree::End(PHCompositeNode *topNode) {
  std::cout << "SInclusiveJetTree::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}



// i/o methods ----------------------------------------------------------------

void SInclusiveJetTree::initializeVariables() {
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
  m_add_particleflow     = true;
  m_add_tracks           = false;
  m_add_EMCal_clusters   = false;
  m_add_HCal_clusters    = false;
  m_jetr                 = 0.4;
  m_jetalgo              = fastjet::antikt_algorithm;
  m_recomb_scheme        = fastjet::pt_scheme;
  m_qualy_plots          = false;
  m_save_dst             = false;
  m_ismc                 = false;
  m_numCst.clear();
  m_jetId.clear();
  m_jetTruId.clear();
  m_jetPt.clear();
  m_jetEta.clear();
  m_jetPhi.clear();
  m_cstZ.clear();
  m_cstDr.clear();
  m_cstJt.clear();
  m_cstEta.clear();
  m_cstPhi.clear();
  m_outFile = new TFile();
  return;
}  // end 'initializeVariables()'



void SInclusiveJetTree::initializeTrees() {
  m_jetTree = new TTree("JetTree", "A tree of jets");
  m_jetTree -> Branch("EvtNumJets",   &m_numJets,       "NumJets/D");
  m_jetTree -> Branch("Parton3_ID",   &m_partonID[0],   "Parton3_ID/I");
  m_jetTree -> Branch("Parton4_ID",   &m_partonID[1],   "Parton4_ID/I");
  m_jetTree -> Branch("Parton3_MomX", &m_partonMomX[0], "Parton3_MomX/D");
  m_jetTree -> Branch("Parton3_MomY", &m_partonMomY[0], "Parton3_MomY/D");
  m_jetTree -> Branch("Parton3_MomZ", &m_partonMomZ[0], "Parton3_MomZ/D");
  m_jetTree -> Branch("Parton4_MomX", &m_partonMomX[1], "Parton4_MomX/D");
  m_jetTree -> Branch("Parton4_MomY", &m_partonMomY[1], "Parton4_MomY/D");
  m_jetTree -> Branch("Parton4_MomZ", &m_partonMomZ[1], "Parton4_MomZ/D");
  m_jetTree -> Branch("JetNumCst",    &m_numCst);
  m_jetTree -> Branch("JetID",        &m_jetId);
  m_jetTree -> Branch("JetTruthID",   &m_jetTruId);
  m_jetTree -> Branch("JetPt",        &m_jetPt);
  m_jetTree -> Branch("JetEta",       &m_jetEta);
  m_jetTree -> Branch("JetPhi",       &m_jetPhi);
  m_jetTree -> Branch("CstZ",         &m_cstZ);
  m_jetTree -> Branch("CstDr",        &m_cstDr);
  m_jetTree -> Branch("CstJt",        &m_cstJt);
  m_jetTree -> Branch("CstEta",       &m_cstEta);
  m_jetTree -> Branch("CstPHi",       &m_cstPhi);
  return;
}  // end 'initializeTrees()'

// end ------------------------------------------------------------------------
