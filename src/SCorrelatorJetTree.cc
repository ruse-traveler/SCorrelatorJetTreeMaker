// 'SCorrelatorJetTree.cc'
// Derek Anderson
// 12.04.202
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)

#define SCORRELATORJETTREE_CC

// user include
#include "SCorrelatorJetTree.h"
#include "SCorrelatorJetTree.io.h"
// f4a includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool includes
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
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
using namespace findNode;

// global constants
static const unsigned long NRange(2);



// ctor/dtor ------------------------------------------------------------------

SCorrelatorJetTree::SCorrelatorJetTree(const string &name, const string &outfile, const bool debug) : SubsysReco(name) {

  // print debug statement
  m_doDebug = debug;
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::SCorrelatorJetTree(const string &name) Calling ctor" << endl;
  }
  m_outfilename = outfile;
  initializeVariables();
  initializeTrees();

}  // end ctor(string, string)



SCorrelatorJetTree::~SCorrelatorJetTree() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::~SCorrelatorJetTree() Calling dtor" << endl;
  }
  delete m_hm;
  delete m_outFile;
  delete m_jetTree;

}  // end dtor



// F4A methods ----------------------------------------------------------------

int SCorrelatorJetTree::Init(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug || (Verbosity() > 5)) {
    cout << "SCorrelatorJetTree::Init(PHCompositeNode *topNode) Initializing..." << endl;
  }

  // intitialize output file
  m_outFile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outFile) {
    cerr << "PANIC: couldn't open SCorrelatorJetTree output file!" << endl;
  }

  // create node for jet-tree
  if (m_save_dst) {
    createJetNode(topNode);
  }

  // create QA histograms
  /* TODO: QA histograms will go here */
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'Init(PHcompositeNode*)'



int SCorrelatorJetTree::process_event(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug || (Verbosity() > 5)) {
    cout << "SCorrelatorJetTree::process_event(PHCompositeNode *topNode) Processing Event" << endl;
  }

  // find jets
  findJets(topNode);
  if (m_ismc) {
    findMcJets(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'process_event(PHCompositeNode*)'



int SCorrelatorJetTree::End(PHCompositeNode *topNode) {

  // print debug statements
  if (m_doDebug || (Verbosity() > 1)) {
    cout << "SCorrelatorJetTree::End(PHCompositeNode *topNode) This is the End..." << endl;
  }

  // save output and close
  m_outFile -> cd();
  m_jetTree -> Write();
  m_outFile -> Write();
  m_outFile -> Close();
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'End(PHcompositeNode*)'



// jet methods ----------------------------------------------------------------

void SCorrelatorJetTree::findJets(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::findJets(PHCompositeNode *topNode) Finding jets..." << endl;
  }

  // declare fastjet objects
  fastjet::JetDefinition     *jetdef = new fastjet::JetDefinition(m_jetalgo, m_jetr, m_recomb_scheme, fastjet::Best);
  vector<fastjet::PseudoJet>  particles;

  // instantiate fastjet map
  map<int, pair<Jet::SRC, int>> fjMap;

  // add constitutents
  const bool doParticleFlow = m_add_particleflow;
  const bool doTracks       = m_add_tracks;
  const bool doCaloClusters = (m_add_EMCal_clusters || m_add_HCal_clusters);
  if (doParticleFlow) addParticleFlow(topNode, particles, fjMap);
  if (doTracks)       addTracks(topNode, particles, fjMap);
  if (doCaloClusters) addClusters(topNode, particles, fjMap);

  // cluster jets
  fastjet::ClusterSequence   jetFinder(particles, *jetdef);
  vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  // prepare vectors for filling
  m_jetNCst.clear();
  m_jetId.clear();
  m_jetTruId.clear();
  m_jetE.clear();
  m_jetPt.clear();
  m_jetEta.clear();
  m_jetPhi.clear();
  m_cstZ.clear();
  m_cstDr.clear();
  m_cstE.clear();
  m_cstJt.clear();
  m_cstEta.clear();
  m_cstPhi.clear();

  // declare vectors for storing constituents
  vector<double> vecCstZ;
  vector<double> vecCstDr;
  vector<double> vecCstE;
  vector<double> vecCstJt;
  vector<double> vecCstEta;
  vector<double> vecCstPhi;
  vecCstZ.clear();
  vecCstDr.clear();
  vecCstE.clear();
  vecCstJt.clear();
  vecCstEta.clear();
  vecCstPhi.clear();

  // fill jet/constituent variables
  unsigned int nJet(0);
  for (unsigned int iJet = 0; iJet < fastjets.size(); ++iJet) {

    // get jet info
    const unsigned int jetNCst  = fastjets[iJet].constituents().size();
    const unsigned int jetID    = iJet;
    const unsigned int jetTruID = 99999;  // FIXME: this will need to be changed to the matched truth jet
    const double       jetPhi   = fastjets[iJet].phi_std();
    const double       jetEta   = fastjets[iJet].pseudorapidity();
    const double       jetE     = fastjets[iJet].E();
    const double       jetPt    = fastjets[iJet].perp();
    const double       jetPx    = fastjets[iJet].px();
    const double       jetPy    = fastjets[iJet].py();
    const double       jetPz    = fastjets[iJet].pz();
    const double       jetP     = sqrt((jetPx * jetPx) + (jetPy * jetPy) + (jetPz * jetPz));

    // clear constituent vectors
    vecCstZ.clear();
    vecCstDr.clear();
    vecCstE.clear();
    vecCstJt.clear();
    vecCstEta.clear();
    vecCstPhi.clear();

    // loop over constituents
    vector<fastjet::PseudoJet> csts = fastjets[iJet].constituents();
    for (unsigned int iCst = 0; iCst < csts.size(); ++iCst) {

      // get constituent info
      const double cstPhi = csts[iCst].phi_std();
      const double cstEta = csts[iCst].pseudorapidity();
      const double cstE   = csts[iCst].E();
      const double cstJt  = csts[iCst].perp();
      const double cstJx  = csts[iCst].px();
      const double cstJy  = csts[iCst].py();
      const double cstJz  = csts[iCst].pz();
      const double cstJ   = ((cstJx * cstJx) + (cstJy * cstJy) + (cstJz * cstJz));
      const double cstZ   = cstJ / jetP;
      const double cstDf  = cstPhi - jetPhi;
      const double cstDh  = cstEta - jetEta;
      const double cstDr  = sqrt((cstDf * cstDf) + (cstDh * cstDh));

      // add csts to vectors
      vecCstZ.push_back(cstZ);
      vecCstDr.push_back(cstDr);
      vecCstE.push_back(cstE);
      vecCstJt.push_back(cstJt);
      vecCstEta.push_back(cstEta);
      vecCstPhi.push_back(cstPhi);
    }  // end constituent loop

    // store jet/cst output
    m_jetNCst.push_back(jetNCst);
    m_jetId.push_back(jetID);
    m_jetTruId.push_back(jetTruID);
    m_jetE.push_back(jetE);
    m_jetPt.push_back(jetPt);
    m_jetEta.push_back(jetEta);
    m_jetPhi.push_back(jetPhi);
    m_cstZ.push_back(vecCstZ);
    m_cstDr.push_back(vecCstDr);
    m_cstE.push_back(vecCstE);
    m_cstJt.push_back(vecCstJt);
    m_cstEta.push_back(vecCstEta);
    m_cstPhi.push_back(vecCstPhi);
    ++nJet;
  }  // end jet loop

  // store evt info
  m_numJets       = nJet;
  m_partonID[0]   = -9999;   // FIXME: store actual value
  m_partonID[1]   = -9999;   // FIXME: store actual value
  m_partonMomX[0] = -9999.;  // FIXME: store actual value
  m_partonMomX[1] = -9999.;  // FIXME: store actual value
  m_partonMomY[0] = -9999.;  // FIXME: store actual value
  m_partonMomY[1] = -9999.;  // FIXME: store actual value
  m_partonMomZ[0] = -9999.;  // FIXME: store actual value
  m_partonMomZ[1] = -9999.;  // FIXME: store actual value
  return;

}  // end 'findJets(PHCompositeNode*)'



void SCorrelatorJetTree::findMcJets(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::findMcJets(PHCompositeNode *topNode) Finding MC jets..." << endl;
  }
  return;

}  // end 'findMcJets(PHCompositeNode*)'



void SCorrelatorJetTree::addParticleFlow(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::addParticleFlow(PHCompositeNode *topNode, vector<PseudoJet>&, map<int, parir<Jet::SRC, int>>&) Adding particle flow elements" << endl;
  }

  // declare pf  objects
  ParticleFlowElementContainer *pflowContainer = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");

  // loop over pf elements
  unsigned int                                iPart     = particles.size();
  ParticleFlowElementContainer::ConstRange    begin_end = pflowContainer -> getParticleFlowElements();
  ParticleFlowElementContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) {

    // get pf element
    ParticleFlowElement *pflow = rtiter -> second;

    // check if good
    const bool isGoodElement = isAcceptableParticleFlow(pflow);
    if(!pflow || !isGoodElement) continue;

    // create pseudojet and add to constituent vector
    const int    pfID = pflow -> get_id();
    const double pfPx = pflow -> get_px();
    const double pfPy = pflow -> get_py();
    const double pfPz = pflow -> get_pz();
    const double pfE  = pflow -> get_e();

    fastjet::PseudoJet fjPartFlow(pfPx, pfPy, pfPz, pfE);
    fjPartFlow.set_user_index(iPart);
    particles.push_back(fjPartFlow);

    // add pf element to fastjet map
    pair<int, pair<Jet::SRC, int>> jetPartFlowPair(iPart, make_pair(Jet::SRC::PARTICLE, pfID));
    fjMap.insert(jetPartFlowPair);
    ++iPart;
  }  // end pf element loop
  return;

}  // end 'addParticleFlow(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



void SCorrelatorJetTree::addTracks(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::addTracks(PHCompositeNode *topNode, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding tracks" << endl;
  }

  // get track map
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE
         << "PANIC: SvtxTrackMap node is missing, can't collect tracks!"
         << endl;
    return;
  }

  // loop over tracks
  unsigned int  iPart = particles.size();
  SvtxTrack    *track = 0;
  for (SvtxTrackMap::Iter iter = trackmap -> begin(); iter != trackmap -> end(); ++iter) {

    // get track
    track = iter -> second;

    // check if good
    const bool isGoodTrack = isAcceptableTrack(track);
    if(!isGoodTrack) continue;

    // create pseudojet and add to constituent vector
    const int    trkID = track -> get_id();
    const double trkPx = track -> get_px();
    const double trkPy = track -> get_py();
    const double trkPz = track -> get_pz();

    fastjet::PseudoJet fjTrack(trkPx, trkPy, trkPz, 0.);  // FIXME: add track energy? (e.g. maybe assume a pion mass...)
    fjTrack.set_user_index(iPart);
    particles.push_back(fjTrack);

    // add track to fastjet map
    pair<int, pair<Jet::SRC, int>> jetTrkPair(iPart, make_pair(Jet::SRC::TRACK, trkID));
    fjMap.insert(jetTrkPair);
    ++iPart;
  }  // end track loop
  return;

}  // end 'addTracks(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



void SCorrelatorJetTree::addClusters(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::addClusters(PHCompositeNode *topNode, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding clusters" << endl;
  }

  // get vertex map
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap) {
    cerr << "SCorrelatorJetTree::getEmcalClusters - Fatal Error - GlobalVertexMap node is missing!\n"
         << "  Please turn on the do_global flag in the main macro in order to reconstruct the global vertex!"
         << endl;
    assert(vertexmap);
    return;
  }
  if (vertexmap -> empty()) {
    cerr << "JetTagging::getEmcalClusters - Fatal Error - GlobalVertexMap node is empty!\n"
         << "  Please turn on the do_global flag in the main macro in order to reconstruct the global vertex!"
         << endl;
    return;
  }

  // grab vertex
  GlobalVertex *vtx = vertexmap -> begin() -> second;
  if (vtx == nullptr) return;

  // add emcal clusters if needed
  int iPart = particles.size();
  if (m_add_EMCal_clusters) {

    // grab em cluster containter
    RawClusterContainer *clustersEMC = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    if (!clustersEMC) {
      cout << PHWHERE
           << "PANIC: EMCal cluster node is missing, can't collect EMCal clusters!"
           << endl;
      return;
    }

    // loop over em clusters
    RawClusterContainer::ConstRange    begin_end_EMC = clustersEMC->getClusters();
    RawClusterContainer::ConstIterator clusIter_EMC;
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC) {

      // grab cluster
      const RawCluster *cluster = clusIter_EMC -> second;

      // construct vertex and get 4-momentum
      const double vX = vtx -> get_x();
      const double vY = vtx -> get_y();
      const double vZ = vtx -> get_z();

      CLHEP::Hep3Vector vertex(vX, vY, vZ);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      // check if good
      const bool isGoodClust = isAcceptableEMCalCluster(E_vec_cluster);
      if(!isGoodClust) continue;

      // create pseudojet and add to constituent vector
      const int    emClustID  = cluster -> get_id();
      const double emClustE   = E_vec_cluster.mag();
      const double emClustPt  = E_vec_cluster.perp();
      const double emClustPhi = E_vec_cluster.getPhi();
      const double emClustPx  = emClustPt * cos(emClustPhi);
      const double emClustPy  = emClustPt * sin(emClustPhi);
      const double emClustPz  = sqrt((emClustE * emClustE) - (emClustPx * emClustPx) - (emClustPy * emClustPy));

      fastjet::PseudoJet fjCluster(emClustPx, emClustPy, emClustPz, emClustE);
      fjCluster.set_user_index(iPart);
      particles.push_back(fjCluster);

      // add em cluster to fastjet map
      pair<int, pair<Jet::SRC, int>> jetEMClustPair(iPart, make_pair(Jet::SRC::CEMC_CLUSTER, emClustID));
      fjMap.insert(jetEMClustPair);
      ++iPart;
    }  // end em cluster loop
  }  // end if (m_add_EMCal_clusters)

  //  add hcal clusters if needed
  if (m_add_HCal_clusters) {

    // grab ih cluster container
    RawClusterContainer *clustersHCALIN = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
    if (!clustersHCALIN) {
      cerr << PHWHERE
           << "PANIC: Inner HCal cluster node is missing, can't collect inner HCal clusters!"
           << endl;
      return;
    }

    // Loop over ih clusters
    RawClusterContainer::ConstRange    begin_end_HCALIN = clustersHCALIN->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALIN;
    for (clusIter_HCALIN = begin_end_HCALIN.first; clusIter_HCALIN != begin_end_HCALIN.second; ++clusIter_HCALIN) {

      // get ih cluster
      const RawCluster *cluster = clusIter_HCALIN -> second;

      // construct vertex and get 4-momentum
      const double vX = vtx -> get_x();
      const double vY = vtx -> get_y();
      const double vZ = vtx -> get_z();

      CLHEP::Hep3Vector vertex(vX, vY, vZ);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      // check if good
      const bool isGoodClust = isAcceptableHCalCluster(E_vec_cluster);
      if (!isGoodClust) continue;

      // create pseudojet and add to constituent vector
      const int    ihClustID  = cluster -> get_id();
      const double ihClustE   = E_vec_cluster.mag();
      const double ihClustPt  = E_vec_cluster.perp();
      const double ihClustPhi = E_vec_cluster.getPhi();
      const double ihClustPx  = ihClustPt * cos(ihClustPhi);
      const double ihClustPy  = ihClustPt * sin(ihClustPhi);
      const double ihClustPz  = sqrt((ihClustE * ihClustE) - (ihClustPx * ihClustPx) - (ihClustPy * ihClustPy));

      fastjet::PseudoJet fjCluster(ihClustPx, ihClustPy, ihClustPz, ihClustE);
      fjCluster.set_user_index(iPart);
      particles.push_back(fjCluster);

      // add ih cluster to fastjet map
      pair<int, pair<Jet::SRC, int>> jetIHClustPair(iPart, make_pair(Jet::SRC::HCALIN_CLUSTER, ihClustID));
      fjMap.insert(jetIHClustPair);
      ++iPart;
    }  // end ih cluster loop

    // grab oh cluster container
    RawClusterContainer *clustersHCALOUT = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
    if (!clustersHCALOUT) {
      cerr << PHWHERE
           << "PANIC: Outer HCal cluster node is missing, can't collect outer HCal clusters!"
           << endl;
      return;
    }

    // loop over oh clusters
    RawClusterContainer::ConstRange    begin_end_HCALOUT = clustersHCALOUT->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALOUT;
    for (clusIter_HCALOUT = begin_end_HCALOUT.first; clusIter_HCALOUT != begin_end_HCALOUT.second; ++clusIter_HCALOUT) {

      // get oh cluster
      const RawCluster *cluster = clusIter_HCALOUT -> second;

      // construct vertex and get 4-momentum
      const double vX = vtx -> get_x();
      const double vY = vtx -> get_y();
      const double vZ = vtx -> get_z();

      CLHEP::Hep3Vector vertex(vX, vY, vZ);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      // check if good
      const bool isGoodClust = isAcceptableHCalCluster(E_vec_cluster);
      if (!isGoodClust) continue;

      // create pseudojet and add to constituent vector
      const int    ohClustID  = cluster -> get_id();
      const double ohClustE   = E_vec_cluster.mag();
      const double ohClustPt  = E_vec_cluster.perp();
      const double ohClustPhi = E_vec_cluster.getPhi();
      const double ohClustPx  = ohClustPt * cos(ohClustPhi);
      const double ohClustPy  = ohClustPt * sin(ohClustPhi);
      const double ohClustPz  = sqrt((ohClustE * ohClustE) - (ohClustPx * ohClustPx) - (ohClustPy * ohClustPy));

      fastjet::PseudoJet fjCluster(ohClustPx, ohClustPy, ohClustPz, ohClustE);
      fjCluster.set_user_index(iPart);
      particles.push_back(fjCluster);

      // add oh cluster to fastjet map
      pair<int, pair<Jet::SRC, int>> jetOHClustPair(iPart, make_pair(Jet::SRC::HCALOUT_CLUSTER, ohClustID));
      fjMap.insert(jetOHClustPair);
      ++iPart;
    }  // end oh cluster loop
  }  // end if (m_add_HCal_clusters)
  return;

}  // end 'addClusters(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



void SCorrelatorJetTree::getTracks(PHCompositeNode *topNode) {
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::getTracks(PHCompositeNode *topNode) Getting tracks" << endl;
  }
  return;
}  // end 'getTracks(PHCompositeNode*)'



// constituent methods --------------------------------------------------------

bool SCorrelatorJetTree::isAcceptableParticleFlow(ParticleFlowElement* pfPart) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::isAcceptableParticleFlow(ParticleFlowElement*) Checking if particle flow element is good" << endl;
  }

  // TODO: explore particle flow cuts
  const double pfEta         = pfPart -> get_eta();
  const bool   isInEtaRange  = ((pfEta > m_particleflow_mineta) && (pfEta < m_particleflow_maxeta));
  const bool   isGoodElement = isInEtaRange;
  return isGoodElement;

}  // end 'isAcceptableParticleFlow(ParticleFlowElement*)'



bool SCorrelatorJetTree::isAcceptableTrack(SvtxTrack *track) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::isAcceptableTrack(SvtxTrack*) Checking if track is good" << endl;
  }

  const double trkPt        = track -> get_pt();
  const double trkEta       = track -> get_eta();
  const bool   isInPtRange  = ((trkPt > m_track_minpt)   && (trkPt < m_track_maxpt));
  const bool   isInEtaRange = ((trkEta > m_track_mineta) && (trkEta < m_track_maxeta));
  const bool   isGoodTrack  = (isInPtRange && isInEtaRange);
  return isGoodTrack;

}  // end 'isAcceptableTrack(SvtxTrack*)'



bool SCorrelatorJetTree::isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::isAcceptableEMCalCluster(CLHEP::Hep3Vector&) Checking if ECal cluster is good" << endl;
  }

  const double clustPt      = E_vec_cluster.perp();
  const double clustEta     = E_vec_cluster.pseudoRapidity();
  const bool   isInPtRange  = ((clustPt > m_EMCal_cluster_minpt)   && (clustPt < m_EMCal_cluster_maxpt));
  const bool   isInEtaRange = ((clustEta > m_EMCal_cluster_mineta) && (clustEta < m_EMCal_cluster_maxeta));
  const bool   isGoodClust  = (isInPtRange && isInEtaRange);
  return isGoodClust;

}  // end 'isAcceptableEMCalCluster(CLHEP::Hep3Vector&)'



bool SCorrelatorJetTree::isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::isAcceptableHCalCluster(CLHEP::Hep3Vector&) Checking if HCal cluster is good" << endl;
  }

  const double clustPt      = E_vec_cluster.perp();
  const double clustEta     = E_vec_cluster.pseudoRapidity();
  const bool   isInPtRange  = ((clustPt > m_HCal_cluster_minpt)   && (clustPt < m_HCal_cluster_maxpt));
  const bool   isInEtaRange = ((clustEta > m_HCal_cluster_mineta) && (clustEta < m_HCal_cluster_maxeta));
  const bool   isGoodClust  = (isInPtRange && isInEtaRange);
  return isGoodClust;

}  // end 'isAcceptableHCalCluster(CLHEP::Hep3Vector&)'



// i/o methods ----------------------------------------------------------------

void SCorrelatorJetTree::initializeVariables() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::initializeVariables() Initializing class members..." << endl;
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
  m_add_particleflow     = true;
  m_add_tracks           = false;
  m_add_EMCal_clusters   = false;
  m_add_HCal_clusters    = false;
  m_jetr                 = 0.4;
  m_jetalgo              = antikt_algorithm;
  m_recomb_scheme        = pt_scheme;
  m_qualy_plots          = false;
  m_save_dst             = false;
  m_ismc                 = false;
  m_numJets              = 0;
  m_partonID[0]          = -9999;
  m_partonID[1]          = -9999;
  m_partonMomX[0]        = -9999.;
  m_partonMomX[1]        = -9999.;
  m_partonMomY[0]        = -9999.;
  m_partonMomY[1]        = -9999.;
  m_partonMomZ[0]        = -9999.;
  m_partonMomZ[1]        = -9999.;
  m_jetNCst.clear();
  m_jetId.clear();
  m_jetTruId.clear();
  m_jetE.clear();
  m_jetPt.clear();
  m_jetEta.clear();
  m_jetPhi.clear();
  m_cstZ.clear();
  m_cstDr.clear();
  m_cstE.clear();
  m_cstJt.clear();
  m_cstEta.clear();
  m_cstPhi.clear();
  m_outFile = new TFile();
  return;

}  // end 'initializeVariables()'



void SCorrelatorJetTree::initializeTrees() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::initializeTrees() Initializing output tree..." << endl;
  }

  // initialize output tree
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
  m_jetTree -> Branch("JetNumCst",    &m_jetNCst);
  m_jetTree -> Branch("JetID",        &m_jetId);
  m_jetTree -> Branch("JetTruthID",   &m_jetTruId);
  m_jetTree -> Branch("JetEnergy",    &m_jetE);
  m_jetTree -> Branch("JetPt",        &m_jetPt);
  m_jetTree -> Branch("JetEta",       &m_jetEta);
  m_jetTree -> Branch("JetPhi",       &m_jetPhi);
  m_jetTree -> Branch("CstZ",         &m_cstZ);
  m_jetTree -> Branch("CstDr",        &m_cstDr);
  m_jetTree -> Branch("CstEnergy",    &m_cstE);
  m_jetTree -> Branch("CstJt",        &m_cstJt);
  m_jetTree -> Branch("CstEta",       &m_cstEta);
  m_jetTree -> Branch("CstPhi",       &m_cstPhi);
  return;

}  // end 'initializeTrees()'



int SCorrelatorJetTree::createJetNode(PHCompositeNode* topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::createJetNode(PHCompositeNode *topNode) Creating jet node..." << endl;
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

  // cant have forward slashes in DST or else you make a subdirectory on save!!!
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

}  // end 'createJetNode(PHCompositeNode*)'



void SCorrelatorJetTree::resetTreeVariables() {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::resetTreeVariables() Resetting tree variables..." << endl;
  }
  m_numJets       = 0;
  m_partonID[0]   = -9999;
  m_partonID[1]   = -9999;
  m_partonMomX[0] = -9999.;
  m_partonMomX[1] = -9999.;
  m_partonMomY[0] = -9999.;
  m_partonMomY[1] = -9999.;
  m_partonMomZ[0] = -9999.;
  m_partonMomZ[1] = -9999.;
  m_jetNCst.clear();
  m_jetId.clear();
  m_jetTruId.clear();
  m_jetE.clear();
  m_jetPt.clear();
  m_jetEta.clear();
  m_jetPhi.clear();
  m_cstZ.clear();
  m_cstDr.clear();
  m_cstE.clear();
  m_cstJt.clear();
  m_cstEta.clear();
  m_cstPhi.clear();
  return;

}  // end 'resetTreeVariables()

// end ------------------------------------------------------------------------
