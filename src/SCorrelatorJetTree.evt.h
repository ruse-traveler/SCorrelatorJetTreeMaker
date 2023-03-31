// ----------------------------------------------------------------------------
// 'SCorrelatorJetTree.evt.h'
// Derek Anderson
// 03.28.2023
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)
// ----------------------------------------------------------------------------

#pragma once

using namespace std;
using namespace findNode;



// event methods --------------------------------------------------------------

void SCorrelatorJetTree::GetEventVariables(PHCompositeNode *topNode) {

  m_recoVtx     = GetRecoVtx(topNode);
  m_recoNumTrks = GetNumTrks(topNode);
  m_recoSumECal = GetSumECalEne(topNode);
  m_recoSumHCal = GetSumHCalEne(topNode);
  if (m_isMC) {
    m_trueNumChrgPars = GetNumChrgPars(topNode);
    m_trueSumPar      = GetSumParEne(topNode);
    m_trueVtx         = GetTrueVtx(topNode);
  }
  return;

}  // end 'GetEventVariables(PHCompositeNode*)'



void SCorrelatorJetTree::GetPartonInfo(PHCompositeNode *topNode) {

  /* TODO will find partons here */
  return;

}  // end 'GetPartonInfo(PHCompositeNode*)'



long SCorrelatorJetTree::GetNumTrks(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetNumTrks(PHCompositeNode*) Calculating no. of tracks..." << endl;
  }

  // loop over tracks
  long          nTrk    = 0;
  SvtxTrack    *track   = 0x0;
  SvtxTrackMap *mapTrks = GetTrackMap(topNode);
  for (SvtxTrackMap::Iter itTrk = mapTrks -> begin(); itTrk != mapTrks -> end(); ++itTrk) {

    // get track
    track = itTrk -> second;
    if (!track) continue;

    // if good, add to count
    const bool isGoodTrack = IsGoodTrack(track);
    if (isGoodTrack) ++nTrk;

  }  // end track loop
  return nTrk;

}  // end 'GetNumTrks(PHCompositeNode*)'



long SCorrelatorJetTree::GetNumChrgPars(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetNumChrgPars(PHCompositeNode*) Calculating no. of charged particles..." << endl;
  }

  // loop over particles
  long             nPar  = 0;
  HepMC::GenEvent *mcEvt = GetMcEvent(topNode);
  for (HepMC::GenEvent::particle_const_iterator itPar = mcEvt -> particles_begin(); itPar != mcEvt -> particles_end(); ++itPar) {

    // check if particle is final state
    const bool isFinalState = ((*itPar) -> status() == 1);
    if (!isFinalState) continue;

    // if good, add to count
    const bool isGoodPar = IsGoodParticle(*itPar);
    if (isGoodPar) ++nPar;

  }  // end particle loop
  return nPar;

}  // end 'GetNumChrgPars(PHCompositeNode*)'



double SCorrelatorJetTree::GetSumECalEne(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetSumECalEne(PHCompositeNode*) Getting sum of ECal energy..." << endl;
  }

  // grab vertex and clusters
  GlobalVertex        *vtx          = GetGlobalVertex(topNode);
  RawClusterContainer *emClustStore = GetClusterStore(topNode, "CLUSTER_CEMC");

  // loop over emcal clusters
  double                             eneECalSum   = 0.;
  RawClusterContainer::ConstRange    emClustRange = emClustStore -> getClusters();
  RawClusterContainer::ConstIterator itEMClust;
  for (itEMClust = emClustRange.first; itEMClust != emClustRange.second; itEMClust++) {

    // grab cluster
    const RawCluster *emClust = itEMClust -> second;
    if (!emClust) continue;

    // construct vertex and get 4-momentum
    const double vX = vtx -> get_x();
    const double vY = vtx -> get_y();
    const double vZ = vtx -> get_z();

    CLHEP::Hep3Vector hepVecVtx     = CLHEP::Hep3Vector(vX, vY, vZ);
    CLHEP::Hep3Vector hepVecEMClust = RawClusterUtility::GetECoreVec(*emClust, hepVecVtx);

    // if good, add to sum
    const bool isGoodECal = IsGoodECal(hepVecEMClust);
    if (isGoodECal) eneECalSum += hepVecEMClust.mag();

  }  // end emcal loop
  return eneECalSum;

}  // end 'GetSumECalEne(PHCompositeNode*)'



double SCorrelatorJetTree::GetSumHCalEne(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetSumHCalEne(PHCompositeNode*) Getting sum of HCal energy..." << endl;
  }

  // grab vertex and clusters
  GlobalVertex        *vtx          = GetGlobalVertex(topNode);
  RawClusterContainer *ihClustStore = GetClusterStore(topNode, "CLUSTER_HCALIN");
  RawClusterContainer *ohClustStore = GetClusterStore(topNode, "CLUSTER_HCALOUT");

  // loop over ihcal clusters
  double                             eneIHCalSum  = 0.;
  RawClusterContainer::ConstRange    ihClustRange = ihClustStore -> getClusters();
  RawClusterContainer::ConstIterator itIHClust;
  for (itIHClust = ihClustRange.first; itIHClust != ihClustRange.second; itIHClust++) {

    // grab cluster
    const RawCluster *ihClust = itIHClust -> second;
    if (!ihClust) continue;

    // construct vertex and get 4-momentum
    const double vX = vtx -> get_x();
    const double vY = vtx -> get_y();
    const double vZ = vtx -> get_z();

    CLHEP::Hep3Vector hepVecVtx     = CLHEP::Hep3Vector(vX, vY, vZ);
    CLHEP::Hep3Vector hepVecIHClust = RawClusterUtility::GetECoreVec(*ihClust, hepVecVtx);

    // if good, add to sum
    const bool isGoodHCal = IsGoodECal(hepVecIHClust);
    if (isGoodHCal) eneIHCalSum += hepVecIHClust.mag();
  }  // end ihcal loop

  // loop over ohcal clusters
  double                             eneOHCalSum  = 0.;
  RawClusterContainer::ConstRange    ohClustRange = ohClustStore -> getClusters();
  RawClusterContainer::ConstIterator itOHClust;
  for (itOHClust = ohClustRange.first; itOHClust != ohClustRange.second; itOHClust++) {

    // grab cluster
    const RawCluster *ohClust = itOHClust -> second;
    if (!ohClust) continue;

    // construct vertex and get 4-momentum
    const double vX = vtx -> get_x();
    const double vY = vtx -> get_y();
    const double vZ = vtx -> get_z();

    CLHEP::Hep3Vector hepVecVtx     = CLHEP::Hep3Vector(vX, vY, vZ);
    CLHEP::Hep3Vector hepVecOHClust = RawClusterUtility::GetECoreVec(*ohClust, hepVecVtx);

    // if good, add to sum
    const bool isGoodHCal = IsGoodECal(hepVecOHClust);
    if (isGoodHCal) eneOHCalSum += hepVecOHClust.mag();
  }  // end ohcal loop

  const double eneHCalSum = eneIHCalSum + eneOHCalSum;
  return eneHCalSum;

}  // end 'GetSumHCalEne(PHCompositeNode*)'



double SCorrelatorJetTree::GetSumParEne(PHCompositeNode *topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetSumParEne(PHComposite*) Calculating sum of particle energy..." << endl;
  }

  // loop over particles
  double           eSumPar = 0.;
  HepMC::GenEvent *mcEvt   = GetMcEvent(topNode);
  for (HepMC::GenEvent::particle_const_iterator itPar = mcEvt -> particles_begin(); itPar != mcEvt -> particles_end(); ++itPar) {

    // check if particle is final state
    const bool isFinalState = ((*itPar) -> status() == 1);
    if (!isFinalState) continue;

    // if good, add to count
    const bool isGoodPar = IsGoodParticle(*itPar, true);
    if (isGoodPar) eSumPar += (*itPar) -> momentum().e();

  }  // end particle loop
  return eSumPar;

}  // end 'GetSumParEne(PHCompositeNode*)'



CLHEP::Hep3Vector SCorrelatorJetTree::GetRecoVtx(PHCompositeNode* topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetRecoVtx(PHComposite*) Getting reconstructed vertex..." << endl;
  }

  const GlobalVertex      *vtx     = GetGlobalVertex(topNode);
  const double             vtxX    = vtx -> get_x();
  const double             vtxY    = vtx -> get_y();
  const double             vtxZ    = vtx -> get_z();
  const CLHEP::Hep3Vector  recoVtx = CLHEP::Hep3Vector(vtxX, vtxY, vtxZ);
  return recoVtx;

}  // end 'GetRecoVtx(PHCompositeNode*)'



CLHEP::Hep3Vector SCorrelatorJetTree::GetTrueVtx(PHCompositeNode* topNode) {

  // print debug statement
  if (m_doDebug) {
    cout << "SCorrelatorJetTree::GetTrueVtx(PHComposite*) Getting truth vertex..." << endl;
  }

  // grab event
  //HepMC::GenEvent *mcEvt = GetMcEvent(topNode);

  // get vertex
  /* TODO vertex grabbing goes here */
  const CLHEP::Hep3Vector trueVtx = CLHEP::Hep3Vector(0., 0., 0.);
  return trueVtx;

}  // end 'GetRecoVtx(PHCompositeNode*)'

// end ------------------------------------------------------------------------
