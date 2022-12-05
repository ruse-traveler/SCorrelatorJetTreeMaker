// 'SInclusiveJetTree.io.h'
// Derek Anderson
// 12.04.202
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)

#pragma once



SInclusiveJetTree::setParticleFlowEtaAcc(double etamin, double etamax) {
  m_particleflow_mineta = etamin;
  m_particleflow_maxeta = etamax;
  return;
}  // end 'setParticleFlowEtaAcc(double, double)'



SInclusiveJetTree::setTrackPtAcc(double ptmin, double ptmax) {
  m_track_minpt = ptmin;
  m_track_maxpt = ptmax;
  return;
}  // end 'setTrackPtAcc(double, double)'



SInclusiveJetTree::setTrackEtaAcc(double etamin, double etamax) {
  m_track_mineta = etamin;
  m_track_maxeta = etamax;
  return;
}  // end 'setTrackEtaAcc(double, double)'



SInclusiveJetTree::setEMCalClusterPtAcc(double ptmin, double ptmax) {
  m_EMCal_cluster_minpt = ptmin;
  m_EMCal_cluster_maxpt = ptmax;
  return;
}  // end 'setEMCalClusterPtAcc(double, double)'



SInclusiveJetTree::setEMCalClusterEtaAcc(double etamin, double etamax) {
  m_EMCal_cluster_mineta = etamin;
  m_EMCal_cluster_maxeta = etamax;
  return;
}  // end 'setEMCalClusterEtaAcc(double, double)'



SInclusiveJetTree::setHCalClusterPtAcc(double ptmin, double ptmax) {
  m_HCal_cluster_minpt = ptmin;
  m_HCal_cluster_maxpt = ptmax;
  return;
}  // end 'setHCalClusterPtAcc(double, double)'



SInclusiveJetTree::setHCalClusterEtaAcc(double etamin, double etamax) {
  m_HCal_cluster_mineta = etamin;
  m_HCal_cluster_maxeta = etamax;
  return;
}  // end 'setHCalClusterEtaAcc(double, double)'



SInclusiveJetTree::setJetAlgo(ALGO jetalgo) {
  switch (jetalgo) {
    case ALGO::ANTIKT:
      m_jetalgo = fastjet::antikt_algorithm;
      break;
    case ALGO::KT:
      m_jetalgo = fastjet::kt_algorithm;
      break;
    case ALGO::CAMBRIDGE:
      m_jetalgo = fastjet::cambridge_algorithm;
      break;
    default:
      m_jetalgo = fastjet::antikt_algorithm;
      break;
  }
  return;
}  // end 'setJetAlgo(ALGO)'



SInclusiveJetTree::setRecombScheme(RECOMB recomb_scheme) {
  switch(recomb_scheme) {
    case RECOMB::E_SCHEME:
      m_recomb_scheme = fastjet::E_scheme;
      break;
    case RECOMB::PT_SCHEME:
      m_recomb_scheme = fastjet::pt_scheme;
      break;
    case RECOMB::PT2_SCHEME:
      m_recomb_scheme = fastjet::pt2_scheme;
      break;
    case RECOMB::ET_SCHEME:
      m_recomb_scheme = fastjet::Et_scheme;
      break;
    case RECOMB::ET2_SCHEME:
      m_recomb_scheme = fastjet::Et2_scheme;
      break;
    default:
      m_recomb_scheme = fastjet::E_scheme;
      break;
  }
  return;
}  // end 'setRecombScheme(RECOMB)'



SInclusiveJetTree::setJetParameters(double r, ALGO jetalgo, RECOMB recomb_scheme) {
  setR(r);
  setJetAlgo(jetalgo);
  setRecombScheme(recomb_scheme);
  return;
}  // end 'setJetParameters(double, ALGO, RECOMB)'

// end ------------------------------------------------------------------------
