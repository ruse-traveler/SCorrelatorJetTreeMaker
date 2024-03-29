// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.evt.h'
// Derek Anderson
// 03.28.2023
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

  // event methods ------------------------------------------------------------

  bool SCorrelatorJetTreeMaker::IsGoodVertex(const CLHEP::Hep3Vector vtx) {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::IsGoodVertex(CLHEP::Hep3Vector) Checking if event is good..." << endl;
    }

    // calculate vr
    const double vr = sqrt((vtx.x() * vtx.x()) + (vtx.y() * vtx.y()));

    // check if event is good
    const bool isInEvtVzRange = ((vtx.z() > m_evtVzRange[0]) && (vtx.z() < m_evtVzRange[1]));
    const bool isInEvtVrRange = ((abs(vr) > m_evtVrRange[0]) && (abs(vr) < m_evtVrRange[1]));
    const bool isGoodVertex   = (isInEvtVzRange && isInEvtVrRange);
    return isGoodVertex;

  }  // end 'IsGoodVertex(CLHEP::Hep3Vector)'



  void SCorrelatorJetTreeMaker::GetEventVariables(PHCompositeNode* topNode) {

    // print debug statement
    if (m_doDebug) {
      cout << "SCorrelatorJetTreeMaker::GetEventVariables(PHCompositeNode*) Grabbing event info..." << endl;
    }

    m_recoVtx     = GetRecoVtx(topNode);
    m_recoNumTrks = GetNumTrks(topNode);
    m_recoSumECal = GetSumECalEne(topNode);
    m_recoSumHCal = GetSumHCalEne(topNode);
    if (m_isMC) {
      m_trueNumChrgPars = GetNumChrgPars(topNode);
      m_trueSumPar      = GetSumParEne(topNode);
    }
    return;

  }  // end 'GetEventVariables(PHCompositeNode*)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
