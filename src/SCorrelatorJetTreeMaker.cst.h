// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.csts.h'
// Derek Anderson
// 01.18.2023
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#pragma once

using namespace std;
using namespace findNode;



namespace SColdQcdCorrelatorAnalysis {

  // constituent methods ------------------------------------------------------

  bool SCorrelatorJetTreeMaker::IsGoodParticle(HepMC::GenParticle* par, const bool ignoreCharge) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodParticle(HepMC::GenParticle*) Checking if MC particle is good..." << endl;
    }

    // check charge if needed
    const bool isJetCharged  = (m_jetType != 1);
    const bool doChargeCheck = (isJetCharged && !ignoreCharge);

    int   parID;
    bool  isGoodCharge;
    float parChrg;
    if (doChargeCheck) {
      parID        = par -> pdg_id();
      parChrg      = GetParticleCharge(parID);
      isGoodCharge = (parChrg != 0.);
    } else {
      isGoodCharge = true;
    }

    const double parEta       = par -> momentum().eta();
    const double parPx        = par -> momentum().px();
    const double parPy        = par -> momentum().py();
    const double parPt        = sqrt((parPx * parPx) + (parPy * parPy));
    const bool   isInPtRange  = ((parPt  > m_parPtRange[0])  && (parPt  < m_parPtRange[1]));
    const bool   isInEtaRange = ((parEta > m_parEtaRange[0]) && (parEta < m_parEtaRange[1]));
    const bool   isGoodPar    = (isGoodCharge && isInPtRange && isInEtaRange);
    return isGoodPar;

  }  // end 'IsGoodParticle(HepMC::GenParticle*, bool)'



  bool SCorrelatorJetTreeMaker::IsGoodTrack(SvtxTrack* track, PHCompositeNode* topNode) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodTrack(SvtxTrack*) Checking if track is good..." << endl;
    }

    // grab track info
    const double trkPt      = track -> get_pt();
    const double trkEta     = track -> get_eta();
    const double trkQual    = track -> get_quality();
    const double trkDeltaPt = GetTrackDeltaPt(track);
    const int    trkNMvtx   = GetNumLayer(track, SUBSYS::MVTX);
    const int    trkNIntt   = GetNumLayer(track, SUBSYS::INTT);
    const int    trkNTpc    = GetNumLayer(track, SUBSYS::TPC);

    // grab track dca
    const auto   trkDca   = GetTrackDcaPair(track, topNode);
    const double trkDcaXY = trkDca.first;
    const double trkDcaZ  = trkDca.second;


    // if above max pt used to fit dca width,
    // use value of fit at max pt
    double ptEvalXY = (trkPt > m_dcaPtFitMaxXY) ? m_dcaPtFitMaxXY : trkPt;
    double ptEvalZ  = (trkPt > m_dcaPtFitMaxZ)  ? m_dcaPtFitMaxZ  : trkPt;

    // check if dca is good
    bool isInDcaRangeXY = false;
    bool isInDcaRangeZ  = false;
    if (m_doDcaSigmaCut) {
      isInDcaRangeXY = (abs(trkDcaXY) < (m_nSigCutXY * (m_fSigDcaXY -> Eval(ptEvalXY))));
      isInDcaRangeZ  = (abs(trkDcaZ)  < (m_nSigCutZ  * (m_fSigDcaZ  -> Eval(ptEvalZ))));
    } else {
      isInDcaRangeXY = ((trkDcaXY > m_trkDcaRangeXY[0]) && (trkDcaXY < m_trkDcaRangeXY[1]));
      isInDcaRangeZ  = ((trkDcaZ  > m_trkDcaRangeZ[0])  && (trkDcaZ  < m_trkDcaRangeZ[1]));
    }  

    // if applying vertex cuts, grab track
    // vertex and check if good
    bool isInVtxRange = true;
    if (m_doVtxCut) {
      CLHEP::Hep3Vector trkVtx = GetTrackVertex(track, topNode);
      isInVtxRange = IsGoodVertex(trkVtx);
    }

    // if using only primary vertex,
    // ignore tracks from other vertices
    if (m_useOnlyPrimVtx) {
      const bool isFromPrimVtx = IsFromPrimaryVtx(track, topNode);
      if (!isFromPrimVtx) {
        isInVtxRange = false;
      }
    }

    // if masking tpc sector boundaries,
    // ignore tracks near boundaries
    bool isGoodPhi = true;
    if (m_maskTpcSectors) {
      isGoodPhi = IsGoodTrackPhi(track);
    }

    // apply cuts
    const bool isSeedGood       = IsGoodTrackSeed(track);
    const bool isInPtRange      = ((trkPt      > m_trkPtRange[0])      && (trkPt      <  m_trkPtRange[1]));
    const bool isInEtaRange     = ((trkEta     > m_trkEtaRange[0])     && (trkEta     <  m_trkEtaRange[1]));
    const bool isInQualRange    = ((trkQual    > m_trkQualRange[0])    && (trkQual    <  m_trkQualRange[1]));
    const bool isInNMvtxRange   = ((trkNMvtx   > m_trkNMvtxRange[0])   && (trkNMvtx   <= m_trkNMvtxRange[1]));
    const bool isInNInttRange   = ((trkNIntt   > m_trkNInttRange[0])   && (trkNIntt   <= m_trkNInttRange[1]));
    const bool isInNTpcRange    = ((trkNTpc    > m_trkNTpcRange[0])    && (trkNTpc    <= m_trkNTpcRange[1]));
    const bool isInDeltaPtRange = ((trkDeltaPt > m_trkDeltaPtRange[0]) && (trkDeltaPt <  m_trkDeltaPtRange[1]));
    const bool isInNumRange     = (isInNMvtxRange && isInNInttRange && isInNTpcRange);
    const bool isInDcaRange     = (isInDcaRangeXY && isInDcaRangeZ);
    const bool isGoodTrack      = (isSeedGood && isGoodPhi && isInPtRange && isInEtaRange && isInQualRange && isInNumRange && isInDcaRange && isInDeltaPtRange && isInVtxRange);
    return isGoodTrack;

  }  // end 'IsGoodTrack(SvtxTrack*)'



  bool SCorrelatorJetTreeMaker::IsGoodFlow(ParticleFlowElement* flow) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodFlow(ParticleFlowElement*) Checking if particle flow element is good..." << endl;
    }

    // TODO: explore particle flow cuts
    const double pfEta        = flow -> get_eta();
    const bool   isInEtaRange = ((pfEta > m_flowEtaRange[0]) && (pfEta < m_flowEtaRange[1]));
    const bool   isGoodFlow   = isInEtaRange;
    return isGoodFlow;

  }  // end 'IsGoodFlow(ParticleFlowElement*)'



  bool SCorrelatorJetTreeMaker::IsGoodECal(CLHEP::Hep3Vector& hepVecECal) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodECal(CLHEP::Hep3Vector&) Checking if ECal cluster is good..." << endl;
    }

    const double clustPt      = hepVecECal.perp();
    const double clustEta     = hepVecECal.pseudoRapidity();
    const bool   isInPtRange  = ((clustPt  > m_ecalPtRange[0])  && (clustPt  < m_ecalPtRange[1]));
    const bool   isInEtaRange = ((clustEta > m_ecalEtaRange[0]) && (clustEta < m_ecalEtaRange[1]));
    const bool   isGoodClust  = (isInPtRange && isInEtaRange);
    return isGoodClust;

  }  // end 'IsGoodECal(CLHEP::Hep3Vector&)'



  bool SCorrelatorJetTreeMaker::IsGoodHCal(CLHEP::Hep3Vector& hepVecHCal) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 1)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodHCal(CLHEP::Hep3Vector&) Checking if HCal cluster is good..." << endl;
    }

    // TODO: explore particle cuts. These should vary with particle charge/species.
    const double clustPt      = hepVecHCal.perp();
    const double clustEta     = hepVecHCal.pseudoRapidity();
    const bool   isInPtRange  = ((clustPt  > m_hcalPtRange[0])  && (clustPt  < m_hcalPtRange[1]));
    const bool   isInEtaRange = ((clustEta > m_hcalEtaRange[0]) && (clustEta < m_hcalEtaRange[1]));
    const bool   isGoodClust  = (isInPtRange && isInEtaRange);
    return isGoodClust;

  }  // end 'IsGoodHCal(CLHEP::Hep3Vector&)'



  bool SCorrelatorJetTreeMaker::IsGoodTrackSeed(SvtxTrack* track) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 2)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodSeedTrack(SvtxTrack*) Checking if track seed is good..." << endl;
    }

    // get track seeds
    TrackSeed* trkSiSeed  = track -> get_silicon_seed();
    TrackSeed* trkTpcSeed = track -> get_tpc_seed();

    // check if one or both seeds are present as needed
    bool isSeedGood = (trkSiSeed && trkTpcSeed);
    if (!m_requireSiSeeds) {
      isSeedGood = (trkSiSeed || trkTpcSeed);
    }
    return isSeedGood;

  }  // end 'IsGoodSeed(SvtxTrack*)'



  bool SCorrelatorJetTreeMaker::IsGoodTrackPhi(SvtxTrack* track, const float phiMaskSize) {

    // print debug statement
    if (m_doDebug && (Verbosity() > 2)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodTrackPhi(SvtxTrack*) Checking if track phi is good..." << endl;
    }

    // TPC sector boundaries:
    //   12 sectors --> ~0.523 rad/sector,
    //   assumed to be symmetric about phi = 0
    // FIXME move to constant in utilities namespace
    const array<float, NTpcSector> phiSectorBoundaries = {
      -2.877,
      -2.354,
      -1.831,
      -1.308,
      -0.785,
      -0.262,
      0.262,
      0.785,
      1.308,
      1.831,
      2.354,
      2.877
    };

    // flag phi as bad if within boundary +- (phiMaskSize / 2)
    const double halfMaskSize = phiMaskSize / 2.;
    const double trkPhi       = track -> get_phi();

    // loop over sector boundaries and check phi
    bool isGoodPhi = true;
    for (const float boundary : phiSectorBoundaries) {
      if ((trkPhi > (boundary - halfMaskSize)) && (trkPhi < (boundary + halfMaskSize))) {
        isGoodPhi = false;
        break;
      }
    }
    return isGoodPhi;

  }  // end 'IsGoodTrackPhi(SvtxTrack*, float)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
