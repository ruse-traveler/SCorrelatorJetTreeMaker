// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMakerConfig.h'
// Derek Anderson
// 03.22.2024
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREEMAKERCONFIG_H
#define SCORRELATORJETTREEMAKERCONFIG_H



namespace SColdQcdCorrelatorAnalysis {

  // SCorrelatorJetTreeMakerConfig definition ---------------------------------

  struct SCorrelatorJetTreeMakerConfig {

    string m_outFileName = "";
    string m_jetTreeName = "";

    bool m_doVtxCut       = false;
    bool m_doQualityPlots = true;
    bool m_requireSiSeeds = true;
    bool m_useOnlyPrimVtx = true;
    bool m_doDcaSigmaCut  = false;
    bool m_maskTpcSectors = false;
    bool m_saveDST        = false;
    bool m_isMC           = true;
    bool m_isEmbed        = false;
    bool m_doDebug        = false;
    bool m_addTracks      = true;
    bool m_addFlow        = false;
    bool m_addECal        = false;
    bool m_addHCal        = false;

    // event acceptance parameters
    // TODO convert most acceptances to pairs/pairs of structs
    double m_evtVzRange[CONST::NRange] = {-10., 10.};
    double m_evtVrRange[CONST::NRange] = {0.0,  0.418};

    // particle acceptance parameters
    double m_parPtRange[CONST::NRange]  = {0.1,  9999.};
    double m_parEtaRange[CONST::NRange] = {-1.1, 1.1};

    // track acceptance parameters
    double m_trkPtRange[CONST::NRange]      = {0.1,  100.};
    double m_trkEtaRange[CONST::NRange]     = {-1.1, 1.1};
    double m_trkQualRange[CONST::NRange]    = {-1.,  10.};
    double m_trkNMvtxRange[CONST::NRange]   = {2.,   100.};
    double m_trkNInttRange[CONST::NRange]   = {1.,   100.};
    double m_trkNTpcRange[CONST::NRange]    = {25.,  100.};
    double m_trkDcaRangeXY[CONST::NRange]   = {-5.,  5.};
    double m_trkDcaRangeZ[CONST::NRange]    = {-5.,  5.};
    double m_trkDeltaPtRange[CONST::NRange] = {0., 0.5};

    // particle flow acceptance parameters
    double m_flowPtRange[CONST::NRange]  = {0.,   9999.};
    double m_flowEtaRange[CONST::NRange] = {-1.1, 1.1};

    // calorimeter acceptance parameters
    double m_ecalPtRange[CONST::NRange]  = {0.,   9999.};
    double m_ecalEtaRange[CONST::NRange] = {-1.1, 1.1};
    double m_hcalPtRange[CONST::NRange]  = {0.,   9999.};
    double m_hcalEtaRange[CONST::NRange] = {-1.1, 1.1};

    // for pt-dependent dca cuts
    TF1*                         m_fSigDcaXY     = NULL;
    TF1*                         m_fSigDcaZ      = NULL;
    double                       m_dcaPtFitMaxXY = 15.;
    double                       m_dcaPtFitMaxZ  = 15.;
    double                       m_nSigCutXY     = 1.;
    double                       m_nSigCutZ      = 1.;
    array<double, CONST::NParam> m_parSigDcaXY   = {1., 1., 1., 1.};
    array<double, CONST::NParam> m_parSigDcaZ    = {1., 1., 1., 1.};

  };  // end SCorrelatorJetTreeMakerConfig

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
