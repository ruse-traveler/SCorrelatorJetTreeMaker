// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMakerConfig.h'
// Derek Anderson
// 03.22.2024
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Initially derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREEMAKERCONFIG_H
#define SCORRELATORJETTREEMAKERCONFIG_H

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // SCorrelatorJetTreeMakerConfig definition ---------------------------------

  struct SCorrelatorJetTreeMakerConfig {

    // system options
    int    verbosity       {0};
    bool   isDebugOn       {false};
    bool   saveDST         {false};
    bool   isSimulation    {true};
    bool   isEmbed         {false};
    string moduleName      {""};
    string outFileName     {""};
    string trueJetTreeName {""};
    string recoJetTreeName {""};

    // system parameters

    // jet options
    bool   addTracks {true};
    bool   addFlow   {false};
    bool   addECal   {false};
    bool   addHCal   {false};
    float  rJet      {0.4};
    string jetAlgo   {"antikt"};
    string jetRecomb {"pt"};
    string jetArea   {"active"};

    // cut options
    bool doVtxCut       {false};
    bool doDcaSigmaCut  {false};
    bool requireSiSeeds {true};
    bool useOnlyPrimVtx {true};
    bool maskTpcSectors {false};

    // constituent cuts
    pair<Types::ParInfo,   Types::ParInfo>   parAccept;
    pair<Types::TrkInfo,   Types::TrkInfo>   trkAccept;
    pair<Types::FlowInfo,  Types::FlowInfo>  flowAccept;
    pair<Types::ClustInfo, Types::ClustInfo> ecalAccept;
    pair<Types::ClustInfo, Types::ClustInfo> hcalAccept;

    // for pt-dependent dca cuts
    pair<float, float> nSigCut;
    pair<float, float> ptFitMax;
    pair<TF1*,  TF1*>  fSigDca;

  };  // end SCorrelatorJetTreeMakerConfig

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
