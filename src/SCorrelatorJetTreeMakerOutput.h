// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMakerOutput.h'
// Derek Anderson
// 03.22.2024
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Initially derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREEMAKEROUTPUT_H
#define SCORRELATORJETTREEMAKEROUTPUT_H

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // SCorrelatorJetTreeMakerOutput definition ---------------------------------

  struct SCorrelatorJetTreeMakerRecoOutput {

    // event-level variables
    Types::RecoInfo evt;

    // jet-level variables
    vector<Types::JetInfo> jets;

    // cst-level variables
    vector<vector<Types::CstInfo>> csts;

    void Reset() {
      evt.Reset();
      jets.clear();
      csts.clear();
      return;
    }  // end 'Reset()'

    void SetTreeAddresses(TTree* tree) {
      tree -> SetBranch("Evt",  "SColdQcdCorrelatorAnalysis::Types::RecoInfo", &evt, 32000, 0);
      tree -> SetBranch("Jets", "std::vector<SColdQcdCorrelatorAnalysis::Types::JetInfo>", &jets, 32000, 0);
      tree -> SetBranch("Csts", "std::vector<std::vector<SColdQcdCorrelatorAnalysis::Types::CstInfo>>", &csts, 32000, 0);
      return;
    }  // end 'SetTreeAddresses(TTree*)'

  };  // end SCorrelatorJetTreeMakerRecoOutput



  struct SCorrelatorJetTreeMakerTruthOutput {

    // event-level variables
    Types::GenInfo evt;

    // jet-level variables
    vector<Types::JetInfo> jets;

    // cst-level variables
    vector<vector<Types::CstInfo>> csts;

    void Reset() {
      evt.Reset();
      jets.clear();
      csts.clear();
      return;
    }  // end 'Reset()'

    void SetTreeAddresses(TTree* tree) {
      tree -> SetBranch("Evt",  "SColdQcdCorrelatorAnalysis::Types::GenInfo", &evt, 32000, 0);
      tree -> SetBranch("Jets", "std::vector<SColdQcdCorrelatorAnalysis::Types::JetInfo>", &jets, 32000, 0);
      tree -> SetBranch("Csts", "std::vector<std::vector<SColdQcdCorrelatorAnalysis::Types::CstInfo>>", &csts, 32000, 0);
      return;
    }  // end 'SetTreeAddresses(TTree*)'

  };  // end SCorrelatorJetTreeMakerTruthOutput



  // SCorrelatorJetTreeMakerLegacyOutput definitions ---------------------------

  struct SCorrelatorJetTreeMakerLegacyRecoOutput {

    // output reco tree event variables
    int    nJets    = numeric_limits<int>::max();
    int    nTrks    = numeric_limits<int>::max();
    double vtxX     = numeric_limits<double>::max();
    double vtxY     = numeric_limits<double>::max();
    double vtxZ     = numeric_limits<double>::max();
    double eSumECal = numeric_limits<double>::max();
    double eSumHCal = numeric_limits<double>::max();

    // output reco tree jet variables
    vector<unsigned long> jetNCst;
    vector<unsigned int>  jetID;
    vector<double>        jetE;
    vector<double>        jetPt;
    vector<double>        jetEta;
    vector<double>        jetPhi;
    vector<double>        jetArea;

    // output reco tree constituent variables
    vector<vector<int>>    cstMatchID;
    vector<vector<double>> cstZ;
    vector<vector<double>> cstDr;
    vector<vector<double>> cstE;
    vector<vector<double>> cstPt;
    vector<vector<double>> cstEta;
    vector<vector<double>> cstPhi;

    void Reset() {
      nJets    = numeric_limits<int>::max();
      nTrks    = numeric_limits<int>::max();
      vtxX     = numeric_limits<double>::max();
      vtxY     = numeric_limits<double>::max();
      vtxZ     = numeric_limits<double>::max();
      eSumECal = numeric_limits<double>::max();
      eSumHCal = numeric_limits<double>::max();
      jetNCst.clear();
      jetID.clear();
      jetE.clear();
      jetPt.clear();
      jetEta.clear();
      jetPhi.clear();
      jetArea.clear();
      cstMatchID.clear();
      cstZ.clear();
      cstDr.clear();
      cstE.clear();
      cstPt.clear();
      cstEta.clear();
      cstPhi.clear();
      return;
    }  // end 'Reset()'

    void SetTreeAddressses(TTree* tree) {
      tree -> Branch("EvtNumJets",    &nJets,    "EvtNumJets/I");
      tree -> Branch("EvtNumTrks",    &nTrks,    "EvtNumTrks/I");
      tree -> Branch("EvtVtxX",       &btxX,     "EvtVtxX/D");
      tree -> Branch("EvtVtxY",       &btxY,     "EvtVtxY/D");
      tree -> Branch("EvtVtxZ",       &btxZ,     "EvtVtxZ/D");
      tree -> Branch("EvtSumECalEne", &eSumECal, "EvtSumECalEne/D");
      tree -> Branch("EvtSumHCalEne", &eSumHCal, "EvtSumHCalEne/D");
      tree -> Branch("JetNumCst",     &jetNCst);
      tree -> Branch("JetID",         &jetID);
      tree -> Branch("JetEnergy",     &jetE);
      tree -> Branch("JetPt",         &jetPt);
      tree -> Branch("JetEta",        &jetEta);
      tree -> Branch("JetPhi",        &jetPhi);
      tree -> Branch("JetArea",       &jetArea);
      tree -> Branch("CstMatchID",    &cstMatchID);
      tree -> Branch("CstZ",          &cstZ);
      tree -> Branch("CstDr",         &cstDr);
      tree -> Branch("CstEnergy",     &cstE);
      tree -> Branch("CstPt",         &cstPt);
      tree -> Branch("CstEta",        &cstEta);
      tree -> Branch("CstPhi",        &cstPhi);
      return;
    }  // end 'SetTreeAddresses(TTree*)'

    void GetTreeMakerOutput(SCorrelatorJetTreeMakerRecoOutput& output) {
      /* TODO fill in */
      return;
    }  // end 'GetTreeMakerOutput(SCorrelatorJetTreeMakerRecoOutput&)'

  };  // end SCorrelatorJetTreeMakerLegacyRecoOutput



  struct SCorrelatorJetTreeMakerLegacyTruthOutput {

    // output truth tree event variables
    int    nJets     = numeric_limits<int>::max();
    int    nChrgPars = numeric_limits<int>::max();
    double eSumPar   = numeric_limits<double>::max();
    double vtxX      = numeric_limits<double>::max();
    double vtxY      = numeric_limits<double>::max();
    double vtxZ      = numeric_limits<double>::max();

    // output truth tree parton variables
    pair<int, int>        partonID = {numeric_limits<int>::max(),    numeric_limits<int>::max()};
    pair<double, double>  partonPX = {numeric_limits<double>::max(), numeric_limits<double>::max()};
    pair<double, double>  partonPY = {numeric_limits<double>::max(), numeric_limits<double>::max()};
    pair<double, double>  partonPZ = {numeric_limits<double>::max(), numeric_limits<double>::max()};

    // output truth tree jet variables
    vector<unsigned long> jetNCst;
    vector<unsigned int>  jetID;
    vector<double>        jetE;
    vector<double>        jetPt;
    vector<double>        jetEta;
    vector<double>        jetPhi;
    vector<double>        jetArea;

    // output truth tree constituent variables
    vector<vector<int>>    cstID;
    vector<vector<int>>    cstEmbedID;
    vector<vector<double>> cstZ;
    vector<vector<double>> cstDr;
    vector<vector<double>> cstE;
    vector<vector<double>> cstPt;
    vector<vector<double>> cstEta;
    vector<vector<double>> cstPhi;

    void Reset() {
      nJets     = numeric_limits<int>::max();
      nChrgPars = numeric_limits<int>::max();
      eSumPar   = numeric_limits<double>::max();
      vtxX      = numeric_limits<double>::max();
      vtxY      = numeric_limits<double>::max();
      vtxZ      = numeric_limits<double>::max();
      partonID  = make_pair(numeric_limits<int>::max(), numeric_limits<int>::max());
      partonPX  = make_pair(numeric_limits<int>::max(), numeric_limits<int>::max());
      partonPY  = make_pair(numeric_limits<int>::max(), numeric_limits<int>::max());
      partonPZ  = make_pair(numeric_limits<int>::max(), numeric_limits<int>::max());
      jetNCst.clear();
      jetID.clear();
      jetE.clear();
      jetPt.clear();
      jetEta.clear();
      jetPhi.clear();
      jetArea.clear();
      cstID.clear();
      cstEmbedID.clear();
      cstZ.clear();
      cstDr.clear();
      cstE.clear();
      cstPt.clear();
      cstEta.clear();
      cstPhi.clear();
      return;
    }  // end 'Reset()'

    void SetTreeAddresses(TTree* tree) {
      tree -> Branch("EvtNumJets",     &nJets,           "EvtNumJets/I");
      tree -> Branch("EvtNumChrgPars", &nChrgPars,       "EvtNumChrgPars/I");
      tree -> Branch("EvtVtxX",        &vtxX,            "EvtVtxX/D");
      tree -> Branch("EvtVtxY",        &vtxY,            "EvtVtxY/D");
      tree -> Branch("EvtVtxZ",        &vtxZ,            "EvtVtxZ/D");
      tree -> Branch("EvtSumParEne",   &eSumPar,         "EvtSumParEne/D");
      tree -> Branch("Parton3_ID",     &partonID.first,  "Parton3_ID/I");
      tree -> Branch("Parton4_ID",     &partonID.second, "Parton4_ID/I");
      tree -> Branch("Parton3_MomX",   &partonPX.first,  "Parton3_MomX/D");
      tree -> Branch("Parton3_MomY",   &partonPY.first,  "Parton3_MomY/D");
      tree -> Branch("Parton3_MomZ",   &partonPZ.first,  "Parton3_MomZ/D");
      tree -> Branch("Parton4_MomX",   &partonPX.second, "Parton4_MomX/D");
      tree -> Branch("Parton4_MomY",   &partonPY.second, "Parton4_MomY/D");
      tree -> Branch("Parton4_MomZ",   &partonPZ.second, "Parton4_MomZ/D");
      tree -> Branch("JetNumCst",      &jetNCst);
      tree -> Branch("JetID",          &jetID);
      tree -> Branch("JetEnergy",      &jetE);
      tree -> Branch("JetPt",          &jetPt);
      tree -> Branch("JetEta",         &jetEta);
      tree -> Branch("JetPhi",         &jetPhi);
      tree -> Branch("JetArea",        &jetArea);
      tree -> Branch("CstID",          &cstID);
      tree -> Branch("CstEmbedID",     &cstEmbedID);
      tree -> Branch("CstZ",           &cstZ);
      tree -> Branch("CstDr",          &cstDr);
      tree -> Branch("CstEnergy",      &cstE);
      tree -> Branch("CstPt",          &cstPt);
      tree -> Branch("CstEta",         &cstEta);
      tree -> Branch("CstPhi",         &cstPhi);
      return;
    }  // end 'SetTreeAddresses(TTree*)'

    void GetTreeMakerOutput(SCorrelatorJetTreeMakerTruthOutput& output) {
      /* TODO fill in */
      return;
    }  // end 'GetTreeMakerOutput(SCorrelatorJetTreeMakerTruthOutput&)'

  };  // end SCorrelatorJetTreeMakerLegacyTruthOutput

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
