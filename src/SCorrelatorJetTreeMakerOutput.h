// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMakerOutput.h'
// Derek Anderson
// 03.22.2024
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREEMAKEROUTPUT_H
#define SCORRELATORJETTREEMAKEROUTPUT_H

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // SCorrelatorJetTreeMakerOutput definition ---------------------------------

  struct SCorrelatorJetTreeMakerTruthOutput {

    /* TODO will go here */

  };  // end SCorrelatorJetTreeMakerTruthOutput



  struct SCorrelatorJetTreeMakerRecoOutput {

  };  // end SCorrelatorJetTreeMakerRecoOutput



  // SCorrelatorJetTreeMakerLegacyOutput definitions ---------------------------

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



    void SetTreeAddresses(TTree* truth) {
      truth -> Branch("EvtNumJets",     &nJets,           "EvtNumJets/I");
      truth -> Branch("EvtNumChrgPars", &nChrgPars,       "EvtNumChrgPars/I");
      truth -> Branch("EvtVtxX",        &vtxX,            "EvtVtxX/D");
      truth -> Branch("EvtVtxY",        &vtxY,            "EvtVtxY/D");
      truth -> Branch("EvtVtxZ",        &vtxZ,            "EvtVtxZ/D");
      truth -> Branch("EvtSumParEne",   &eSumPar,         "EvtSumParEne/D");
      truth -> Branch("Parton3_ID",     &partonID.first,  "Parton3_ID/I");
      truth -> Branch("Parton4_ID",     &partonID.second, "Parton4_ID/I");
      truth -> Branch("Parton3_MomX",   &partonPX.first,  "Parton3_MomX/D");
      truth -> Branch("Parton3_MomY",   &partonPY.first,  "Parton3_MomY/D");
      truth -> Branch("Parton3_MomZ",   &partonPZ.first,  "Parton3_MomZ/D");
      truth -> Branch("Parton4_MomX",   &partonPX.second, "Parton4_MomX/D");
      truth -> Branch("Parton4_MomY",   &partonPY.second, "Parton4_MomY/D");
      truth -> Branch("Parton4_MomZ",   &partonPZ.second, "Parton4_MomZ/D");
      truth -> Branch("JetNumCst",      &jetNCst);
      truth -> Branch("JetID",          &jetID);
      truth -> Branch("JetEnergy",      &jetE);
      truth -> Branch("JetPt",          &jetPt);
      truth -> Branch("JetEta",         &jetEta);
      truth -> Branch("JetPhi",         &jetPhi);
      truth -> Branch("JetArea",        &jetArea);
      truth -> Branch("CstID",          &cstID);
      truth -> Branch("CstEmbedID",     &cstEmbedID);
      truth -> Branch("CstZ",           &cstZ);
      truth -> Branch("CstDr",          &cstDr);
      truth -> Branch("CstEnergy",      &cstE);
      truth -> Branch("CstPt",          &cstPt);
      truth -> Branch("CstEta",         &cstEta);
      truth -> Branch("CstPhi",         &cstPhi);
      return;
    }  // end 'SetTreeAddresses(TTree*)'



    void GetTreeMakerOutput(SCorrelatorJetTreeMakerTruthOutput& output) {
      /* TODO fill in */
      return;
    }  // end 'GetTreeMakerOutput(SCorrelatorJetTreeMakerTruthOutput&)'

  };  // end SCorrelatorJetTreeMakerLegacyTruthOutput



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



    void SetRecoTreeAddressses(TTree* reco) {
      reco -> Branch("EvtNumJets",    &nJets,    "EvtNumJets/I");
      reco -> Branch("EvtNumTrks",    &nTrks,    "EvtNumTrks/I");
      reco -> Branch("EvtVtxX",       &btxX,     "EvtVtxX/D");
      reco -> Branch("EvtVtxY",       &btxY,     "EvtVtxY/D");
      reco -> Branch("EvtVtxZ",       &btxZ,     "EvtVtxZ/D");
      reco -> Branch("EvtSumECalEne", &eSumECal, "EvtSumECalEne/D");
      reco -> Branch("EvtSumHCalEne", &eSumHCal, "EvtSumHCalEne/D");
      reco -> Branch("JetNumCst",     &jetNCst);
      reco -> Branch("JetID",         &jetID);
      reco -> Branch("JetEnergy",     &jetE);
      reco -> Branch("JetPt",         &jetPt);
      reco -> Branch("JetEta",        &jetEta);
      reco -> Branch("JetPhi",        &jetPhi);
      reco -> Branch("JetArea",       &jetArea);
      reco -> Branch("CstMatchID",    &cstMatchID);
      reco -> Branch("CstZ",          &cstZ);
      reco -> Branch("CstDr",         &cstDr);
      reco -> Branch("CstEnergy",     &cstE);
      reco -> Branch("CstPt",         &cstPt);
      reco -> Branch("CstEta",        &cstEta);
      reco -> Branch("CstPhi",        &cstPhi);
      return;
    }  // end 'SetRecoTreeAddresses(TTree*)'



    void GetTreeMakerOutput(SCorrelatorJetTreeMakerRecoOutput& output) {
      /* TODO fill in */
      return;
    }  // end 'GetTreeMakerOutput(SCorrelatorJetTreeMakerRecoOutput&)'

  };  // end SCorrelatorJetTreeMakerLegacyOutput

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
