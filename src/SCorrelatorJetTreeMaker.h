// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.h'
// Derek Anderson
// 12.04.2022
//
// A module to produce a tree of jets for the sPHENIX
// Cold QCD Energy-Energy Correlator analysis.
//
// Initially derived from code by Antonio Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREEMAKER_H
#define SCORRELATORJETTREEMAKER_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// c++ utilities
#include <map>
#include <array>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <cstdlib>
#include <utility>
// root libraries
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TDirectory.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
// hepmc libraries
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
// f4a utilities
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool libraries
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
// truth utilities
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
// jet utilities
#include <jetbase/Jet.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetMapv1.h>
#include <jetbase/FastJetAlgo.h>
// calo utilities
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>
// trackng utilities
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
// particle flow utilities
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>
// vtx utilities
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Tools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Types.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Constants.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Interfaces.h"
// analysis definitions
#include "SCorrelatorJetTreeMakerConfig.h"
#include "SCorrelatorJetTreeMakerOutput.h"

#pragma GCC diagnostic pop

using namespace std;
using namespace fastjet;
using namespace findNode;



namespace SColdQcdCorrelatorAnalysis {

  // SCorrelatorJetTreeMaker definition ---------------------------------------

  class SCorrelatorJetTreeMaker : public SubsysReco {

    public:

      // public enums
      enum ALGO {
        ANTIKT    = 0,
        KT        = 1,
        CAMBRIDGE = 2
      };
      enum RECOMB {
        E_SCHEME   = 0,
        PT_SCHEME  = 1,
        PT2_SCHEME = 2,
        ET_SCHEME  = 3,
        ET2_SCHEME = 4
      };

      // ctor/dtor
      SCorrelatorJetTreeMaker(const string& name = "SCorrelatorJetTreeMaker", const string& outFile = "correlator_jet_tree.root", const bool isMC = false, const bool isEmbed = false, const bool debug = false);
      ~SCorrelatorJetTreeMaker() override;

      // F4A methods
      int Init(PHCompositeNode*)          override;
      int process_event(PHCompositeNode*) override;
      int End(PHCompositeNode*)           override;

    private:

      // constants
      enum CONST {
        NPart      = 2,
        NComp      = 3,
        NParam     = 4,
        NRange     = 2,
        NMoment    = 2,
        NInfoQA    = 9,
        NJetType   = 2,
        NCstType   = 5,
        NObjType   = 9,
        NDirectory = 6,
        NMvtxLayer = 3,
        NInttLayer = 8,
        NTpcLayer  = 48,
        NTpcSector = 12
      };

      // qa info & tracking subsystems
      enum SUBSYS   {MVTX, INTT, TPC};
      enum CST_TYPE {PART_CST, TRACK_CST, FLOW_CST, ECAL_CST, HCAL_CST};
      enum OBJECT   {TRACK, ECLUST, HCLUST, FLOW, PART, TJET, RJET, TCST, RCST};
      enum INFO     {PT, ETA, PHI, ENE, QUAL, DCAXY, DCAZ, DELTAPT, NTPC};

      // event methods (*.evt.h)
      bool IsGoodVertex(const CLHEP::Hep3Vector vtx);
      void GetEventVariables(PHCompositeNode* topNode);

      // jet methods (*.jet.h)
      void FindTrueJets(PHCompositeNode* topNode);
      void FindRecoJets(PHCompositeNode* topNode);
      void AddParticles(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddTracks(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddFlow(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddECal(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddHCal(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);

      // constituent methods (*.cst.h)
      bool IsGoodParticle(HepMC::GenParticle* par, const bool ignoreCharge = false);
      bool IsGoodTrack(SvtxTrack* track, PHCompositeNode* topNode);
      bool IsGoodFlow(ParticleFlowElement* flow);
      bool IsGoodECal(CLHEP::Hep3Vector& hepVecECal);
      bool IsGoodHCal(CLHEP::Hep3Vector& hepVecHCal);
      bool IsGoodTrackSeed(SvtxTrack* track);
      bool IsGoodTrackPhi(SvtxTrack* track, const float phiMaskSize = 0.01);  // FIXME make user configurable

      // system methods (*.sys.h)
      void InitVariables();
      void InitHists();
      void InitTrees();
      void InitFuncs();
      void InitEvals(PHCompositeNode* topNode);
      void FillTrueTree();
      void FillRecoTree();
      void SaveOutput();
      void ResetVariables();
      void DetermineEvtsToGrab(PHCompositeNode* topNode);
      int  CreateJetNode(PHCompositeNode* topNode);

      // F4A/utility members
      Fun4AllHistoManager* m_histMan   = NULL;
      SvtxEvalStack*       m_evalStack = NULL;
      SvtxTrackEval*       m_trackEval = NULL;

      // io members
      TFile*    m_outFile     = NULL;
      TTree*    m_recoTree    = NULL;
      TTree*    m_trueTree    = NULL;
      JetMapv1* m_recoJetMap  = NULL;
      JetMapv1* m_trueJetMap  = NULL;

      // system members
      vector<int>   m_vecEvtsToGrab;
      map<int, int> m_mapCstToEmbedID;

      // jet parameters
      double               m_jetR         = 0.4;
      uint32_t             m_jetType      = 0;
      JetAlgorithm         m_jetAlgo      = antikt_algorithm;
      JetDefinition*       m_trueJetDef   = NULL;
      JetDefinition*       m_recoJetDef   = NULL;
      ClusterSequence*     m_trueClust    = NULL;
      ClusterSequence*     m_recoClust    = NULL;
      RecombinationScheme  m_recombScheme = pt_scheme;

      // event, jet members
      long long         m_partonID[CONST::NPart];
      CLHEP::Hep3Vector m_partonMom[CONST::NPart];
      CLHEP::Hep3Vector m_trueVtx;
      CLHEP::Hep3Vector m_recoVtx;
      vector<PseudoJet> m_trueJets;
      vector<PseudoJet> m_recoJets;

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
