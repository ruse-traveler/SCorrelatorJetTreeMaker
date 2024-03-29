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
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <utility>
// root libraries
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
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
#include <scorrelatorutilities/Tools.h>
#include <scorrelatorutilities/Types.h>
#include <scorrelatorutilities/Constants.h>
#include <scorrelatorutilities/Interfaces.h>
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

      // ctor/dtor
      SCorrelatorJetTreeMaker(const string& name = "SCorrelatorJetTreeMaker", const bool debug = false);
      SCorrelatorJetTreeMaker(SCorrelatorJetTreeMakerConfig& config);
      ~SCorrelatorJetTreeMaker() override;

      // F4A methods
      int Init(PHCompositeNode*)          override;
      int process_event(PHCompositeNode*) override;
      int End(PHCompositeNode*)           override;

      // setters
      void SetConfig(const SCorrelatorJetTreeMakerConfig& config) {m_config = config;}

      // getters
      SCorrelatorJetTreeMakerConfig GetConfig() const {return m_config;}

    private:

      // io members
      TFile*    m_outFile     = NULL;
      TTree*    m_recoTree    = NULL;
      TTree*    m_trueTree    = NULL;
      JetMapv1* m_recoJetMap  = NULL;
      JetMapv1* m_trueJetMap  = NULL;

      // track evaluator members
      SvtxEvalStack* m_evalStack = NULL;
      SvtxTrackEval* m_trackEval = NULL;

      // system members
      vector<int>   m_vecEvtsToGrab;
      map<int, int> m_mapCstToEmbedID;

      // jet parameters
      double               m_jetR         = 0.4;
      uint32_t             m_jetType      = 0;
      unique_ptr<JetDefinition>   m_trueJetDef;
      unique_ptr<JetDefinition>   m_recoJetDef;
      unique_ptr<ClusterSequence> m_trueClust    = NULL;
      ClusterSequence*     m_recoClust    = NULL;
      RecombinationScheme  m_recombScheme = pt_scheme;

      // event, jet members
      long long         m_partonID[CONST::NPart];
      CLHEP::Hep3Vector m_partonMom[CONST::NPart];
      CLHEP::Hep3Vector m_trueVtx;
      CLHEP::Hep3Vector m_recoVtx;
      vector<PseudoJet> m_trueJets;
      vector<PseudoJet> m_recoJets;

      // event methods (*.evt.h)
      bool IsGoodVertex(const CLHEP::Hep3Vector vtx);

      // jet methods (*.jet.h)
      void FindTrueJets(PHCompositeNode* topNode);
      void FindRecoJets(PHCompositeNode* topNode);
      void AddParticles(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddTracks(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddFlow(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddECal(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);
      void AddHCal(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap);

      // constituent methods (*.cst.h)
      bool IsGoodParticle(Types::ParInfo& par);
      bool IsGoodTrack(Types::TrkInfo& trk);
      bool IsGoodFlow(Types::FlowInfo& flow);
      bool IsGoodClust(Types::ClustInfo& clust, const int subsys);

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

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
