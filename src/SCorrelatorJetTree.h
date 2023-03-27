// ----------------------------------------------------------------------------
// 'SCorrelatorJetTree.h'
// Derek Anderson
// 12.04.2022
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)
// ----------------------------------------------------------------------------

#ifndef SCORRELATORJETTREE_H
#define SCORRELATORJETTREE_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// standard c include
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
// f4a include
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool includes
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
// g4 includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4jets/Jet.h>
#include <g4jets/Jetv1.h>
#include <g4jets/JetMap.h>
#include <g4jets/JetMapv1.h>
#include <g4jets/FastJetAlgo.h>
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
// tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
// calo includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>
// particle flow includes
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>
// fastjet includes
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
// hepmc includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
// root includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TDirectory.h>

#pragma GCC diagnostic pop

using namespace std;
using namespace fastjet;
using namespace findNode;

// forward declarations
class TH1;
class TFile;
class TTree;
class PHG4Particle;
class PHCompositeNode;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class PHG4TruthInfoContainer;
class Fun4AllHistoManager;
class RawClusterContainer;
class RawCluster;
class GlobalVertex;
class SvtxTrackMap;
class JetRecoEval;
class SvtxTrackEval;
class SvtxTrack;
class ParticleFlowElement;

// global constants
static const size_t NPart(2);
static const size_t NComp(3);
static const size_t NRange(2);
static const size_t NMoment(2);
static const size_t NInfoQA(4);
static const size_t NJetType(2);
static const size_t NCstType(4);
static const size_t NObjType(9);
static const size_t NDirectory(NObjType - 3);
static const double MassPion(0.140);



// SCorrelatorJetTree definition ----------------------------------------------

class SCorrelatorJetTree : public SubsysReco {

  public:

    // enums
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
    enum OBJECT {
      TRACK  = 0,
      ECLUST = 1,
      HCLUST = 2,
      FLOW   = 3,
      PART   = 4,
      TJET   = 5,
      RJET   = 6,
      TCST   = 7,
      RCST   = 8
    };
    enum CST_TYPE {
      TRACK_CST = 0,
      CALO_CST  = 1,
      FLOW_CST  = 2,
      PART_CST  = 3
    };
    enum INFO {
      PT  = 0,
      ETA = 1,
      PHI = 2,
      ENE = 3
    };

    // ctor/dtor
    SCorrelatorJetTree(const string &name = "SCorrelatorJetTree", const string &outfile = "correlator_jet_tree.root", const bool isMC = false, const bool debug = false);
    ~SCorrelatorJetTree() override;

    // F4A methods
    int Init(PHCompositeNode *)          override;
    int process_event(PHCompositeNode *) override;
    int End(PHCompositeNode *)           override;

    // particle flow setters (*.io.h)
    void SetParticleFlowMinEta(double etamin) {m_particleflow_mineta = etamin;}
    void SetParticleFlowMaxEta(double etamax) {m_particleflow_maxeta = etamax;}
    void SetParticleFlowEtaAcc(double etamin, double etamax);
    // particle flow getters
    double GetParticleFlowMinEta() {return m_particleflow_mineta;}
    double GetParticleFlowMaxEta() {return m_particleflow_maxeta;}

    // track setters (*.io.h)
    void SetTrackMinPt(double ptmin)   {m_track_minpt = ptmin;}
    void SetTrackMaxPt(double ptmax)   {m_track_maxpt = ptmax;}
    void SetTrackMinEta(double etamin) {m_track_mineta = etamin;}
    void SetTrackMaxEta(double etamax) {m_track_maxeta = etamax;}
    void SetTrackPtAcc(double ptmin, double ptmax);
    void SetTrackEtaAcc(double etamin, double etamax);
    // track getters
    double GetTrackMinPt()  {return m_track_minpt;}
    double GetTrackMaxPt()  {return m_track_maxpt;}
    double GetTrackMinEta() {return m_track_mineta;}
    double GetTrackMaxEta() {return m_track_maxeta;}

    // emcal setters (*.io.h)
    void SetEMCalClusterMinPt(double ptmin)   {m_EMCal_cluster_minpt = ptmin;}
    void SetEMCalClusterMaxPt(double ptmax)   {m_EMCal_cluster_maxpt = ptmax;}
    void SetEMCalClusterMinEta(double etamin) {m_EMCal_cluster_mineta = etamin;}
    void SetEMCalClusterMaxEta(double etamax) {m_EMCal_cluster_maxeta = etamax;}
    void SetEMCalClusterPtAcc(double ptmin, double ptmax);
    void SetEMCalClusterEtaAcc(double etamin, double etamax);
    // emcal getters
    double GetEMCalClusterMinPt()  {return m_EMCal_cluster_minpt;}
    double GetEMCalClusterMaxPt()  {return m_EMCal_cluster_maxpt;}
    double GetEMCalClusterMinEta() {return m_EMCal_cluster_mineta;}
    double GetEMCalClusterMaxEta() {return m_EMCal_cluster_maxeta;}

    // hcal setters (*.io.h)
    void SetHCalClusterMinPt(double ptmin)   {m_HCal_cluster_minpt = ptmin;}
    void SetHCalClusterMaxPt(double ptmax)   {m_HCal_cluster_maxpt = ptmax;}
    void SetHCalClusterMinEta(double etamin) {m_HCal_cluster_mineta = etamin;}
    void SetHCalClusterMaxEta(double etamax) {m_HCal_cluster_maxeta = etamax;}
    void SetHCalClusterPtAcc(double ptmin, double ptmax);
    void SetHCalClusterEtaAcc(double etamin, double etamax);
    // hcal getters
    double GetHCalClusterMinPt()  {return m_HCal_cluster_minpt;}
    double GetHCalClusterMaxPt()  {return m_HCal_cluster_maxpt;}
    double GetHCalClusterMinEta() {return m_HCal_cluster_mineta;}
    double GetHCalClusterMaxEta() {return m_HCal_cluster_maxeta;}

    // particle setters
    void SetParticleMinPt(double ptmin)   {m_MC_particle_minpt = ptmin;}
    void SetParticleMaxPt(double ptmax)   {m_MC_particle_maxpt = ptmax;}
    void SetParticleMinEta(double etamin) {m_MC_particle_mineta = etamin;}
    void SetParticleMaxEta(double etamax) {m_MC_particle_maxeta = etamax;}
    void SetParticlePtAcc(double ptmin, double ptmax);
    void SetParticleEtaAcc(double etamin, double etamx);

    // constituent setters (*.io.h)
    void SetAddParticleFlow(bool b)  {m_add_particleflow = b;}
    void SetAddTracks(bool b)        {m_add_tracks = b;}
    void SetAddEMCalClusters(bool b) {m_add_EMCal_clusters = b;}
    void SetAddHCalClusters(bool b)  {m_add_HCal_clusters = b;}
    // constituent getters
    bool GetAddParticleFlow()  {return m_add_particleflow;}
    bool GetAddTracks()        {return m_add_tracks;}
    bool GetAddEMCalClusters() {return m_add_EMCal_clusters;}
    bool GetAddHCalClusters()  {return m_add_HCal_clusters;}

    // jet setters (*.io.h)
    void SetR(double r)             {m_jetr    = r;}
    void SetType(unsigned int type) {m_jetType = type;}
    void SetJetAlgo(ALGO jetalgo);
    void SetRecombScheme(RECOMB recomb_scheme);
    void SetJetParameters(double r, unsigned int type, ALGO jetalgo, RECOMB recomb_scheme);
    // jet getters
    double                       GetR()            {return m_jetr;}
    unsigned int                 GetType()         {return m_jetType;}
    JetAlgorithm        GetJetAlgo()      {return m_jetalgo;}
    RecombinationScheme GetRecombScheme() {return m_recomb_scheme;}

    // i/o setters
    void SetDoQualityPlots(bool q)     {m_doQualityPlots = q;}
    void SetJetContainerName(string n) {m_jetcontainer_name = n;}
    void SetSaveDST(bool s)            {m_save_dst = s;}
    void SetIsMC(bool b)               {m_ismc = b;}
    void SetSaveDSTMC(bool s)          {m_save_truth_dst = s;}
    // i/o getters
    bool   GetDoQualityPlots()   {return m_doQualityPlots;}
    string GetJetContainerName() {return m_jetcontainer_name;}
    bool   GetSaveDST()          {return m_save_dst;}
    bool   GetIsMC()             {return m_ismc;}
    bool   GetSaveDSTMC()        {return m_save_truth_dst;}

  private:

    // event methods (*.evt.h)
    void   FindPartons(PHCompositeNode *topNode);
    long   GetNumTrks(PHCompositeNode *topNode);
    long   GetNumChrgPars(PHCompositeNode *topNode);
    double GetSumECalEne(PHCompositeNode *topNode);
    double GetSumHCalEne(PHCompositeNode *topNode);
    double GetSumNeutralEne(PHCompositeNode *topNode);
    // jet methods (*.jet.h)
    void FindJets(PHCompositeNode *topNode);
    void FindMcJets(PHCompositeNode *topNode);
    void MatchJets();
    void AddParticleFlow(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap);
    void AddTracks(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap);
    void AddClusters(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap);
    void AddParticles(PHCompositeNode *topNode, vector<PseudoJet> &particles, map<int, pair<Jet::SRC, int>> &fjMap);
    bool IsJetGoodMatch();
    // constituent methods (*.cst.h)
    bool  IsGoodParticleFlow(ParticleFlowElement *pfPart);
    bool  IsGoodTrack(SvtxTrack *track);
    bool  IsGoodEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
    bool  IsGoodHCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
    bool  IsGoodParticle(HepMC::GenParticle *part);
    bool  IsCstGoodMatch();
    float GetParticleCharge(const int pid);
    // system methods (*.sys.h)
    void InitVariables();
    void InitHists();
    void InitTrees();
    void FillTrees();
    void SaveOutput();
    void ResetTreeVariables();
    int  CreateJetNode(PHCompositeNode* topNode);

    // F4A histogram manager
    Fun4AllHistoManager *m_hm;

    // particle flow variables
    double m_particleflow_mineta;
    double m_particleflow_maxeta;
    // track variables
    double m_track_minpt;
    double m_track_maxpt;
    double m_track_mineta;
    double m_track_maxeta;
    // emcal variables
    double m_EMCal_cluster_minpt;
    double m_EMCal_cluster_maxpt;
    double m_EMCal_cluster_mineta;
    double m_EMCal_cluster_maxeta;
    // hcal variables
    double m_HCal_cluster_minpt;
    double m_HCal_cluster_maxpt;
    double m_HCal_cluster_mineta;
    double m_HCal_cluster_maxeta;
    // particle variables
    double m_MC_particle_minpt;
    double m_MC_particle_maxpt;
    double m_MC_particle_mineta;
    double m_MC_particle_maxeta;

    // constituent parameters
    bool m_add_particleflow;
    bool m_add_tracks;
    bool m_add_EMCal_clusters;
    bool m_add_HCal_clusters;
    // jet parameters
    double                m_jetr;
    unsigned int          m_jetType;
    JetAlgorithm          m_jetalgo;
    RecombinationScheme   m_recomb_scheme;
    JetMapv1             *m_jetMap;
    JetMapv1             *m_truth_jetMap;
    // i/o parameters
    string m_outfilename;
    string m_jetcontainer_name;
    bool   m_doQualityPlots;
    bool   m_save_dst;
    bool   m_save_truth_dst;
    bool   m_ismc;
    bool   m_doDebug;

    // output file & trees
    TFile *m_outFile;
    TTree *m_recoTree;
    TTree *m_trueTree;
    TTree *m_matchTree;

    // QA histograms
    TH1D *m_hJetArea[NJetType];
    TH1D *m_hJetNumCst[NJetType];
    TH1D *m_hNumObject[NObjType];
    TH1D *m_hSumCstEne[NCstType];
    TH1D *m_hObjectQA[NObjType][NInfoQA];
    TH1D *m_hNumCstAccept[NCstType][NMoment];

    // for output
    vector<PseudoJet>         m_recoJets;
    vector<PseudoJet>         m_trueJets;
    vector<vector<PseudoJet>> m_recoCsts;
    vector<vector<PseudoJet>> m_trueCsts;

    // output reco event variables
    unsigned long m_recoNumJets           = 0;
    long long     m_recoPartonID[NPart]   = {-9999,  -9999};
    double        m_recoPartonMomX[NPart] = {-9999., -9999.};
    double        m_recoPartonMomY[NPart] = {-9999., -9999.};
    double        m_recoPartonMomZ[NPart] = {-9999., -9999.};
    // output reco jet variables
    vector<unsigned long> m_recoJetNCst;
    vector<unsigned int>  m_recoJetId;
    vector<unsigned int>  m_recoJetTruId;
    vector<double>        m_recoJetE;
    vector<double>        m_recoJetPt;
    vector<double>        m_recoJetEta;
    vector<double>        m_recoJetPhi;
    vector<double>        m_recoJetArea;
    // output reco constituent variables
    vector<vector<double>> m_recoCstZ;
    vector<vector<double>> m_recoCstDr;
    vector<vector<double>> m_recoCstE;
    vector<vector<double>> m_recoCstJt;
    vector<vector<double>> m_recoCstEta;
    vector<vector<double>> m_recoCstPhi;

    // output truth event variables
    unsigned long m_trueNumJets           = 0;
    long long     m_truePartonID[NPart]   = {-9999,  -9999};
    double        m_truePartonMomX[NPart] = {-9999., -9999.};
    double        m_truePartonMomY[NPart] = {-9999., -9999.};
    double        m_truePartonMomZ[NPart] = {-9999., -9999.};
    // output truth jet variables
    vector<unsigned long> m_trueJetNCst;
    vector<unsigned int>  m_trueJetId;
    vector<unsigned int>  m_trueJetTruId;
    vector<double>        m_trueJetE;
    vector<double>        m_trueJetPt;
    vector<double>        m_trueJetEta;
    vector<double>        m_trueJetPhi;
    vector<double>        m_trueJetArea;
    // output truth constituent variables
    vector<vector<double>> m_trueCstZ;
    vector<vector<double>> m_trueCstDr;
    vector<vector<double>> m_trueCstE;
    vector<vector<double>> m_trueCstJt;
    vector<vector<double>> m_trueCstEta;
    vector<vector<double>> m_trueCstPhi;

    // output match event variables
    /* will go here */

};

#endif

// end ------------------------------------------------------------------------
