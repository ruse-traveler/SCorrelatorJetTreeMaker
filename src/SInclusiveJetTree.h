// 'SInclusiveJetTree.cc'
// Derek Anderson
// 12.04.202
//
// Class to construct a tree of
// jets from a specified set of
// events.
//
// Derived from code by Antonio
// Silva (thanks!!)

#ifndef SINCLUSIVEJETTREE_H
#define SINCLUSIVEJETTREE_H

// f4a include
#include <fun4all/SubsysReco.h>
// phool includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>
// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
// misc includes
#include <HepMC/GenEvent.h>
#include <g4jets/Jetv1.h>
#include <g4jets/JetMapv1.h>
// standard c include
#include <string>
#include <vector>

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
static const unsigned long NPart(2);
static const unsigned long NComp(3);



class SInclusiveJetTree : public SubsysReco {

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

    // ctor/dtor
    SInclusiveJetTree(const std::string &name = "SInclusiveJetTree", const std::string &outfile = "inclusive_jet_tree.root");
    ~SInclusiveJetTree() override;

    // F4A methods
    int Init(PHCompositeNode *topNode)          override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode)           override;

    // particle flow setters
    void setParticleFlowMinEta(double etamin) {m_particleflow_mineta = etamin;}
    void setParticleFlowMaxEta(double etamax) {m_particleflow_maxeta = etamax;}
    void setParticleFlowEtaAcc(double etamin, double etamax);
    // particle flow getters
    double getParticleFlowMinEta() {return m_particleflow_mineta;}
    double getParticleFlowMaxEta() {return m_particleflow_maxeta;}

    // track setters
    void setTrackMinPt(double ptmin)   {m_track_minpt = ptmin;}
    void setTrackMaxPt(double ptmax)   {m_track_maxpt = ptmax;}
    void setTrackMinEta(double etamin) {m_track_mineta = etamin;}
    void setTrackMaxEta(double etamax) {m_track_maxeta = etamax;}
    void setTrackPtAcc(double ptmin, double ptmax);
    void setTrackEtaAcc(double etamin, double etamax);
    // track getters
    double getTrackMinPt()  {return m_track_minpt;}
    double getTrackMaxPt()  {return m_track_maxpt;}
    double getTrackMinEta() {return m_track_mineta;}
    double getTrackMaxEta() {return m_track_maxeta;}

    // emcal setters
    void setEMCalClusterMinPt(double ptmin)   {m_EMCal_cluster_minpt = ptmin;}
    void setEMCalClusterMaxPt(double ptmax)   {m_EMCal_cluster_maxpt = ptmax;}
    void setEMCalClusterMinEta(double etamin) {m_EMCal_cluster_mineta = etamin;}
    void setEMCalClusterMaxEta(double etamax) {m_EMCal_cluster_maxeta = etamax;}
    void setEMCalClusterPtAcc(double ptmin, double ptmax);
    void setEMCalClusterEtaAcc(double etamin, double etamax);
    // emcal getters
    double getEMCalClusterMinPt()  {return m_EMCal_cluster_minpt;}
    double getEMCalClusterMaxPt()  {return m_EMCal_cluster_maxpt;}
    double getEMCalClusterMinEta() {return m_EMCal_cluster_mineta;}
    double getEMCalClusterMaxEta() {return m_EMCal_cluster_maxeta;}

    // hcal setters
    void setHCalClusterMinPt(double ptmin)   {m_HCal_cluster_minpt = ptmin;}
    void setHCalClusterMaxPt(double ptmax)   {m_HCal_cluster_maxpt = ptmax;}
    void setHCalClusterMinEta(double etamin) {m_HCal_cluster_mineta = etamin;}
    void setHCalClusterMaxEta(double etamax) {m_HCal_cluster_maxeta = etamax;}
    void setHCalClusterPtAcc(double ptmin, double ptmax);
    void setHCalClusterEtaAcc(double etamin, double etamax);
    // hcal getters
    double getHCalClusterMinPt()  {return m_HCal_cluster_minpt;}
    double getHCalClusterMaxPt()  {return m_HCal_cluster_maxpt;}
    double getHCalClusterMinEta() {return m_HCal_cluster_mineta;}
    double getHCalClusterMaxEta() {return m_HCal_cluster_maxeta;}

    // constituent setters
    void setAddParticleFlow(bool b)  {m_add_particleflow = b;}
    void setAddTracks(bool b)        {m_add_tracks = b;}
    void setAddEMCalClusters(bool b) {m_add_EMCal_clusters = b;}
    void setAddHCalClusters(bool b)  {m_add_HCal_clusters = b;}
    // constituent getters
    bool getAddParticleFlow()  {return m_add_particleflow;}
    bool getAddTracks()        {return m_add_tracks;}
    bool getAddEMCalClusters() {return m_add_EMCal_clusters;}
    bool getAddHCalClusters()  {return m_add_HCal_clusters;}

    // jet setters
    void setR(double r) {m_jetr = r;}
    void setJetAlgo(ALGO jetalgo);
    void setRecombScheme(RECOMB recomb_scheme);
    void setJetParameters(double r, ALGO jetalgo, RECOMB recomb_scheme);
    // jet getters
    double                       getR()            {return m_jetr;}
    fastjet::JetAlgorithm        getJetAlgo()      {return m_jetalgo;}
    fastjet::RecombinationScheme getRecombScheme() {return m_recomb_scheme;}

    // i/o setters
    void setMakeQualityPlots(bool q)        {m_qualy_plots = q;}
    void setJetContainerName(std::string n) {m_jetcontainer_name = n;}
    void setSaveDST(bool s)                 {m_save_dst = s;}
    void setIsMC(bool b)                    {m_ismc = b;}
    void setSaveDSTMC(bool s)               {m_save_truth_dst = s;}
    // i/o getters
    bool        getMakeQualityPlots() {return m_qualy_plots;}
    std::string getJetContainerName() {return m_jetcontainer_name;}
    bool        getSaveDST()          {return m_save_dst;}
    bool        getIsMC()             {return m_ismc;}
    bool        getSaveDSTMC()        {return m_save_truth_dst;}

  private:

    // jet methods
    void findJets(PHCompositeNode *topNode);
    void addParticleFlow(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
    void addTracks(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
    void addClusters(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
    void getTracks(PHCompositeNode *topNode);
    void findNonRecMC(PHCompositeNode *topNode);
    void doMCLoop(PHCompositeNode *topNode);
    // constituent methods
    bool isAcceptableParticleFlow(ParticleFlowElement* pfPart);
    bool isAcceptableTrack(SvtxTrack *track);
    bool isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
    bool isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
    // io methods
    void initializeVariables();
    void initializeTrees();
    int  createJetNode(PHCompositeNode* topNode);
    void resetTreeVariables();

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

    // constituent parameters
    bool m_add_particleflow;
    bool m_add_tracks;
    bool m_add_EMCal_clusters;
    bool m_add_HCal_clusters;
    // jet parameters
    double                        m_jetr;
    fastjet::JetAlgorithm         m_jetalgo;
    fastjet::RecombinationScheme  m_recomb_scheme;
    JetMapv1                     *m_taggedJetMap;
    JetMapv1                     *m_truth_taggedJetMap;
    // i/o parameters
    std::string  m_outfilename;
    int          m_tag_pdg;
    bool         m_qualy_plots;
    bool         m_save_dst;
    bool         m_save_truth_dst;
    bool         m_ismc;
    unsigned int m_jet_id(0);
    unsigned int m_truth_jet_id(0);

    // output file & tree
    TFile *m_outFile;
    TTree *m_jetTree;

    // output event variables
    unsigned long m_numJets           = 0;
    long long     m_partonID[NPart]   = {-9999,  -9999};
    double        m_partonMomX[NPart] = {-9999., -9999.};
    double        m_partonMomY[NPart] = {-9999., -9999.};
    double        m_partonMomZ[NPart] = {-9999., -9999.};
    // output jet variables
    unsigned long m_numCts = 0;
    double        m_jetPt  = -9999.;
    double        m_jetEta = -9999.;
    double        m_jetPhi = -9999.;
    // output constituent variables
    double        m_cstZ   = -9999.;
    double        m_cstDr  = -9999.;
    double        m_cstJt  = -9999.;
    double        m_cstEta = -9999.;
    double        m_cstPhi = -9999.;

};

#endif

// end ------------------------------------------------------------------------
