// ----------------------------------------------------------------------------
// 'SCorrelatorJetTreeMaker.csts.h'
// Derek Anderson
// 01.18.2023
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

  // analysis methods ---------------------------------------------------------

  void SCorrelatorJetTreeMaker::MakeJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::MakeJets(PHCompositeNode*) Making jets" << endl;
    }

    // prepare variables for jet finding
    /* TODO will go here */

    MakeRecoJets(topNode);
    if (m_config.isSimulation) {
      MakeTrueJets(topNode);
    }
    return;

  }  // end 'MakeJets(PHCompositeNode*)'



  void SCorrelatorJetTreeMaker::MakeRecoJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 2)) {
      cout << "SCorrelatorJetTreeMaker::MakeRecoJets(PHCompositeNode*) Making reco jets..." << endl;
    }

    // declare constituent vector and jet map
    //   - TODO make class members
    vector<fastjet::PseudoJet>    particles;
    map<int, pair<Jet::SRC, int>> fjMap;

    // add constitutents
    //   - FIXME tie added consituents to jet type
    if (m_config.addTracks) AddTracks(topNode, particles, fjMap);
    if (m_config.addFlow)   AddFlow(topNode, particles, fjMap);
    if (m_config.addECal)   AddClusts(topNode, {Const::Subsys::EMCal}, particles, fjMap);
    if (m_config.addHCal)   AddClusts(topNode, {Const::Subsys::IHCal, Const::Subsys::OHCal}, particles, fjMap);

    // set reco jet definition
    //   - FIXME move to initialization
    m_recoJetDef = make_unique<JetDefinition>(
      Const::MapStringOntoFJAlgo()[ m_config.jetAlgo ],
      m_config.rJet,
      Const::MapStringOntoFJRecomb()[ m_config.recombScheme ],
      fastjet::Best
    );

    // cluster jets
    m_recoClust = make_unique<ClusterSequence>(particles, *m_recoJetDef);
    m_recoJets  = m_recoClust -> inclusive_jets();
    return;

  }  // end 'MakeRecoJets(PHCompositeNode*)'



  void SCorrelatorJetTreeMaker::MakeTrueJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 2)) {
      cout << "SCorrelatorJetTreeMaker::MakeTrueJets(PHCompositeNode*) Making truth jets..." << endl;
    }

    // declare constituent vector and jet map
    vector<fastjet::PseudoJet>     particles;
    map<int, pair<Jet::SRC, int>>  fjMapMC;

    // add constituents
    AddParticles(topNode, particles, fjMapMC);

    // set truth jet definition
    //   - FIXME move to initialization
    m_trueJetDef = make_unique<JetDefinition>(
      Const::MapStringOntoFJAlgo()[ m_config.jetAlgo ],
      m_config.rJet,
      Const::MapStringOntoFJRecomb()[ m_config.recombScheme ],
      fastjet::Best
    );

    // run clustering, grab jets, and return
    m_trueClust = make_unique<ClusterSequence>(particles, *m_trueJetDef);
    m_trueJets  = m_trueClust -> inclusive_jets();
    return;

  }  // end 'MakeTrueJets(PHCompositeNode*)'



  void SCorrelatorJetTreeMaker::AddTracks(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 3)) {
      cout << "SCorrelatorJetTreeMaker::AddTracks(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding tracks to fastjet input" << endl;
    }

    // loop over tracks
    int64_t       iCst    = particles.size();
    SvtxTrackMap* mapTrks = Interfaces::GetTrackMap(topNode);
    for (
      SvtxTrackMap::Iter itTrk = mapTrks -> begin();
      itTrk != mapTrks -> end();
      ++itTrk
    ) {

      // get track
      SvtxTrack* track = itTrk -> second;
      if (!track) continue;

      // grab info
      Types::TrkInfo info(track, topNode);

      // check if good
      const bool isGoodTrack = IsGoodTrack(info, topNode);
      if (!isGoodTrack) continue;

      // grab barcode of matching particle
      int matchID;
      if (m_config.isSimulation) {
        matchID = Tools::GetMatchID(track, m_svtxEval);
      } else {
        matchID = -1;
      }

      // make pseudojet
      fastjet::PseudoJet pseudojet(
        info.GetPX(),
        info.GetPY(),
        info.GetPZ(),
        info.GetEne()
      );
      fjTrack.set_user_index(matchID);

      // add to lists
      //   - TODO can fold into a templated function
      //   - FIXME probably can get rid of declaration
      pair<int, pair<Jet::SRC, int>> jetTrkPair(iCst, make_pair(Jet::SRC::TRACK, info.GetID()));
      particles.push_back(paseudojet);
      fjMap.insert(jetTrkPair);
      ++iCst;

    }  // end track loop
    return;

  }  // end 'AddTracks(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



  void SCorrelatorJetTreeMaker::AddFlow(PHCompositeNode* topNode, vector<PseudoJet>& particles, map<int, pair<Jet::SRC, int>>& fjMap) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 3)) {
      cout << "SCorrelatorJetTreeMaker::AddFlow(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding particle flow elements to fastjet input" << endl;
    }

    // abort if jets should be charged
    if (m_config.jetType == Const::JetType::Charged) {
      cerr << "SCorrelatorJetTreeMaker::AddFlow: PANIC: trying to add particle flow elements to charged jets! Aborting" << endl;
      assert(m_config.jetType != Const::JetType::Charged);
    }

    // loop over pf elements
    int64_t                       iCst      = particles.size();
    ParticleFlowElementContainer* flowStore = Interfaces::GetFlowStore(topNode);
    for (
      ParticleFlowElementContainer::ConstIterator itFlow = flowStore -> getParticleFlowElements().first;
      itFlow != flowStore -> getParticleFlowElements().second;
      ++itFlow
    ) {

      // get pf element
      ParticleFlowElement* flow = itFlow -> second;
      if (!flow) continue;

      // grab info
      Types::FlowInfo info(flow);

      // check if good
      const bool isGoodFlow = IsGoodFlow(info);
      if (!isGoodFlow) continue;

      // make pseudojet
      fastjet::PseudoJet pseudojet(
        info.GetPX(),
        info.GetPY(),
        info.GetPZ(),
        info.GetEne()
      );
      pseudojet.set_user_index(iCst);

      // add to lists
      //   - TODO can fold into a templated function
      //   - FIXME probably can do without declaration
      pair<int, pair<Jet::SRC, int>> jetPartFlowPair(iCst, make_pair(Jet::SRC::PARTICLE, info.GetID()));
      particles.push_back(pseudojet);
      fjMap.insert(jetPartFlowPair);
      ++iCst;

    }  // end pf element loop
    return;

  }  // end 'AddFlow(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



  void SCorrelatorJetTreeMaker::AddClusts(
    PHCompositeNode* topNode,
    vector<Const::Subsys> vecSubsysToAdd,
    vector<PseudoJet>& particles,
    map<int, pair<Jet::SRC, int>>& fjMap
  ) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 3)) {
      cout << "SCorrelatorJetTreeMaker::AddClusters(PHCompositeNode*, vector<Const::Subsys>, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding calo clusters to fastjet input" << endl;
    }

    // abort if jets should be charged
    if (m_config.jetType == Const::JetType::Charged) {
      cerr << "SCorrelatorJetTreeMaker::AddClusts: ABORT: trying to add calo clusters to charged jets! Aborting" << endl;
      assert(m_config.jetType != Const::JetType::Charged);
    }

    // loop over subsystems to add
    int64_t iCst = particles.size();
    for (Const::Subsys subsys : vecSubsysToAdd) {

      // loop over clusters
      RawClusterContainer::ConstRange clusters = Interfaces::GetClusters(topNode, Const::MapIndexOnotoNode()[ subsys ]);
      for (
        RawClusterContainer::ConstIterator itClust = clusters.first;
        itClust != clusters.second;
        ++itClust
      ) { 

        // get cluster
        RawCluster* cluster = itClust -> second;
        if (!cluster) continue;

        // get primary reconstructed vtx
        ROOT::Math::XYZVector vtx = Interfaces::GetRecoVtx(topNode);

        // grab info
        Types::ClustInfo info(cluster, vtx, subsys);

        // check if good
        const bool isGoodClust = IsGoodClust(info);
        if (!isGoodClust) continue;

        // make pseudojet
        fastjet::PseudoJet pseudojet(
          info.GetPX(),
          info.GetPY(),
          info.GetPZ(),
          info.GetEne()
        );
        pseudojet.set_user_index(iCst);

        // add to lists
        //   - TODO can fold into a templated function
        //   - FIXME probably can get rid of declaration
        pair<int, pair<Jet::SRC, int>> jetClustPair(
          iCst,
          make_pair(
            Const::MapIndexOntoSrc()[ subsys ],
            info.GetID()
          )
        );
        particles.push_back(pseudojet);
        fjMap.insert(jetClustPair);
        ++iCst;

      }  // end cluster loop
    }  // end subsystem loop
    return;

  }  // end 'AddClusts(PHCompositeNode* topNode, vector<Const::Subsys>, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'



  void SCorrelatorJetTreeMaker::AddParticles(
    PHCompositeNode* topNode,
    vector<PseudoJet>& particles,
    map<int, pair<Jet::SRC, int>>& fjMap
  ) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 3)) {
      cout << "SCorrelatorJetTreeMaker::AddParticles(PHComposite*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&) Adding MC particles..." << endl;
    }

    // loop over relevant subevents
    int64_t iCst = particles.size();
    for (const int evtToGrab : m_vecEvtsToGrab) {

      // grab subevent
      HepMC::GenEvent* mcEvt = Interfaces::GetGenEvent(topNode, evtToGrab);

      // grab embedding ID
      const int embedID = GetEmbedID(topNode, evtToGrab);

      // loop over particles in subevent
      for (
        HepMC::GenEvent::particle_const_iterator itPar = mcEvt -> particles_begin();
        itPar != mcEvt -> particles_end();
        ++itPar
      ) {

        // grab info
        Types::ParInfo info(*itPar, embedID);
        if (!info.IsFinalState()) continue;

        // check if particle is good
        const bool isGoodPar = IsGoodParticle(info);
        if (!isGoodPar) continue;

        // map barcode onto relevant embeddingID
        m_mapCstToEmbedID[parID] = embedID;

        // create pseudojet & add to constituent vector
        fastjet::PseudoJet pseudojet(
          info.GetPX(),
          info.GetPY(),
          info.GetPZ(),
          info.GetEne()
        );
        pseudojet.set_user_index(info.GetBarcode());

        // add to lists
        //   - TODO can fold into a templated function
        //   - FIXME probably can get rid of declaration
        pair<int, pair<Jet::SRC, int>> jetParPair(iCst, make_pair(Jet::SRC::PARTICLE, info.GetBarcode()));
        particles.push_back(pseudojet);
        fjMap.insert(jetPartPair);
        ++iCst;

      }  // end particle loop
    }  // end subevent loop
    return;

  }  // end 'AddParticles(PHCompositeNode*, vector<PseudoJet>&, map<int, pair<Jet::SRC, int>>&)'


  bool IsGoodParticle(Types::ParInfo& info) {

  }

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
