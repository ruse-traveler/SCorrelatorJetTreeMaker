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

  void SCorrelatorJetTreeMaker::GetEventVariables(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::GetEventVariables(PHCompositeNode*) Grabbing event info..." << endl;
    }

    // get indices of relevant subevents
    m_vecEvtsToGrab = Tools::GrabSubevents(topNode, subEvtOpt, m_config.isEmbed);

    // set event info
    m_recoOutput.evt.SetInfo(topNode);
    if (m_config.isSimulation) {
      m_trueOutput.evt.SetInfo(topNode, m_config.isEmbed, m_vecEvtsToGrab);
    }
    return;

  }  // end 'GetEventVariables(PHCompositeNode*)'



  void SCorrelatorJetTreeMaker::MakeJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 1)) {
      cout << "SCorrelatorJetTreeMaker::MakeJets(PHCompositeNode*) Making jets" << endl;
    }

    // prepare variables for jet finding
    ResetJetVariables();

    // do jet finding
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

    // cluster jets
    ClusterSequence clustering(particles, *m_recoJetDef);
    m_recoJets = clustering.inclusive_jets();
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

    // run clustering, grab jets, and return
    ClusterSequence clustering(particles, *m_trueJetDef);
    m_trueJets = clustering.inclusive_jets();
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
      const bool isGoodTrack = IsGoodTrack(info, track, topNode);
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
        const bool isGoodClust = IsGoodCluster(info);
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



  bool SCorrelatorJetTreeMaker::IsGoodTrack(Types::TrkInfo& info, SvtxTrack* track, PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 4)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodTrack(Types:TrkInfo&, SvtxTrack*, PHCompositeNode*) Checking if track is good..." << endl;
    }

    // if needed, check if dca is in pt-dependent range
    bool isInDcaSigma = true;
    if (m_config.doDcaSigCut) {
      isInDcaSigma = info.IsInSigmaDcaCut(
        m_config.nSigCut,
        m_config.ptFitMax,
        m_config.fSigDca
      );
    }

    // if needed, check if track vertex is in acceptance
    const bool isInVtxRange = m_config.doVtxCut ? IsGoodVtx((info.GetVX(), info.GetVY(), info.GetVZ())) : true;

    // if needed, check if track is from primary vertex,
    const bool isFromPrimVtx = m_config.useOnlyPrimVtx ? info.IsFromPrimaryVtx(topNode) : true;

    // check seed, if in acceptance, and return overall goodness
    const bool isSeedGood = Tools::IsGoodTrackSeed(track, m_config.requireSiSeed);
    const bool isInAccept = info.IsInAcceptance(m_config.trkAccept);
    return (isInDcaSigma && isFromPrimVtx && isSeedGood && isInAccept);

  }  // end 'IsGoodTrack(Types::TrkInfo&, SvtxTrack*, PHCompositeNode*)'



  bool SCorrelatorJetTreeMaker::IsGoodFlow(Types::FlowInfo& info) {

    // print debug statement
    if (m_config.isebugOn && (m_config.verbosity > 4)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodFlow(Types::FlowInfo&) Checking if particle flow element is good..." << endl;
    }

    const bool isInAccept = info.IsInAcceptance(m_config.flowAccept);
    return isInAccept;

  }  // end 'IsGoodFlow(Types::FlowInfo&)'



  bool IsGoodCluster(Types::ClustInfo& info, Const::Subsys subsys) {

    // print debug statement
    if (m_dconfig.isDebugOn && (m_config.verbosity > 4)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodCluster(types::ClustInfo&, Const::subsys) Checking if cluster is good..." << endl;
    }

    // check if in acceptance
    bool isInAccept;
    switch (subsys) {

      case Const::Subsys::EMCal:
        isInAccept = info.IsInAcceptance(m_config.ecalAccept);
        break;

      case Const::Subsys::IHCal:
        [[fallthrough]];

      case Const::Subsys::OHCal:
        isInAccept = info.IsInAcceptance(m_config.hcalAccept);
        break;

      default:
        isInAccept = true;
        break;
    }
    return isInAccept;

  }  // end 'IsGoodCluster(Types::ClustInfo&, Const::Subsys)'



  bool IsGoodParticle(Types::ParInfo& info) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 4)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodParticle(Types::ParInfo&) Checking if particle is good..." << endl;
    }

    // check charge if needed
    //   - TODO allow for different truth/reco jet types
    bool isGoodCharge;
    switch (m_config.jetType) {

      case Const::JetType::Charged:
        isGoodCharge = (info.GetCharge() != 0.);
        break;

      case Const::JetType::Neutral:
        isGoodCharge = (info.GetCharge() == 0.);
        break;

      case Const::JetType::Full:
        [[fallthrough]];

      default:
        isGoodCharge = true;
        break;
    }

    // check if in acceptance and return
    const bool isInAccept = info.IsInAcceptance(m_config.parAccept);
    return (isGoodCharge && isInAccept);

  }  // end 'IsGoodParticle(Types::ParInfo&)'



  bool SCorrelatorJetTreeMaker::IsGoodVertex(const ROOT::Math::XYZVector vtx) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SCorrelatorJetTreeMaker::IsGoodVertex(ROOT::Math::XYZVector) Checking if vertex is good..." << endl;
    }

    // check if vertex is in acceptance
    const bool isInVrAccept = ((vtx.Rho() >= m_config.vrAccept.first) && (vtx.Rho() <= m_config.vzAccept.second));
    const bool isInVzAccept = ((vtx.Z()   >= m_config.vzAccept.frist) && (vtx.Z()   <= m_config.vzAccept.second));
    return (isInVrAccept && isInVzAccept);

  }  // end 'IsGoodVertex(const ROOT::Math::XYZVector)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
