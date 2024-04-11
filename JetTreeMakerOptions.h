// ----------------------------------------------------------------------------
// 'JetTreeMakerOptions.h'
// Derek Anderson
// 04.11.2024
//
// Options for the SCorrelatorJetTreeMaker module
// ----------------------------------------------------------------------------

#ifndef CORRELATORQAMAKEROPTIONS_H
#define CORRELATORQAMAKEROPTIONS_H

// c++ utilities
#include <string>
#include <utility>
// analysis libraries
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Types.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Constants.h"
#include "/sphenix/user/danderson/install/include/scorrelatorjettreemaker/SCorrelatorJetTreeMakerConfig.h"

// make common namespaces implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis;



namespace JetTreeMakerOptions {

  // track & particle flow parameters
  const bool   runTracking(false);
  const bool   doTruthTableReco(false);
  const double nSigma(1.5);

  // topo cluster parameters
  const double showerR(0.025);
  const double noiseLevels[NTopoPar]   = {0.0025, 0.006, 0.03};
  const double significance[NTopoPar]  = {4.0,    2.0,   0.0};
  const double localMinE[NTopoPar]     = {1.0,    2.0,   0.5};
  const bool   enableHCal[NTopoClusts] = {false, true};
  const bool   enableECal[NTopoClusts] = {true, false};
  const bool   doSplit(true);
  const bool   allowCorners(true);

  // jet tree general parameters
  const bool isMC(true);
  const bool isEmbed(true);
  const bool doDebug(false);
  const bool saveDst(true);
  const bool doVtxCut(false);
  const bool doQuality(true);
  const bool requireSiSeeds(true);
  const bool useOnlyPrimVtx(true);
  const bool doDcaSigmaCut(false);
  const bool maskTpcSectors(false);
  const bool addTracks(true);
  const bool addECal(false);
  const bool addHCal(false);
  const bool addParticleFlow(false);

  // jet tree jet parameters
  const double       jetRes  = 0.4;
  const unsigned int jetType = 0;
  const auto         jetAlgo = SCorrelatorJetTreeMaker::ALGO::ANTIKT;
  const auto         jetReco = SCorrelatorJetTreeMaker::RECOMB::PT_SCHEME;

  // event acceptance
  const pair<double, double> vzEvtRange = {-10., 10.};
  const pair<double, double> vrEvtRange = {0.0,  0.418};

  // particle acceptance
  const pair<double, double> ptParRange  = {0.,   9999.};
  const pair<double, double> etaParRange = {-1.1, 1.1};

  // track acceptance
  const pair<double, double> ptTrackRange      = {0.2,  100.};
  const pair<double, double> etaTrackRange     = {-1.1, 1.1};
  const pair<double, double> qualTrackRange    = {0.,   10.};
  const pair<double, double> nMvtxTrackRange   = {2.,   100.};
  const pair<double, double> nInttTrackRange   = {1.,   100.};
  const pair<double, double> nTpcTrackRange    = {24.,  100.};
  const pair<double, double> dcaTrackRangeXY   = {-5.,  5.};
  const pair<double, double> dcaTrackRangeZ    = {-5.,  5.};
  const pair<double, double> deltaPtTrackRange = {0., 0.5};

  // for pt dependent dca cuts
  const pair<double, double> dcaPtFitMax      = {15., 15.};
  const pair<double, double> nDcaSigmaTrack   = {3.,  3.};
  const vector<double>       dcaSigmaParamsXY = {-0.0095, 0.091, -0.029};
  const vector<double>       dcaSigmaParamsZ  = {1.73,    26.1,  -9.45};

  // particle flow acceptance
  const pair<double, double> ptFlowRange  = {0.2,  9999.};
  const pair<double, double> etaFlowRange = {-1.1, 1.1};

  // calo acceptance
  const pair<double, double> ptECalRange  = {0.3,  9999.};
  const pair<double, double> etaECalRange = {-1.1, 1.1};
  const pair<double, double> ptHCalRange  = {0.3,  9999.};
  const pair<double, double> etaHCalRange = {-1.1, 1.1};



  // bundle acceptances into pairs --------------------------------------------

  pair<Types::TrkInfo, Types::TrkInfo> GetTrkAccept() {

    // create maximal range
    pair<Types::TrkInfo, Types::TrkInfo> trkAccept = {
      Types::TrkInfo(Const::Init::Minimize),
      Types::TrkInfo(Const::Init::Maximize)
    };

    // set specific bounds

  }  // end 'GetTrkAccept()'




  // set up configuration -----------------------------------------------------

  SCorrelatorJetTreeMakerConfig GetConfig(const int verbosity, const string outFile) {

    SCorrelatorJetTreeMakerConfig cfg {

    };
    return cfg;

  }  // end 'GetConfig()'

}  // end JetTreeMakerOptions namespace

#endif

// end ------------------------------------------------------------------------

