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

// user include
#include "SInclusiveJetTree.h"
// f4a includes
#include <fun4all/Fun4AllReturnCodes.h>
// phool includes
#include <phool/PHCompositeNode.h>




//____________________________________________________________________________..
SInclusiveJetTree::SInclusiveJetTree(const std::string &name):
 SubsysReco(name)
{
  std::cout << "SInclusiveJetTree::SInclusiveJetTree(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
SInclusiveJetTree::~SInclusiveJetTree()
{
  std::cout << "SInclusiveJetTree::~SInclusiveJetTree() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int SInclusiveJetTree::Init(PHCompositeNode *topNode)
{
  std::cout << "SInclusiveJetTree::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::InitRun(PHCompositeNode *topNode)
{
  std::cout << "SInclusiveJetTree::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::process_event(PHCompositeNode *topNode)
{
  std::cout << "SInclusiveJetTree::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "SInclusiveJetTree::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::EndRun(const int runnumber)
{
  std::cout << "SInclusiveJetTree::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::End(PHCompositeNode *topNode)
{
  std::cout << "SInclusiveJetTree::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SInclusiveJetTree::Reset(PHCompositeNode *topNode)
{
 std::cout << "SInclusiveJetTree::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void SInclusiveJetTree::Print(const std::string &what) const
{
  std::cout << "SInclusiveJetTree::Print(const std::string &what) const Printing info for " << what << std::endl;
}

// end ------------------------------------------------------------------------
