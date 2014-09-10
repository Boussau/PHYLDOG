/*
 * Copyright or © or Copr. Centre National de la Recherche Scientifique
 * contributor : Bastien Boussau (2009-2013)
 * 
 * bastien.boussau@univ-lyon1.fr
 * 
 * This software is a computer program whose purpose is to simultaneously build 
 * gene and species trees when gene families have undergone duplications and 
 * losses. It can analyze thousands of gene families in dozens of genomes 
 * simultaneously, and was presented in an article in Genome Research. Trees and 
 * parameters are estimated in the maximum likelihood framework, by maximizing 
 * theprobability of alignments given the species tree, the gene trees and the 
 * parameters of duplication and loss.
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

#include "COALGeneTreeLikelihood.h"
#include "GenericTreeExplorationAlgorithms.h"

#include <iostream>

using namespace bpp;




/******************************************************************************/

/**
 * @brief Build a new COALGeneTreeLikelihood object.
 *
 * @param params The parameters to parse.	 
 * @param spTree The species tree
 * @throw Exception in an error occured.
 */
COALGeneTreeLikelihood::COALGeneTreeLikelihood(std::string file , 
                                               map<string, string> params, 
                                               TreeTemplate<Node> & spTree) 
throw (exception) : GeneTreeLikelihood(file, params, spTree)
{
  
  size_t speciesTreeNodeNumber = spTree.getNumberOfNodes();
  for (int i=0; i< speciesTreeNodeNumber ; i++) 
  {
    //Coalescent model:
    coalBl_.push_back(1.0);
  }  
  
  //coalCounts_: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
  std::vector< std::vector< std::vector< unsigned int > > > coalCounts_2;
  std::vector< std::vector<unsigned int> > coalCounts_3;
  std::vector< unsigned int > coalCounts_4;
  for (unsigned int j = 0 ; j < 2 ; j++ ) {
    coalCounts_4.push_back(0);
  }
  for (unsigned int j = 0 ; j < spTree_->getNumberOfNodes() ; j++ ) {
    coalCounts_3.push_back(coalCounts_4);
  }
  for (unsigned int j = 0 ; j < 3 ; j++ ) {
    coalCounts_2.push_back(coalCounts_3);
  }
  for (unsigned int j = 0 ; j < rootedTree_->getNumberOfNodes() ; j++ ) {
    coalCounts_.push_back(coalCounts_2);
  }
  tentativeCoalCounts_ = coalCounts_;
  for (unsigned int i = 0 ; i < coalBl_.size() ; i++) {
    num12Lineages_.push_back(0);
    num22Lineages_.push_back(0);
  }
  computeNumLineagesFromCoalCounts ();
  
}



// /******************************************************************************/
// COALGeneTreeLikelihood::COALGeneTreeLikelihood(
//   const Tree & tree,
//   SubstitutionModel * model,
//   DiscreteDistribution * rDist,
//   TreeTemplate<Node> & spTree,  
//   TreeTemplate<Node> & rootedTree, 
//   TreeTemplate<Node> & geneTreeWithSpNames,
//   const std::map <std::string, std::string> seqSp,
//   std::map <std::string,int> spId,
//   std::vector < std::vector < std::vector < std::vector<unsigned int> > > > coalCounts,
//   std::vector < double > coalBl,
//   int speciesIdLimitForRootPosition,
//   int heuristicsLevel,
//   int & MLindex, 
//   bool checkRooted,
//   bool verbose,
//   bool rootOptimization, 
//   bool considerSequenceLikelihood, 
//   unsigned int sprLimit)
// throw (Exception):
// GeneTreeLikelihood(tree,
//                    model,
//                    rDist,
//                    spTree,  
//                    rootedTree, 
//                    geneTreeWithSpNames,
//                    seqSp,
//                    spId,
//                    speciesIdLimitForRootPosition,
//                    heuristicsLevel,
//                    MLindex, 
//                    checkRooted,
//                    verbose,
//                    rootOptimization, 
//                    considerSequenceLikelihood, 
//                    sprLimit)
// {
//   coalBl_ = coalBl;
//   coalCounts_ = coalCounts;
//   tentativeCoalCounts_ = coalCounts;
//   for (unsigned int i = 0 ; i < coalBl_.size() ; i++) {
//     num12Lineages_.push_back(0);
//     num22Lineages_.push_back(0);
//   }
//   computeNumLineagesFromCoalCounts ();
// }

/******************************************************************************/

COALGeneTreeLikelihood::COALGeneTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  TreeTemplate<Node> & spTree,  
  TreeTemplate<Node> & rootedTree,  
  TreeTemplate<Node> & geneTreeWithSpNames,
  const std::map <std::string, std::string> seqSp,
  std::map <std::string,int> spId,
  std::vector < std::vector < std::vector < std::vector < unsigned int > > > > coalCounts,
  std::vector < double > coalBl,
  int speciesIdLimitForRootPosition,  
  int heuristicsLevel,
  int & MLindex, 
  std::map <std::string, std::string > params,
  bool checkRooted,
  bool verbose, 
  bool rootOptimization, 
  bool considerSequenceLikelihood, 
  unsigned int sprLimit)
throw (Exception):
GeneTreeLikelihood(tree,
                   data,
                   model,
                   rDist,
                   spTree,  
                   rootedTree, 
                   geneTreeWithSpNames,
                   seqSp,
                   spId,
                   speciesIdLimitForRootPosition,
                   heuristicsLevel,
                   MLindex, 
		   params,
                   checkRooted,
                   verbose,
                   rootOptimization, 
                   considerSequenceLikelihood, 
                   sprLimit)
{
  coalBl_ = coalBl;
  coalCounts_ = coalCounts;
  tentativeCoalCounts_ = coalCounts;
  for (unsigned int i = 0 ; i < coalBl_.size() ; i++) {
    num12Lineages_.push_back(0);
    num22Lineages_.push_back(0);
  }
  
  computeNumLineagesFromCoalCounts ();
}

/******************************************************************************/

COALGeneTreeLikelihood::COALGeneTreeLikelihood(const COALGeneTreeLikelihood & lik):
GeneTreeLikelihood(lik)
{
  coalBl_ = lik.coalBl_;
  coalCounts_ = lik.coalCounts_;
  tentativeCoalCounts_ = lik.tentativeCoalCounts_;
  num12Lineages_ = lik.num12Lineages_;
  num22Lineages_ = lik.num22Lineages_;
  
}

/******************************************************************************/

COALGeneTreeLikelihood & COALGeneTreeLikelihood::operator=(const COALGeneTreeLikelihood & lik)
{
  GeneTreeLikelihood::operator=(lik);
  coalBl_ = lik.coalBl_;
  coalCounts_ = lik.coalCounts_;
  tentativeCoalCounts_ = lik.tentativeCoalCounts_;
  num12Lineages_ = lik.num12Lineages_;
  num22Lineages_ = lik.num22Lineages_;
  return *this;
}




/******************************************************************************/

COALGeneTreeLikelihood::~COALGeneTreeLikelihood()
{
  if (levaluator_) delete levaluator_;
  if (spTree_) delete spTree_;
  if (rootedTree_) delete rootedTree_;
  if (geneTreeWithSpNames_) delete geneTreeWithSpNames_;
}



/******************************************************************************/

void COALGeneTreeLikelihood::initParameters()
{
  // std::cout << "in initParameters"<<std::endl;
  if (considerSequenceLikelihood_) {
    levaluator_->initialize();
  }
  if (heuristicsLevel_>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    MPI::COMM_WORLD.Abort(1);
    exit(-1);
  }
  else {
    /*std::cout << "BEFORE: "<<std::endl;
     *		for (unsigned int i = 0 ; i < coalCounts_[0][0].size() ; i++) {
     *			if (spTree_->getNode(i)->isLeaf() ) {
     *				std::cout << "Leaf i: "<<i<<" coalCounts_[0][0][i][0]: "<< coalCounts_[0][0][i][0] <<" coalCounts_[0][0][i][1] " << coalCounts_[0][0][i][1] << std::endl;
  }
  }*/
    scenarioLikelihood_ = findMLCoalReconciliationDR (spTree_, rootedTree_, 
                                                      seqSp_, spId_, 
                                                      coalBl_, 
                                                      MLindex_, 
                                                      coalCounts_, 
                                                      nodesToTryInNNISearch_);
    computeNumLineagesFromCoalCounts ();
    
  }
  MLindex_ = -1;
}



/******************************************************************************/

void COALGeneTreeLikelihood::resetMLindex() {
  MLindex_ = -1;
}




/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double COALGeneTreeLikelihood::getLogLikelihood() const
{
  cout << "(COAL)sequence likelihood=" <<  levaluator_->getLogLikelihood() << "/nscenario likelihood=" << scenarioLikelihood_ << endl;
  double ll = 0;
  if (considerSequenceLikelihood_) {
    ll = levaluator_->getLogLikelihood() + scenarioLikelihood_;
  }
  else {
    ll = scenarioLikelihood_;
  }
  return ll;
}
/******************************************************************************/
//returns -loglikelihood
//As a reminder: _sequenceLikelihood < 0, scenarioLikelihood_ < 0, getValue >0
double COALGeneTreeLikelihood::getValue() const  
throw (Exception)
{
  {
    if(!levaluator_->isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
    if (considerSequenceLikelihood_) {
      return (- levaluator_->getLogLikelihood() - scenarioLikelihood_);
    }
    else {
      return (- scenarioLikelihood_);
    }
  }
}



/******************************************************************************/
// /*
// 
// void COALGeneTreeLikelihood::fireParameterChanged(const ParameterList & params)
// {
//     if (considerSequenceLikelihood_) {
//         levaluator_->getNniLk()->applyParameters();
//         levaluator_->getNniLk()->matchParametersValues(params);
//     }
//     //If we need to update the reconciliation likelihood
//     if (optimizeReconciliationLikelihood_) {
//         computeReconciliationLikelihood();
//     }
//     
//     
// }*/


/******************************************************************************/
// /*
// 
// void COALGeneTreeLikelihood::computeSequenceLikelihood()
// {
//     if ( considerSequenceLikelihood_ && (optimizeSequenceLikelihood_==true) ) {
//         levaluator_->getNniLk()->computeTreeLikelihood ();
//     }
// }
// */



/******************************************************************************/

void COALGeneTreeLikelihood::computeReconciliationLikelihood()
{
  if (heuristicsLevel_>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    MPI::COMM_WORLD.Abort(1);
    exit(-1);
  }
  else {
    //Compute the COAL likelihood
    scenarioLikelihood_ = findMLCoalReconciliationDR (spTree_, rootedTree_, 
                                                      seqSp_, spId_, 
                                                      coalBl_, 
                                                      tentativeMLindex_, 
                                                      tentativeCoalCounts_, 
                                                      tentativeNodesToTryInNNISearch_); 
    MLindex_ = tentativeMLindex_;
    coalCounts_ = tentativeCoalCounts_;
    nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
    computeNumLineagesFromCoalCounts ();
    
  }
}



// /******************************************************************************/
// 
// 
// void COALGeneTreeLikelihood::computeTreeLikelihood()
// {
//     if (considerSequenceLikelihood_ )
//     {
//         computeSequenceLikelihood();
//     }
//     computeReconciliationLikelihood();  
// }

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.
 * 
 */
double COALGeneTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  //If the NNI is around a branch where a duplication was found, 
  //or if we just try all branches because the starting gene trees are parsimonious in
  //numbers of DL.
  if ((nodesToTryInNNISearch_.count(nodeId)==1) /*|| DLStartingGeneTree_*/) {
    TreeTemplate<Node> * treeForNNI = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
    
    tentativeMLindex_ = MLindex_;
    tentativeCoalCounts_ = coalCounts_;
    
    tentativeNodesToTryInNNISearch_.clear();
    
    //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
    //We use the same procedure as in doNNI !
    const Node * son    = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->getNode(nodeId);
    
    
    if(!son->hasFather()) throw NodeException("COALGeneTreeLikelihood::testNNI(). Node 'son' must not be the root node.", nodeId);
    const Node * parent = son->getFather();
    
    if(!parent->hasFather()) throw NodeException("COALGeneTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent->getId());
    const Node * grandFather = parent->getFather();
    
    //From here: Bifurcation assumed.
    //In case of multifurcation, an arbitrary uncle is chosen.
    //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
    unsigned int parentPosition = grandFather->getSonPosition(parent);
    
    Node * sonForNNI    = treeForNNI->getNode(nodeId);
    Node * parentForNNI = sonForNNI->getFather();
    Node * grandFatherForNNI = parentForNNI->getFather();
    Node * uncleForNNI = grandFatherForNNI->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
    parentForNNI->removeSon(sonForNNI);
    grandFatherForNNI->removeSon(uncleForNNI);
    parentForNNI->addSon(uncleForNNI);
    grandFatherForNNI->addSon(sonForNNI);
    
    //Now we root the tree sent to findMLReconciliation as in rootedTree_
    int id = treeForNNI->getRootNode()->getId();
    if(TreeTemplateTools::hasNodeWithId(*(rootedTree_->getRootNode()->getSon(0)),id)) {
      treeForNNI->newOutGroup(rootedTree_->getRootNode()->getSon(1)->getId());
    }
    else {
      treeForNNI->newOutGroup(rootedTree_->getRootNode()->getSon(0)->getId());
    }
    double candidateScenarioLk = 0;
    totalIterations_ = totalIterations_+1;
    
    double logL = getLogLikelihood();
    double bestlogL = logL;
    
    
    //Compute the COAL likelihood
    candidateScenarioLk =  findMLCoalReconciliationDR (spTree_, treeForNNI, 
                                                       seqSp_, spId_, 
                                                       coalBl_, 
                                                       tentativeMLindex_, 
                                                       tentativeCoalCounts_, 
                                                       tentativeNodesToTryInNNISearch_, false); 
    
    
    if (considerSequenceLikelihood_ ) 
    {
      if  (candidateScenarioLk >  scenarioLikelihood_) 
      { //If it is worth computing the sequence likelihood
        //Retrieving arrays of interest:
        
        levaluator_->setAlternativeTree(treeForNNI);             
        double newLkMinusOldLk = levaluator_->getAlternativeLogLikelihood() -   getSequenceLikelihood();
        
        double tot = - candidateScenarioLk + scenarioLikelihood_ + newLkMinusOldLk;
        
        
        if ( tot < 0) {
          tentativeScenarioLikelihood_=candidateScenarioLk;
        }
        return tot;                
      }
      else 
      {
        tentativeMLindex_ = -1;
        return 1;
      }
    }
    else 
    {
      if  (candidateScenarioLk >  scenarioLikelihood_) 
      {
        tentativeScenarioLikelihood_= candidateScenarioLk;
        return (- candidateScenarioLk + scenarioLikelihood_);
      }
      else 
      {
        tentativeMLindex_ = -1;
        return 1;
      }
    }
  }
  else {
    tentativeMLindex_ = -1;
    return 1;
  }
}

/*******************************************************************************/

void COALGeneTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  levaluator_->acceptAlternativeTree();
  
  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  
  MLindex_ = tentativeMLindex_;
  
  nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
  
  scenarioLikelihood_ = tentativeScenarioLikelihood_;// + _brLikFunction->getValue();
  
  //Now we need to update rootedTree_
  TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
  
  //First we root this temporary tree as in rootedTree_ (same lines as in testNNI)
  int id = tree->getRootNode()->getId();
  if(TreeTemplateTools::hasNodeWithId(*(rootedTree_->getRootNode()->getSon(0)),id)) {
    tree->newOutGroup(rootedTree_->getRootNode()->getSon(1)->getId());
  }
  else {
    tree->newOutGroup(rootedTree_->getRootNode()->getSon(0)->getId());
  }
  //Then we root this tree according to MLindex
  tree->newOutGroup(MLindex_);
  //We update rootedTree_
  if (rootedTree_) delete rootedTree_;
  rootedTree_ = tree->clone();
  delete tree;
  //we need to update the sequence likelihood
  if (considerSequenceLikelihood_ ) 
  {
    OptimizeSequenceLikelihood(true);
  }
  OptimizeReconciliationLikelihood(true);
}

/*******************************************************************************/

std::vector < std::vector < std::vector< std::vector<unsigned int> > > > COALGeneTreeLikelihood::getCoalCounts() const{
  return coalCounts_;
}



void COALGeneTreeLikelihood::computeNumLineagesFromCoalCounts () {  
  for (unsigned int i = 0 ; i < coalCounts_[0][0].size() ; i++) {
    if ( coalCounts_[0][0][i][0]  == 1 && coalCounts_[0][0][i][1] == 2 ) {
      num12Lineages_[i] = 1; 
      num22Lineages_[i] = 0; 
    }
    else  if ( coalCounts_[0][0][i][0]  == coalCounts_[0][0][i][1]  && coalCounts_[0][0][i][1] != 1 ){
      num12Lineages_[i] = 0; 
      num22Lineages_[i] = 1; 
    }
    else {
      num12Lineages_[i] = 0; 
      num22Lineages_[i] = 0; 
    }	
    
  }   
  
  return;
}



std::vector< unsigned int > COALGeneTreeLikelihood::getNum12Lineages() const {
  return num12Lineages_;
}

std::vector< unsigned int > COALGeneTreeLikelihood::getNum22Lineages() const {
  return num22Lineages_;
}


/*******************************************************************************/

std::vector <double> COALGeneTreeLikelihood::getCoalBranchLengths() const{
  return coalBl_;
}

/*******************************************************************************/

void COALGeneTreeLikelihood::setCoalBranchLengths (std::vector < double > coalBl)
{    
  coalBl_ = coalBl;
}

/*******************************************************************************/

int COALGeneTreeLikelihood::getRootNodeindex(){
  return MLindex_;
}

/*******************************************************************************/
/*
 * void COALGeneTreeLikelihood::resetSequenceLikelihood(){
 * _sequenceLikelihood = UNLIKELY;
 * }
 */
/*******************************************************************************/

double COALGeneTreeLikelihood::getSequenceLikelihood() const {
  //TODO différence entre getValue et getLogLikelihood
  //a priori question de signe
  //dispensable ?
  return levaluator_->getLogLikelihood();
}
/*******************************************************************************/

void COALGeneTreeLikelihood::initialize() {
  levaluator_->initialize();
  return;
}



/*******************************************************************************/
void COALGeneTreeLikelihood::print () const {
  std::cout << "Species tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getSpTree(), true)<<std::endl;
  std::cout << "Gene family rooted tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
  std::cout << "Gene family tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (*(levaluator_->getTree()), true)<<std::endl;
  std::cout << "Branch lengths"<<std::endl;
  VectorTools::print(getCoalBranchLengths());
  std::cout << "Root index"<<std::endl;
  std::cout << MLindex_ <<std::endl;
}




/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * To do all this as fast as possible, we optimize opnly a few branch lengths on the SPR tree, 
 * and we use a simple recursion for that.
 ************************************************************************/
void COALGeneTreeLikelihood::refineGeneTreeSPRsFast (map<string, string> params) {
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    //std::cout << "WE DONT DO SPRS"<<std::endl;
    computeReconciliationLikelihood();
    return;
  }
  // std::cout << "WE DO SPRS"<<std::endl;
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  // if (getLogLikelihood()==UNLIKELY) 
  computeReconciliationLikelihood();
  //    ReconciliationTreeLikelihood * bestTreeLogLk = this->clone();
  double logL = getLogLikelihood();
  
  double bestlogL = logL;
  double candidateScenarioLk ;
  double bestSequenceLogL = getSequenceLikelihood();
  double bestScenarioLk = getScenarioLikelihood();
  unsigned int numIterationsWithoutImprovement = 0;
  breadthFirstreNumber (*rootedTree_);
  
  string parentDup;
  string nodeDup;
  string numLoss = "0";
  
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
  
  
  while (numIterationsWithoutImprovement < rootedTree_->getNumberOfNodes() - 2)
  {
    
    annotateGeneTreeWithDuplicationEvents (*spTree_, 
                                           *rootedTree_, 
                                           rootedTree_->getRootNode(), 
                                           seqSp_, spId_); 
    
    for (int nodeForSPR=rootedTree_->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
    {
      Node * n = rootedTree_->getNode(nodeForSPR);
      if (n->hasBranchProperty("L")) {
        numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
      }
      if ( numLoss != "0"  ) {
        
        buildVectorOfRegraftingNodesCoalGeneTree(*spTree_, *rootedTree_, nodeForSPR, sprLimit_, nodeIdsToRegraft);
        
        betterTree = false;
        for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
        {
          if (treeForSPR) 
          {
            delete treeForSPR;
            treeForSPR = 0;
          }
          treeForSPR = rootedTree_->clone();
          
          nodesToUpdate = makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false, true);
          
          //Compute the COAL likelihood
          candidateScenarioLk =  findMLCoalReconciliationDR (spTree_, treeForSPR, 
                                                             seqSp_, spId_, 
                                                             coalBl_, 
                                                             tentativeMLindex_, 
                                                             tentativeCoalCounts_, 
                                                             tentativeNodesToTryInNNISearch_, false); 
          
          if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the scenario likelihood is not bad
          {
            
            if (computeSequenceLikelihoodForSPR) {
              levaluator_->setAlternativeTree ( treeForSPR );
              logL = candidateScenarioLk - levaluator_->getAlternativeLogLikelihood();
            }
            else {
              logL = candidateScenarioLk - bestSequenceLogL;
            }
          }
          else { 
            logL =logL - 10;
          }
          //If the candidate tree has a DL + sequence Lk better than the current best
          if (logL - 0.01 > bestlogL) 
          {
            levaluator_->acceptAlternativeTree();
            std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
            betterTree = true;
            bestlogL = logL;
            bestScenarioLk = candidateScenarioLk;
            if ( computeSequenceLikelihoodForSPR ) {
              bestSequenceLogL = levaluator_->getAlternativeLogLikelihood();
            }
            if ( bestTree ) {
              delete bestTree;
              bestTree = 0;
            }                     
            bestTree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
            
            //Rooting bestTree as in TreeForSPR:
            vector<Node*> rlkNodes = bestTree->getNodes();
            //  std::cout << "BEFORE NEWOUTGROUP: "<<TreeTemplateTools::treeToParenthesis(*bestTree, true)<< std::endl;
            for (unsigned int j = 0 ; j < rlkNodes.size() ; j++) {
              if (rlkNodes[j]->hasNodeProperty("outgroupNode")) {
                if (bestTree->getRootNode() == rlkNodes[j]) {
                  if (j < rlkNodes.size()-1) 
                  {
                    bestTree->rootAt(rlkNodes[rlkNodes.size()-1]);   
                  }
                  else {
                    bestTree->rootAt(rlkNodes[rlkNodes.size()-2]);
                  }
                }
                bestTree->newOutGroup( rlkNodes[j] );
                break;
              }
            }
          }
        }
        if (betterTree) //If, among all the SPRs tried, a better tree has been found 
        {
          
          logL = bestlogL; 
          numIterationsWithoutImprovement = 0;
          if (treeForSPR) 
          {
            delete treeForSPR;
            treeForSPR = 0;
          }
          if (rootedTree_) 
          {
            delete rootedTree_;
            rootedTree_ = 0;
          }
          rootedTree_ = bestTree->clone();
          breadthFirstreNumberAndResetProperties (*rootedTree_);
          
          if (bestTree) {
            delete bestTree;
            bestTree = 0;
          }
        }
        else 
        {
          numIterationsWithoutImprovement++;
        }
        
        if (treeForSPR) 
        {
          delete treeForSPR;
          treeForSPR = 0;
        }
      }
      else {
        numIterationsWithoutImprovement++;
      }
    }
  }
  
  scenarioLikelihood_ = bestScenarioLk;
  //final branch lengths optimization
  
  rootedTree_->resetNodesId();
  
  //One more reconciliation, to update the "_num*Lineages" vectors.
  computeReconciliationLikelihood();
  
  Nhx *nhx = new Nhx();
  annotateGeneTreeWithDuplicationEvents (*spTree_, 
                                         *rootedTree_, 
                                         rootedTree_->getRootNode(), 
                                         seqSp_, spId_); 
  cout << "Reconciled tree: "<<endl;
  nhx->write(*rootedTree_, cout);
  
  if (bestTree) {
    delete bestTree;
    bestTree = 0;
  }
  if (nhx) delete nhx;
  
}




/************************************************************************
 * Tries all NNIs, and accepts NNIs that improve the likelihood as soon as
 * they have been tried.
 ************************************************************************/
void COALGeneTreeLikelihood::refineGeneTreeNNIs(map<string, string> params, unsigned int verbose ) {
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do NNIs
    computeReconciliationLikelihood();
    return;
  }
  bool test = true;
  do
  { 
    TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
    vector<Node *> nodes = tree->getNodes();
    
    vector<Node *> nodesSub = nodes;
    for(unsigned int i = nodesSub.size(); i>0; i--)
    {
      // !!! must not reach i==0 because of unsigned int
      if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.      
      else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node. 
    }
    
    // Test all NNIs:
    test = false;
    for(unsigned int i = 0; !test && i < nodesSub.size(); i++)
    {
      Node* node = nodesSub[i];
      double diff = testNNI(node->getId());
      if(verbose >= 3)
      {
        ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
        + " at " + TextTools::toString(node->getFather()->getId()),
                                        TextTools::toString(diff));
      }
      
      if(diff < 0.)
      { //Good NNI found...
        if(verbose >= 2)
        {
          ApplicationTools::displayResult("   Swapping node " + TextTools::toString(node->getId())
          + " at " + TextTools::toString(node->getFather()->getId()),
                                          TextTools::toString(diff));
        }
        doNNI(node->getId());
        test = true;
        
        if(verbose >= 1)
          ApplicationTools::displayResult("   Current value", TextTools::toString(getValue(), 10));
      }
    }
  }
  while(test);
}


/************************************************************************
 * Tells if the gene family is single copy (1 gene per sp)
 ************************************************************************/
bool COALGeneTreeLikelihood::isSingleCopy() {
  vector<string> names = TreeTemplateTools::getLeavesNames(*(geneTreeWithSpNames_->getRootNode() ) );
  vector<string> uniqueNames = VectorTools::unique(names);
  if (uniqueNames.size() == names.size() ) {
    return true;
  }
  return false;
}





