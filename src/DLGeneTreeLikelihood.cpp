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


#include "DLGeneTreeLikelihood.h"
#include "GenericTreeExplorationAlgorithms.h"

using namespace bpp;

#define FAST 0


  /******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(std::string file, 
					     std::map<std::string, std::string> params, 
					     TreeTemplate<Node> & spTree ) 
  throw (exception) : GeneTreeLikelihood(file, params, spTree)
  {
    size_t speciesTreeNodeNumber = spTree.getNumberOfNodes();
    for (int i=0; i< speciesTreeNodeNumber ; i++) 
    {
        //We fill the vectors with 0.1s until they are the right sizes.
        //DL model:
        lossExpectedNumbers_.push_back(0.1);
        duplicationExpectedNumbers_.push_back(0.11);
        num0Lineages_.push_back(0);
        num1Lineages_.push_back(0);
        num2Lineages_.push_back(0);
    }
  tentativeNum0Lineages_ =num0Lineages_;
  tentativeNum1Lineages_ =num1Lineages_; 
  tentativeNum2Lineages_ =num2Lineages_;
  }

/******************************************************************************/
DLGeneTreeLikelihood::DLGeneTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  TreeTemplate<Node> & spTree,
  TreeTemplate<Node> & rootedTree,
  TreeTemplate<Node> & geneTreeWithSpNames,
  const std::map <std::string, std::string> seqSp,
  std::map <std::string,int> spId,
  //std::vector <int> & lossNumbers, 
  std::vector <double> & lossProbabilities, 
  //std::vector <int> & duplicationNumbers, 
  std::vector <double> & duplicationProbabilities,
  //std::vector <int> & branchNumbers,
  std::vector <int> & num0Lineages,
  std::vector <int> & num1Lineages,
  std::vector <int> & num2Lineages, 
  int speciesIdLimitForRootPosition,
  int heuristicsLevel,
  int & MLindex, 
  bool checkRooted,
  bool verbose, 
  bool rootOptimization, 
  bool considerSequenceLikelihood, 
  bool DLStartingGeneTree, 
  unsigned int sprLimit)
throw (Exception):
GeneTreeLikelihood(tree,
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
		   checkRooted,
		   verbose,
		   rootOptimization, 
		   considerSequenceLikelihood, 
		   sprLimit)

{
  lossExpectedNumbers_ = lossProbabilities;
  duplicationExpectedNumbers_ = duplicationProbabilities;
  num0Lineages_=num0Lineages;
  num1Lineages_=num1Lineages;
  num2Lineages_=num2Lineages;
  tentativeNum0Lineages_ =num0Lineages;
  tentativeNum1Lineages_ =num1Lineages; 
  tentativeNum2Lineages_ =num2Lineages;
  tentativeMLindex_ = MLindex;
  DLStartingGeneTree_ = DLStartingGeneTree;
}

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  TreeTemplate<Node> & spTree,
  TreeTemplate<Node> & rootedTree,
  TreeTemplate<Node> & geneTreeWithSpNames,
  const std::map <std::string, std::string> seqSp,
  std::map <std::string,int> spId,
  std::vector <double> & lossProbabilities, 
  std::vector <double> & duplicationProbabilities, 
  std::vector <int> & num0Lineages,
  std::vector <int> & num1Lineages,
  std::vector <int> & num2Lineages, 
  int speciesIdLimitForRootPosition, 
  int heuristicsLevel,
  int & MLindex, 
  bool checkRooted,
  bool verbose,
  bool rootOptimization, 
  bool considerSequenceLikelihood, 
  bool DLStartingGeneTree,
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
		   checkRooted,
		   verbose,
		   rootOptimization, 
		   considerSequenceLikelihood, 
		   sprLimit)
{
  lossExpectedNumbers_ = lossProbabilities;
  duplicationExpectedNumbers_ = duplicationProbabilities; 
  num0Lineages_=num0Lineages;
  num1Lineages_=num1Lineages;
  num2Lineages_=num2Lineages;
  tentativeNum0Lineages_ =num0Lineages;
  tentativeNum1Lineages_ =num1Lineages; 
  tentativeNum2Lineages_ =num2Lineages;
  DLStartingGeneTree_ = DLStartingGeneTree;
}

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(const DLGeneTreeLikelihood & lik):
GeneTreeLikelihood(lik)
{
  lossExpectedNumbers_ = lik.lossExpectedNumbers_;
  duplicationExpectedNumbers_ = lik.duplicationExpectedNumbers_; 
  num0Lineages_=lik.num0Lineages_;
  num1Lineages_=lik.num1Lineages_;
  num2Lineages_=lik.num2Lineages_;
  tentativeNum0Lineages_ =lik.tentativeNum0Lineages_;
  tentativeNum1Lineages_ =lik.tentativeNum1Lineages_;
  tentativeNum2Lineages_ =lik.tentativeNum2Lineages_;
  DLStartingGeneTree_ = lik.DLStartingGeneTree_;
}

/******************************************************************************/

DLGeneTreeLikelihood & DLGeneTreeLikelihood::operator=(const DLGeneTreeLikelihood & lik)
{
  GeneTreeLikelihood::operator=(lik);
  lossExpectedNumbers_ = lik.lossExpectedNumbers_;
  duplicationExpectedNumbers_ = lik.duplicationExpectedNumbers_;
  num0Lineages_=lik.num0Lineages_;
  num1Lineages_=lik.num1Lineages_;
  num2Lineages_=lik.num2Lineages_;
  tentativeNum0Lineages_ =lik.tentativeNum0Lineages_;
  tentativeNum1Lineages_ =lik.tentativeNum1Lineages_;
  tentativeNum2Lineages_ =lik.tentativeNum2Lineages_;
  DLStartingGeneTree_ = lik.DLStartingGeneTree_;
  return *this;
}




/******************************************************************************/

DLGeneTreeLikelihood::~DLGeneTreeLikelihood()
{
  if (nniLk_) delete nniLk_;
  if (spTree_) delete spTree_;
  if (rootedTree_) delete rootedTree_;
  if (geneTreeWithSpNames_) delete geneTreeWithSpNames_;
}



/******************************************************************************/

void DLGeneTreeLikelihood::initParameters()
{
  //     std::cout << "in initParameters"<<std::endl;
  if (considerSequenceLikelihood_) {
    nniLk_->initParameters();
  }
  if (heuristicsLevel_>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    exit(-1);
    //    scenarioLikelihood_ = findMLReconciliation (&spTree_, &rootedTree_, seqSp_, lossNumbers_, lossExpectedNumbers_, duplicationNumbers_, duplicationExpectedNumbers_, MLindex_, branchNumbers_, _speciesIdLimitForRootPosition_, heuristicsLevel_, num0Lineages_, num1Lineages_, num2Lineages_, nodesToTryInNNISearch_); 
  }
  else {
    //  std::cout << "in initParameters 2"<<std::endl;
    
    scenarioLikelihood_ = findMLReconciliationDR (spTree_, rootedTree_, 
						  seqSp_, spId_, lossExpectedNumbers_, 
						  duplicationExpectedNumbers_, MLindex_, 
						  num0Lineages_, num1Lineages_,
						  num2Lineages_, nodesToTryInNNISearch_); 
    //Rooting bestTree as in TreeForSPR:
    //std::cout << "in initParameters 3"<<std::endl;
    
    vector<Node*> nodes = rootedTree_->getNodes();
    //  std::cout << "in initParameters 4"<<std::endl;
    
    for (unsigned int j = 0 ; j < nodes.size() ; j++) {
      // std::cout << "in initParameters 5"<<std::endl;
      
      if (nodes[j]->hasNodeProperty("outgroupNode")) {
	//    std::cout << "in initParameters 6"<<std::endl;
	
	if (rootedTree_->getRootNode() == nodes[j]) {
	  if (j < nodes.size()-1) 
	  {
	    rootedTree_->rootAt(nodes[nodes.size()-1]);   
	  }
	  else {
	    rootedTree_->rootAt(nodes[nodes.size()-2]);
	  }
	};
	rootedTree_->newOutGroup( nodes[j] );
	//    std::cout << "FOUND"<<std::endl;
	break;
      }
    }
  }
  MLindex_ = -1;
  // std::cout << "in ReconciliationTreeLikelihood::initParameters : num0Lineages_, num1Lineages_, num2Lineages_ : "<< TextTools::toString(VectorTools::sum(num0Lineages_))<<" "<<TextTools::toString(VectorTools::sum(num1Lineages_))<<" "<< TextTools::toString(VectorTools::sum(num2Lineages_))<<std::endl;
  
  //    std::cout <<"INITIAL scenarioLikelihood_ "<<scenarioLikelihood_ <<" nniLk_->getLogLikelihood(): "<< nniLk_->getLogLikelihood()<<std::endl;
}



/******************************************************************************/

void DLGeneTreeLikelihood::resetMLindex() {
  MLindex_ = -1;
}




/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double DLGeneTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  // _sequenceLikelihood=UNLIKELY;
  //  if (_sequenceLikelihood ==UNLIKELY) {
  
  //TEST
  /*
   *     if ((_sequenceLikelihood == UNLIKELY)||(optimizeSequenceLikelihood_==true)) {
   *     Vdouble * lik = & likelihoodData_->getRootRateSiteLikelihoodArray();
   *     const std::vector<unsigned int> * w = & likelihoodData_->getWeights();
   *     for(unsigned int i = 0; i < nbDistinctSites_; i++)
   *     {
   *     ll += (* w)[i] * log((* lik)[i]);
   *     //std::cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << std::endl;
}
_sequenceLikelihood = ll;
// std::cout <<"getLogLikelihood _sequenceLikelihood "<< _sequenceLikelihood<<std::endl;
}
*/
  //END TEST
  
  // ll = getValue();
  
  /* else {
   *     std::cout <<"NOT RECOMPUTING IN GETLOGLIKELIHOOD : "<< _sequenceLikelihood<<std::endl; 
}*/
  //Now we add the scenario likelihood
  //  std::cout << "COMPUTING LOGLIKELIHOOD :" << ll  << " SCENARIO LOGLIKELIHOOD :" << scenarioLikelihood_ <<" TOTAL : "<< ll + scenarioLikelihood_ <<std::endl;
  
  //TEST 16 02 2010
  // DRHomogeneousTreeLikelihood::getLogLikelihood();
  //  setMinuslogLikelihood_ (_sequenceLikelihood);
  if (considerSequenceLikelihood_) {
    ll = levaluator_->getLogLikelihood() + scenarioLikelihood_;
  }
  else {
    ll = scenarioLikelihood_;
  }
  // ll = _sequenceLikelihood ;
  return ll;
}
/******************************************************************************/
//returns -loglikelihood
//As a reminder: _sequenceLikelihood < 0, scenarioLikelihood_ < 0, getValue >0
double DLGeneTreeLikelihood::getValue() const  
throw (Exception)
{
  {
    if(!levaluator_->getNniLk()->isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
    // return (-getLogLikelihood());
    //TEST 16 02 2010
    // std::cout<<"\t\t\t_sequenceLikelihood: "<<_sequenceLikelihood<< " scenarioLikelihood_: "<<scenarioLikelihood_<<std::endl;
    if (considerSequenceLikelihood_) {
      std::cout << "CONSIDERING" <<std::endl;
      return (- levaluator_->getLogLikelihood() - scenarioLikelihood_);
    }
    else {
            std::cout << "NOT CONSIDERING" <<std::endl;

      return (- scenarioLikelihood_);
    }
    //return (minusLogLik_ - scenarioLikelihood_);
    //return (-minusLogLik_);
  }
}



/******************************************************************************/


void DLGeneTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  if (considerSequenceLikelihood_) {
    levaluator_->getNniLk()->applyParameters();
    levaluator_->getNniLk()->matchParametersValues(params);
  }
  //If we need to update the reconciliation likelihood
  if (optimizeReconciliationLikelihood_) {
    computeReconciliationLikelihood();
  }
  
  
}


/******************************************************************************/


void DLGeneTreeLikelihood::computeSequenceLikelihood()
{
  if ( considerSequenceLikelihood_ && (optimizeSequenceLikelihood_==true) ) {
    levaluator_->getNniLk()->computeTreeLikelihood ();
    /*
     *        nniLk_->computeSubtreeLikelihoodPostfix(dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getRootNode());
     *        nniLk_->computeSubtreeLikelihoodPrefix(nniLk_->getTree().getRootNode());
     *        nniLk_->computeRootLikelihood();*/
  }
}




/******************************************************************************/

void DLGeneTreeLikelihood::computeReconciliationLikelihood()
{
  resetLossesAndDuplications(*spTree_, /*lossNumbers_, */lossExpectedNumbers_, /*duplicationNumbers_, */duplicationExpectedNumbers_);
  if (heuristicsLevel_>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    exit(-1);
    //    scenarioLikelihood_ = findMLReconciliation (&spTree_, &rootedTree_, seqSp_, lossNumbers_, lossExpectedNumbers_, duplicationNumbers_, duplicationExpectedNumbers_, MLindex_, branchNumbers_, _speciesIdLimitForRootPosition_, heuristicsLevel_, num0Lineages_, num1Lineages_, num2Lineages_, nodesToTryInNNISearch_); 
  }
  else {
    //    scenarioLikelihood_ = findMLReconciliationDR (&spTree_, &rootedTree_, seqSp_, spId_, lossExpectedNumbers_, duplicationExpectedNumbers_, MLindex_, num0Lineages_, num1Lineages_, num2Lineages_, nodesToTryInNNISearch_);
    scenarioLikelihood_ = findMLReconciliationDR (spTree_, rootedTree_, 
						  seqSp_, spId_, 
						  lossExpectedNumbers_, duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
    MLindex_ = tentativeMLindex_;
    num0Lineages_ = tentativeNum0Lineages_;
    num1Lineages_ = tentativeNum1Lineages_;
    num2Lineages_ = tentativeNum2Lineages_;
    nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
  }
}



/******************************************************************************/


void DLGeneTreeLikelihood::computeTreeLikelihood()
{
  if (considerSequenceLikelihood_ )
  {
    computeSequenceLikelihood();
  }
  computeReconciliationLikelihood();  
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.
 * 
 */
double DLGeneTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  //int nodeId = son->getId();
  //  std::cout<<"IN TESTNNI, nodeId "<< nodeId << "nodesToTryInNNISearch_.size() "<< nodesToTryInNNISearch_.size()<<std::endl;
  // std::cout << "before "<<TreeTemplateTools::treeToParenthesis (*tree_, true)<<std::endl;
  //If the NNI is around a branch where a duplication was found, 
  //or if we just try all branches because the starting gene trees are parsimonious in
  //numbers of DL.
  if ((nodesToTryInNNISearch_.count(nodeId)==1) /*|| DLStartingGeneTree_*/) {
    TreeTemplate<Node> * treeForNNI = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->clone();
    
    tentativeMLindex_ = MLindex_;
    /* tentativeLossNumbers_ = lossNumbers_;
     *         tentativeDuplicationNumbers_ = duplicationNumbers_;
     *         tentativeBranchNumbers_ = branchNumbers_;*/
    tentativeNum0Lineages_ = num0Lineages_;
    tentativeNum1Lineages_ = num1Lineages_;
    tentativeNum2Lineages_ =num2Lineages_;
    tentativeNodesToTryInNNISearch_.clear();
    
    //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
    //We use the same procedure as in doNNI !
    const Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->getNode(nodeId);
    
    
    if(!son->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'son' must not be the root node.", nodeId);
    const Node * parent = son->getFather();
    
    if(!parent->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent->getId());
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
    double ScenarioMLValue = 0;
    totalIterations_ = totalIterations_+1;
    
    /* //If we want to optimize the root or if we are at the first try
     *         if (rootOptimization_){
     *         if (heuristicsLevel_>0) {
     *         std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
     *         exit(-1);
     *         // ScenarioMLValue =  findMLReconciliation (&spTree_, treeForNNI, seqSp_, tentativeLossNumbers_, lossExpectedNumbers_, tentativeDuplicationNumbers_, duplicationExpectedNumbers_, tentativeMLindex_, tentativeBranchNumbers_, _speciesIdLimitForRootPosition_, heuristicsLevel_, tentativeNum0Lineages_, tentativeNum1Lineages_, tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
  }
  else {
    ScenarioMLValue =  findMLReconciliationDR (&spTree_, treeForNNI, seqSp_, spId_, lossExpectedNumbers_, duplicationExpectedNumbers_, tentativeMLindex_, tentativeNum0Lineages_, tentativeNum1Lineages_, tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
  }
  }
  else {
    if (heuristicsLevel_>0) {
      std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
      exit(-1);
      // ScenarioMLValue =  findMLReconciliation (&spTree_, treeForNNI, seqSp_, tentativeLossNumbers_, lossExpectedNumbers_, tentativeDuplicationNumbers_, duplicationExpectedNumbers_, tentativeMLindex_, tentativeBranchNumbers_, _speciesIdLimitForRootPosition_, heuristicsLevel_, tentativeNum0Lineages_, tentativeNum1Lineages_, tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
  }
  else {
    ScenarioMLValue =  findMLReconciliationDR (&spTree_, treeForNNI, seqSp_, spId_, lossExpectedNumbers_, duplicationExpectedNumbers_, tentativeMLindex_, tentativeNum0Lineages_, tentativeNum1Lineages_, tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
  }
  
  
  
  
  }*/
    ScenarioMLValue =  findMLReconciliationDR (spTree_, treeForNNI/*&rootedTree_*/, 
					       seqSp_, spId_, 
					       lossExpectedNumbers_, duplicationExpectedNumbers_, 
					       tentativeMLindex_, 
					       tentativeNum0Lineages_, tentativeNum1Lineages_, 
					       tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_, false); //false so that _tentativeNum*Lineages are not updated
    
    
    if (treeForNNI) delete treeForNNI;
    // std::cout<<"???WORTH computing the sequence likelihood "<< ScenarioMLValue<< " "<< scenarioLikelihood_<<std::endl;
    if (considerSequenceLikelihood_ ) 
    {
      if  (ScenarioMLValue >  scenarioLikelihood_) 
      { //If it is worth computing the sequence likelihood
	//Retrieving arrays of interest:
	// std::cout << "before nniLk_->testNNI(nodeId);"<<std::endl;
	double newLkMinusOldLk = levaluator_->getNniLk()->testNNI(nodeId);
	// std::cout << "after nniLk_->testNNI(nodeId);"<<std::endl;
	
	
	/*
	 *                const DRASDRTreeLikelihoodNodeData * parentData = & nniLk_->getLikelihoodData()->getNodeData(parent->getId());
	 *                const VVVdouble                    * sonArray   = & parentData->getLikelihoodArrayForNeighbor(son->getId());
	 *                std::vector<const Node *> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
	 *                unsigned int nbParentNeighbors = parentNeighbors.size();
	 *                std::vector<const VVVdouble *> parentArrays(nbParentNeighbors);
	 *                std::vector<const VVVdouble *> parentTProbs(nbParentNeighbors);
	 *                for(unsigned int k = 0; k < nbParentNeighbors; k++)
	 *                {
	 *                    const Node * n = parentNeighbors[k]; // This neighbor
	 *                    parentArrays[k] = & parentData->getLikelihoodArrayForNeighbor(n->getId()); 
	 *                    parentTProbs[k] = & pxy_[n->getId()];
      }
      const DRASDRTreeLikelihoodNodeData * grandFatherData = & nniLk_->getLikelihoodData()->getNodeData(grandFather->getId());
      const VVVdouble                    * uncleArray      = & grandFatherData->getLikelihoodArrayForNeighbor(uncle->getId()); 
      std::vector<const Node *> grandFatherNeighbors = TreeTemplateTools::getRemainingNeighbors(grandFather, parent, uncle);
      unsigned int nbGrandFatherNeighbors = grandFatherNeighbors.size();
      
      std::vector<const VVVdouble *> grandFatherArrays;
      std::vector<const VVVdouble *> grandFatherTProbs;
      for(unsigned int k = 0; k < nbGrandFatherNeighbors; k++)
      {
      const Node * n = grandFatherNeighbors[k]; // This neighbor
      if(grandFather->getFather() == NULL || n != grandFather->getFather())
      {
      grandFatherArrays.push_back(& grandFatherData->getLikelihoodArrayForNeighbor(n->getId())); 
      grandFatherTProbs.push_back(& pxy_[n->getId()]);
      }
      }
      //Compute array 1: grand father array
      VVVdouble array1 = *sonArray;
      resetLikelihoodArray(array1);
      grandFatherArrays.push_back(sonArray);
      grandFatherTProbs.push_back(& pxy_[son->getId()]);
      if(grandFather->hasFather())
      {
      computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, & grandFatherData->getLikelihoodArrayForNeighbor(grandFather->getFather()->getId()), & pxy_[grandFather->getId()], array1, nbGrandFatherNeighbors, nbDistinctSites_, nbClasses_, nbStates_, false); 
      
      }  
      else
      {
      computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, array1, nbGrandFatherNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false); 
      
      //This is the root node, we have to account for the ancestral frequencies:
      for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
      for(unsigned int j = 0; j < nbClasses_; j++)
      {
      for(unsigned int x = 0; x < nbStates_; x++)
	array1[i][j][x] *= rootFreqs_[x];
      }
      }
      }
      
      
      //Compute array 2: parent array
      VVVdouble array2 = *uncleArray;
      resetLikelihoodArray(array2);
      parentArrays.push_back(uncleArray);
      parentTProbs.push_back(& pxy_[uncle->getId()]);
      computeLikelihoodFromArrays(parentArrays, parentTProbs, array2, nbParentNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false); 
      //Initialize BranchLikelihood:
      brLikFunction_->initModel(nniLk->getModel(), nniLk->getRateDistribution());
      brLikFunction_->initLikelihoods(&array1, &array2);
      ParameterList parameters; 
      unsigned int pos = 0;
      while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
      if(pos == nodes_.size()) throw Exception("NNIHomogeneousTreeLikelihood::testNNI. Invalid node id.");
      Parameter brLen = getParameter("BrLen" + TextTools::toString(pos));
      brLen.setName("BrLen");
      parameters.addParameter(brLen);
      brLikFunction_->setParameters(parameters);
      //Re-estimate branch length:
      brentOptimizer_->setMessageHandler(NULL);
      brentOptimizer_->setFunction(brLikFunction_);
      brentOptimizer_->getStopCondition()->setTolerance(0.1);
      brentOptimizer_->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
      brentOptimizer_->init(parameters); 
      brentOptimizer_->optimize(); 
      brLenNNIValues_[nodeId] = brentOptimizer_->getParameters().getParameter("BrLen").getValue();
      brLikFunction_->resetLikelihoods(); //Array1 and Array2 will be destroyed after this function call.
      //We should not keep pointers towards them...
      */
	
	//Return the resulting likelihood:
	//   double temp = getValue() ; 
	
	//std::cout<<"temp "<< temp<< "; brLikFunction_->getValue()"<< brLikFunction_->getValue()<<std::endl;
	
	//      double tot = brLikFunction_->getValue() - ScenarioMLValue - temp; // -newsequencelogLk - (newscenariologLk) - (-currenttotallogLk); if <0, worth doing 
	//     if (tot<0) 
	double tot = - ScenarioMLValue + scenarioLikelihood_ + newLkMinusOldLk;
	// if (newLkMinusOldLk<0) 
	if (tot < 0)
	{
	  tentativeScenarioLikelihood_=ScenarioMLValue;
	}
	return tot;
	// std::cout << "after "<<TreeTemplateTools::treeToParenthesis (tree_, true)<<std::endl;
	
      }
      else 
      {
	tentativeMLindex_ = -1;
	return 1;
      }
    }
    else 
    {
      if  (ScenarioMLValue >  scenarioLikelihood_) 
      {
	tentativeScenarioLikelihood_=ScenarioMLValue;
	return (- ScenarioMLValue + scenarioLikelihood_);
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

void DLGeneTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  //std::cout<<"\t\t\tIN DONNI "<< std::endl;
  //std::cout << TreeTemplateTools::treeToParenthesis(*tree_, true) << std::endl;
  //Perform the topological move, the likelihood array will have to be recomputed...
  
  levaluator_->getNniLk()->doNNI(nodeId);
  
  /*
   *    Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getNode(nodeId);
   *    if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'son' must not be the root node.", nodeId);
   *    Node * parent = son->getFather();
   *    if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'parent' must not be the root node.", parent->getId());
   *    Node * grandFather = parent->getFather();
   *    //From here: Bifurcation assumed.
   *    //In case of multifurcation, an arbitrary uncle is chosen.
   *    //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
   *    unsigned int parentPosition = grandFather->getSonPosition(parent);
   *    Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
   *    //Swap nodes:
   *    parent->removeSon(son);
   *    grandFather->removeSon(uncle);
   *    parent->addSon(uncle);
   *    grandFather->addSon(son);
   *    unsigned int pos = 0;
   *    while(pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
   *    if(pos == nodes_.size()) throw Exception("ReconciliationTreeLikelihood::doNNI. Invalid node id.");
   *    
   *    //Julien Proposition :
   *    
   *    std::string name = "BrLen"+ TextTools::toString(pos);
   * 
   *    
   *    if (considerSequenceLikelihood_ ) 
   *    {
   *        if(brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
   *        {
   *            double length = brLenNNIValues_[nodeId];
   *            brLenParameters_.setParameterValue(name, length);
   *            getParameter_(name).setValue(length);
   *            parent->setDistanceToFather(length);
}
else std::cerr << "ERROR, branch not found: " << nodeId << std::endl;
try { brLenNNIParams_.addParameter(brLenParameters_.getParameter(name)); }
catch(ParameterException & ex)
{
std::cerr << "DEBUG:" << std::endl;
brLenNNIParams_.printParameters(std::cerr);
std::cerr << "DEBUG:" << name << std::endl;
}
brLenNNIParams_[brLenNNIParams_.size()-1].removeConstraint();
}
else 
{
double length = 0.1;
// brLenParameters_.setParameterValue(name, length);
// getParameter_(name).setValue(length);
parent->setDistanceToFather(length);
}

*/
  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  
  MLindex_ = tentativeMLindex_;
  
  nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
  
  // std::cout <<"tentativeScenarioLikelihood_: "<<tentativeScenarioLikelihood_<<std::endl;
  scenarioLikelihood_ = tentativeScenarioLikelihood_;// + _brLikFunction->getValue();
  
  //Now we need to update rootedTree_
  TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->clone();
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

std::vector <int> DLGeneTreeLikelihood::getDuplicationNumbers(){
  return duplicationNumbers_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getLossNumbers(){
  return lossNumbers_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getBranchNumbers(){
  return branchNumbers_;
}


/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get0LineagesNumbers() const {
  return num0Lineages_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get1LineagesNumbers() const {
  return num1Lineages_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get2LineagesNumbers() const {
  return num2Lineages_;
}

/*******************************************************************************/

void DLGeneTreeLikelihood::setExpectedNumbers(std::vector <double> duplicationProbabilities, std::vector <double> lossProbabilities){
  
  lossExpectedNumbers_ = lossProbabilities;
  duplicationExpectedNumbers_ = duplicationProbabilities; 
}

/*******************************************************************************/

int DLGeneTreeLikelihood::getRootNodeindex(){
  return MLindex_;
}

/*******************************************************************************/
/*
 * void DLGeneTreeLikelihood::resetSequenceLikelihood(){
 *    _sequenceLikelihood = UNLIKELY;
 * }
 */
/*******************************************************************************/

double DLGeneTreeLikelihood::getSequenceLikelihood() {
  return levaluator_->getLogLikelihood(); 
}
/*******************************************************************************/

void DLGeneTreeLikelihood::initialize() {
  levaluator_->getNniLk()->initialize();
  return;
}



/*******************************************************************************/
void DLGeneTreeLikelihood::print () const {
  std::cout << "Species tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getSpTree(), true)<<std::endl;
  std::cout << "Gene family rooted tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
  std::cout << "Gene family tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (levaluator_->getTree(), true)<<std::endl;
  std::cout << "0 lineage numbers"<<std::endl;
  VectorTools::print(get0LineagesNumbers());
  std::cout << "1 lineage numbers"<<std::endl;
  VectorTools::print(get1LineagesNumbers());
  std::cout << "2 lineages numbers"<<std::endl;
  VectorTools::print(get2LineagesNumbers());
  std::cout << "Expected numbers of losses"<<std::endl;
  VectorTools::print(lossExpectedNumbers_);
  std::cout << "Expected numbers of duplications"<<std::endl;
  VectorTools::print(duplicationExpectedNumbers_);
  std::cout << "Root index"<<std::endl;
  std::cout << MLindex_ <<std::endl;
  
  
}





/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * To do all this as fast as possible, we optimize only a few branch lengths on the SPR tree, 
 * and we use a simple recursion for that.
 * WORKS.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRsFast2 (map<string, string> params) {
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    //std::cout << "WE DONT DO SPRS"<<std::endl;
    computeReconciliationLikelihood();
    return;
  }
  //   std::cout << "WE DO SPRS; seq lk: "<< getSequenceLikelihood()  <<std::endl;
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  // if (getLogLikelihood()==UNLIKELY) 
  computeReconciliationLikelihood();
  double logL = getLogLikelihood();
  
  double bestlogL = logL;
  double candidateScenarioLk ;
  double bestSequenceLogL = getSequenceLikelihood();
  double bestScenarioLk = getScenarioLikelihood();
  std::cout << "LOGL: "<<logL << "ScenarioLK: "<< bestScenarioLk <<"; sequenceLK: "<<getSequenceLikelihood() << std::endl;
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
	
	buildVectorOfRegraftingNodesGeneTree(*spTree_, *rootedTree_, nodeForSPR, sprLimit_, nodeIdsToRegraft);
	
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
	  
	  //Compute the DL likelihood
	  candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
							 seqSp_, spId_, 
						  lossExpectedNumbers_, 
						  duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, 
						  tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, 
						  tentativeNodesToTryInNNISearch_, false); 
	  
	  if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the DL likelihood is not bad
	  {
	    
	    if (computeSequenceLikelihoodForSPR) {
	      
	      levaluator_->setAlternativeTree(treeForSPR);
            
	      
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
	    if (computeSequenceLikelihoodForSPR) {
	      bestSequenceLogL = levaluator_->getAlternativeLogLikelihood();

	    }
	    if (bestTree) {
	      delete bestTree;
	      bestTree = 0;
	    }
	    
	    bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->clone();
	    //Rooting bestTree as in TreeForSPR:
	    vector<Node*> rlkNodes = bestTree->getNodes();
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
		};
		bestTree->newOutGroup( rlkNodes[j] );
		//  std::cout << "FOUND"<<std::endl;
		break;
	      }
	    }
	    writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
	    
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
	  //  breadthFirstreNumber (*rootedTree_);
	  breadthFirstreNumberAndResetProperties (*rootedTree_);
	  
	  if (bestTree) {
	    delete bestTree;
	    bestTree = 0;
	  }
	  
	  
	}
	else 
	{
	  // logL = bestlogL;  
	  numIterationsWithoutImprovement++;
	  //  std::cout <<"\t\t\tSPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement << "; Total number  of iterations: "<< index << std::endl;
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
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * To do all this as fast as possible, we compute the sequence likelihood only for the most promising SPRs
 * given their reconciliation score.
 * We optimize all branch lengths on the SPR tree.
 * We use a simple recursion for that.
 * WORKS.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRsFast3 (map<string, string> params) {
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    //std::cout << "WE DONT DO SPRS"<<std::endl;
    computeReconciliationLikelihood();
    return;
  }
  //   std::cout << "WE DO SPRS; seq lk: "<< getSequenceLikelihood()  <<std::endl;
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  // if (getLogLikelihood()==UNLIKELY) 
  computeReconciliationLikelihood();
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
  std::map < double, TreeTemplate<Node> * >  treesToOptimizeSeqLk ;
  std::map < double, TreeTemplate<Node> * >::reverse_iterator it;
  double bestCurrentCandidateScenarioLk;
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
  
  
  while (numIterationsWithoutImprovement < rootedTree_->getNumberOfNodes() - 2)
  {
    
    annotateGeneTreeWithDuplicationEvents (*spTree_, 
					   *rootedTree_, 
					   rootedTree_->getRootNode(), 
					   seqSp_, spId_); 
    
    for (unsigned int nodeForSPR=rootedTree_->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
    {
      Node * n = rootedTree_->getNode(nodeForSPR);
      if (n->hasBranchProperty("L")) {
	numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
      }
      if ( numLoss != "0"  ) {
	treesToOptimizeSeqLk.clear();
	buildVectorOfRegraftingNodesGeneTree(*spTree_, *rootedTree_, nodeForSPR, sprLimit_, nodeIdsToRegraft);
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
	  
	  //Compute the DL likelihood
	  candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
							 seqSp_, spId_, 
						  lossExpectedNumbers_, 
						  duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, 
						  tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, 
						  tentativeNodesToTryInNNISearch_, false); 
	  while (treesToOptimizeSeqLk.count(candidateScenarioLk) > 0 )
	  {
	    candidateScenarioLk = candidateScenarioLk + NumConstants::SMALL();
	  }
	  treesToOptimizeSeqLk[candidateScenarioLk] = treeForSPR->clone();
	}
	if (treesToOptimizeSeqLk.size() > 0)
	  bestCurrentCandidateScenarioLk = (* (treesToOptimizeSeqLk.rbegin() ) ).first;
	for ( it=treesToOptimizeSeqLk.rbegin() ; it != treesToOptimizeSeqLk.rend(); it++ ) {
	  candidateScenarioLk = (*it).first;
	  if ( ( candidateScenarioLk > bestScenarioLk) && ( candidateScenarioLk < bestCurrentCandidateScenarioLk + 10 * NumConstants::SMALL()) )//If the likelihood is better than the current best, and one of the best ones to test
	  {
	    //We compute the sequence likelihood
	    if (computeSequenceLikelihoodForSPR) {
	      if (treeForSPR) 
	      {
		delete treeForSPR;
		treeForSPR = 0;
	      }
	      treeForSPR = (*it).second->clone();
	      
              levaluator_->setAlternativeTree(treeForSPR);
              
	      logL = candidateScenarioLk - levaluator_->getAlternativeLogLikelihood();
	      
	      
	      //If the candidate tree has a DL + sequence Lk better than the current best
	      if (logL - 0.01 > bestlogL) 
	      {
                // levaluator: since it the best tree, setting it as the current one
                levaluator_->acceptAlternativeTree();
                
		std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
		betterTree = true;
		bestlogL = logL;
		bestScenarioLk = candidateScenarioLk;
		if (computeSequenceLikelihoodForSPR) {
		  bestSequenceLogL = levaluator_->getLogLikelihood();
		}
		if (bestTree) {
		  delete bestTree;
		  bestTree = 0;
		}
		bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->clone();
		//Rooting bestTree as in TreeForSPR:
		vector<Node*> rlkNodes = bestTree->getNodes();
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
		    };
		    bestTree->newOutGroup( rlkNodes[j] );
		    //  std::cout << "FOUND"<<std::endl;
		    break;
		  }
		}
		writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
		break; //If we have found one topology better than the current one for seqlk+scenlk
	      }
	      else {
		//   copyContentsFrom(*bestTreeLogLk);
		//  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
	      }
	    }
	    else {
	      logL = candidateScenarioLk - bestSequenceLogL;
	    }
	    
	  }
	  else {
	    break; //No more likelihoods better than the current best
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
	  scenarioLikelihood_ = bestScenarioLk;
	  //  breadthFirstreNumber (*rootedTree_);
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
      else
      {
	numIterationsWithoutImprovement++;
      }
    }
  }
  
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
 * Implements the moves by Muffato and Crollius. Goes through the tree
 * while computing the number of genes ng and the number of species ns downstream, 
 * and in cases where a duplication is predicted but ng << 2 * ns, 
 * rearranges the tree to resolve the duplication.
 * First classical double-recursive tree traversals are used to find the root of the tree,
 * then another post-order tree traversal is done to compute ng and ns, 
 * and at the same time produce rearranged trees. Then the DL likelihood is computed
 * for all rearranged trees, and the best ones have their sequence likelihood computed.
 * 
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeMuffato (map<string, string> params) {
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    std::cout << "WE DONT DO SPRS"<<std::endl;
    computeReconciliationLikelihood();
    return;
  }
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree = false;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  
  computeReconciliationLikelihood();
  
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
  std::map < double, TreeTemplate<Node> * >  treesToOptimizeSeqLk ;
  std::map < double, TreeTemplate<Node> * >::reverse_iterator it;
  //double bestCurrentCandidateScenarioLk;
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
  
  
  while (1)
  {
    std::cout << "Starting Muffato optimization loop"<<std::endl;
    
    annotateGeneTreeWithScoredDuplicationEvents (*spTree_, 
						 *rootedTree_, 
						 rootedTree_->getRootNode(), 
						 seqSp_, spId_); 
    
    double editionThreshold = ApplicationTools::getDoubleParameter("muffato.edition.threshold", params, 0.3, "", false, false);
    
    treeForSPR = rootedTree_->clone();
    
    editDuplicationNodesMuffato(*spTree_, *treeForSPR, treeForSPR->getRootNode(), editionThreshold);
    //    std::cout <<"Here 13"<<std::endl;
    std::cout <<"Edited tree after Muffato move: \n"<< TreeTemplateTools::treeToParenthesis(*treeForSPR, false, "Score") << "\n" << std::endl;
    
    //Compute the DL likelihood
    candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
						   seqSp_, spId_, 
						   lossExpectedNumbers_, 
						   duplicationExpectedNumbers_, 
						   tentativeMLindex_, 
						   tentativeNum0Lineages_, 
						   tentativeNum1Lineages_, 
						   tentativeNum2Lineages_, 
						   tentativeNodesToTryInNNISearch_, false); 
    std::cout <<"Candidate scenario lk for Muffato rearranged tree: "<< candidateScenarioLk <<std::endl;
    
    
    //We compute the sequence likelihood
    if (computeSequenceLikelihoodForSPR) {
      
      unsigned int numEval = 0;
      
      levaluator_->setAlternativeTree(treeForSPR);
      
      logL = candidateScenarioLk - levaluator_->getAlternativeLogLikelihood();
      
      //If the candidate tree has a DL + sequence Lk better than the current best
      if (logL - 0.1 > bestlogL) 
      {
        // then we accept the tree
        levaluator_->acceptAlternativeTree();
        
	std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
	betterTree = true;
	bestlogL = logL;
	bestScenarioLk = candidateScenarioLk;
	if (computeSequenceLikelihoodForSPR) {
	  bestSequenceLogL = levaluator_->getAlternativeLogLikelihood();
	}
	
	if (bestTree) {
	  delete bestTree;
	  bestTree = 0;
	}
	bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(levaluator_->getTree()))->clone();
	//Rooting bestTree as in TreeForSPR:
	vector<Node*> rlkNodes = bestTree->getNodes();
	
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
	    };
	    bestTree->newOutGroup( rlkNodes[j] );
	    break;
	  }
	}
	writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
	
      }
      else {
	break;
	
      }
    }
    else {
      
      logL = candidateScenarioLk - bestSequenceLogL;
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
      scenarioLikelihood_ = bestScenarioLk;
      breadthFirstreNumberAndResetProperties (*rootedTree_);
      
      if (bestTree) {
	delete bestTree;
	bestTree = 0;
      }
      
    }
    else {
      break;
    }
  }
  
  rootedTree_->resetNodesId();
  
  
  //One more reconciliation, to update the "_num*Lineages" vectors.
  computeReconciliationLikelihood();
  
  Nhx *nhx = new Nhx();
  annotateGeneTreeWithDuplicationEvents (*spTree_, 
					 *rootedTree_, 
					 rootedTree_->getRootNode(), 
					 seqSp_, spId_); 
  cout << "Muffato reconciled tree: "<<endl;
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
void DLGeneTreeLikelihood::refineGeneTreeNNIs(map<string, string> params, unsigned int verbose ) {
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do NNIs
    computeReconciliationLikelihood();
    return;
  }
  bool test = true;
  do
  { 
    TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
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
	nniLk_->topologyChangeTested(   TopologyChangeEvent ());
	nniLk_->topologyChangeSuccessful(   TopologyChangeEvent ());
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
bool DLGeneTreeLikelihood::isSingleCopy() {
  vector<string> names = TreeTemplateTools::getLeavesNames(*(geneTreeWithSpNames_->getRootNode() ) );
  vector<string> uniqueNames = VectorTools::unique(names);
  if (uniqueNames.size() == names.size() ) {
    return true;
  }
  return false;
}





