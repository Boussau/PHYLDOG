//
// File: ReconciliationTreeLikelihood.cpp
// Created by: Bastien Boussau from NNIHomogeneousTreeLikelihood.cpp by Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "ReconciliationTreeLikelihood.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/AutoParameter.h>

// From the STL:
#include <iostream>
//using namespace std;
using namespace bpp;

/******************************************************************************/
ReconciliationTreeLikelihood::ReconciliationTreeLikelihood(
                                                           const Tree & tree,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
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
                                                           bool rootOptimization)
throw (Exception):
  NNIHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose), _spTree(0),_rootedTree(0),_seqSp(seqSp), _spId(spId)
{
  _spTree = spTree.clone();
  _rootedTree = rootedTree.clone();

  //_lossNumbers = lossNumbers;
  _lossProbabilities = lossProbabilities;
  //_duplicationNumbers = duplicationNumbers;
  _duplicationProbabilities = duplicationProbabilities;
  //_branchNumbers = branchNumbers;
  _num0Lineages=num0Lineages;
  _num1Lineages=num1Lineages;
  _num2Lineages=num2Lineages;
  _scenarioLikelihood = UNLIKELY;
  _sequenceLikelihood = UNLIKELY;
  _MLindex = MLindex;
  _rootOptimization = rootOptimization; 
  /*_tentativeDuplicationNumbers = duplicationNumbers;
  _tentativeLossNumbers = lossNumbers; 
  _tentativeBranchNumbers = branchNumbers;*/
  _tentativeNum0Lineages =num0Lineages;
  _tentativeNum1Lineages =num1Lineages; 
  _tentativeNum2Lineages =num2Lineages;
  _tentativeMLindex = MLindex;
  _totalIterations = 0;
  _counter = 0;
  _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
  _heuristicsLevel = heuristicsLevel;
  _optimizeSequenceLikelihood = true;
  _optimizeReconciliationLikelihood = true;
  // _listOfPreviousRoots = new std::vector <int> ();
}

/******************************************************************************/

ReconciliationTreeLikelihood::ReconciliationTreeLikelihood(
                                                           const Tree & tree,
                                                           const SiteContainer & data,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
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
                                                           bool rootOptimization)
throw (Exception):
NNIHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose), 
_spTree(0), _rootedTree(0), _seqSp (seqSp), _spId(spId)
{
  _spTree = spTree.clone();
  _rootedTree = rootedTree.clone();
  _lossProbabilities = lossProbabilities;
  _duplicationProbabilities = duplicationProbabilities; 
  _num0Lineages=num0Lineages;
  _num1Lineages=num1Lineages;
  _num2Lineages=num2Lineages;
  _scenarioLikelihood = UNLIKELY;
  _sequenceLikelihood = UNLIKELY;
  _MLindex = MLindex;
  _rootOptimization = rootOptimization; 
  _tentativeNum0Lineages =num0Lineages;
  _tentativeNum1Lineages =num1Lineages; 
  _tentativeNum2Lineages =num2Lineages;
  _tentativeMLindex = MLindex;
  _totalIterations = 0; 
  _counter = 0;
  _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
  _heuristicsLevel = heuristicsLevel;
  _optimizeSequenceLikelihood = true;
  _optimizeReconciliationLikelihood = true;
}

/******************************************************************************/

ReconciliationTreeLikelihood::ReconciliationTreeLikelihood(const ReconciliationTreeLikelihood & lik):
  NNIHomogeneousTreeLikelihood(lik) , _spTree(0), _rootedTree(0), _seqSp (lik._seqSp), _spId(lik._spId)
{
  _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone()) ;
  _rootedTree = dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone()) ;
  _lossProbabilities = lik._lossProbabilities;
  _duplicationProbabilities = lik._duplicationProbabilities; 
  _num0Lineages=lik._num0Lineages;
  _num1Lineages=lik._num1Lineages;
  _num2Lineages=lik._num2Lineages;
  _scenarioLikelihood = lik._scenarioLikelihood;
  _sequenceLikelihood = lik._sequenceLikelihood;
  _MLindex = lik._MLindex;
  _rootOptimization = lik._rootOptimization; 
  _tentativeNum0Lineages =lik._tentativeNum0Lineages;
  _tentativeNum1Lineages =lik._tentativeNum1Lineages;
  _tentativeNum2Lineages =lik._tentativeNum2Lineages;
  _tentativeMLindex = lik._MLindex;
  _totalIterations = lik._totalIterations;
  _counter = lik._counter;
  _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
  _heuristicsLevel = lik._heuristicsLevel;
  _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
  _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
  _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
  _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;

}

/******************************************************************************/

ReconciliationTreeLikelihood & ReconciliationTreeLikelihood::operator=(const ReconciliationTreeLikelihood & lik)
{
  NNIHomogeneousTreeLikelihood::operator=(lik);
  if (_spTree) delete _spTree;
  _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone());
  if (_rootedTree) delete _rootedTree;
  _rootedTree= dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone());
  _spId = lik._spId;
  _lossProbabilities = lik._lossProbabilities;
  _duplicationProbabilities = lik._duplicationProbabilities;
  _num0Lineages=lik._num0Lineages;
  _num1Lineages=lik._num1Lineages;
  _num2Lineages=lik._num2Lineages;
  _scenarioLikelihood = lik._scenarioLikelihood;
  _sequenceLikelihood = lik._sequenceLikelihood;
  _MLindex = lik._MLindex;
  _rootOptimization = lik._rootOptimization;
  _tentativeNum0Lineages =lik._tentativeNum0Lineages;
  _tentativeNum1Lineages =lik._tentativeNum1Lineages;
  _tentativeNum2Lineages =lik._tentativeNum2Lineages;
  _tentativeMLindex = lik._MLindex;
  _totalIterations = lik._totalIterations;
  _counter = lik._counter;
  _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
  _heuristicsLevel = lik._heuristicsLevel;
  _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
  _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
  _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
  _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;

  return *this;
}




/******************************************************************************/

ReconciliationTreeLikelihood::~ReconciliationTreeLikelihood()
{
    if (_spTree) delete _spTree;
    if (_rootedTree) delete _rootedTree; 
}

/******************************************************************************/

/*ReconciliationTreeLikelihood::~ReconciliationTreeLikelihood()
{
  NNIHomogeneousTreeLikelihood::~NNIHomogeneousTreeLikelihood();
  }*/

void ReconciliationTreeLikelihood::initParameters()
{
 // std::cout << "in initParameters"<<std::endl;
  NNIHomogeneousTreeLikelihood::initParameters();

  if (_heuristicsLevel>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    exit(-1);
//    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
  }
  else {
    _scenarioLikelihood = findMLReconciliationDR (_spTree, _rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _MLindex, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
  }
  _MLindex = -1;
  // std::cout << "in ReconciliationTreeLikelihood::initParameters : _num0Lineages, _num1Lineages, _num2Lineages : "<< TextTools::toString(VectorTools::sum(_num0Lineages))<<" "<<TextTools::toString(VectorTools::sum(_num1Lineages))<<" "<< TextTools::toString(VectorTools::sum(_num2Lineages))<<std::endl;
  
  // std::cout <<"INITIAL _scenarioLikelihood "<<_scenarioLikelihood<<std::endl;
}



/******************************************************************************/

void ReconciliationTreeLikelihood::resetMLindex() {
  _MLindex = -1;
}



/******************************************************************************/

/*
double ReconciliationTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  l*=exp(_scenarioLikelihood);
  return l;
}
*/
/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double ReconciliationTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
 // _sequenceLikelihood=UNLIKELY;
//  if (_sequenceLikelihood ==UNLIKELY) {

  //TEST
 /*
  if ((_sequenceLikelihood == UNLIKELY)||(_optimizeSequenceLikelihood==true)) {
    Vdouble * lik = & likelihoodData_->getRootRateSiteLikelihoodArray();
    const std::vector<unsigned int> * w = & likelihoodData_->getWeights();
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
        ll += (* w)[i] * log((* lik)[i]);
        //std::cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << std::endl;
      }
    _sequenceLikelihood = ll;
   // std::cout <<"getLogLikelihood _sequenceLikelihood "<< _sequenceLikelihood<<std::endl;
  }
  */
  //END TEST
  
 // ll = getValue();

  /* else {
    std::cout <<"NOT RECOMPUTING IN GETLOGLIKELIHOOD : "<< _sequenceLikelihood<<std::endl; 
  }*/
  //Now we add the scenario likelihood
  //  std::cout << "COMPUTING LOGLIKELIHOOD :" << ll  << " SCENARIO LOGLIKELIHOOD :" << _scenarioLikelihood <<" TOTAL : "<< ll + _scenarioLikelihood <<std::endl;

  //TEST 16 02 2010
 // DRHomogeneousTreeLikelihood::getLogLikelihood();
//  setMinuslogLikelihood_ (_sequenceLikelihood);
  ll = _sequenceLikelihood + _scenarioLikelihood;
 // ll = _sequenceLikelihood ;
  return ll;
}
/******************************************************************************/
//returns -loglikelihood

double ReconciliationTreeLikelihood::getValue() const  
throw (Exception)
{
   {
    if(!isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
   // return (-getLogLikelihood());
     //TEST 16 02 2010
    // std::cout<<"\t\t\t_sequenceLikelihood: "<<_sequenceLikelihood<< " _scenarioLikelihood: "<<_scenarioLikelihood<<std::endl;
     return (- _sequenceLikelihood - _scenarioLikelihood);
     //return (minusLogLik_ - _scenarioLikelihood);
     //return (-minusLogLik_);
  }
}



/******************************************************************************/


void ReconciliationTreeLikelihood::fireParameterChanged(const ParameterList & params)
{

  applyParameters();
  
  if(rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
     || model_->getParameters().getCommonParametersWith(params).size() > 0)
    {
      //Rate parameter changed, need to recompute all probs:
      computeAllTransitionProbabilities();
    }
  else if(params.size() > 0)
    {
      //We may save some computations:
      for(unsigned int i = 0; i < params.size(); i++)
        {
          std::string s = params[i].getName();
          if(s.substr(0,5) == "BrLen")
            {
              //Branch length parameter:
              computeTransitionProbabilitiesForNode(nodes_[TextTools::to<unsigned int>(s.substr(5))]);
            }
        }
    }
  
  computeSequenceLikelihood();

 // computeTreeLikelihood();
  if(computeFirstOrderDerivatives_)
    {
      computeTreeDLikelihoods();  
    }
  if(computeSecondOrderDerivatives_)
    {
      computeTreeD2Likelihoods();
    }
  double ll =0.0;
  //minusLogLik_ = - getLogLikelihood();
  Vdouble * lik = & getLikelihoodData()->getRootRateSiteLikelihoodArray();
  const std::vector<unsigned int> * w = & getLikelihoodData()->getWeights();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      ll += (* w)[i] * log((* lik)[i]);
      //std::cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << std::endl;
    }
  _sequenceLikelihood = ll;

  minusLogLik_ = - _sequenceLikelihood ;   
  //setMinuslogLikelihood_(- _sequenceLikelihood);
  //If we need to update the reconciliation likelihood
  if (_optimizeReconciliationLikelihood) {
    computeReconciliationLikelihood();
  }

  
}


/******************************************************************************/


void ReconciliationTreeLikelihood::computeSequenceLikelihood()
{
  if ((_sequenceLikelihood == UNLIKELY)||(_optimizeSequenceLikelihood==true)) {
    computeSubtreeLikelihoodPostfix(tree_->getRootNode());
    computeSubtreeLikelihoodPrefix(tree_->getRootNode());
    computeRootLikelihood();
    

  }

  
}
  



/******************************************************************************/

void ReconciliationTreeLikelihood::computeReconciliationLikelihood()
{
  
  resetLossesAndDuplications(*_spTree, /*_lossNumbers, */_lossProbabilities, /*_duplicationNumbers, */_duplicationProbabilities);
  if (_heuristicsLevel>0) {
    std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
    exit(-1);
    //    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
  }
  else {

    //    _scenarioLikelihood = findMLReconciliationDR (&_spTree, &_rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _MLindex, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
    _scenarioLikelihood = findMLReconciliationDR (_spTree, _rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
    _MLindex = _tentativeMLindex;
    _num0Lineages = _tentativeNum0Lineages;
    _num1Lineages = _tentativeNum1Lineages;
    _num2Lineages = _tentativeNum2Lineages;
    _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;


  }
  

}



/******************************************************************************/


void ReconciliationTreeLikelihood::computeTreeLikelihood()
{

  computeSequenceLikelihood();
  computeReconciliationLikelihood();  
  
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.

*/
double ReconciliationTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  //int nodeId = son->getId();
//  std::cout<<"IN TESTNNI, nodeId "<< nodeId << "_nodesToTryInNNISearch.size() "<< _nodesToTryInNNISearch.size()<<std::endl;
 // std::cout << "before "<<TreeTools::treeToParenthesis (*tree_, true)<<std::endl;
  if (_nodesToTryInNNISearch.count(nodeId)==1) {
   // std::cout << "before "<<TreeTools::treeToParenthesis (tree_, true)<<std::endl;
   // std::cout<<"Node To Try"<<std::endl;
    TreeTemplate<Node> * treeForNNI = tree_->clone();
    
    _tentativeMLindex = _MLindex;
   /* _tentativeLossNumbers = _lossNumbers;
    _tentativeDuplicationNumbers = _duplicationNumbers;
    _tentativeBranchNumbers = _branchNumbers;*/
    _tentativeNum0Lineages = _num0Lineages;
    _tentativeNum1Lineages = _num1Lineages;
    _tentativeNum2Lineages =_num2Lineages;
    _tentativeNodesToTryInNNISearch.clear();
    
    //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
    //We use the same procedure as in doNNI !
        const Node * son    = tree_->getNode(nodeId);


    if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", nodeId);
    const Node * parent = son->getFather();

    if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent->getId());
    const Node * grandFather = parent->getFather();

    //From here: Bifurcation assumed.
    //In case of multifurcation, an arbitrary uncle is chosen.
    //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
    unsigned int parentPosition = grandFather->getSonPosition(parent);
    const Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);

    Node * sonForNNI    = treeForNNI->getNode(nodeId);
    Node * parentForNNI = sonForNNI->getFather();
    Node * grandFatherForNNI = parentForNNI->getFather();
    Node * uncleForNNI = grandFatherForNNI->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
    parentForNNI->removeSon(sonForNNI);
    grandFatherForNNI->removeSon(uncleForNNI);
    parentForNNI->addSon(uncleForNNI);
    grandFatherForNNI->addSon(sonForNNI);

    //Now we root the tree sent to findMLReconciliation as in _rootedTree
    int id = treeForNNI->getRootNode()->getId();
    if(TreeTemplateTools::hasNodeWithId(*(_rootedTree->getRootNode()->getSon(0)),id)) {
      treeForNNI->newOutGroup(_rootedTree->getRootNode()->getSon(1)->getId());
    }
    else {
      treeForNNI->newOutGroup(_rootedTree->getRootNode()->getSon(0)->getId());
    }
    double ScenarioMLValue = 0;
    _totalIterations = _totalIterations+1;

    /* //If we want to optimize the root or if we are at the first try
    if (_rootOptimization){
      if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        exit(-1);
       // ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      else {
        ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
    }
    else {
      if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        exit(-1);
       // ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      else {
        ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      
      
      
      
    }*/
    ScenarioMLValue =  findMLReconciliationDR (_spTree, treeForNNI/*&_rootedTree*/, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 

    
    delete treeForNNI;
  //  std::cout<<"???WORTH computing the sequence likelihood "<< ScenarioMLValue<< " "<< _scenarioLikelihood<<std::endl;

    if (ScenarioMLValue >  _scenarioLikelihood) { //If it is worth computing the sequence likelihood
      //Retrieving arrays of interest:

      const DRASDRTreeLikelihoodNodeData * parentData = & getLikelihoodData()->getNodeData(parent->getId());
      const VVVdouble                    * sonArray   = & parentData->getLikelihoodArrayForNeighbor(son->getId());
      std::vector<const Node *> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
      unsigned int nbParentNeighbors = parentNeighbors.size();
      std::vector<const VVVdouble *> parentArrays(nbParentNeighbors);
      std::vector<const VVVdouble *> parentTProbs(nbParentNeighbors);
      for(unsigned int k = 0; k < nbParentNeighbors; k++)
        {
          const Node * n = parentNeighbors[k]; // This neighbor
          parentArrays[k] = & parentData->getLikelihoodArrayForNeighbor(n->getId()); 
          parentTProbs[k] = & pxy_[n->getId()];
        }
      
      const DRASDRTreeLikelihoodNodeData * grandFatherData = & getLikelihoodData()->getNodeData(grandFather->getId());
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
      brLikFunction_->initModel(model_, rateDistribution_);
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
            
      
      //Return the resulting likelihood:
      double temp = getValue() ; 
      
      //std::cout<<"temp "<< temp<< "; brLikFunction_->getValue()"<< brLikFunction_->getValue()<<std::endl;
      
      double tot = brLikFunction_->getValue() - ScenarioMLValue - temp; //positive + positive + negative; if <0, worth doing 
      if (tot<0) {
        _tentativeScenarioLikelihood=ScenarioMLValue;
      }
      return tot;
     // std::cout << "after "<<TreeTools::treeToParenthesis (tree_, true)<<std::endl;

    }
    else {
      _tentativeMLindex = -1;
      return 1;
    }
    
  }
  else {
    _tentativeMLindex = -1;
    return 1;
  }
}
    
/*******************************************************************************/

void ReconciliationTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
 // std::cout<<"\t\t\tIN DONNI "<< std::endl;

  //Perform the topological move, the likelihood array will have to be recomputed...
  Node * son    = tree_->getNode(nodeId);
  if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'son' must not be the root node.", nodeId);
  Node * parent = son->getFather();
  if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'parent' must not be the root node.", parent->getId());
  Node * grandFather = parent->getFather();
  //From here: Bifurcation assumed.
  //In case of multifurcation, an arbitrary uncle is chosen.
  //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  unsigned int parentPosition = grandFather->getSonPosition(parent);
  Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
  //Swap nodes:
  parent->removeSon(son);
  grandFather->removeSon(uncle);
  parent->addSon(uncle);
  grandFather->addSon(son);
  unsigned int pos = 0;
  while(pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
  if(pos == nodes_.size()) throw Exception("ReconciliationTreeLikelihood::doNNI. Invalid node id.");
  
  //Julien Proposition :
  
  std::string name = "BrLen"+ TextTools::toString(pos);
  /*  for(unsigned int i = 0; i < nodes_.size(); i++)
    if(nodes_[i]->getId() == parent->getId())
      {
      name += TextTools::toString(i);
      break;
      }
  */



  if(brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
  {
    double length = brLenNNIValues_[nodeId];
    brLenParameters_.setParameterValue(name, length);
    /*  Parameter p = getParameter(name);
    p.setValue(length);
    parent->setDistanceToFather(length);*/
    getParameter_(name).setValue(length);
    parent->setDistanceToFather(length);
  }
  else std::cerr << "ERROR, branch not found: " << nodeId << std::endl;
  try { brLenNNIParams_.addParameter(brLenParameters_.getParameter(name)); }
  catch(ParameterException & ex)
  {
    std::cerr << "DEBUG:" << std::endl;
    brLenNNIParams_.printParameters(std::cerr);
    std::cerr << "DEBUG:" << name << std::endl;
  }


  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  brLenNNIParams_[brLenNNIParams_.size()-1].removeConstraint();
  _MLindex = _tentativeMLindex;
  _duplicationNumbers = _tentativeDuplicationNumbers;
  _lossNumbers = _tentativeLossNumbers;
  _branchNumbers = _tentativeBranchNumbers;

  _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
  

  _scenarioLikelihood = _tentativeScenarioLikelihood;// + _brLikFunction->getValue();

  //Now we need to update _rootedTree
  TreeTemplate<Node> * tree = tree_->clone();
  //First we root this temporary tree as in _rootedTree (same lines as in testNNI)
  int id = tree->getRootNode()->getId();
  if(TreeTemplateTools::hasNodeWithId(*(_rootedTree->getRootNode()->getSon(0)),id)) {
    tree->newOutGroup(_rootedTree->getRootNode()->getSon(1)->getId());
  }
  else {
    tree->newOutGroup(_rootedTree->getRootNode()->getSon(0)->getId());
  }
  //Then we root this tree according to MLindex
  tree->newOutGroup(_MLindex);
  //We update _rootedTree
  if (_rootedTree) delete _rootedTree;
  _rootedTree = tree->clone();
  delete tree;
  //we need to update the sequence likelihood
  OptimizeSequenceLikelihood(true);
  OptimizeReconciliationLikelihood(true);

}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::getDuplicationNumbers(){
  return _duplicationNumbers;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::getLossNumbers(){
return _lossNumbers;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::getBranchNumbers(){
return _branchNumbers;
}


/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::get0LineagesNumbers() const {
  return _num0Lineages;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::get1LineagesNumbers() const {
  return _num1Lineages;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::get2LineagesNumbers() const {
  return _num2Lineages;
}

/*******************************************************************************/

void ReconciliationTreeLikelihood::setProbabilities(std::vector <double> duplicationProbabilities, std::vector <double> lossProbabilities){
  
  _lossProbabilities = lossProbabilities;
  _duplicationProbabilities = duplicationProbabilities; 
}

/*******************************************************************************/

int ReconciliationTreeLikelihood::getRootNodeindex(){
  return _MLindex;
}

/*******************************************************************************/

void ReconciliationTreeLikelihood::resetSequenceLikelihood(){
   _sequenceLikelihood = UNLIKELY;
}

/*******************************************************************************/

double ReconciliationTreeLikelihood::getSequenceLikelihood() {
  return _sequenceLikelihood; 
}

/*******************************************************************************/
void ReconciliationTreeLikelihood::print () const {
  std::cout << "Species tree:"<<std::endl;
  std::cout << TreeTools::treeToParenthesis (getSpTree(), true)<<std::endl;
  std::cout << "Gene family rooted tree:"<<std::endl;
  std::cout << TreeTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
  std::cout << "Gene family tree:"<<std::endl;
  std::cout << TreeTools::treeToParenthesis (getTree(), true)<<std::endl;
  std::cout << "0 lineage numbers"<<std::endl;
  VectorTools::print(get0LineagesNumbers());
  std::cout << "1 lineage numbers"<<std::endl;
  VectorTools::print(get1LineagesNumbers());
  std::cout << "2 lineages numbers"<<std::endl;
  VectorTools::print(get2LineagesNumbers());
  std::cout << "Expected numbers of losses"<<std::endl;
  VectorTools::print(_lossProbabilities);
  std::cout << "Expected numbers of duplications"<<std::endl;
  VectorTools::print(_duplicationProbabilities);
  std::cout << "Root index"<<std::endl;
  std::cout << _MLindex <<std::endl;


}




