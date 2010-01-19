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
                                                           const TreeTemplate<Node> & tree,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
                                                           const std::map <std::string, std::string> seqSp,
                                                           std::map <std::string,int> spId,
                                                           std::vector <int> & lossNumbers, 
                                                           std::vector <double> & lossProbabilities, 
                                                           std::vector <int> & duplicationNumbers, 
                                                           std::vector <double> & duplicationProbabilities,
                                                           std::vector <int> & branchNumbers,
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
  NNIHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose), _spTree(spTree),_rootedTree(rootedTree),_seqSp(seqSp), _spId(spId)
{
  _lossNumbers = lossNumbers;
  _lossProbabilities = lossProbabilities;
  _duplicationNumbers = duplicationNumbers;
  _duplicationProbabilities = duplicationProbabilities;
  _branchNumbers = branchNumbers;
  _num0Lineages=num0Lineages;
  _num1Lineages=num1Lineages;
  _num2Lineages=num2Lineages;
  _scenarioLikelihood = UNLIKELY;
  _sequenceLikelihood = UNLIKELY;
  _MLindex = MLindex;
  _rootOptimization = rootOptimization; 
  _tentativeDuplicationNumbers = duplicationNumbers;
  _tentativeLossNumbers = lossNumbers; 
  _tentativeBranchNumbers = branchNumbers;
  _tentativeNum0Lineages =num0Lineages;
  _tentativeNum1Lineages =num1Lineages; 
  _tentativeNum2Lineages =num2Lineages;
  _tentativeMLindex = MLindex;
  _totalIterations = 0;
  _counter = 0;
  _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
  _heuristicsLevel = heuristicsLevel;
  _optimizeSequenceLikelihood = true;
  // _listOfPreviousRoots = new std::vector <int> ();
}

/******************************************************************************/

ReconciliationTreeLikelihood::ReconciliationTreeLikelihood(
                                                           const TreeTemplate<Node> & tree,
                                                           const SiteContainer & data,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
                                                           const std::map <std::string, std::string> seqSp,
                                                           std::map <std::string,int> spId,
                                                           std::vector <int> & lossNumbers, 
                                                           std::vector <double> & lossProbabilities, 
                                                           std::vector <int> & duplicationNumbers, 
                                                           std::vector <double> & duplicationProbabilities, 
                                                           std::vector <int> & branchNumbers, 
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
NNIHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose), _spTree(spTree),_rootedTree(rootedTree), _seqSp (seqSp), _spId(spId)
{
  _lossNumbers = lossNumbers;
  _lossProbabilities = lossProbabilities;
  _duplicationNumbers = duplicationNumbers;
  _duplicationProbabilities = duplicationProbabilities; 
  _branchNumbers = branchNumbers;
  _num0Lineages=num0Lineages;
  _num1Lineages=num1Lineages;
  _num2Lineages=num2Lineages;
  _scenarioLikelihood = UNLIKELY;
  _sequenceLikelihood = UNLIKELY;
  _MLindex = MLindex;
  _rootOptimization = rootOptimization; 
  _tentativeDuplicationNumbers = duplicationNumbers;
  _tentativeLossNumbers = lossNumbers; 
  _tentativeBranchNumbers = branchNumbers; 
  _tentativeNum0Lineages =num0Lineages;
  _tentativeNum1Lineages =num1Lineages; 
  _tentativeNum2Lineages =num2Lineages;
  _tentativeMLindex = MLindex;
  _totalIterations = 0; 
  _counter = 0;
  _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
  _heuristicsLevel = heuristicsLevel;
  _optimizeSequenceLikelihood = true;
  // _listOfPreviousRoots = new std::vector <int> ();
}

/******************************************************************************/

ReconciliationTreeLikelihood::ReconciliationTreeLikelihood(const ReconciliationTreeLikelihood & lik):
  NNIHomogeneousTreeLikelihood(lik), _spTree(lik._spTree), _rootedTree(lik._rootedTree), _seqSp (lik._seqSp), _spId(lik._spId)
{
  _lossNumbers = lik._lossNumbers;
  _duplicationNumbers = lik._duplicationNumbers;
  _lossProbabilities = lik._lossProbabilities;
  _duplicationProbabilities = lik._duplicationProbabilities; 
  _branchNumbers = lik._branchNumbers;
  _num0Lineages=lik._num0Lineages;
  _num1Lineages=lik._num1Lineages;
  _num2Lineages=lik._num2Lineages;
  _scenarioLikelihood = lik._scenarioLikelihood;
  _sequenceLikelihood = lik._sequenceLikelihood;
  _MLindex = lik._MLindex;
  _rootOptimization = lik._rootOptimization; 
  _tentativeDuplicationNumbers = lik._duplicationNumbers;
  _tentativeLossNumbers = lik._lossNumbers; 
  _tentativeBranchNumbers = lik._branchNumbers;
  _tentativeNum0Lineages =lik._tentativeNum0Lineages;
  _tentativeNum1Lineages =lik._tentativeNum1Lineages;
  _tentativeNum2Lineages =lik._tentativeNum2Lineages;
  _tentativeMLindex = lik._MLindex;
  _totalIterations = lik._totalIterations;
  _counter = lik._counter;
 // _listOfPreviousRoots = lik._listOfPreviousRoots;
  _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
  _heuristicsLevel = lik._heuristicsLevel;
  _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
  _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
  _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;

}

/******************************************************************************/

ReconciliationTreeLikelihood & ReconciliationTreeLikelihood::operator=(const ReconciliationTreeLikelihood & lik)
{
  NNIHomogeneousTreeLikelihood::operator=(lik);
  _spTree = lik._spTree;
  _rootedTree=lik._rootedTree;
  //_seqSp = lik._seqSp; No need for that a priori
  _spId = lik._spId;
  _lossNumbers = lik._lossNumbers;
  _duplicationNumbers = lik._duplicationNumbers;
  _lossProbabilities = lik._lossProbabilities;
  _duplicationProbabilities = lik._duplicationProbabilities;
  _branchNumbers = lik._branchNumbers;  
  _num0Lineages=lik._num0Lineages;
  _num1Lineages=lik._num1Lineages;
  _num2Lineages=lik._num2Lineages;
  _scenarioLikelihood = lik._scenarioLikelihood;
  _sequenceLikelihood = lik._sequenceLikelihood;
  _MLindex = lik._MLindex;
  _rootOptimization = lik._rootOptimization;
  _tentativeDuplicationNumbers = lik._duplicationNumbers;
  _tentativeLossNumbers = lik._lossNumbers; 
  _tentativeBranchNumbers = lik._branchNumbers; 
  _tentativeNum0Lineages =lik._tentativeNum0Lineages;
  _tentativeNum1Lineages =lik._tentativeNum1Lineages;
  _tentativeNum2Lineages =lik._tentativeNum2Lineages;
  _tentativeMLindex = lik._MLindex;
  _totalIterations = lik._totalIterations;
  _counter = lik._counter;
//  _listOfPreviousRoots = lik._listOfPreviousRoots;
  _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
  _heuristicsLevel = lik._heuristicsLevel;
  _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
  _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
  _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
  return *this;
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
 //   std::cout <<"CLOCKBEFOREfindMLReconciliation "<<clock()<<std::endl;
    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
 //   std::cout <<"CLOCKAFTERfindMLReconciliation "<<clock()<<std::endl;
  }
  else {
 //   std::cout <<"CLOCKBEFOREfindMLReconciliationDR "<<clock()<<std::endl;
    _scenarioLikelihood = findMLReconciliationDR (&_spTree, &_rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _MLindex, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
 //   std::cout <<"CLOCKAFTERfindMLReconciliationDR "<<clock()<<std::endl;
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


double ReconciliationTreeLikelihood::getLikelihood() const
{
//  std::cout <<"IN getLikelihood"<<std::endl;
  double l = 1.;
 // if (_sequenceLikelihood ==UNLIKELY) {
  //TEST
/*
  Vdouble * lik = & likelihoodData_->getRootRateSiteLikelihoodArray();
    const std::vector<unsigned int> * w = & likelihoodData_->getWeights();
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
        l *= std::pow((*lik)[i], (int)(* w)[i]);
      }
      _sequenceLikelihood = log(l);
 */
  //END TEST
  //  std::cout <<"_sequenceLikelihood "<< _sequenceLikelihood<<std::endl;
  // }
 /* else {
    std::cout <<"NOT RECOMPUTING IN GETLIKELIHOOD"<<std::endl; 
  }*/
//  l = exp( getValue());
  
  l*=exp(_scenarioLikelihood);
  return l;
}

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

  ll = _sequenceLikelihood + _scenarioLikelihood;
 // ll = _sequenceLikelihood ;
  return ll;
}
/******************************************************************************/

double ReconciliationTreeLikelihood::getValue() const  
throw (Exception)
{
   {
    if(!isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
   // return (-getLogLikelihood());
     return (minusLogLik_ - _scenarioLikelihood);
     //return (-minusLogLik_);
  }
}



/******************************************************************************/


void ReconciliationTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
 // std::cout << "!!!!!!!!!!!!!!!!!!!!!!IN FIRE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
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
  Vdouble * lik = & likelihoodData_->getRootRateSiteLikelihoodArray();
  const std::vector<unsigned int> * w = & likelihoodData_->getWeights();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      ll += (* w)[i] * log((* lik)[i]);
      //std::cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << std::endl;
    }
  _sequenceLikelihood = ll;
  minusLogLik_ = -_sequenceLikelihood;
/*cout.precision(10); 
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!IN FIRE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< minusLogLik_<<std::endl;
  */
}


/******************************************************************************/


void ReconciliationTreeLikelihood::computeSequenceLikelihood()
{
  if ((_sequenceLikelihood == UNLIKELY)||(_optimizeSequenceLikelihood==true)) {
    computeSubtreeLikelihoodPostfix(tree_->getRootNode());
    computeSubtreeLikelihoodPrefix(tree_->getRootNode());
    computeRootLikelihood();
    
    //The following line is to update _sequenceLikelihood
    /* resetSequenceLikelihood();
     getLogLikelihood();*/
    // std::cout <<"IN computeTreeLikelihood : "<<_sequenceLikelihood<<std::endl;
  }
  // std::cout <<"After computeRootLikelihood : "<<getValue()<<std::endl;
  
}
  



/******************************************************************************/


void ReconciliationTreeLikelihood::computeTreeLikelihood()
{
  //  std::cout <<"IN computeTreeLikelihood : "<<_sequenceLikelihood<<std::endl;
  if ((_sequenceLikelihood == UNLIKELY)||(_optimizeSequenceLikelihood==true)) {
    computeSubtreeLikelihoodPostfix(tree_->getRootNode());
    computeSubtreeLikelihoodPrefix(tree_->getRootNode());
    computeRootLikelihood();
    
    //The following line is to update _sequenceLikelihood
/* resetSequenceLikelihood();
    getLogLikelihood();*/
  
  }
  // std::cout <<"After computeRootLikelihood : "<<getValue()<<std::endl;
 
  
  resetLossesAndDuplications(_spTree, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities);
  if (_heuristicsLevel>0) {
    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
  }
  else {
       _scenarioLikelihood = findMLReconciliationDR (&_spTree, &_rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _MLindex, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
  }
  
  
  
//  _scenarioLikelihood = findMLReconciliation(&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch);
  // std::cout <<"After findMLReconciliation : "<<_scenarioLikelihood<<std::endl;
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.

*/
double ReconciliationTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  // _tree->resetNodesId();
  /*
    std::vector <int> nodeIds=_tree->getNodesId ();
    for (int i =0; i<nodeIds.size() ; i++ ) {
    std::cout << "nodeIds[i] "<<nodeIds[i]<<std::endl;
    }
    
    
    for (int i =0; i<84 ; i++ ) {
    std::cout << "Test : "<<i<<std::endl;
    Parameter brLen = *_parameters.getParameter("BrLen" + TextTools::toString(i));
    std::cout <<"test passed"<<std::endl;
    }
  */
/*  bool toTry=false;
  for (int i=0; i<_nodesToTryInNNISearch.size(); i++) {
    if (nodeId==_nodesToTryInNNISearch[i]) {
      toTry=true;
      break;
    }
  }*/
  
  if (_nodesToTryInNNISearch.count(nodeId)==1) {
  
  /*  if (_optimizeSequenceLikelihood) {
      OptimizeSequenceLikelihood(false);
    }
    */
    
    TreeTemplate<Node> * treeForNNI = tree_->clone();
    
    _tentativeMLindex = _MLindex;
    _tentativeLossNumbers = _lossNumbers;
    _tentativeDuplicationNumbers = _duplicationNumbers;
    _tentativeBranchNumbers = _branchNumbers;
    _tentativeNum0Lineages = _num0Lineages;
    _tentativeNum1Lineages = _num1Lineages;
    _tentativeNum2Lineages =_num2Lineages;
    _tentativeNodesToTryInNNISearch.clear();
    
    //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
    //We use the same procedure as in doNNI !
    const Node * son    = tree_->getNode(nodeId);
    
    if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
    const Node * parent = son->getFather();
    
    if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
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
    if(TreeTemplateTools::hasNodeWithId(*(_rootedTree.getRootNode()->getSon(0)),id)) {
      treeForNNI->newOutGroup(_rootedTree.getRootNode()->getSon(1)->getId());
    }
    else {
      treeForNNI->newOutGroup(_rootedTree.getRootNode()->getSon(0)->getId());
    }
    double ScenarioMLValue = 0;
    _totalIterations = _totalIterations+1;
    //If we want to optimize the root or if we are at the first try
    if (_rootOptimization){
      //PROVISORY
      //ScenarioMLValue = findMLReconciliationSafe(&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers,_duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers);
      // ScenarioMLValue = findMLReconciliation(&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers,_duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages,_tentativeNum2Lineages);  
      
      
      if (_heuristicsLevel>0) {
        ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI/*&_rootedTree*/, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      else {
        ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI/*&_rootedTree*/, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      
      
      
      
    }
    else {
      //   ScenarioMLValue = findMLReconciliation(&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers,_duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages);
      if (_heuristicsLevel>0) {
        ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI/*&_rootedTree*/, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      else {
        ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI/*&_rootedTree*/, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
      }
      
      
      
      
    }
    
    //double ScenarioMLValue = findMLReconciliation (&_spTree, treeForNNI, _seqSp, _lossProbabilities, _duplicationProbabilities);
    delete treeForNNI;
   // std::cout <<"returned ScenarioMLValue : "<<ScenarioMLValue<< " against _scenarioLikelihood : "<<_scenarioLikelihood<<std::endl;
    if (ScenarioMLValue >  _scenarioLikelihood) { //If it is worth computing the sequence likelihood
      //getting the nodes of interest in the unrooted tree
    //  std::cout << "NNI giving a more likely scenario: we compute the sequence likelihood" <<std::endl;

      
      
      //Retrieving arrays of interest:
      const DRASDRTreeLikelihoodNodeData * parentData = & likelihoodData_->getNodeData(parent->getId());
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
      
      const DRASDRTreeLikelihoodNodeData * grandFatherData = & likelihoodData_->getNodeData(grandFather->getId());
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
          //computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, array1, nbGrandFatherNeighbors + 1, _nbDistinctSites, _nbClasses, _nbStates, false); 
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
      //	std::cout <<"HERER8"<<std::endl;
      //Initialize BranchLikelihood:
      brLikFunction_->initModel(model_, rateDistribution_);
      brLikFunction_->initLikelihoods(&array1, &array2);
      ParameterList parameters;	
      
      std::string name = "BrLen" + TextTools::toString(parent->getId()); 
      
      //Julien Proposition :
      
      std::string p = "BrLen";
      for(unsigned int i = 0; i < nodes_.size(); i++)
        if(nodes_[i]->getId() == parent->getId())
          {
            p += TextTools::toString(i);
            break;
          }
      Parameter brLen = getParameter(p);
      
      
      //Parameter brLen = *_parameters.getParameter("BrLen" + TextTools::toString(parent->getId()));
      
      
      
      brLen.setName("BrLen"); 
      parameters.addParameter(brLen);
      brLikFunction_->setParameters(parameters);
      //Re-estimate branch length:
     // std::cout << " Setting Message handler To NULL"<<std::endl;
      brentOptimizer_->setMessageHandler(NULL);
      brentOptimizer_->setFunction(brLikFunction_);
      brentOptimizer_->getStopCondition()->setTolerance(0.1);
      brentOptimizer_->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
      brentOptimizer_->init(parameters); 
      brentOptimizer_->setMessageHandler(NULL);
      brentOptimizer_->optimize(); 
      
      brLenNNIValues_[nodeId] = brentOptimizer_->getParameters().getParameter("BrLen").getValue();
      brLikFunction_->resetLikelihoods(); //Array1 and Array2 will be destroyed after this function call.
      //We should not keep pointers towards them...
      
    //  std::cout <<"_brLenNNIValues[nodeId] "<<_brLenNNIValues[nodeId]<<std::endl;
      
      
      //Return the resulting likelihood:
      double temp = getValue() ; 
      
     // std::cout << " TestNNI : getValue () : " << temp ;
      
      /* std::cout << " TestNNI : getValue () : " << temp ;
       std::cout << " TestNNI : ScenarioMLValue : " << ScenarioMLValue ;
       std::cout << " TestNNI : _brLikFunction->getValue() : " << _brLikFunction->getValue() <<std::endl;*/
      double tot = brLikFunction_->getValue() - ScenarioMLValue - temp;
      if (tot<0) {
        _tentativeScenarioLikelihood=ScenarioMLValue;
        // std::cout << "NNI done, better scenario with "<< VectorTools::sum(_tentativeLossNumbers) <<" losses and "<<VectorTools::sum(_tentativeDuplicationNumbers)<<" duplications"<<std::endl;
        //    std::cout << "BETTER TOTAL LIKELIHOOD"<<std::endl;
      }
      /*  if (tot<0) {
       _MLindex = MLindex;  
       _lossNumbers = lossNumbers;
       _duplicationNumbers = duplicationNumbers;
       _branchNumbers = branchNumbers;
       
       }*/
      //std::cout << "OUT 1 testNNI "<<std::endl;
      return tot;
    }
    else {
         // std::cout << "NNI giving a less likely scenario: we do not compute the sequence likelihood" <<std::endl;
      _tentativeMLindex = -1;
     // std::cout << "OUT 2 testNNI "<<std::endl;
      return getValue();
    }
    
  }
  else {
    _tentativeMLindex = -1;
   // std::cout << "OUT 3 testNNI "<<std::endl;
    return getValue();
  }
}
    
/*******************************************************************************/

void ReconciliationTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
// std::cout << "\n\n/////////////////////////////////DO NNI////////////////////////////////////\n\n"<<std::endl;	
  //Perform the topological move, the likelihood array will have to be recomputed...
  Node * son    = tree_->getNode(nodeId);
	if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'son' must not be the root node.", son);
  Node * parent = son->getFather();
	if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'parent' must not be the root node.", parent);
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
	//The two following functions produce the same error : _seqSp becomes empty !
	//	std::cout <<"_seqSp.size()"<<_seqSp.size()<<std::endl; 
	//	std::cout <<"in do NNI BEFORE : _scenarioLikelihood: " <<_scenarioLikelihood <<std::endl;
	double temp = getValue() ;
	/*	std::cout <<"Before findMLReconciliationSafe : MLindex"<<_MLindex <<std::endl;
	  _scenarioLikelihood = findMLReconciliationSafe(&_spTree, _tree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers);
	std::cout <<"After findMLReconciliationSafe : MLindex"<<_MLindex <<std::endl;*/
	
	//_scenarioLikelihood = findMLReconciliation (&_spTree, _tree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, MLindex);
	//	std::cout <<"in do NNI AFTER : _scenarioLikelihood: " <<_scenarioLikelihood <<std::endl;
//	temp = getValue() ;

//Julien Proposition :
  
  std::string name = "BrLen";
  for(unsigned int i = 0; i < nodes_.size(); i++)
    if(nodes_[i]->getId() == parent->getId())
      {
	name += TextTools::toString(i);
	break;
      }

  //  Parameter brLen = *_parameters.getParameter("BrLen" + TextTools::toString(parent->getId()));
  
  // Parameter brLen = *_parameters.getParameter(p);

 
  if(brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
  {
    double length = brLenNNIValues_[nodeId];
    brLenParameters_.setParameterValue(name, length);
    Parameter p = getParameter(name);
    p.setValue(length);
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


  /*
    std::string name = "BrLen" + TextTools::toString(parent->getId());
    if(_brLenNNIValues.find(nodeId) != _brLenNNIValues.end())
    {
    double length = _brLenNNIValues[nodeId];
    _brLenParameters.setParameterValue(name, length);
    Parameter * p = _parameters.getParameter(name);
    if(p) p->setValue(length);
    parent->setDistanceToFather(length);
    }
    else std::cerr << "ERROR, branch not found: " << nodeId << std::endl;
    try { _brLenNNIParams.addParameter(*_brLenParameters.getParameter(name)); }
    catch(ParameterException & ex)
    {
    std::cerr << "DEBUG:" << std::endl;
    _brLenNNIParams.printParameters(cerr);
    std::cerr << "DEBUG:" << name << std::endl;
    }
  */


  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  brLenNNIParams_[brLenNNIParams_.size()-1].removeConstraint();
//  (brLenNNIParams_.rbegin())->removeConstraint();
  // std::cout << "before affectations"<<std::endl;


  // _tree->resetNodesId(); //Fait planter ????
  /*
    if (_MLindex == _tentativeMLindex) {
    std::cout << "!!!!!!!!!!!! SAME !!!!!!!!!!!!!!!!!"<<std::endl;
    _counter +=1;
    std::cout << "COUNTER"<<_counter<<std::endl;
    }
    else {
    std::cout << "Different : "<< _MLindex<< " != "<<_tentativeMLindex <<std::endl;
    if (!VectorTools::contains(_listOfPreviousRoots, _tentativeMLindex)) {
    std::cout <<"Previously not contained !"<<std::endl;
    _listOfPreviousRoots.push_back(_tentativeMLindex);
    if (!_tree->getNode(_tentativeMLindex)->isLeaf()){
    if (!VectorTools::contains(_listOfPreviousRoots, _tree->getNode(_tentativeMLindex)->getSon(0)->getId())) {
    _listOfPreviousRoots.push_back(_tree->getNode(_tentativeMLindex)->getSon(0)->getId());
    }
    if (!VectorTools::contains(_listOfPreviousRoots, _tree->getNode(_tentativeMLindex)->getSon(1)->getId())) {
    _listOfPreviousRoots.push_back(_tree->getNode(_tentativeMLindex)->getSon(1)->getId());
    }
    }
    if (_tree->getNode(_tentativeMLindex)->hasFather()){
    _listOfPreviousRoots.push_back(_tree->getNode(_tentativeMLindex)->getFather()->getId());
    }
    std::cout <<"Lots Added !"<<std::endl;
    _counter = 0;
    }
    else {
    _counter +=1;
    }
    std::cout <<"AFTER : ";
    for (int i = 0; i<_listOfPreviousRoots.size();i++){
    std::cout << _listOfPreviousRoots[i]<< " ";
    }
    std::cout <<" The End"<<std::endl;
    //    _counter = 0;
    _MLindex = _tentativeMLindex;
    }*/
  _MLindex = _tentativeMLindex;
  // std::cout << "IN DO_NNI _MLindex : "<< _MLindex << "_tentativeMLindex: " << _tentativeMLindex<<std::endl;
  _duplicationNumbers = _tentativeDuplicationNumbers;
  _lossNumbers = _tentativeLossNumbers;
  _branchNumbers = _tentativeBranchNumbers;

  _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
  


  // _tree->resetNodesId();
  /* if (_tree->isRooted()) {
     _tree->unroot();
     }std::cout <<"AFTER STUFF 2"<<std::endl;
     _tree->resetNodesId();
     _tree->newOutGroup(_MLindex); 
     _tree->resetNodesId();*/

  _scenarioLikelihood = _tentativeScenarioLikelihood;// + _brLikFunction->getValue();

  //Now we need to update _rootedTree
  TreeTemplate<Node> * tree = tree_->clone();
  //First we root this temporary tree as in _rootedTree (same lines as in testNNI)
  int id = tree->getRootNode()->getId();
  if(TreeTemplateTools::hasNodeWithId(*(_rootedTree.getRootNode()->getSon(0)),id)) {
    tree->newOutGroup(_rootedTree.getRootNode()->getSon(1)->getId());
  }
  else {
    tree->newOutGroup(_rootedTree.getRootNode()->getSon(0)->getId());
  }
  //Then we root this tree according to MLindex
  tree->newOutGroup(_MLindex);
  //We update _rootedTree
  _rootedTree = *(tree->clone());
  //we need to update the sequence likelihood
  OptimizeSequenceLikelihood(true);
  //resetSequenceLikelihood();//_sequenceLikelihood = UNLIKELY;
/*  std::cout << "before treeToParenthesis"<<std::endl;
  std::cout << TreeTools::treeToParenthesis (*tree, true)<<std::endl;
  std::cout << "END doNNI"<<std::endl;*/
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

std::vector <int> ReconciliationTreeLikelihood::get0LineagesNumbers(){
  return _num0Lineages;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::get1LineagesNumbers(){
  return _num1Lineages;
}

/*******************************************************************************/

std::vector <int> ReconciliationTreeLikelihood::get2LineagesNumbers(){
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

