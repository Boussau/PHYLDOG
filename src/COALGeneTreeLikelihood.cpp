//
// File: COALGeneTreeLikelihood.cpp
// Created by: Bastien Boussau
// Created on: Tue Oct 04 18:14:51 2011
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


#include "COALGeneTreeLikelihood.h"
#include "GenericTreeExplorationAlgorithms.h"

// From Utils:
//#include <Bpp/Text/TextTools.h>
//#include <Bpp/App/ApplicationTools.h>


// From the STL:
//#include <iostream>


//using namespace std;
using namespace bpp;

/******************************************************************************/
COALGeneTreeLikelihood::COALGeneTreeLikelihood(
                       const Tree & tree,
                       SubstitutionModel * model,
                       DiscreteDistribution * rDist,
                       TreeTemplate<Node> & spTree,  
                       TreeTemplate<Node> & rootedTree, 
                       TreeTemplate<Node> & geneTreeWithSpNames,
                       const std::map <std::string, std::string> seqSp,
                       std::map <std::string,int> spId,
                       std::vector < std::vector < std::vector < std::vector<unsigned int> > > > coalCounts,
                       std::vector < double > coalBl,
                       int speciesIdLimitForRootPosition,
                       int heuristicsLevel,
                       int & MLindex, 
                       bool checkRooted,
                       bool verbose,
                       bool rootOptimization, 
                       bool considerSequenceLikelihood, 
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
    _coalBl = coalBl;
    _coalCounts = coalCounts;
    _tentativeCoalCounts = coalCounts;
    for (unsigned int i = 0 ; i < _coalBl.size() ; i++) {
        num12Lineages_.push_back(0);
        num22Lineages_.push_back(0);
    }
    computeNumLineagesFromCoalCounts ();
}

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
                   checkRooted,
                   verbose,
                   rootOptimization, 
                   considerSequenceLikelihood, 
                   sprLimit)
{
    _coalBl = coalBl;
    _coalCounts = coalCounts;
    _tentativeCoalCounts = coalCounts;
    for (unsigned int i = 0 ; i < _coalBl.size() ; i++) {
        num12Lineages_.push_back(0);
        num22Lineages_.push_back(0);
    }

    computeNumLineagesFromCoalCounts ();
}

/******************************************************************************/

COALGeneTreeLikelihood::COALGeneTreeLikelihood(const COALGeneTreeLikelihood & lik):
GeneTreeLikelihood(lik)
{
    _coalBl = lik._coalBl;
    _coalCounts = lik._coalCounts;
    _tentativeCoalCounts = lik._tentativeCoalCounts;
    num12Lineages_ = lik.num12Lineages_;
    num22Lineages_ = lik.num22Lineages_;

}

/******************************************************************************/

COALGeneTreeLikelihood & COALGeneTreeLikelihood::operator=(const COALGeneTreeLikelihood & lik)
{
    GeneTreeLikelihood::operator=(lik);
    _coalBl = lik._coalBl;
    _coalCounts = lik._coalCounts;
    _tentativeCoalCounts = lik._tentativeCoalCounts;
    num12Lineages_ = lik.num12Lineages_;
    num22Lineages_ = lik.num22Lineages_;
    return *this;
}




/******************************************************************************/

COALGeneTreeLikelihood::~COALGeneTreeLikelihood()
{
    if (nniLk_) delete nniLk_;
    if (_spTree) delete _spTree;
    if (_rootedTree) delete _rootedTree;
    if (_geneTreeWithSpNames) delete _geneTreeWithSpNames;
}



/******************************************************************************/

void COALGeneTreeLikelihood::initParameters()
{
    // std::cout << "in initParameters"<<std::endl;
    if (_considerSequenceLikelihood) {
        nniLk_->initParameters();
    }
    if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
        exit(-1);
     }
    else {
		/*std::cout << "BEFORE: "<<std::endl;
		for (unsigned int i = 0 ; i < _coalCounts[0][0].size() ; i++) {
			if (_spTree->getNode(i)->isLeaf() ) {
				std::cout << "Leaf i: "<<i<<" _coalCounts[0][0][i][0]: "<< _coalCounts[0][0][i][0] <<" _coalCounts[0][0][i][1] " << _coalCounts[0][0][i][1] << std::endl;
			}
		}*/
        _scenarioLikelihood = findMLCoalReconciliationDR (_spTree, _rootedTree, 
                                    _seqSp, _spId, 
                                    _coalBl, 
                                    _MLindex, 
                                    _coalCounts, 
                                    _nodesToTryInNNISearch);
		computeNumLineagesFromCoalCounts ();

    }
    _MLindex = -1;
}



/******************************************************************************/

void COALGeneTreeLikelihood::resetMLindex() {
    _MLindex = -1;
}




/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double COALGeneTreeLikelihood::getLogLikelihood() const
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
    if (_considerSequenceLikelihood) {
        ll = nniLk_->getLogLikelihood() + _scenarioLikelihood;
    }
    else {
        ll = _scenarioLikelihood;
    }
    // ll = _sequenceLikelihood ;
    return ll;
}
/******************************************************************************/
//returns -loglikelihood
//As a reminder: _sequenceLikelihood < 0, _scenarioLikelihood < 0, getValue >0
double COALGeneTreeLikelihood::getValue() const  
throw (Exception)
{
    {
        if(!nniLk_->isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
        // return (-getLogLikelihood());
        //TEST 16 02 2010
        // std::cout<<"\t\t\t_sequenceLikelihood: "<<_sequenceLikelihood<< " _scenarioLikelihood: "<<_scenarioLikelihood<<std::endl;
        if (_considerSequenceLikelihood) {
            return (- nniLk_->getLogLikelihood() - _scenarioLikelihood);
        }
        else {
            return (- _scenarioLikelihood);
        }
        //return (minusLogLik_ - _scenarioLikelihood);
        //return (-minusLogLik_);
    }
}



/******************************************************************************/


void COALGeneTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
    if (_considerSequenceLikelihood) {
        nniLk_->applyParameters();
        nniLk_->matchParametersValues(params);
    }
    //If we need to update the reconciliation likelihood
    if (_optimizeReconciliationLikelihood) {
        computeReconciliationLikelihood();
    }
    
    
}


/******************************************************************************/


void COALGeneTreeLikelihood::computeSequenceLikelihood()
{
    if ( _considerSequenceLikelihood && (_optimizeSequenceLikelihood==true) ) {
        nniLk_->computeTreeLikelihood ();
        /*
         nniLk_->computeSubtreeLikelihoodPostfix(dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getRootNode());
         nniLk_->computeSubtreeLikelihoodPrefix(nniLk_->getTree().getRootNode());
         nniLk_->computeRootLikelihood();*/
    }
}




/******************************************************************************/

void COALGeneTreeLikelihood::computeReconciliationLikelihood()
{
    if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
        exit(-1);
    }
    else {
        //Compute the COAL likelihood
		/*std::cout << "BEFORE: "<<std::endl;
		for (unsigned int i = 0 ; i < _coalCounts[0][0].size() ; i++) {
			if (_spTree->getNode(i)->isLeaf() ) {
				std::cout << "Leaf i: "<<i<<" _coalCounts[0][0][i][0]: "<< _coalCounts[0][0][i][0] <<" _coalCounts[0][0][i][1] " << _coalCounts[0][0][i][1] << std::endl;
			}
		}*/

        _scenarioLikelihood = findMLCoalReconciliationDR (_spTree, _rootedTree, 
                                                           _seqSp, _spId, 
                                                           _coalBl, 
                                                           _tentativeMLindex, 
                                                           _tentativeCoalCounts, 
                                                           _tentativeNodesToTryInNNISearch); 
        _MLindex = _tentativeMLindex;
        _coalCounts = _tentativeCoalCounts;
        _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
		computeNumLineagesFromCoalCounts ();

    }
}



/******************************************************************************/


void COALGeneTreeLikelihood::computeTreeLikelihood()
{
    if (_considerSequenceLikelihood )
    {
        computeSequenceLikelihood();
    }
    computeReconciliationLikelihood();  
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.
 
 */
double COALGeneTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
    //int nodeId = son->getId();
    //  std::cout<<"IN TESTNNI, nodeId "<< nodeId << "_nodesToTryInNNISearch.size() "<< _nodesToTryInNNISearch.size()<<std::endl;
    // std::cout << "before "<<TreeTemplateTools::treeToParenthesis (*tree_, true)<<std::endl;
    //If the NNI is around a branch where a duplication was found, 
    //or if we just try all branches because the starting gene trees are parsimonious in
    //numbers of DL.
    if ((_nodesToTryInNNISearch.count(nodeId)==1) /*|| _DLStartingGeneTree*/) {
        TreeTemplate<Node> * treeForNNI = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
        
        _tentativeMLindex = _MLindex;
        /* _tentativeLossNumbers = _lossNumbers;
         _tentativeDuplicationNumbers = _duplicationNumbers;
         _tentativeBranchNumbers = _branchNumbers;*/
        _tentativeCoalCounts = _coalCounts;

        _tentativeNodesToTryInNNISearch.clear();
        
        //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
        //We use the same procedure as in doNNI !
        const Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getNode(nodeId);
        
        
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
        
         //Compute the COAL likelihood
	
        ScenarioMLValue =  findMLCoalReconciliationDR (_spTree, treeForNNI, 
                                                           _seqSp, _spId, 
                                                           _coalBl, 
                                                           _tentativeMLindex, 
                                                           _tentativeCoalCounts, 
                                                           _tentativeNodesToTryInNNISearch, false); 

        
        
        if (treeForNNI) delete treeForNNI;
         //std::cout<<"???WORTH computing the sequence likelihood "<< ScenarioMLValue<< " "<< _scenarioLikelihood<<std::endl;
        if (_considerSequenceLikelihood ) 
        {
            if  (ScenarioMLValue >  _scenarioLikelihood) 
            { //If it is worth computing the sequence likelihood
                //Retrieving arrays of interest:
                // std::cout << "before nniLk_->testNNI(nodeId);"<<std::endl;
                double newLkMinusOldLk = nniLk_->testNNI(nodeId);
//                 std::cout << "after nniLk_->testNNI(nodeId); "<< newLkMinusOldLk <<std::endl;

                double tot = - ScenarioMLValue + _scenarioLikelihood + newLkMinusOldLk;
		//		std::cout << "TOT: "<<tot<<std::endl;
                // if (newLkMinusOldLk<0) 
                if (tot < 0)
                {
                    _tentativeScenarioLikelihood=ScenarioMLValue;
                }
                return tot;
                // std::cout << "after "<<TreeTemplateTools::treeToParenthesis (tree_, true)<<std::endl;
                
            }
            else 
            {
                _tentativeMLindex = -1;
                return 1;
            }
        }
        else 
        {
            if  (ScenarioMLValue >  _scenarioLikelihood) 
            {
                _tentativeScenarioLikelihood=ScenarioMLValue;
                return (- ScenarioMLValue + _scenarioLikelihood);
            }
            else 
            {
                _tentativeMLindex = -1;
                return 1;
            }
        }
    }
    else {
        _tentativeMLindex = -1;
        return 1;
    }
}

/*******************************************************************************/

void COALGeneTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  //  std::cout<<"\t\t\tIN DONNI "<< std::endl;
    //std::cout << TreeTemplateTools::treeToParenthesis(*tree_, true) << std::endl;
    //Perform the topological move, the likelihood array will have to be recomputed...
    
    nniLk_->doNNI(nodeId);
    
    /*
     Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getNode(nodeId);
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
     
     
     if (_considerSequenceLikelihood ) 
     {
     if(brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
     {
     double length = brLenNNIValues_[nodeId];
     brLenParameters_.setParameterValue(name, length);
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
    
    _MLindex = _tentativeMLindex;
    
    _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
    
    // std::cout <<"_tentativeScenarioLikelihood: "<<_tentativeScenarioLikelihood<<std::endl;
    _scenarioLikelihood = _tentativeScenarioLikelihood;// + _brLikFunction->getValue();
    
    //Now we need to update _rootedTree
    TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
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
    if (_considerSequenceLikelihood ) 
    {
        OptimizeSequenceLikelihood(true);
    }
    OptimizeReconciliationLikelihood(true);
}

/*******************************************************************************/

std::vector < std::vector < std::vector< std::vector<unsigned int> > > > COALGeneTreeLikelihood::getCoalCounts() const{
    return _coalCounts;
}



void COALGeneTreeLikelihood::computeNumLineagesFromCoalCounts () {  
//	unsigned int id = _rootedTree->getRootNode()->getId();
	//std::cout << "id: "<< id <<std::endl;

	//List all keys supposed to be species
	/*std::cout << "spIds: "<< std::endl;
	for (std::map<std::string, int >::iterator it2 = _spId.begin(); it2 != _spId.end(); it2++) {
		std::cout << it2->second << " ";
	}
	std::cout << std::endl;*/

   for (unsigned int i = 0 ; i < _coalCounts[0][0].size() ; i++) {
       /* if ( _coalCounts[0][0][i][0]  == 1 && _coalCounts[0][0][i][1] == 2) {
            num12Lineages_[i] = 1; 
            num22Lineages_[i] = 0; 
        }
        else if (_coalCounts[0][0][i][0]  == 2 && _coalCounts[0][0][i][1] == 2) {
            num12Lineages_[i] = 0; 
            num22Lineages_[i] = 1; 
        }
		else {
			num12Lineages_[i] = 0; 
            num22Lineages_[i] = 0; 
		}*/
	   if ( _coalCounts[0][0][i][0]  == 1 && _coalCounts[0][0][i][1] == 2 ) {
		   num12Lineages_[i] = 1; 
		   num22Lineages_[i] = 0; 
	   }
	   else  if ( _coalCounts[0][0][i][0]  == _coalCounts[0][0][i][1]  && _coalCounts[0][0][i][1] != 1 ){
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
    return _coalBl;
}

/*******************************************************************************/

void COALGeneTreeLikelihood::setCoalBranchLengths (std::vector < double > coalBl)
{    
    _coalBl = coalBl;
}

/*******************************************************************************/

int COALGeneTreeLikelihood::getRootNodeindex(){
    return _MLindex;
}

/*******************************************************************************/
/*
 void COALGeneTreeLikelihood::resetSequenceLikelihood(){
 _sequenceLikelihood = UNLIKELY;
 }
 */
/*******************************************************************************/

double COALGeneTreeLikelihood::getSequenceLikelihood() {
    return nniLk_->getValue(); 
}
/*******************************************************************************/

void COALGeneTreeLikelihood::initialize() {
    nniLk_->initialize();
    return;
}



/*******************************************************************************/
void COALGeneTreeLikelihood::print () const {
    std::cout << "Species tree:"<<std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis (getSpTree(), true)<<std::endl;
    std::cout << "Gene family rooted tree:"<<std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
    std::cout << "Gene family tree:"<<std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis (nniLk_->getTree(), true)<<std::endl;
    /*std::cout << "Counts of coalescence events"<<std::endl;
    VectorTools::print(getCoalCounts());*/
    std::cout << "Branch lengths"<<std::endl;
    VectorTools::print(getCoalBranchLengths());
    std::cout << "Root index"<<std::endl;
    std::cout << _MLindex <<std::endl;
}


/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 ************************************************************************/
void COALGeneTreeLikelihood::refineGeneTreeSPRs(map<string, string> params) {
    /*  std::cout <<"\t\t\tStarting MLSearch : current tree : "<< std::endl;
     std::cout<< TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;*/
    if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
        //We don't do SPRs
        //std::cout << "WE DONT DO SPRS"<<std::endl;
        computeReconciliationLikelihood();
        return;
    }
    // std::cout << "WE DO SPRS"<<std::endl;
    
    std::vector <int> nodeIdsToRegraft;
    bool betterTree;
    TreeTemplate<Node> * treeForSPR = 0;
    TreeTemplate<Node> * bestTree = 0;
    // if (getLogLikelihood()==UNLIKELY) 
    computeReconciliationLikelihood();
    //    ReconciliationTreeLikelihood * bestTreeLogLk = this->clone();
    double logL = getLogLikelihood();
    //std::cout << "logL before mapping: " <<getSequenceLikelihood()<<std::endl;
    //initial perturbation of branch lengths, to use suboptimal ones for speed
    /*ParameterList bls = nniLk_->getBranchLengthsParameters () ;
     for (unsigned int i  = 0 ; i < bls.size()/3 ; i++) {
     bls[i].setValue(0.1);
     }
     nniLk_->matchParametersValues(bls);*/
    //  std::cout << "logL before mapping 2: " <<getSequenceLikelihood()<<std::endl;    
    //TEST 14 03 2012, if it works, it's inefficient
    
    NNIHomogeneousTreeLikelihood * drlk = 0;
    drlk = nniLk_->clone();
    // optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (nniLk_), 0.1, params);
    
    /*
     int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
     bool backupOpt = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
     params[ std::string("optimization.max_number_f_eval")] = 100;
     params[ std::string("optimization.topology")] = "false";
     PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
     params[ std::string("optimization.max_number_f_eval")] = backup;
     params[ std::string("optimization.topology")] = backupOpt;*/
    
    
    auto_ptr<BackupListener> backupListener;
    //unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", params, 1, "", true, false);
    double tolerance = 0.1;
    unsigned int tlEvalMax = 1000000;
    OutputStream* messageHandler = 0 ; 
    OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (drlk), 
                                                       drlk->getParameters(), backupListener.get(), 
                                                       tolerance, tlEvalMax, messageHandler, messageHandler, 0);
    
    
    
    delete nniLk_;
    nniLk_ = drlk->clone();
    
    delete drlk;
    drlk = 0;
    
    /*  optimizeBLForSPRs(dynamic_cast<DRTreeLikelihood*> (nniLk_),
     0.1, params);*/
    
    
    
    // std::cout<< "Unrooted starting tree: "<<TreeTemplateTools::treeToParenthesis(nniLk_->getTree(), true)<< std::endl;
    /*
     bls = nniLk_->getBranchLengthsParameters () ;
     for (unsigned int i  = 0 ; i < bls.size() ; i++) {
     if (_rootedTree->getNode(i)->hasFather()) {
     _rootedTree->getNode(i)->setDistanceToFather(bls[i].getValue());
     }
     }*/
    
    logL = getLogLikelihood();
    //  std::cout << "logL after mapping: " <<getSequenceLikelihood()<<std::endl;
    
    double bestlogL = logL;
    double candidateScenarioLk ;
    double bestSequenceLogL = getSequenceLikelihood();
    double bestScenarioLk = getScenarioLikelihood();
    // std::cout << "LOGL: "<<logL << "ScenarioLK: "<< bestScenarioLk <<"; sequenceLK: "<<getSequenceLikelihood() << std::endl;
    unsigned int numIterationsWithoutImprovement = 0;
    NNIHomogeneousTreeLikelihood * bestDrlk = 0;
    breadthFirstreNumber (*_rootedTree);
    
    
    // DRHomogeneousTreeLikelihood drlk;
    //  std::cout<< "Starting tree: "<<TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    
    string parentDup;
    string nodeDup;
    string numLoss = "0";
    
    bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    
    
    while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    {
        
        annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                               *_rootedTree, 
                                               _rootedTree->getRootNode(), 
                                               _seqSp, _spId); 
        
        for (int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
        {
            Node * n = _rootedTree->getNode(nodeForSPR);
            /*            if (n->getFather()->hasBranchProperty("Ev")) {
             parentDup = (dynamic_cast<const BppString *>(n->getFather()->getBranchProperty("Ev")))->toSTL() ;
             }*/
            if (n->hasBranchProperty("L")) {
                numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
            }
            /* else {
             parentDup = "F" ;}*/
            //   std::cout << "parentDup : "<< parentDup <<std::endl;
            /*  if ( _rootedTree->getNode(nodeForSPR)->hasNodeProperty("D")) { 
             nodeDup = (dynamic_cast<const BppString *>(_rootedTree->getNode(nodeForSPR)->getNodeProperty("D")))->toSTL() ; }
             else {
             nodeDup = "F" ; }*/
            //if ( parentDup == "D"  /*|| nodeDup == "Y"*/) {
            if ( numLoss != "0"  /*|| nodeDup == "Y"*/) {
                
                //CHANGE12 11 2011 buildVectorOfRegraftingNodesLimitedDistance(*_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                buildVectorOfRegraftingNodesCoalGeneTree(*_spTree, *_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                
                /*              
                 if (_rootedTree->getNode(nodeForSPR)->isLeaf()) {
                 std::cout << "nodeForSPR: "<< _rootedTree->getNode(nodeForSPR)->getName()<<" ; "<< nodeForSPR << "; Node ids to regraft: " <<std::endl;
                 }
                 else {
                 std::cout << "nodeForSPR: "<< nodeForSPR << "; Node ids to regraft: " <<std::endl;
                 }
                 VectorTools::print(nodeIdsToRegraft);
                 */
                betterTree = false;
                for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
                {
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    treeForSPR = _rootedTree->clone();
                    // treeForSPR->getRootNode()->getSon(0)->setName("outgroupNode");
                    
                    makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false);
                    
                    //breadthFirstreNumber (*treeForSPR);
                    
                    //Compute the COAL likelihood
                    candidateScenarioLk =  findMLCoalReconciliationDR (_spTree, treeForSPR, 
                                                                   _seqSp, _spId, 
                                                                   _coalBl, 
                                                                   _tentativeMLindex, 
                                                                   _tentativeCoalCounts, 
                                                                   _tentativeNodesToTryInNNISearch, false); 
                    // std::cout<< "candidate tree: "<< TreeTemplateTools::treeToParenthesis(*treeForSPR, true)<< std::endl;
                    
                    if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the DL likelihood is not bad
                    {
                        
                        if (computeSequenceLikelihoodForSPR) {
                            
                            // std::cout << "good candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                            if (drlk) {
                                delete drlk;
                                drlk = 0;
                            }
                            // std::cout << "COMPUTING SEQLK: " << TreeTemplateTools::treeToParenthesis(*treeForSPR, true) <<std::endl;
                            drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, 
                                                                      *(nniLk_->getData()), 
                                                                      nniLk_->getSubstitutionModel(), 
                                                                      nniLk_->getRateDistribution(), 
                                                                      true, false);
                            //  drlk  = new DRTreeLikelihood (*treeForSPR, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                            
                            
                            //drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                            /*          ParameterList bls = nniLk_->getBranchLengthsParameters () ;
                             for (unsigned int i  = 0 ; i < bls.size()/3 ; i++) {
                             bls[i].setValue(0.1);
                             }
                             drlk->matchParametersValues(bls);
                             */
                            drlk->initialize();
                            /*std::cout << "TEST getValue: "<< drlk->getValue() <<std::endl;
                             std::cout<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            
                            //  std::cout << "Good SPR; Sequence lk before optimization: "<< -drlk->getValue() << std::endl;
                            //  _rootedTree = treeForSPR->clone();
                            // makeSPR(*tree_, nodeForSPR, nodeIdsToRegraft[i]);
                            // tree_->unroot();
                            //  tree_ = treeForSPR->clone();
                            
                            //Compute sequence likelihood and improve branch lengths on the gene tree.
                            /*  std::cout<< TreeTemplateTools::treeToParenthesis(*tree_, true)<< std::endl;
                             optimizeBLMapping(this,
                             0.1);
                             */
                            // tree_->unroot();
                            //fireParameterChanged(getParameters());
                            //Compute sequence likelihood and improve branch lengths on the gene tree.
                            //  std::cout<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;
                            
                            //  params[ std::string("optimization.topology")] = "false";
                            
                            /*TEMP18012012 */
                            /*optimizeBLMappingForSPRs does not work for some families...
                             optimizeBLMappingForSPRs(dynamic_cast<DRTreeLikelihood*> (drlk),
                             0.1, params);*/
                            /*                            optimizeBLForSPRs(dynamic_cast<DRTreeLikelihood*> (drlk),
                             0.1, params);
                             */  
                            /* std::cout << "before norm opt: " <<std::endl;
                             std::cout<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            
                            //TEST 13 03 2012
                            OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (drlk), drlk->getParameters(), backupListener.get(), tolerance, tlEvalMax, messageHandler, messageHandler, 0);
                            
                            /*
                             int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
                             bool backupOpt = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
                             params[ std::string("optimization.max_number_f_eval")] = 100;
                             params[ std::string("optimization.topology")] = "false";
                             PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
                             params[ std::string("optimization.max_number_f_eval")] = backup;
                             params[ std::string("optimization.topology")] = backupOpt;
                             */
                            
                            
                            /* std::cout << "TEST getValue2: "<< drlk->getValue() <<std::endl;
                             std::cout<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            
                            /*std::cout << "Good SPR; Sequence lk after optimization: "<< -drlk->getValue() << std::endl;
                             std::cout<<"Good SPR; optimized tree: "<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            logL = candidateScenarioLk - drlk->getValue();
                        }
                        else {
                            logL = candidateScenarioLk - bestSequenceLogL;
                        }
                        
                    }
                    else { 
                        // std::cout << "bad candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                        logL =logL - 10;
                    }
                    //If the candidate tree has a DL + sequence Lk better than the current best
                    if (logL - 0.01 > bestlogL) 
                    {
                        betterTree = true;
                        bestlogL =logL;
                        bestScenarioLk = candidateScenarioLk;
                        if (computeSequenceLikelihoodForSPR) {
                            bestSequenceLogL = drlk->getValue();
                            
                            if (bestDrlk) {
                                delete bestDrlk;
                                bestDrlk = 0;
                            }
                            
                            bestDrlk = drlk->clone();
                        }
                        if (bestTree) {
                            delete bestTree;
                            bestTree = 0;
                        }
                        
                        bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->clone();
                        //Rooting bestTree as in TreeForSPR:
                        vector<Node*> drlkNodes = bestTree->getNodes();
                        //  std::cout << "BEFORE NEWOUTGROUP: "<<TreeTemplateTools::treeToParenthesis(*bestTree, true)<< std::endl;
                        for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
                            if (drlkNodes[j]->hasNodeProperty("outgroupNode")) {
                                if (bestTree->getRootNode() == drlkNodes[j]) {
                                    if (j < drlkNodes.size()-1) 
                                    {
                                        bestTree->rootAt(drlkNodes[drlkNodes.size()-1]);   
                                    }
                                    else {
                                        bestTree->rootAt(drlkNodes[drlkNodes.size()-2]);
                                    }
                                };
                                bestTree->newOutGroup( drlkNodes[j] );
                                //  std::cout << "FOUND"<<std::endl;
                                break;
                            }
                        }
                        //   std::cout << "AFTER NEWOUTGROUP: "<<TreeTemplateTools::treeToParenthesis(*bestTree, true)<< std::endl;
                        
                        /*
                         bestTree = treeForSPR->clone();
                         //                        bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();                                                                             
                         
                         //Putting optimized branch lengths onto bestTree                                                                                                                                           
                         vector<int> drlkNodes = bestTree->getNodesId();
                         for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
                         if ((dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather()))
                         bestTree->getNode(drlkNodes[j])->setDistanceToFather(dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->getNode(drlkNodes[j])->getDistanceToFather());
                         }
                         */
                        
                        
                        //                        bestTree = drlk->getTree().clone();
                        /*   bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->clone();
                         int id = treeForSPR->getRootNode()->getSon(0)->getId();
                         bestTree->newOutGroup(id);*/
                        
                        /*  std::cout << "\t\t\tSPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
                         std::cout << "\t\t\tTotal likelihood: "<<logL <<"; Reconciliation likelihood: "<< bestScenarioLk << ", Sequence likelihood: "<< bestSequenceLogL <<", for tree: "<< std::endl;*/
                        /*TreeTemplateTools::treeToParenthesis(bestTreeLogLk->getRootedTree()) << */
                    }
                    else {
                        //   copyContentsFrom(*bestTreeLogLk);
                        //  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
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
                    if (_rootedTree) 
                    {
                        delete _rootedTree;
                        _rootedTree = 0;
                    }
                    _rootedTree = bestTree->clone();
                    breadthFirstreNumber (*_rootedTree);
                    
                    if (bestTree) {
                        delete bestTree;
                        bestTree = 0;
                    }
                    
                    /* bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->clone();
                     vector<int> drlkNodes = bestTree->getNodesId();
                     
                     for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
                     if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
                     _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
                     }*/
                    
                    /*    std::cout <<"\t\t\tSPRs: Improvement! : "<<numIterationsWithoutImprovement<< std::endl;
                     std::cout << "\t\t\tNew total Likelihood value "<<logL<< std::endl;
                     std::cout <<"\t\t\tNumber of gene trees tried : "<<index<< std::endl;*/
                    /*  if (bestTree) {
                     delete bestTree;
                     bestTree = 0;
                     }*/
                    if (computeSequenceLikelihoodForSPR) {                        
                        if (nniLk_) {
                            delete nniLk_;
                            nniLk_ = 0;
                        }
                        nniLk_ = bestDrlk->clone();                        
                    }
                    /* std::cout << "NOT ROOTED: "<<std::endl;
                     std::cout<< TreeTemplateTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;
                     std::cout << "ROOTED: "<<std::endl;
                     std::cout<< TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
                     */
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
                if (bestTree) {
                    delete bestTree;
                    bestTree = 0;
                }
                if (bestDrlk) {
                    delete bestDrlk;
                    bestDrlk = 0;
                }
                if (drlk) {
                    delete drlk;
                    drlk = 0;
                }
            }
            else {
                numIterationsWithoutImprovement++;
            }
        }
    }
    
    _scenarioLikelihood = bestScenarioLk;
    //final branch lengths optimization
    
    _rootedTree->resetNodesId();
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    /*Do we need this extra lk optimization?
     bestTree = _rootedTree->clone();
     
     drlk = new NNIHomogeneousTreeLikelihood (*bestTree, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
     drlk->initialize();
     
     PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
     
     if (nniLk_) {
     delete nniLk_;
     nniLk_ = 0;
     }
     
     nniLk_ = drlk->clone();
     if (bestTree) {
     delete bestTree;
     bestTree = 0;
     }
     bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
     
     vector<int> drlkNodes = bestTree->getNodesId();
     for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
     if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
     _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
     }
     */
    //One more reconciliation, to update the "_num*Lineages" vectors.
    computeReconciliationLikelihood();
    
    Nhx *nhx = new Nhx();
    annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                           *_rootedTree, 
                                           _rootedTree->getRootNode(), 
                                           _seqSp, _spId); 
    cout << "Reconciled tree: "<<endl;
    nhx->write(*_rootedTree, cout);
    /*  cout << "unrooted tree: "<<endl;
     nhx->write(*bestTree, cout);*/
    //    std::cout<< TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    if (nhx) delete nhx;
    
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
    
    FastRHomogeneousTreeLikelihood * rlk = 0;
    rlk = new FastRHomogeneousTreeLikelihood (nniLk_->getTree(), 
                                              *(nniLk_->getData()), 
                                              nniLk_->getSubstitutionModel(), 
                                              nniLk_->getRateDistribution(), 
                                              true, false);
    
    rlk->initialize();
    rlk->initializeLikelihoodData();
    
    
    auto_ptr<BackupListener> backupListener;
    //unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", params, 1, "", true, false);
    double tolerance = 0.1;
    unsigned int tlEvalMax = 1000000;
    OutputStream* messageHandler = 0 ; 
    OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (rlk), 
                                                       rlk->getParameters(), backupListener.get(), 
                                                       tolerance, tlEvalMax, messageHandler, messageHandler, 0);
    
    
    
    delete nniLk_;
    nniLk_ = new NNIHomogeneousTreeLikelihood (rlk->getTree(), 
                                               *(rlk->getData()), 
                                               rlk->getSubstitutionModel(), 
                                               rlk->getRateDistribution(), 
                                               true, false);
    nniLk_->initialize();
    
    //  nniLk_->reInit();//To update the likelihood data vectors
    delete rlk;
    rlk = 0;
    
    
    logL = getLogLikelihood();
    //  std::cout << "logL after mapping: " <<getSequenceLikelihood()<<std::endl;
    
    double bestlogL = logL;
    double candidateScenarioLk ;
    double bestSequenceLogL = getSequenceLikelihood();
    double bestScenarioLk = getScenarioLikelihood();
    // std::cout << "LOGL: "<<logL << "ScenarioLK: "<< bestScenarioLk <<"; sequenceLK: "<<getSequenceLikelihood() << std::endl;
    unsigned int numIterationsWithoutImprovement = 0;
    FastRHomogeneousTreeLikelihood * bestRlk = 0;
    breadthFirstreNumber (*_rootedTree);
    
    
    // DRHomogeneousTreeLikelihood drlk;
    //  std::cout<< "Starting tree: "<<TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    
    string parentDup;
    string nodeDup;
    string numLoss = "0";
    
    bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    
    
    while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    {
        
        annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                               *_rootedTree, 
                                               _rootedTree->getRootNode(), 
                                               _seqSp, _spId); 
        
        for (int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
        {
            Node * n = _rootedTree->getNode(nodeForSPR);
            if (n->hasBranchProperty("L")) {
                numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
            }
            if ( numLoss != "0"  ) {
                
                buildVectorOfRegraftingNodesCoalGeneTree(*_spTree, *_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                
                betterTree = false;
                for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
                {
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    treeForSPR = _rootedTree->clone();
                    
                    nodesToUpdate = makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false, true);
                    
                    //Compute the COAL likelihood
                    candidateScenarioLk =  findMLCoalReconciliationDR (_spTree, treeForSPR, 
                                                                       _seqSp, _spId, 
                                                                       _coalBl, 
                                                                       _tentativeMLindex, 
                                                                       _tentativeCoalCounts, 
                                                                       _tentativeNodesToTryInNNISearch, false); 
                    
                    if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the scenario likelihood is not bad
                    {
                        
                        if (computeSequenceLikelihoodForSPR) {
                            
                            // std::cout << "good candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                            if (rlk) {
                                delete rlk;
                                rlk = 0;
                            }
                            //  std::cout << "COMPUTING SEQLK: " << TreeTemplateTools::treeToParenthesis(*treeForSPR, true) <<std::endl;
                            ParameterList pl ;//= rlk->getBranchLengthsParameters();
                            std::auto_ptr<Constraint> brLenConstraint;
//                            brLenConstraint.reset(new IncludingInterval(0.000001, 10000));
                            brLenConstraint.reset(new IntervalConstraint(0.000001, 10000, true, true));
                            //Only optimizes branch lengths likely to have changed because of the SPR
                            for (unsigned int k = 0 ; k < nodesToUpdate.size() ; k++) {
                                for (unsigned int j = 0 ; j < nodesToUpdate[k]->getNumberOfSons() ; j++) {
                                    if (! (VectorTools::contains (nodesToUpdate, nodesToUpdate[k]->getSon(j) ) ) ) {
                                        nodesToUpdate[k]->getSon(j)->setNodeProperty("toComp", BppString("N"));
                                    }
                                }
                                if (nodesToUpdate[k]->hasDistanceToFather()) {
                                    brLenConstraint->clone();
                                    if (nodesToUpdate[k]->getDistanceToFather() < 0.000001) {
                                        nodesToUpdate[k]->setDistanceToFather(0.000001);
                                    }
                                    pl.addParameter(Parameter("BrLen" + TextTools::toString(nodesToUpdate[k]->getId()), nodesToUpdate[k]->getDistanceToFather(), brLenConstraint->clone(), true));
                                    
                                }
                            }
                            
                            rlk  = new FastRHomogeneousTreeLikelihood (*treeForSPR, 
                                                                       *(nniLk_->getData()), 
                                                                       nniLk_->getSubstitutionModel(), 
                                                                       nniLk_->getRateDistribution(), 
                                                                       true, false);
                            
                            rlk->initialize();
                            rlk->initializeLikelihoodData();
                            
                            
                            OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (rlk), 
                                                                               pl, backupListener.get(), 
                                                                               tolerance, tlEvalMax, messageHandler, messageHandler, 0);
                            
                            logL = candidateScenarioLk - rlk->getValue();
                        }
                        else {
                            logL = candidateScenarioLk - bestSequenceLogL;
                        }
                        
                    }
                    else { 
                        // std::cout << "bad candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                        logL =logL - 10;
                    }
                    //If the candidate tree has a DL + sequence Lk better than the current best
                    if (logL - 0.01 > bestlogL) 
                    {
                        betterTree = true;
                        bestlogL =logL;
                        bestScenarioLk = candidateScenarioLk;
                        if (computeSequenceLikelihoodForSPR) {
                            bestSequenceLogL = rlk->getValue();
                            
                            if (bestRlk) {
                                delete bestRlk;
                                bestRlk = 0;
                            }
                            
                            bestRlk = rlk->clone();
                        }
                        if (bestTree) {
                            delete bestTree;
                            bestTree = 0;
                        }
                        
                        bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(bestRlk->getTree()))->clone();
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
                                //  std::cout << "FOUND"<<std::endl;
                                break;
                            }
                        }
                        
                    }
                    else {
                        //   copyContentsFrom(*bestTreeLogLk);
                        //  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
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
                    if (_rootedTree) 
                    {
                        delete _rootedTree;
                        _rootedTree = 0;
                    }
                    _rootedTree = bestTree->clone();
                    //  breadthFirstreNumber (*_rootedTree);
                    breadthFirstreNumberAndResetProperties (*_rootedTree);
                    
                    if (bestTree) {
                        delete bestTree;
                        bestTree = 0;
                    }
                    
                    
                    if (computeSequenceLikelihoodForSPR) {                        
                        if (nniLk_) {
                            delete nniLk_;
                            nniLk_ = 0;
                        }
                        nniLk_ = new NNIHomogeneousTreeLikelihood (bestRlk->getTree(), 
                                                                   *(bestRlk->getData()), 
                                                                   bestRlk->getSubstitutionModel(), 
                                                                   bestRlk->getRateDistribution(), 
                                                                   true, false);
                        nniLk_->initialize();
                        // nniLk_->fireParameterChanged(nniLk_->getParameters()); //To update the likelihood data vectors
                        
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
                if (bestTree) {
                    delete bestTree;
                    bestTree = 0;
                }
                if (bestRlk) {
                    delete bestRlk;
                    bestRlk = 0;
                }
                if (rlk) {
                    delete rlk;
                    rlk = 0;
                }
            }
            else {
                numIterationsWithoutImprovement++;
            }
        }
    }
    
    _scenarioLikelihood = bestScenarioLk;
    //final branch lengths optimization
    
    _rootedTree->resetNodesId();
    if (rlk) {
        delete rlk;
        rlk = 0;
    }
    
    //One more reconciliation, to update the "_num*Lineages" vectors.
    computeReconciliationLikelihood();
    
    Nhx *nhx = new Nhx();
    annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                           *_rootedTree, 
                                           _rootedTree->getRootNode(), 
                                           _seqSp, _spId); 
    cout << "Reconciled tree: "<<endl;
    nhx->write(*_rootedTree, cout);
    
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    if (rlk) {
        delete rlk;
        rlk = 0;
    }
    if (nhx) delete nhx;
    
}




/************************************************************************
 * Tries all possible (and potentially valuable) SPRs and computes the DL likelihood.
 * For the three best ones, if the DL lk is better than the current best, computes the sequence lk.
 * Repeat until no SPR can improve the gene tree. 
 * This algorithm is slower and less efficient than the other one.
 ************************************************************************/
void COALGeneTreeLikelihood::refineGeneTreeSPRs2(map<string, string> params) { 
    /*  std::cout <<"\t\t\tStarting MLSearch : current tree : "<< std::endl;
     std::cout<< TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;*/
    if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
        //We don't do SPRs
        //  std::cout << "WE DONT DO SPRS"<<std::endl;
        computeReconciliationLikelihood();
        return;
    }
    //  std::cout << "WE DO SPRS"<<std::endl;
    
    std::vector <int> nodeIdsToRegraft;
    bool betterTree = false;
    TreeTemplate<Node> * treeForSPR = 0;
    TreeTemplate<Node> * bestTree = 0;
    // if (getLogLikelihood()==UNLIKELY) 
    computeReconciliationLikelihood();
    //    ReconciliationTreeLikelihood * bestTreeLogLk = this->clone();
    double logL = getLogLikelihood();
    
    std::map <double, std::vector<int> > allDLLKsAndIndices;
    std::vector< TreeTemplate<Node>* > allTreesTested;
    unsigned int index = 0;
    
    optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (nniLk_),
                             0.1, params);
    logL = getLogLikelihood();
    
    double bestlogL = logL;
    double candidateScenarioLk ;
    double bestSequenceLogL = getSequenceLikelihood();
    double bestScenarioLk = getScenarioLikelihood();
    // std::cout << "bestLogL: "<< bestlogL << "bestScenarioLK: "<< bestScenarioLk <<"; bestSequenceLK: "<<bestSequenceLogL << std::endl;
    unsigned int numIterationsWithoutImprovement = 0;
    NNIHomogeneousTreeLikelihood * drlk = 0;
    NNIHomogeneousTreeLikelihood * bestDrlk = 0;
    breadthFirstreNumber (*_rootedTree);
    
    
    // DRHomogeneousTreeLikelihood drlk;
    //  std::cout<< "Starting tree: "<<TreeTemplateTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    
    string parentDup;
    string nodeDup;
    string numLoss = "0";
    
    bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    
    
    while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    {
        betterTree=false;
        allDLLKsAndIndices.clear();
        allTreesTested.clear();        
        index = 0;
        annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                               *_rootedTree, 
                                               _rootedTree->getRootNode(), 
                                               _seqSp, _spId); 
        for (unsigned int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
        {
            Node * n = _rootedTree->getNode(nodeForSPR);
            /*            if (n->getFather()->hasBranchProperty("Ev")) {
             parentDup = (dynamic_cast<const BppString *>(n->getFather()->getBranchProperty("Ev")))->toSTL() ;
             }*/
            if (n->hasBranchProperty("L")) {
                numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
            }
            if ( numLoss != "0" ) {
                buildVectorOfRegraftingNodesCoalGeneTree(*_spTree, *_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
                {
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    treeForSPR = _rootedTree->clone();
                    
                    makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false);
                    
                    //breadthFirstreNumber (*treeForSPR);
                    
                    //Compute the COAL likelihood
                    candidateScenarioLk =  findMLCoalReconciliationDR (_spTree, treeForSPR, 
                                                                       _seqSp, _spId, 
                                                                       _coalBl, 
                                                                       _tentativeMLindex, 
                                                                       _tentativeCoalCounts, 
                                                                       _tentativeNodesToTryInNNISearch, false); 
                    if (allDLLKsAndIndices.find(- candidateScenarioLk) != allDLLKsAndIndices.end() ) {
                        allDLLKsAndIndices[- candidateScenarioLk].push_back(index);
                    }
                    else {
                        allDLLKsAndIndices[- candidateScenarioLk] = vector<int> ();
                        allDLLKsAndIndices[- candidateScenarioLk].push_back(index);
                    }
                    allTreesTested.push_back (treeForSPR->clone());
                    index = index +1;
                    // logL =logL - 10;
                }
            } // End if ( numLoss != "0" )
            else {
                numIterationsWithoutImprovement++;
            }
        } //End for (int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--)
        
        unsigned int numberOfTopologiesToTry = 3;
        
        vector <int> treeIndicesToTry;
        vector <double> scenarioLksToTry;
        for(std::map<double, std::vector<int> >::iterator it = allDLLKsAndIndices.begin(); it != allDLLKsAndIndices.end(); it++){
            // std::cout << "Tested scenario Lk: "<< it->first <<" Found for "<<it->second.size()<<" topologies"<<std::endl;
            if (- it->first > bestScenarioLk) {
                for (unsigned int i = 0 ; i < it->second.size(); i++) {
                    scenarioLksToTry.push_back(- it->first);
                    treeIndicesToTry.push_back(it->second[i]);
                }
                if (treeIndicesToTry.size() >= numberOfTopologiesToTry) break;
            }
        }
        // std::cout << "treeIndicesToTry.size(): " << treeIndicesToTry.size() <<std::endl;
        /////For the top k candidate topologies, we compute the loglk
        for (unsigned int k = 0 ; k < treeIndicesToTry.size() ; k++)
        {
            if (treeForSPR) 
            {
                delete treeForSPR;
                treeForSPR = 0;
            }
            treeForSPR = allTreesTested[treeIndicesToTry[k]]->clone();
            candidateScenarioLk = scenarioLksToTry[k];
            
            
            if (computeSequenceLikelihoodForSPR) {
                if (drlk) {
                    delete drlk;
                    drlk = 0;
                }
                // std::cout << "COMPUTING SEQLK: " << TreeTemplateTools::treeToParenthesis(*treeForSPR, true) <<std::endl;
                drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                drlk->initialize();
                optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (drlk),
                                         0.1, params);
                logL = candidateScenarioLk - drlk->getValue();
            }
            else {
                logL = candidateScenarioLk - bestSequenceLogL;
            }
            //If the candidate tree has a DL + sequence Lk better than the current best
            if (logL - 0.01 > bestlogL) 
            {
                betterTree = true;
                bestlogL =logL;
                bestScenarioLk = candidateScenarioLk;
                if (computeSequenceLikelihoodForSPR) {
                    bestSequenceLogL = drlk->getValue();
                    
                    if (bestDrlk) {
                        delete bestDrlk;
                        bestDrlk = 0;
                    }
                    
                    bestDrlk = drlk->clone();
                }
                if (bestTree) {
                    delete bestTree;
                    bestTree = 0;
                }
                
                bestTree = treeForSPR->clone();
            }
            else {
                //   copyContentsFrom(*bestTreeLogLk);
                //  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
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
            if (_rootedTree) 
            {
                delete _rootedTree;
                _rootedTree = 0;
            }
            _rootedTree = bestTree->clone();
            breadthFirstreNumber (*_rootedTree);
            
            if (bestTree) {
                delete bestTree;
                bestTree = 0;
            }
            if (computeSequenceLikelihoodForSPR) {                        
                if (nniLk_) {
                    delete nniLk_;
                    nniLk_ = 0;
                }
                nniLk_ = bestDrlk->clone();                        
            }
        }
        else 
        {
            numIterationsWithoutImprovement = _rootedTree->getNumberOfNodes() - 2;
        }
        
        if (treeForSPR) 
        {
            delete treeForSPR;
            treeForSPR = 0;
        }
        if (bestTree) {
            delete bestTree;
            bestTree = 0;
        }
        if (bestDrlk) {
            delete bestDrlk;
            bestDrlk = 0;
        }
        if (drlk) {
            delete drlk;
            drlk = 0;
        }
    }// End while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    
    _scenarioLikelihood = bestScenarioLk;
    //final branch lengths optimization
    
    _rootedTree->resetNodesId();
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    
    bestTree = _rootedTree->clone();
    
    drlk = new NNIHomogeneousTreeLikelihood (*bestTree, 
                                             *(nniLk_->getData()), 
                                             nniLk_->getSubstitutionModel(), 
                                             nniLk_->getRateDistribution(), 
                                             true, false);
    drlk->initialize();
    
    int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
    bool backupOpt = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
    params[ std::string("optimization.max_number_f_eval")] = 100;
    params[ std::string("optimization.topology")] = "false";
    PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
    params[ std::string("optimization.max_number_f_eval")] = backup;
    params[ std::string("optimization.topology")] = backupOpt;
    
    if (nniLk_) {
        delete nniLk_;
        nniLk_ = 0;
    }
    
    nniLk_ = drlk->clone();
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
    
    vector<int> drlkNodes = bestTree->getNodesId();
    for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
        if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
            _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
    }
    
    //One more reconciliation, to update the "_num*Lineages" vectors.
    
    computeReconciliationLikelihood();
    
    Nhx *nhx = new Nhx();
    annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                           *_rootedTree, 
                                           _rootedTree->getRootNode(), 
                                           _seqSp, _spId); 
    cout << "Reconciled tree: "<<endl;
    nhx->write(*_rootedTree, cout);
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    if (drlk) {
        delete drlk;
        drlk = 0;
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
                nniLk_->topologyChangeTested( 	TopologyChangeEvent ());
                nniLk_->topologyChangeSuccessful( 	TopologyChangeEvent ());
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
    vector<string> names = TreeTemplateTools::getLeavesNames(*(_geneTreeWithSpNames->getRootNode() ) );
    vector<string> uniqueNames = VectorTools::unique(names);
    if (uniqueNames.size() == names.size() ) {
        return true;
    }
    return false;
}




