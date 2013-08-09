//
//  GeneTreeLikelihood.cpp
//  phyldog
//
//  Created by Bastien Boussau on 15/04/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "GeneTreeLikelihood.h"

using namespace bpp;


GeneTreeLikelihood::GeneTreeLikelihood() {
    _totalIterations = 0;
    _counter = 0;
}


/**
 * @brief Build a new ReconciliationTreeLikelihood object.
 *
 * @param tree The tree to use.
 * @param model The substitution model to use.
 * @param rDist The rate across sites distribution to use.
 * @param spTree The species tree
 * @param rootedTree rooted version of the gene tree
 * @param seqSp link between sequence and species names
 * @param spId link between species name and species ID
 * @param lossNumbers vector to store loss numbers per branch
 * @param lossProbabilities vector to store expected numbers of losses per branch
 * @param duplicationNumbers vector to store duplication numbers per branch
 * @param duplicationProbabilities vector to store expected numbers of duplications per branch
 * @param branchNumbers vector to store branch numbers in the tree
 * @param num0Lineages vectors to store numbers of branches ending with a loss
 * @param num1Lineages vectors to store numbers of branches ending with 1 gene
 * @param num2Lineages vectors to store numbers of branches ending with 2 genes
 * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
 * @param heuristicsLevel type of heuristics used
 * @param MLindex ML rooting position
 * @param checkRooted Tell if we have to check for the tree to be unrooted.
 * If true, any rooted tree will be unrooted before likelihood computation.
 * @param verbose Should I display some info?
 * @throw Exception in an error occured.
 */
GeneTreeLikelihood::GeneTreeLikelihood(
                   const Tree & tree,
                   SubstitutionModel * model,
                   DiscreteDistribution * rDist,
                   TreeTemplate<Node> & spTree,  
                   TreeTemplate<Node> & rootedTree, 
                   TreeTemplate<Node> & geneTreeWithSpNames,
                   const std::map <std::string, std::string> seqSp,
                   std::map <std::string,int> spId,
                   int speciesIdLimitForRootPosition,
                   int heuristicsLevel,
                   int & MLindex, 
                   bool checkRooted,
                   bool verbose,
                   bool rootOptimization, 
                   bool considerSequenceLikelihood, 
                   unsigned int sprLimit)
throw (Exception):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp(seqSp), _spId(spId)
{
    
    nniLk_ = new NNIHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose); 
    _spTree = spTree.clone();
    _rootedTree = rootedTree.clone();
    _geneTreeWithSpNames = geneTreeWithSpNames.clone();
    _scenarioLikelihood = UNLIKELY;
    // _sequenceLikelihood = UNLIKELY;
    _MLindex = MLindex;
    _rootOptimization = rootOptimization; 
    _tentativeMLindex = MLindex;
    _totalIterations = 0;
    _counter = 0;
    _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
    _heuristicsLevel = heuristicsLevel;
    _optimizeSequenceLikelihood = true;
    _optimizeReconciliationLikelihood = true;
    _considerSequenceLikelihood = considerSequenceLikelihood;
    sprLimit_ = sprLimit;
    // _listOfPreviousRoots = new std::vector <int> ();
}


/**
 * @brief Build a new ReconciliationTreeLikelihood object.
 *
 * @param tree The tree to use.
 * @param data Sequences to use.
 * @param model The substitution model to use.
 * @param rDist The rate across sites distribution to use.
 * @param spTree The species tree
 * @param rootedTree rooted version of the gene tree
 * @param seqSp link between sequence and species names
 * @param spId link between species name and species ID
 * @param lossNumbers vector to store loss numbers per branch
 * @param lossProbabilities vector to store expected numbers of losses per branch
 * @param duplicationNumbers vector to store duplication numbers per branch
 * @param duplicationProbabilities vector to store expected numbers of duplications per branch
 * @param branchNumbers vector to store branch numbers in the tree
 * @param num0Lineages vectors to store numbers of branches ending with a loss
 * @param num1Lineages vectors to store numbers of branches ending with 1 gene
 * @param num2Lineages vectors to store numbers of branches ending with 2 genes
 * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
 * @param heuristicsLevel type of heuristics used
 * @param MLindex ML rooting position     
 * @param checkRooted Tell if we have to check for the tree to be unrooted.
 * If true, any rooted tree will be unrooted before likelihood computation.
 * @param verbose Should I display some info?
 * @throw Exception in an error occured.
 */

GeneTreeLikelihood::GeneTreeLikelihood(
                   const Tree & tree,
                   const SiteContainer & data,
                   SubstitutionModel * model,
                   DiscreteDistribution * rDist,
                   TreeTemplate<Node> & spTree,  
                   TreeTemplate<Node> & rootedTree,  
                   TreeTemplate<Node> & geneTreeWithSpNames,
                   const std::map <std::string, std::string> seqSp,
                   std::map <std::string,int> spId,
                   int speciesIdLimitForRootPosition,  
                   int heuristicsLevel,
                   int & MLindex, 
                   bool checkRooted,
                   bool verbose, 
                   bool rootOptimization, 
                   bool considerSequenceLikelihood, 
                   unsigned int sprLimit)
throw (Exception):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp (seqSp), _spId(spId)
{
    nniLk_ = new NNIHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose);
    _spTree = spTree.clone();
    _rootedTree = rootedTree.clone();
    _geneTreeWithSpNames = geneTreeWithSpNames.clone();
    _scenarioLikelihood = UNLIKELY;
    _MLindex = MLindex;
    _rootOptimization = rootOptimization; 
    _tentativeMLindex = MLindex;
    _totalIterations = 0; 
    _counter = 0;
    _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
    _heuristicsLevel = heuristicsLevel;
    _optimizeSequenceLikelihood = true;
    _optimizeReconciliationLikelihood = true;
    _considerSequenceLikelihood = considerSequenceLikelihood;
    sprLimit_ = sprLimit;
}


/**
 * @brief Copy constructor.
 */ 
GeneTreeLikelihood::GeneTreeLikelihood(const GeneTreeLikelihood & lik):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp (lik._seqSp), _spId(lik._spId)
{
    nniLk_ = lik.nniLk_->clone(); 
    _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone()) ;
    _rootedTree = dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone()) ;
    _geneTreeWithSpNames = dynamic_cast<TreeTemplate<Node> *> (lik._geneTreeWithSpNames->clone()) ;
    _scenarioLikelihood = lik._scenarioLikelihood;
    _MLindex = lik._MLindex;
    _rootOptimization = lik._rootOptimization; 
    _tentativeMLindex = lik._MLindex;
    _totalIterations = lik._totalIterations;
    _counter = lik._counter;
    _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
    _heuristicsLevel = lik._heuristicsLevel;
    _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
    _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
    _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
    _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;
    _considerSequenceLikelihood = lik._considerSequenceLikelihood;
    sprLimit_ = lik.sprLimit_;
}

GeneTreeLikelihood & GeneTreeLikelihood::operator=(const GeneTreeLikelihood & lik)
{
    if (nniLk_) delete nniLk_;
    nniLk_ = lik.nniLk_->clone(); 
    if (_spTree) delete _spTree;
    _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone());
    if (_rootedTree) delete _rootedTree;
    _rootedTree= dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone());
    if (_geneTreeWithSpNames) delete _geneTreeWithSpNames;
    _geneTreeWithSpNames = dynamic_cast<TreeTemplate<Node> *> (lik._geneTreeWithSpNames->clone()) ;
    _spId = lik._spId;
    _scenarioLikelihood = lik._scenarioLikelihood;
    _MLindex = lik._MLindex;
    _rootOptimization = lik._rootOptimization;
    _tentativeMLindex = lik._MLindex;
    _totalIterations = lik._totalIterations;
    _counter = lik._counter;
    _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
    _heuristicsLevel = lik._heuristicsLevel;
    _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
    _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
    _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
    _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;
    _considerSequenceLikelihood = lik._considerSequenceLikelihood;
    sprLimit_ = lik.sprLimit_;
    return *this;
}




