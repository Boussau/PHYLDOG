/*
 *  SpeciesTreeLikelihood.h
 *  
 *
 *  Created by boussau on 18/02/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */

#ifndef _SPECIESTREELIKELIHOOD_H_
#define _SPECIESTREELIKELIHOOD_H_

// From NumCalc:
#include <NumCalc/Parametrizable.h>
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/Parameter.h>
#include <NumCalc/AbstractParametrizable.h>
#include <NumCalc/Functions.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/PhylogeneticsApplicationTools.h>


//From the BOOST library 
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>

#include "ReconciliationTools.h"
#include "SpeciesTreeExploration.h"


namespace mpi = boost::mpi;


namespace bpp
{
  
	class SpeciesTreeLikelihood :
  public Function, 
  public AbstractParametrizable
  {
  private:
  //For MPI communication:
  mpi::communicator world_;
  int server_;
  int size_;
  //Parameter list
  std::map<std::string, std::string> params_;
  //Information related to the clients
  std::vector <int> numbersOfGenesPerClient_;
  int assignedNumberOfGenes_;
  std::vector<std::string> assignedFilenames_;
  std::vector< std::vector<std::string> > listOfOptionsPerClient_;
  //Vectors of expected numbers of events per branch
  std::vector<double> duplicationExpectedNumbers_;
  std::vector<double> lossExpectedNumbers_;
  std::vector<double> backupDuplicationExpectedNumbers_;
  std::vector<double> backupLossExpectedNumbers_;
  //Vectors of counts of times 0, 1, 2+ lineages were found at the end of
  //a species tree branch
  std::vector <int> num0Lineages_;
  std::vector <int> num1Lineages_; 
  std::vector <int> num2Lineages_;
  std::vector <int> bestNum0Lineages_; 
  std::vector <int> bestNum1Lineages_; 
  std::vector <int> bestNum2Lineages_; 
  std::vector <std::vector<int> > allNum0Lineages_;
  std::vector <std::vector<int> > allNum1Lineages_;
  std::vector <std::vector<int> > allNum2Lineages_;
  //Species tree, and its string representation
  TreeTemplate<Node> * tree_;
  TreeTemplate<Node> * bestTree_;
  TreeTemplate<Node> * currentTree_;
  std::string currentSpeciesTree_;
  //Index of the tree under study, and index of the most likely tree
  int index_;
  int bestIndex_;
  //Whether we should rearrange the species tree
  bool optimizeSpeciesTreeTopology_; 
  //Whether we should stop tree search
  bool stop_;
  //logLk and best logLk
  double logL_;
  double bestlogL_;
  //Whether we should rearrange the gene trees
  bool rearrange_;
  //Number of iterations of the search algorithm without improvement
  int numIterationsWithoutImprovement_;
  //How far can we regraft subtrees when doing a spr
  int sprLimit_;
  //string giving the kind of optimization to perform
  std::string branchExpectedNumbersOptimization_;
  //Map giving the expected percentage of genes in genomes missing 
  //because of poor sequencing
  std::map <std::string, int> genomeMissing_;
  //Node number in the species tree
  int speciesTreeNodeNumber_;
  //vectors to keep NNI and root likelihoods, for making aLRTs.
  std::vector <double > NNILks_;
  std::vector <double > rootLks_;
  //Time limit: the program has to stop before this limit (in hours)
  int timeLimit_;
  //When the program stops, it knows at what step of the algorithm it is
  int currentStep_;
  //Suffix to add to all files output by the program
  std::string suffix_;

  public:

  //Simple constructor
  SpeciesTreeLikelihood(const mpi::communicator& world, 
                        int & server,
                        int & size, 
                        std::map<std::string, std::string> & params) :
  AbstractParametrizable(""),
  world_(world), server_(server), size_(size), params_(params)
  {
  parseOptions();
  Parameter p("coefDup", 1, &Parameter::R_PLUS_STAR);
  addParameter_(p);
  Parameter p2("coefLoss", 1, &Parameter::R_PLUS_STAR);
  addParameter_(p2);
  }
  
  
  
  //Constructor
  /*
  SpeciesTreeLikelihood(const mpi::communicator& world, 
                        int server, 
                        TreeTemplate<Node> *tree, 
                        int &index, 
                        int &bestIndex,  
                        bool &stop, 
                        double &logL, 
                        double &bestlogL, 
                        std::vector<int> &num0Lineages, 
                        std::vector<int> &num1Lineages, 
                        std::vector<int> &num2Lineages, 
                        std::vector<int> &bestNum0Lineages, 
                        std::vector<int> &bestNum1Lineages, 
                        std::vector<int> &bestNum2Lineages, 
                        std::vector< std::vector<int> > &allNum0Lineages, 
                        std::vector< std::vector<int> > &allNum1Lineages, 
                        std::vector< std::vector<int> > &allNum2Lineages, 
                        std::vector<double> &lossExpectedNumbers, 
                        std::vector<double> &duplicationExpectedNumbers, 
                        bool rearrange, int &numIterationsWithoutImprovement, 
                        std::string & branchProbaOptimization, 
                        std::map < std::string, int> genomeMissing) :
  AbstractParametrizable(""),
  world_(world), server_(server), tree_(tree), index_(index),
  bestIndex_(bestIndex), stop_(stop),
  logL_(logL), bestlogL_(bestlogL),
  num0Lineages_(num0Lineages), num1Lineages_(num1Lineages), 
  num2Lineages_(num2Lineages), 
  allNum0Lineages_(allNum0Lineages), allNum1Lineages_(allNum1Lineages), 
  allNum2Lineages_(allNum2Lineages),
  lossExpectedNumbers_(lossExpectedNumbers), 
  duplicationExpectedNumbers_(duplicationExpectedNumbers),
  backupLossExpectedNumbers_(lossProbabilities), 
  backupDuplicationExpectedNumbers_(duplicationProbabilities),
  rearrange_(rearrange), 
  numIterationsWithoutImprovement_(numIterationsWithoutImprovement), 
  branchExpectedNumbersOptimization_(branchProbaOptimization), 
  genomeMissing_(genomeMissing)
  {
  Parameter p("coefDup", 1, &Parameter::R_PLUS_STAR);
  addParameter_(p);
  Parameter p2("coefLoss", 1, &Parameter::R_PLUS_STAR);
  addParameter_(p2);
  speciesTreeNodeNumber_ = tree->getNumberOfNodes();
  NNILks_ (2*tree->getNumberOfLeaves()-2, NumConstants::VERY_BIG);
  }
  */
  
  //Copy constructor
  SpeciesTreeLikelihood(const SpeciesTreeLikelihood& stl) :
  AbstractParametrizable(stl),
  world_(stl.world_), server_(stl.server_), tree_(stl.tree_), index_(stl.index_),
  bestIndex_(stl.bestIndex_), stop_(stl.stop_),
  logL_(stl.logL_), bestlogL_(stl.bestlogL_),
  num0Lineages_(stl.num0Lineages_), num1Lineages_(stl.num1Lineages_), 
  num2Lineages_(stl.num2Lineages_), 
  allNum0Lineages_(stl.allNum0Lineages_), allNum1Lineages_(stl.allNum1Lineages_), 
  allNum2Lineages_(stl.allNum2Lineages_),
  lossExpectedNumbers_(stl.lossExpectedNumbers_), 
  duplicationExpectedNumbers_(stl.duplicationExpectedNumbers_),
  backupDuplicationExpectedNumbers_(stl.backupDuplicationExpectedNumbers_),
  backupLossExpectedNumbers_(stl.backupLossExpectedNumbers_), 
  rearrange_(stl.rearrange_), 
  numIterationsWithoutImprovement_(stl.numIterationsWithoutImprovement_), 
  branchExpectedNumbersOptimization_(stl.branchExpectedNumbersOptimization_), 
  genomeMissing_(stl.genomeMissing_), 
  speciesTreeNodeNumber_(stl.speciesTreeNodeNumber_), NNILks_(stl.NNILks_),
  rootLks_(stl.rootLks_), timeLimit_(stl.timeLimit_), currentStep_(stl.currentStep_),
  suffix_(stl.suffix_)
  {}
  
  
  //= operator
  SpeciesTreeLikelihood& operator=(const SpeciesTreeLikelihood& stl)
  {
  AbstractParametrizable::operator=(stl);
  world_ = stl.world_;
  server_ = stl.server_;
  tree_ = stl.tree_;
  index_ = stl.index_;
  bestIndex_ = stl.bestIndex_;
  stop_ = stl.stop_;
  logL_ = stl.logL_;
  bestlogL_ = stl.bestlogL_;
  num0Lineages_ = stl.num0Lineages_;
  num1Lineages_ = stl.num1Lineages_;
  num2Lineages_ = stl.num2Lineages_;
  allNum0Lineages_ = stl.allNum0Lineages_;
  allNum1Lineages_ = stl.allNum1Lineages_;
  allNum2Lineages_ = stl.allNum2Lineages_;
  lossExpectedNumbers_ = stl.lossExpectedNumbers_;
  duplicationExpectedNumbers_ =stl.duplicationExpectedNumbers_;
  backupLossExpectedNumbers_ = stl.backupLossExpectedNumbers_;
  backupDuplicationExpectedNumbers_ =stl.backupDuplicationExpectedNumbers_;
  rearrange_ = stl.rearrange_;
  numIterationsWithoutImprovement_ = stl.numIterationsWithoutImprovement_;
  branchExpectedNumbersOptimization_ = stl.branchExpectedNumbersOptimization_;
  genomeMissing_ = stl.genomeMissing_;
  speciesTreeNodeNumber_ = stl.speciesTreeNodeNumber_;
  NNILks_ = stl.NNILks_;
  rootLks_ = stl.rootLks_;
  timeLimit_ = stl.timeLimit_;
  currentStep_ = stl.currentStep_;
  suffix_ = stl.suffix_;
  return *this;
  }
  
  //Destructor
  virtual ~SpeciesTreeLikelihood() 
  {			
    delete currentTree_;
    delete bestTree_;
    delete tree_;
  }
  
  //Clone function
  SpeciesTreeLikelihood* clone() const { return new SpeciesTreeLikelihood(*this); }
  
  public:
  //Function useful to optimize the parameters of a species tree through Bio++ optimization functions.
  void setParameters(const ParameterList &parameters)
  throw (ParameterNotFoundException, ConstraintException)
  {
  setParametersValues(parameters);
  }
  
  //Get the logL of the species tree
  double getValue() const throw (Exception) { return logL_; }
  
  //If a parameter has changed
  void fireParameterChanged(const ParameterList & parameters)
  {
  updateDuplicationAndLossExpectedNumbers();
  computeLogLikelihood();
  }
  
  //Updates the parameters of the DL process.
  void updateDuplicationAndLossExpectedNumbers();
  
  //Initializes various fields in the species tree
  void initialize();
  
  //Does a ML search for the best species tree
  void MLsearch();
  
  protected:
  //Computes the loglk of the species tree
  void computeLogLikelihood();
  //Parses the options and builds the SpeciesTreeLikelihoodObject
  void parseOptions();
  };
  
  
  
  
  
  
  
  
  
}


#endif //_SPECIESTREELIKELIHOOD_H_
