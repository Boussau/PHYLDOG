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

    std::vector<double> duplicationProbabilities_;
    std::vector<double> lossProbabilities_;
    std::vector<double> backupDuplicationProbabilities_;
    std::vector<double> backupLossProbabilities_;
    std::vector <int> num0Lineages_;
    std::vector <int> num1Lineages_; 
    std::vector <int> num2Lineages_;
    TreeTemplate<Node> * tree_;
    std::string currentSpeciesTree_;
    int index_;
    int bestIndex_;
    bool stop_;
    double logL_;
    double bestlogL_;
    std::vector <int> bestNum0Lineages_; 
    std::vector <int> bestNum1Lineages_; 
    std::vector <int> bestNum2Lineages_; 
    std::vector <std::vector<int> > allNum0Lineages_;
		std::vector <std::vector<int> > allNum1Lineages_;
		std::vector <std::vector<int> > allNum2Lineages_;
    bool rearrange_;
    int numIterationsWithoutImprovement_;
    std::string branchProbaOptimization_;
    std::map <std::string, int> genomeMissing_;
    
    
  public:
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
                          std::vector<double> &lossProbabilities, 
                          std::vector<double> &duplicationProbabilities, 
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
    lossProbabilities_(lossProbabilities), 
    duplicationProbabilities_(duplicationProbabilities),
    backupLossProbabilities_(lossProbabilities), 
    backupDuplicationProbabilities_(duplicationProbabilities),
    rearrange_(rearrange), 
    numIterationsWithoutImprovement_(numIterationsWithoutImprovement), 
    branchProbaOptimization_(branchProbaOptimization), 
    genomeMissing_(genomeMissing)
    {
      Parameter p("coefDup", 1, &Parameter::R_PLUS_STAR);
      addParameter_(p);
      Parameter p2("coefLoss", 1, &Parameter::R_PLUS_STAR);
      addParameter_(p2);
    }
    
    SpeciesTreeLikelihood(const SpeciesTreeLikelihood& stl) :
    AbstractParametrizable(stl),
    world_(stl.world_), server_(stl.server_), tree_(stl.tree_), index_(stl.index_),
    bestIndex_(stl.bestIndex_), stop_(stl.stop_),
    logL_(stl.logL_), bestlogL_(stl.bestlogL_),
    num0Lineages_(stl.num0Lineages_), num1Lineages_(stl.num1Lineages_), 
    num2Lineages_(stl.num2Lineages_), 
    allNum0Lineages_(stl.allNum0Lineages_), allNum1Lineages_(stl.allNum1Lineages_), 
    allNum2Lineages_(stl.allNum2Lineages_),
    lossProbabilities_(stl.lossProbabilities_), 
    duplicationProbabilities_(stl.duplicationProbabilities_),
    backupDuplicationProbabilities_(stl.backupDuplicationProbabilities_),
    backupLossProbabilities_(stl.backupLossProbabilities_), 
    rearrange_(stl.rearrange_), 
    numIterationsWithoutImprovement_(stl.numIterationsWithoutImprovement_), 
    branchProbaOptimization_(stl.branchProbaOptimization_), 
    genomeMissing_(stl.genomeMissing_)
    {}
    
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
      lossProbabilities_ = stl.lossProbabilities_;
      duplicationProbabilities_ =stl.duplicationProbabilities_;
      backupLossProbabilities_ = stl.backupLossProbabilities_;
      backupDuplicationProbabilities_ =stl.backupDuplicationProbabilities_;
      rearrange_ = stl.rearrange_;
      numIterationsWithoutImprovement_ = stl.numIterationsWithoutImprovement_;
      branchProbaOptimization_ = stl.branchProbaOptimization_;
      genomeMissing_ = stl.genomeMissing_;
      
      return *this;
    }
    
    virtual ~SpeciesTreeLikelihood() {}
    
    SpeciesTreeLikelihood* clone() const { return new SpeciesTreeLikelihood(*this); }
    
  public:
    void setParameters(const ParameterList &parameters)
    throw (ParameterNotFoundException, ConstraintException)
    {
      setParametersValues(parameters);
    }
    
    double getValue() const throw (Exception) { return logL_; }
    
    void fireParameterChanged(const ParameterList & parameters)
    {
      updateDuplicationAndLossRates();
      computeLogLikelihood();
    }
    
    void updateDuplicationAndLossRates();
    
  protected:
    void computeLogLikelihood();
  
  };









}


#endif //_SPECIESTREELIKELIHOOD_H_
