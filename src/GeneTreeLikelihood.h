/*
 * Copyright or © or Copr. Centre National de la Recherche Scientifique
 * contributor : Bastien Boussau (2009-2013)
 * 
 * bastien.boussau@univ-lyon1.fr
 * 
 * This software is a bioinformatics computer program whose purpose is to
 * simultaneously build gene and species trees when gene families have
 * undergone duplications and losses. It can analyze thousands of gene
 * families in dozens of genomes simultaneously, and was presented in
 * an article in Genome Research. Trees and parameters are estimated
 * in the maximum likelihood framework, by maximizing the probability
 * of alignments given the species tree, the gene trees and the parameters
 * of duplication and loss.
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


#ifndef GeneTreeLikelihood_h
#define GeneTreeLikelihood_h

#include <exception>

#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Nhx.h>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/AutoParameter.h>

#include "FastRHomogeneousTreeLikelihood.h"
#include "ReconciliationTools.h"
#include "GeneTreeAlgorithms.h"
#include "LikelihoodEvaluator.h"

#include <Bpp/Text/StringTokenizer.h>

#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>


namespace mpi = boost::mpi;

/**
 * @brief This class adds support for reconciliation to a species tree to the NNIHomogeneousTreeLikelihood class.
 */
class GeneTreeLikelihood
{
protected:
  /***
   * The new implementation of a likelihood estimator
   */
  LikelihoodEvaluator * levaluator_;
  
  //  bpp::TreeTemplate<bpp::Node> * _tree;
  bpp::TreeTemplate<bpp::Node> * spTree_;
  bpp::TreeTemplate<bpp::Node> * rootedTree_;
  bpp::TreeTemplate<bpp::Node> * geneTreeWithSpNames_;
  std::map <std::string, std::string> seqSp_; //link between sequence and species
  std::map <std::string, int> spId_;
  std::set <int> nodesToTryInNNISearch_;
  double scenarioLikelihood_;
  //  mutable double _sequenceLikelihood;
  int MLindex_;
  bool rootOptimization_;
  mutable std::set <int> tentativeNodesToTryInNNISearch_;
  mutable int tentativeMLindex_;
  mutable double tentativeScenarioLikelihood_;
  mutable int totalIterations_;
  mutable int counter_;
  mutable std::vector <int> listOfPreviousRoots_;
  int _speciesIdLimitForRootPosition_;
  mutable bool optimizeSequenceLikelihood_;
  mutable bool optimizeReconciliationLikelihood_;
  mutable bool considerSequenceLikelihood_;
  //unsigned int sprLimit_;
  std::map <std::string, std::string > params_;
  unsigned int sprLimitGeneTree_;
  double timeLimit_;
  double elapsedTime_;
  
public:
  
  GeneTreeLikelihood(); 
  
  /**
   * @brief Build a new DLGeneTreeLikelihood object.
   *
   * @param params The parameters to parse.
   * @param spTree The species tree
   * @throw Exception if an error occured.
   */
  GeneTreeLikelihood(std::string file , map<string, string> params, bpp::TreeTemplate<bpp::Node> & spTree) throw (exception);
  
  
  
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
   * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
   * @param MLindex ML rooting position     
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @throw Exception if an error occured.
   */
  GeneTreeLikelihood(
    const Tree & tree,
    const SiteContainer & data,
    SubstitutionModel * model,
    DiscreteDistribution * rDist,
    bpp::TreeTemplate<bpp::Node> & spTree,
    bpp::TreeTemplate<bpp::Node> & rootedTree,
    bpp::TreeTemplate<bpp::Node> & geneTreeWithSpNames,
    const std::map <std::string, std::string> seqSp,
    std::map <std::string,int> spId,
    int speciesIdLimitForRootPosition,
    int & MLindex,
    std::map <std::string, std::string > params,
    bool checkRooted = true,
    bool verbose = false,
    bool rootOptimization = false,
    bool considerSequenceLikelihood = true,
    unsigned int sprLimitGeneTree = 2)
  throw (Exception);
  
  
  // unload likelihood evaluator
  void unload();
  
  /**
   * @brief Copy constructor.
   */ 
  GeneTreeLikelihood(const GeneTreeLikelihood & lik);
  
  GeneTreeLikelihood & operator=(const GeneTreeLikelihood & lik);
  
  virtual ~GeneTreeLikelihood() {};
  
  
  
  #ifndef NO_VIRTUAL_COV
  GeneTreeLikelihood*
  #else
  bpp::Clonable*
  #endif
  clone() const { return new GeneTreeLikelihood(*this); }
  
  double getScenarioLikelihood() const throw (Exception) { return scenarioLikelihood_; }
  
  void setSpTree(bpp::TreeTemplate<bpp::Node> & spTree) { if (spTree_) delete spTree_; spTree_ = spTree.clone(); }
  
  void setSpId(std::map <std::string, int> & spId) {spId_ = spId;}
  
  //   ParameterList getParameters() {return nniLk_->getParameters();}
  
  bpp::TreeTemplate<bpp::Node> & getSpTree() const {return *spTree_;}
  
  bpp::TreeTemplate<bpp::Node> & getRootedTree() const {return *rootedTree_;}
  
  bpp::TreeTemplate<bpp::Node> & getGeneTreeWithSpNames() const {return *geneTreeWithSpNames_;}
  
  std::map <std::string, std::string> getSeqSp() {return seqSp_;}
  
  void OptimizeSequenceLikelihood(bool yesOrNo) const  {
    optimizeSequenceLikelihood_ = yesOrNo;
  }
  
  void OptimizeReconciliationLikelihood(bool yesOrNo) const {
    optimizeReconciliationLikelihood_ = yesOrNo;
  }
  
  LikelihoodEvaluator* getSequenceLikelihoodObject() const {
    return levaluator_;
  }
  
  unsigned int getSprLimitGeneTree() const {
    return sprLimitGeneTree_; 
  }
  
  
  bool isInitialized() {
    return levaluator_->isInitialized();
  }
  unsigned int seqsToRemove();
  
  void setGeneTree(bpp::TreeTemplate<bpp::Node>* tree, bpp::TreeTemplate<bpp::Node>* rootedTree) ;
  
  std::map <std::string, std::string > getParams () {
    return params_;
  }
  
  std::string getLikelihoodMethod () {
    return ApplicationTools::getStringParameter("likelihood.evaluator", params_, "PLL");;
  }
  
  
};

#endif
