/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a bioinformatics computer program whose purpose is to
simultaneously build gene and species trees when gene families have
undergone duplications and losses. It can analyze thousands of gene
families in dozens of genomes simultaneously, and was presented in
an article in Genome Research. Trees and parameters are estimated
in the maximum likelihood framework, by maximizing theprobability
of alignments given the species tree, the gene trees and the parameters
of duplication and loss.

This software is governed by the CeCILL license under French law and
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


#ifndef LikelihoodEvaluator_hpp
#define LikelihoodEvaluator_hpp

#include<string>

#include<Bpp/Phyl/Node.h>
#include<Bpp/Phyl/TreeTemplate.h>
#include<Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>


extern "C" {
#include <pll/pll.h>
}



class LikelihoodEvaluator:
public bpp::Clonable
{
  
private:
  
  
  
  std::string treeFile;
  std::string alignmentFile;
  
  
  /** @name Data and commands for PLL
  *  PLL specific objects and routine calls.
  */
  ///@{
  
  /**
  * PLL attributes. Example:
  * attr.rateHetModel     = PLL_GAMMA;
    attr.fastScaling      = PLL_FALSE;
    attr.saveMemory       = PLL_FALSE;
    attr.useRecom         = PLL_FALSE;
    attr.randomNumberSeed = 0xDEADBEEF;
  */
  pllInstanceAttr attr_PLL;
  
  /**
  The PLL tree.
  */
  pllInstance * tr_PLL;
  
  /**
  The PLL alignment.
  */
  pllAlignmentData * alignmentData_PLL;
  
  /**
  The PLL newick data.
  */
  pllNewickTree * newick_PLL;
  
  /**
  The PLL partitions list.
  */
  partitionList * partitions_PLL;
  
  /**
  The PLL partition information.
  */  
  pllQueue * partitionInfo_PLL;
    
  
  /**
  Loads the PLL alignment
  */
  loadPLLalignment(char* path);
  
  
  /**
  Loads the PLL tree.
  */
  loadPLLnewick(char* path);
  
  /**
  Initializing partitions for PLL.
  */
  loadPLLpartitions(char* path);
  
  /**
  Initializing PLL tree: tr_PLL.
  */
  initializePLLtree();
  
  /**
  Updating the PLL tree with PLL newick.
  */
  updatePLLtreeWithPLLnewick();
  
  
  
  ///@}
  
  
  /** @name Translation data
  *  Data used for conversion purposes BPP <-> PLL
  */
  ///@{
  
  /**
  Alignment strictly formatted for PLL
  */
  std::string fastaForPLL;
  
  /**
  Newick strictly formatted for PLL
  */
  std::string newickForPLL;
  
  /**
  Defines the real sequences names to simplified ones for PLL
  */
  std::map<std::string,std::string> realToStrict;
  
  /**
  Defines simplified names for PLL to real sequence ones
  */
  std::map<std::string,std::string> strictToReal;

  ///@}
  
  
  /** @name Original BPP Data
  *  BPP data are used as reference in this wrapper
  */
  ///@{
  
  /**
  Alignment strictly formatted for PLL
  */
  bpp::TreeTemplate<bpp::Node> * tree;
  
  /**
  Newick strictly formatted for PLL
  */
  std::string newickForPLL;
  
  /**
  Defines the real sequences names to simplified ones for PLL
  */
  std::map<std::string,std::string> realToStrict;
  
  /**
  Defines simplified names for PLL to real sequence ones
  */
  std::map<std::string,std::string> strictToReal;

  ///@}
  
  
  /** @name Likelihood management with BPP
  *  BPP methods using NNIHomogeneousTreeLikelihood
  */
  ///@{
  
  bpp::NNIHomogeneousTreeLikelihood * nniLk;
  
  
  /**
   * @brief initialize the NNIHomogeneousTreeLikelihood object for BPP management
   */
  initialize_BPP_nniLk();
  
  /**
   * @brief substitution model managed by BPP
   */
  bpp::SubstitutionModel * substitutionModel;
  
  /**
   * @brief discrete distribution managed by BPP
   * GeneTreeLikelihood correspondance: *rDist* The rate across sites distribution to use
   */
  bpp::DiscreteDistribution * rateDistribution;
  
  /**
   * @brief discrete distribution managed by BPP
   */
  bool mustUnrootTrees;
  
  
  
  ///@}
  
  
  bool verbose;
  
  
  
public:
  
  enum LikelihoodMethod{PLL,BPP};
  
  
  
  /** @name Constructors
  */
  ///@{
  
  /**
  * @brief nnilk like contructor number 1
  */
  LikelihoodEvaluator(bpp::TreeTemplate<bpp::Node> * tree, bpp::SubstitutionModel* model, bpp::DiscreteDistribution * rateDistribution, bool mustUnrootTrees, bool verbose=false);
  
  
  /**
  * @brief nnilk like contructor number 2
  */
  LikelihoodEvaluator(bpp::TreeTemplate<bpp::Node> * tree, bpp::SiteContainer data, bpp::SubstitutionModel* model, bpp::DiscreteDistribution * rateDistribution, bool mustUnrootTrees, bool verbose=false);
  
  /**
  * @brief empty contructor
  */
  LikelihoodEvaluator();
  
  /**
  * @brief copy from another object
  */
  LikelihoodEvaluator(LikelihoodEvaluator const &leval);
  
  ///@}
  
  /**
  * @brief returns the Log Likelihood computed with the default method
  * @return double log likelihood
  */
  double getLogLikelihood();
  
  /**
  * @brief returns the Log Likelihood computed with specified method
  * @param likelihood the likelihood method
  * @return double log likelihood
  */
  double getLogLikelihood(LikelihoodEvaluator::LikelihoodMethod likelihoodMethod);
  
  /**
  * @brief returns the params of the nnilk member
  * @return parameters
  */
  bpp::ParameterList getParameters();
  
  
  LikelihoodEvaluator* clone();
  
  /**
  * @brief Temporary way to directly access the associated nniLk
  * @return a pointer to the associated nniLk of this object
  */
  bpp::NNIHomogeneousTreeLikelihood * getnniLk();
  
  
};

#else

class LikelihoodEvaluator;

#endif