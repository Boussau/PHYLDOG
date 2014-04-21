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


/** TODO
 * 
 * - déraciner les arbres proposés.
 */

/**
 * A little explanation about how this class works. It is a replacement to
 * NNIHomogeneousTreeLikelihood.
 * First step: one constructs the object with the parameters (tree and
 * alignment path).
 * Second step: one can modify the tree and/or the alignment accordingly
 * to one’s criterias, then the estimator has to be initialized.
 * 
 * Once initialized, there always is a computed likelihood for the main tree
 * and if available, for the alternative tree.
 * 
 */


public bpp::Clonable
{
  
  
public:
  
  enum LikelihoodMethod{PLL,BPP};
  
  
private:
  
  LikelihoodMethod method;
    
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
  PLL_loadAlignment(char* path);
  
  
  /**
  Loads the PLL tree.
  */
  PLL_loadNewick(char* path);
  
  /**
  Initializing partitions for PLL.
  */
  PLL_loadPartitions(char* path);
  
  /**
  Initializing PLL tree: tr_PLL.
  */
  PLL_initializeTree();
  
  /**
  Updating the PLL tree with PLL newick.
  */
  updatePLLtreeWithPLLnewick();
  
  
  
  ///@}
  
  
  /** @name Translation data
  *  Data used for conversion purposes BPP <-> PLL
  * 
  * Here, strict data is a copy of data, designed for pll: just alphanum chars in names
  */
  ///@{
  
  /**
  Main tree, strict version
  */
  bpp::TreeTemplate<bpp::Node>* strictTree;
  
  /**
  Alternative tree, strict version
  */
  bpp::TreeTemplate<bpp::Node>* strictAlternativeTree;
  
  /**
  Alignment strictly formatted for PLL
  */
  VectorSiteContainer* fastaForPLL;
  
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
  

  /**
  * Loads the names of the sequences and fill the strict vector
  */
  void loadStrictNamesFromAlignment();

  /**
  * Assign the strict name to the leaves of the tree
  * @param targetTree modify the tree in place and set the leaves names to their
  * strict version.
  */
  void convertTreeToStrict(bpp::TreeTemplate< bpp::Node >& targetTree);
  
  /**
  * Assign the strict name to the leaves of the alternative tree
  */
  void updateStrictTree();
  
  
  /**
  * Assign the strict name to the leaves of the alternative tree
  */
  void updateStrictAlternativeTree();
  
  /**
  * Initialize PLL with the right data
  */
  void initialize_PLL();
  
  /**
   * Writes alignment files for PLL: the alignment and the partition file
   * @param prefix the family name, so many alignment/partitions files can
   * 			co-exist in the same working directory.
   */
  void writeAlignmentFilesForPLL(std::string prefix);
  
  /**
   * Transform the tree to a clean string, send it to PLL and get the LLK
   * then destroy temporary objects.
   * @param prefix a BPP tree
   */
  double optimizeBranchLengthAndReturnLogLikelihood(bpp::TreeTemplate<bpp::Node>);

  

  ///@}
  
  
  /** @name Internal BPP Data
  *  BPP data are used as reference in this wrapper
  */
  ///@{
  
  /**
  The BPP tree
  */
  bpp::TreeTemplate<bpp::Node> * tree;
  
  /**
  An alternative tree to be tested
  */
  bpp::TreeTemplate<bpp::Node> * alternativeTree;
  
  /**
  Map of parameters string -> string, BPP style
  */
  std::map<std::string, std::string> params;
  
  /**
  last computed likelihood whatever the method was
  */
  double lastLikelihood;
  
  /**
  last computed likelihood for the alternative tree whatever the method was
  */
  double alternativeLikelihood;
  
  ///@}
  
  
  /** @name Likelihood management with BPP
  *  BPP methods using NNIHomogeneousTreeLikelihood
  */
  ///@{
  /**
    * @brief NniLk of the tree
    */
  bpp::NNIHomogeneousTreeLikelihood * nniLk;
  
  /**
  * @brief NniLk of the alternative tree
  */
  bpp::NNIHomogeneousTreeLikelihood * nniLkAlternative;
  
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
  
  /**
   * @brief sites of the sequences used by bpp
   */
  bpp::VectorSiteContainer * sites;
  
  /**
   * @brief alphabet of the sequences used by bpp
   */
  bpp::Alphabet * alphabet;
  
  
  
  ///@}
  
  
  bool verbose;
  
  bool initialized;
  
  
  
public:
  
  enum LikelihoodMethod{PLL,BPP};
  
  /** @name Accessors
  * let the outside acces to our data in order to perform filtering operations
  */
  ///@{
  
  
  /**
   * @brief alphabet of the sequences used by bpp
   */
  bpp::Alphabet* getAlphabet();
  
  /**
  * @brief sites of the sequences used by bpp
  */
  bpp::VectorSiteContainer* getSites();
  
  /**
   * @brief discrete distribution managed by BPP
   * GeneTreeLikelihood correspondance: *rDist* The rate across sites distribution to use
   */
  bpp::DiscreteDistribution * getRateDistribution();
  
  /**
   * @brief substitution model managed by BPP
   */
  bpp::SubstitutionModel * getSubstitutionModel();
  
  /**
   * @brief tree, managed by BPP
   */
  bpp::TreeTemplate<N> * getTree();
  
  ///@}
  
  
  /** @name Alternative Tree
  * managing and testing the alternative tree
  */
  ///@{
  
  /**
   * @brief alternative tree, managed by BPP
   */
  bpp::TreeTemplate<N> * getAlternativeTree();
  
  /**
  * @brief set the alternative tree to a new one
  */
  void setAlternativeTree(bpp::TreeTemplate<N>* newAlternative);
  

  /**
  * @brief get the likelihood of the alternative tree
  */
  double getAlternativeLogLikelihood();
  
  /**
  * @brief replace the main tree by the alternative one
  */
  void acceptAlternativeTree();
  

  
  
   
  
  
  ///@}
  
  /** @name Constructors, Destructor and validation functions
  */
  ///@{
  
  /**
  * @brief nnilk like contructor number 1
  */
  LikelihoodEvaluator(bpp::TreeTemplate<bpp::Node> * tree, bpp::SubstitutionModel* model, bpp::DiscreteDistribution * rateDistribution, bool mustUnrootTrees, bool verbose=false);
  
  
  /**
  * @brief nnilk like contructor number 2
  */
  LikelihoodEvaluator(bpp::TreeTemplate<bpp::Node> * tree, bpp::SiteContainer* data, bpp::SubstitutionModel* model, bpp::DiscreteDistribution * rateDistribution, bool mustUnrootTrees, bool verbose=false);
  
  /**
  * @brief empty contructor
  */
  LikelihoodEvaluator();
  
  /**
  * @brief default constructor with params
  * @param params map of parameters (string -> string), BPP style
  */
  LikelihoodEvaluator(std::map<std::string, std::string> params);
  
  
  /**
  * @brief copy from another object
  */
  LikelihoodEvaluator(LikelihoodEvaluator const &leval);
  
  /**
  * @brief destuctor
  */
  ~LikelihoodEvaluator();
  
  /**
  @brief this method is triggered once the tree and the aligment
  * have been (optionaly) filtered by the extrenal accessors
  * Once executed, it is not possible to modify the alignment.
  */
  void initialize();
  
  
  /**
  @brief has initializeEvaluator() been launched?
  @return the boolean value
  */
  bool isInitialized();
  
  
  
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
  

  
  LikelihoodEvaluator* clone();
  
  
  
  
};

#else

class LikelihoodEvaluator;

#endif