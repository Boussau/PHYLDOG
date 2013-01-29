/*
 *  GeneTreeAlgorithms.h
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 29/06/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

#ifndef _GENETREEALGORITHMS_H_
#define _GENETREEALGORITHMS_H_


// From SeqLib:
/*#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
*/

// From PhylLib:
/*#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
*/

//#include <Bpp/Phyl/Io/Nhx.h>
//#include <Bpp/Phyl/Mapping.all>
#include <Bpp/Phyl/Likelihood/DRTreeLikelihoodTools.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>

//#include <Bpp/Phyl/OptimizationTools.h>


// From NumCalc:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Function.all>

// From Utils:
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/Text/KeyvalTools.h>




#include "ReconciliationTools.h"
#include "GenericTreeExplorationAlgorithms.h"

/**************************************************************************
 * This function creates a sequence tree from a species tree and a std::map 
 * containing the link between the species and their sequences.
 **************************************************************************/
TreeTemplate<Node> * buildARandomSequenceTreeFromASpeciesTree (std::map <std::string, 
                                                               std::deque<std::string> > & spSeqs, 
                                                               TreeTemplate<Node> & tree, 
                                                               std::map <std::string, std::string> & spSelectedSeq);
/**************************************************************************
 * This function creates a sequence tree from a species tree and the std::map 
 * containing the link between the species and the putative orthologous sequence.
 **************************************************************************/
TreeTemplate<Node> * buildASequenceTreeFromASpeciesTreeAndCorrespondanceMap (TreeTemplate<Node> & tree, 
                                                                             std::map <std::string, 
                                                                             std::string> & spSelectedSeq);
/******************************************************************************/
// This function refines branch lengths of a gene tree.
/******************************************************************************/
double refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                             TreeTemplate<Node>  *& unrootedGeneTree, 
                                                             VectorSiteContainer * sites, 
                                                             SubstitutionModel* model, 
                                                             DiscreteDistribution* rDist, 
                                                             string file, Alphabet *alphabet, bool mapping=false);
/******************************************************************************/
// This function maps substitutions in a gene tree.
/******************************************************************************/

/*vector< vector<unsigned int> > getCountsPerBranch(
                                                  DRTreeLikelihood& drtl,
                                                  const vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg,
                                                  SubstitutionCount *count,
                                                  bool stationarity = true,
                                                  double threshold = -1);*/
/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/

void optimizeBLMapping(
                       DRTreeLikelihood* tl,
                       double precision);

/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/
void optimizeBLMappingForSPRs(
                              DRTreeLikelihood* tl,
                              double precision, map<string, string> params);

/******************************************************************************/
// This function optimizes branch lengths in a gene tree without substitution mapping
/******************************************************************************/
void optimizeBLForSPRs(
                       DRTreeLikelihood* tl,
                       double precision, map<string, string> params);


/******************************************************************************/
// This function builds a bionj tree
/******************************************************************************/
TreeTemplate<Node>  * buildBioNJTree (std::map<std::string, std::string> & params, 
                                      SiteContainer * sites, 
                                      SubstitutionModel* model, 
                                      DiscreteDistribution* rDist, 
                                      Alphabet *alphabet);
/******************************************************************************/
// This function refines a gene tree topology and branch lengths using the PhyML 
// algorithm.
/******************************************************************************/
void refineGeneTreeUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                TreeTemplate<Node>  *& unrootedGeneTree, 
                                                VectorSiteContainer * sites, 
                                                SubstitutionModel* model, 
                                                DiscreteDistribution* rDist, 
                                                string file, 
                                                Alphabet *alphabet);


/*
string parenthesisWithSpeciesNamesToGeneTree (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp ) {  
  //,                                             std::vector<string> &geneNames) {
  
}
*/

/**************************************************************************
 * This function produces a string version of a gene tree, 
 * with gene names replaced by species names. 
 **************************************************************************/
string geneTreeToParenthesisWithSpeciesNames (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp );

/**************************************************************************
 * This function produces a gene tree from a string version in which 
 * gene names have been changed to include species names. 
 **************************************************************************/
TreeTemplate<Node> * parenthesisPlusSpeciesNamesToGeneTree (string geneTreeStr) ;
  

/**************************************************************************
 * This function produces a string version of a gene tree, 
 * with gene names changed to include species names. 
 **************************************************************************/
string geneTreeToParenthesisPlusSpeciesNames (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp );


/**************************************************************************
 * This function produces a gene tree with leaves annotated with species names.
 **************************************************************************/
void annotateGeneTreeWithSpeciesNames (TreeTemplate<Node> * geneTree,
                                       std::map<std::string, std::string > seqSp ) ;

/**************************************************************************
 * This function optimizes a gene tree based on the reconciliation score only.
 * It uses SPRs and NNIs, and calls findMLReconciliationDR to compute the likelihood.
 **************************************************************************/
double refineGeneTreeDLOnly (TreeTemplate<Node> * spTree, 
                             TreeTemplate<Node> *& geneTree, 
                             std::map<std::string, std::string > seqSp,
                             std::map<std::string, int > spID,
                             std::vector< double> &lossExpectedNumbers, 
                             std::vector < double> &duplicationExpectedNumbers, 
                             int & MLindex, 
                             std::vector <int> &num0lineages, 
                             std::vector <int> &num1lineages, 
                             std::vector <int> &num2lineages, 
                             std::set <int> &nodesToTryInNNISearch);

/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void buildVectorOfRegraftingNodesGeneTree(TreeTemplate<Node> &spTree, 
                                          TreeTemplate<Node> &tree, 
                                          int nodeForSPR, 
                                          int distance, 
                                          std::vector <int> & nodeIdsToRegraft) ;

/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void getAllCandidateBranchingPointsFromSpeciesID (TreeTemplate<Node> &tree, 
                                                  std::vector <std::string> spIds, 
                                                  std::vector <int> & allNodeIds) ;

/**************************************************************************
 * This recursive function returns node ids with a given species id.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIds(Node * node, string spId, std::vector <int> & allNodeIds);

/**************************************************************************
 * This recursive function returns node ids with a given species id, upstream from Node node.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIdsUpstream(Node * node, string spId, std::vector <int> & allNodeIds);


/**************************************************************************
 * This recursive function returns node ids sons of nodes with a given species id.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void getSonsOfNodesWithSimilarSpeciesIds(Node * node, string spId, std::vector <int> & allNodeIds) ;

/**************************************************************************
 * This function returns a vector of branching points for gene tree SPR in the coalescent framework.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void buildVectorOfRegraftingNodesCoalGeneTree(TreeTemplate<Node> &spTree, 
											  TreeTemplate<Node> &tree, 
											  int nodeForSPR, 
											  int distance, 
											  std::vector <int> & nodeIdsToRegraft) ;

/**
 * @brief Optimize branch lengths parameters of a TreeLikelihood function.
 *
 * Uses Newton's method.
 *
 * A condition over function values is used as a stop condition for the algorithm.
 *
 * @see NewtonBrentMetaOptimizer
 *
 * @param tl             A pointer toward the TreeLikelihood object to optimize.
 * @param parameters     The list of parameters to optimize. The intersection of branch length parameters and the input set will be used. Use tl->getBranchLengthsParameters() in order to estimate all branch length parameters.
 * @param target         Current maximum likelihood value. If early rounds of optimization suggest we are not going to improve upon it, stop optimization.
 * @param listener       A pointer toward an optimization listener, if needed.
 * @param tolerance      The tolerance to use in the algorithm.
 * @param tlEvalMax      The maximum number of function evaluations.
 * @param messageHandler The massage handler.
 * @param profiler       The profiler.
 * @param verbose        The verbose level.
 * @param optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
 * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
 * @throw Exception any exception thrown by the Optimizer.
 */
unsigned int optimizeBranchLengthsParameters(
													DiscreteRatesAcrossSitesTreeLikelihood* tl,
													const ParameterList& parameters,
													double target,
													OptimizationListener* listener     = 0,
													double tolerance                   = 0.000001,
													unsigned int tlEvalMax             = 1000000,
													OutputStream* messageHandler       = ApplicationTools::message,
													OutputStream* profiler             = ApplicationTools::message,
													unsigned int verbose               = 1,
													const std::string& optMethodDeriv  = OptimizationTools::OPTIMIZATION_NEWTON)
throw (Exception);


/**************************************************************************                                                                                                                                                                                                                                                                                               *
  **************************************************************************/

void editDuplicationNodesMuffato(TreeTemplate<Node> & spTree, 
				 TreeTemplate<Node> & geneTree,
				 Node * node,
				 double editionThreshold) ;
void recoverSAndSpPresentInSubtree ( TreeTemplate<Node> & spTree, Node * node ) ;

Node * removeNodesWithDegree1 ( Node * node) ;


#endif //_GENETREEALGORITHMS_H_


