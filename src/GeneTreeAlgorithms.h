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
 * in the maximum likelihood framework, by maximizing theprobability
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

#ifndef _GENETREEALGORITHMS_H_
#define _GENETREEALGORITHMS_H_

// From PhylLib:
#include <Bpp/Phyl/Likelihood/DRTreeLikelihoodTools.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/DRTreeLikelihood.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Protein/JCprot.h>


// From NumCalc:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>


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


//const double DIST = 0.1;


/**************************************************************************
 * This function creates a sequence tree from a species tree and a std::map 
 * containing the link between the species and their sequences.
 **************************************************************************/
bpp::TreeTemplate<bpp::Node> * buildARandomSequenceTreeFromASpeciesTree
(
 std::map <std::string,
 std::deque<std::string> > & spSeqs, 
 bpp::TreeTemplate<bpp::Node> & tree, 
 std::map <std::string, std::string> & spSelectedSeq
);
/**************************************************************************
 * This function creates a sequence tree from a species tree and the std::map 
 * containing the link between the species and the putative orthologous sequence.
 **************************************************************************/
bpp::TreeTemplate<bpp::Node> * buildASequenceTreeFromASpeciesTreeAndCorrespondanceMap
(
 bpp::TreeTemplate<bpp::Node> & tree, 
 std::map <std::string, 
 std::string> & spSelectedSeq
);
/******************************************************************************/
// This function refines branch lengths of a gene tree.
/******************************************************************************/
double refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly
(
 std::map<std::string, std::string> & params, 
 bpp::TreeTemplate<bpp::Node>  *& unrootedGeneTree, 
 bpp::SiteContainer * sites, 
 bpp::SubstitutionModel* model, 
 bpp::DiscreteDistribution* rDist, 
 std::string file, bpp::Alphabet *alphabet, bool mapping=false
);

/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/

void optimizeBLMapping
(
 bpp::DRTreeLikelihood* tl,
 double precision
);

/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/
void optimizeBLMappingForSPRs
(
 bpp::DRTreeLikelihood* tl,
 double precision, std::map<std::string, std::string> params
);

/******************************************************************************/
// This function optimizes branch lengths in a gene tree without substitution mapping
/******************************************************************************/
void optimizeBLForSPRs
(
 bpp::DRTreeLikelihood* tl,
 double precision, std::map<std::string, std::string> params
);


/******************************************************************************/
// This function builds a bionj tree
/******************************************************************************/
bpp::TreeTemplate<bpp::Node>  * buildBioNJTree
(
 std::map<std::string, std::string> & params, 
 bpp::SiteContainer * sites, 
 bpp::SubstitutionModel* model, 
 bpp::DiscreteDistribution* rDist, 
 bpp::Alphabet *alphabet
);
/******************************************************************************/
// This function refines a gene tree topology and branch lengths using the PhyML 
// algorithm.
/******************************************************************************/
void refineGeneTreeUsingSequenceLikelihoodOnly
(
 std::map<std::string, std::string> & params, 
 bpp::TreeTemplate<bpp::Node>  *& unrootedGeneTree, 
 bpp::VectorSiteContainer * sites, 
 bpp::SubstitutionModel* model, 
 bpp::DiscreteDistribution* rDist, 
 std::string file, 
 bpp::Alphabet *alphabet
);


/**************************************************************************
 * This function produces a std::string version of a gene tree, 
 * with gene names replaced by species names. 
 **************************************************************************/
std::string geneTreeToParenthesisWithSpeciesNames
(
 bpp::TreeTemplate<bpp::Node> * geneTree,
 std::map<std::string, std::string > seqSp
);

/**************************************************************************
 * This function produces a gene tree from a std::string version in which 
 * gene names have been changed to include species names. 
 **************************************************************************/
bpp::TreeTemplate<bpp::Node> * parenthesisPlusSpeciesNamesToGeneTree(std::string geneTreeStr);


/**************************************************************************
 * This function produces a std::string version of a gene tree, 
 * with gene names changed to include species names. 
 **************************************************************************/
std::string geneTreeToParenthesisPlusSpeciesNames
(
  bpp::TreeTemplate<bpp::Node> * geneTree,
 std::map<std::string, std::string > seqSp
);


/**************************************************************************
 * This function produces a gene tree with leaves annotated with species names.
 **************************************************************************/
void annotateGeneTreeWithSpeciesNames
(
  bpp::TreeTemplate<bpp::Node> * geneTree,
 std::map<std::string, std::string > seqSp
) ;

/**************************************************************************
 * This function optimizes a gene tree based on the reconciliation score only.
 * It uses SPRs and NNIs, and calls findMLReconciliationDR to compute the likelihood.
 **************************************************************************/
double refineGeneTreeDLOnly
(
  bpp::TreeTemplate<bpp::Node> * spTree, 
 bpp::TreeTemplate<bpp::Node> *& geneTree, 
 std::map<std::string, std::string > seqSp,
 std::map<std::string, int > spID,
 std::vector< double> &lossExpectedNumbers, 
 std::vector < double> &duplicationExpectedNumbers, 
 int & MLindex, 
 std::vector <int> &num0lineages, 
 std::vector <int> &num1lineages, 
 std::vector <int> &num2lineages, 
 std::set <int> &nodesToTryInNNISearch
);

/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void buildVectorOfRegraftingNodesGeneTree
(
  bpp::TreeTemplate<bpp::Node> &spTree, 
 bpp::TreeTemplate<bpp::Node> &tree, 
 int nodeForSPR, 
 int distance, 
 std::vector <int> & nodeIdsToRegraft
) ;

/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void getAllCandidateBranchingPointsFromSpeciesID
(
  bpp::TreeTemplate<bpp::Node> &tree, 
 std::vector <std::string> spIds, 
 std::vector <int> & allNodeIds
) ;

/**************************************************************************
 * This recursive function returns node ids with a given species id.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIds(bpp::Node * node, std::string spId, std::vector <int> & allNodeIds);

/**************************************************************************
 * This recursive function returns node ids with a given species id, upstream from bpp::Node node.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIdsUpstream(bpp::Node * node, std::string spId, std::vector <int> & allNodeIds);


/**************************************************************************
 * This recursive function returns node ids sons of nodes with a given species id.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void getSonsOfNodesWithSimilarSpeciesIds(bpp::Node * node, std::string spId, std::vector <int> & allNodeIds) ;

/**************************************************************************
 * This function returns a vector of branching points for gene tree SPR in the coalescent framework.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void buildVectorOfRegraftingNodesCoalGeneTree
(
  bpp::TreeTemplate<bpp::Node> &spTree, 
 bpp::TreeTemplate<bpp::Node> &tree, 
 int nodeForSPR, 
 int distance, 
 std::vector <int> & nodeIdsToRegraft
) ;

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
unsigned int optimizeBranchLengthsParameters
(
  bpp::DiscreteRatesAcrossSitesTreeLikelihood* tl,
 const bpp::ParameterList& parameters,
 double target,
 bpp::OptimizationListener* listener     = 0,
 double tolerance                   = 0.000001,
 unsigned int tlEvalMax             = 1000000,
 bpp::OutputStream* messageHandler       = bpp::ApplicationTools::message,
 bpp::OutputStream* profiler             = bpp::ApplicationTools::message,
 unsigned int verbose               = 1,
 const std::string& optMethodDeriv  = bpp::OptimizationTools::OPTIMIZATION_NEWTON
)
throw (bpp::Exception);


/**************************************************************************                                                                                                                                                                                                                                                                                               *
 **************************************************************************/

bool editDuplicationNodesMuffato(bpp::TreeTemplate<bpp::Node> & spTree,
                 bpp::TreeTemplate<bpp::Node> & geneTree,
                 double editionThreshold) ;


void editDuplicationNodesMuffato(bpp::TreeTemplate<bpp::Node> & spTree,
                 bpp::TreeTemplate<bpp::Node> & geneTree,
                 bpp::Node * node,
                 double editionThreshold, 
                 bool& edited) ;



void editDuplicationNodesMuffato
(
  bpp::TreeTemplate<bpp::Node> & spTree, 
 bpp::TreeTemplate<bpp::Node> & geneTree,
 bpp::Node * node,
 double editionThreshold) ;
 void recoverSAndSpPresentInSubtree ( bpp::TreeTemplate<bpp::Node> & spTree, bpp::Node * node ) ;
 
 bpp::Node * removeNodesWithDegree1 ( bpp::Node * node
 
 ) ;
 
 
 #endif //_GENETREEALGORITHMS_H_
 
 
 
