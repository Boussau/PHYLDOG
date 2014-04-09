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



extern "C" {
#include <pll/pll.h>
}

#include<iostream>

#include <Bpp/Seq/Container/SiteContainerTools.h>


#include "LikelihoodEvaluator.h"
#include "ReconciliationTools.h"




using namespace std;
using namespace bpp;


LikelihoodEvaluator::initializePLLtree(){
  
  // PLL
  /* Set the PLL instance attributes */
  pllInstanceAttr attr;
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_TRUE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = 8;            /* This only affects the pthreads version */
  
  tr_PLL = pllCreateInstance (&attr);
}


LikelihoodEvaluator::LikelihoodEvaluator(map<string, string> params):
  params(params)
{
 
  // loading data, imported from GeneTreeLikelihood (Bastien)
  
  std::vector <std::string> spNames;
  bool cont = false;
  alphabet =  getAlphabetFromOptions(params, cont);

  if (!cont)
    throw(Exception("Unable to load this family"));
 
  sites = getSequencesFromOptions(params, alphabet, cont);

  if (!cont)
    throw(Exception("Unable to load this family"));
 
  substitutionModel = getModelFromOptions(params, alphabet, sites, cont);

  if (!cont)
    throw(Exception("Unable to load this family"));
    
  if (substitutionModel->getName() != "RE08")
    SiteContainerTools::changeGapsToUnknownCharacters(*sites, cont);
 
  if (!cont)
    throw(Exception("Unable to load this family"));
 
  rateDistribution = getRateDistributionFromOptions(params, substitutionModel, cont);
  
  try 
  {
    tree = getTreeFromOptions(params, alphabet, sites, substitutionModel, rateDistribution);
  }
  catch (std::exception& e)
  {
    std::cout << e.what() <<"; Unable to get a proper gene tree for family "<<file<<"; avoiding this family."<<std::endl;
    cont=false;
  }
    
}

LikelihoodEvaluator::LikelihoodEvaluator(LikelihoodEvaluator &levaluator){
  nniLk = levaluator.nniLk->clone();
  tree = levaluator.tree->clone();
}

LikelihoodEvaluator::loadPLLalignment(char* path)
{
  /* Parse a PHYLIP/FASTA file */
  pllAlignmentData * alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, path);
  if (!alignmentData)
  {
    throw Exception("PLL: Error while parsing " + path);
  }
}


LikelihoodEvaluator::loadPLLnewick(char* path)
{
  newick_PLL = pllNewickParseFile(path);
  if (!newick_PLL)
  {
    throw Exception("PLL: Error while parsing newick file");
  }
  if (!pllValidateNewick (newick_PLL))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    throw Exception("PLL: Invalid phylogenetic tree.");
  }
}

LikelihoodEvaluator::loadPLLpartitions(char* path)
{
  /* Parse the partitions file into a partition queue structure */
  partitionInfo_PLL = pllPartitionParse (path);
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (partitionInfo_PLL, alignmentData_PLL))
  {
    throw Exception("Error: Partitions do not cover all sites.");
  }
  
  /* Commit the partitions and build a partitions structure */
  partitions_PLL = pllPartitionsCommit (partitionInfo_PLL, alignmentData_PLL);
  
  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo_PLL);
  
  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (alignmentData_PLL, partitions_PLL);
}

LikelihoodEvaluator::updatePLLtreeWithPLLnewick()
{
  pllTreeInitTopologyNewick (tr_PLL, newick_PLL, PLL_FALSE);
    
  // cout << "PLL: Connect the alignment and partition structure with the tree structure" << std::endl ;
  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (tr_PLL, alignmentData_PLL, partitions_PLL, PLL_DEEP_COPY))
  {
    throw Exception("PLL: Incompatible tree/alignment combination.");
  }
}


LikelihoodEvaluator::initialize_BPP_nniLk()
{
  nniLk = new NNIHomogeneousTreeLikelihood(tree, sites, substitutionModel, rateDistribution, mustUnrootTrees, verbose);
//   with *sites*, found in GeneTreeLikelihood
//    new NNIHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, levaluator_->getSubstitutionModel(), levaluator_->getRateDistribution(), true, true); 
  nniLk->initParameters();
}


LikelihoodEvaluator* LikelihoodEvaluator::clone()
{
  return(new LikelihoodEvaluator(this));
}


bpp::ParameterList LikelihoodEvaluator::getParameters()
{
  return(this->nniLk->getParameters());
}


NNIHomogeneousTreeLikelihood * LikelihoodEvaluator::getNniLk()
{
  return(this->nniLk);
}


LikelihoodEvaluator::~LikelihoodEvaluator()
{
  delete nniLk;
}

NNIHomogeneousTreeLikelihood* LikelihoodEvaluator::initialize()
{
  //for the moment, BPP only
  initialize_BPP_nniLk();
}


void LikelihoodEvaluator::setAlternativeTree(TreeTemplate* newAlternative)
{
  if(alternativeTree != 00)
    delete alternativeTree;
  alternativeTree = newAlternative->clone();
}

void LikelihoodEvaluator::acceptAlternativeTree()
{
  delete tree;
  tree = alternativeTree;
  delete nniLk;
  nniLk = nniLkAlternative;
}





TreeTemplate* LikelihoodEvaluator::getAlternativeTree()
{
  return alternativeTree;
}

bool LikelihoodEvaluator::isInitialized()
{
  return initialized;
}

Alphabet* LikelihoodEvaluator::getAlphabet()
{
  return alphabet;
}

double LikelihoodEvaluator::getAlternativeLogLikelihood()
{
  return lastAlternativeLikelihood;
}

double LikelihoodEvaluator::getLogLikelihood()
{
  return lastLikelihood;
}

DiscreteDistribution* LikelihoodEvaluator::getRateDistribution()
{
  return rateDistribution;
}

VectorSiteContainer* LikelihoodEvaluator::getSites()
{
  return sites;
}

bpp::TreeTemplate<N>* LikelihoodEvaluator::getTree()
{
  return tree;
}

