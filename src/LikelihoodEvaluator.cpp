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

/*extern "C" {
#include <pll/pllInternal.h>
}*/

#include <iostream>
#include <string>
#include <fstream>
#include <boost/graph/graph_traits.hpp>

#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplateTools.h>


#include "LikelihoodEvaluator.h"
#include "ReconciliationTools.h"




using namespace std;
using namespace bpp;


void LikelihoodEvaluator::PLL_initializePLLInstance(){
  
  /* Set the PLL instance attributes */
  PLL_attributes.rateHetModel     = PLL_GAMMA;
  PLL_attributes.fastScaling      = PLL_TRUE;
  PLL_attributes.saveMemory       = PLL_FALSE;
  PLL_attributes.useRecom         = PLL_FALSE;
  PLL_attributes.randomNumberSeed = 0xDEADBEEF;
  PLL_attributes.numberOfThreads  = 8;            /* This only affects the pthreads version */
  
  PLL_instance = pllCreateInstance (&PLL_attributes);
}


LikelihoodEvaluator::LikelihoodEvaluator(map<string, string> params):
  params(params), alternativeTree(00), initialized(false)
{
  loadDataFromParams();
  tolerance_ = 0.5;
}

void LikelihoodEvaluator::loadDataFromParams(){
  
  // set the name of this evaluator
  istringstream tempName(ApplicationTools::getStringParameter("input.sequence.file",params,"rnd"));
  while(std::getline(tempName,name,'/'))
    ;
  tempName.str(name);
  tempName.clear();
  std::getline(tempName,name,'.');
  if(name.size() == 0)
    name = "unnamed";
  
  cout << "@@Instanciating a Likelihood Evaluator named " << name << endl;
    
  string methodString = ApplicationTools::getStringParameter("likelihood.evaluator",params,"PLL");
  method = (methodString == "PLL"? PLL:BPP);
  
  // loading data, imported from GeneTreeLikelihood (Bastien)
  
  std::vector <std::string> spNames;
  bool cont = false;
  alphabet =  getAlphabetFromOptions(params, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Alphabet"));
 
  sites = getSequencesFromOptions(params, alphabet, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Sequences"));
 
  substitutionModel = getModelFromOptions(params, alphabet, sites, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Model"));
   
 
  rateDistribution = getRateDistributionFromOptions(params, substitutionModel, cont);

  
  //For the moment, we have not filtered sequences, 
  //so we don't want to spend time computing a tree that we may discard soon.
  //Instead we just create a random tree.
  std::vector<std::string> leafNames = sites->getSequencesNames();
  tree = TreeTemplateTools::getRandomTree(leafNames, false);

/*  
  try 
  {
    bool cont = true;
    tree = getTreeFromOptions(params, alphabet, sites, substitutionModel, rateDistribution, cont);
  }
  catch (std::exception& e)
  {
    std::cout << e.what() <<"; Unable to get a proper gene tree for family <<file<< avoiding this family."<<std::endl;
    cont=false;
  }*/
}


// TODO: is clonable?
// LikelihoodEvaluator::LikelihoodEvaluator(LikelihoodEvaluator &levaluator){
//   nniLk = levaluator.nniLk->clone();
//   tree = levaluator.tree->clone();
// }

void LikelihoodEvaluator::PLL_loadAlignment(string path)
{
  /* Parse a PHYLIP/FASTA file */
  PLL_alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, path.c_str());
  if (!PLL_alignmentData)
  {
    throw Exception("PLL: Error while parsing " + path);
  }
}


void LikelihoodEvaluator::PLL_loadNewick_fromFile(string path)
{
  PLL_newick = pllNewickParseFile(path.c_str());
  if (!PLL_newick)
  {
    throw Exception("PLL: Error while parsing newick file");
  }
  if (!pllValidateNewick (PLL_newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    throw Exception("PLL: Invalid phylogenetic tree.");
  }
}

void LikelihoodEvaluator::PLL_loadNewick_fromString(string newick)
{
  PLL_newick = pllNewickParseString (newick.c_str());
  if (!PLL_newick)
  {
    throw Exception("PLL: Error while parsing newick string: " + newick);
  }
  if (!pllValidateNewick (PLL_newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    throw Exception("PLL: Invalid phylogenetic tree.");
  }
}


void LikelihoodEvaluator::PLL_loadPartitions(string path)
{
  /* Parse the partitions file into a partition queue structure */
  PLL_partitionInfo = pllPartitionParse (path.c_str());
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (PLL_partitionInfo, PLL_alignmentData))
  {
    throw Exception("PLL: Partitions do not cover all sites.");
  }
  
  /* Commit the partitions and build a partitions structure */
  PLL_partitions = pllPartitionsCommit (PLL_partitionInfo, PLL_alignmentData);
  
  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&PLL_partitionInfo);
  
  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (PLL_alignmentData, PLL_partitions);
}

void LikelihoodEvaluator::PLL_connectTreeAndAlignment()
{
  pllTreeInitTopologyNewick (PLL_instance, PLL_newick, PLL_FALSE);
    
  // cout << "PLL: Connect the alignment and partition structure with the tree structure" << std::endl ;
  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (PLL_instance, PLL_alignmentData, PLL_partitions, PLL_DEEP_COPY))
  {
    throw Exception("PLL: Incompatible tree/alignment combination.");
  }
}


void LikelihoodEvaluator::initialize_BPP_nniLk()
{
  nniLk = new NNIHomogeneousTreeLikelihood(*tree, *sites, substitutionModel, rateDistribution, mustUnrootTrees, verbose);
  
  nniLk->initParameters();
  nniLk->initialize();
  logLikelihood = nniLk->getValue();
}

void LikelihoodEvaluator::initialize_PLL()
{

  // #1 PREPARING
  // must have the strict names loaded
  loadStrictNamesFromAlignment_forPLL();
  writeAlignmentFilesForPLL();
  
  // preparing the tree
  
  // PLL process
  alpha_ = 1.0;
/*  baseFreq_[0]=0.25;
  baseFreq_[1]=0.25;
  baseFreq_[2]=0.25;
  baseFreq_[3]=0.25;
  subsMatrix_[0] = 1/6;
  subsMatrix_[1] = 1/6;  
  subsMatrix_[2] = 1/6;
  subsMatrix_[3] = 1/6;
  subsMatrix_[4] = 1/6;
  subsMatrix_[5] = 1/6;*/

  PLL_initializePLLInstance();
  PLL_loadAlignment(fileNamePrefix + "alignment.fasta");
  PLL_loadPartitions(fileNamePrefix + "partition.txt");
  
  logLikelihood = PLL_optimizeBranchLengthsAndParameters(&tree);
  
  //logLikelihood = PLL_evaluate(&tree);
  
}

void LikelihoodEvaluator::setTree(TreeTemplate<Node> * newTree)
{
  if(!isInitialized())
    tree = newTree->clone();
  else
    throw Exception("This evaluator has already be initialized. It is forbidden to modify it now.");
}



double LikelihoodEvaluator::PLL_evaluate(TreeTemplate<Node>** treeToEvaluate)
{
  
  //TODO debug remove
  Newick debugTree;
  stringstream debugSS;
  debugTree.write(**treeToEvaluate,debugSS);
//  cout << "tree to evaluate:\n" << debugSS.str() << endl;
  
  // preparing the tree
  TreeTemplate<Node>* treeForPLL = (*treeToEvaluate)->clone();
  convertTreeToStrict(treeForPLL);
  
  // getting the root
  bool wasRooted = (treeForPLL->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted)
  {
    Node* root = treeForPLL->getRootNode();
    Node* son1 = root->getSon(0);
    Node* son2 = root->getSon(1);
    vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
    vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
    leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
    leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
    treeForPLL->unroot();
  }
  
  
  Newick newickForPll;
  stringstream newickStingForPll;
  newickForPll.write(*treeForPLL,newickStingForPll);
  PLL_loadNewick_fromString(newickStingForPll.str());
  delete treeForPLL;
  
  // processing by PLL
  PLL_connectTreeAndAlignment();
//  pllSetFixedAlpha(alpha_, 0, PLL_partitions, PLL_instance);
//  pllSetFixedBaseFrequencies(baseFreq_, 4, 0, PLL_partitions, PLL_instance);
  //pllSetFixedSubstitutionMatrix(subsMatrix_, 6, 0, PLL_partitions, PLL_instance);
  pllInitModel(PLL_instance, PLL_partitions, PLL_alignmentData);
 // pllOptimizeBranchLengths (PLL_instance, PLL_partitions, 64);
 // pllOptimizeModelParameters(PLL_instance, PLL_partitions, 0.1);

/*std::cout << "HEHEHEHE " << std::endl;
std::cout << "tolerance_ "<< tolerance_ << std::endl;
std::cout << "HEHEHEHE 2" << std::endl;*/
  pllOptimizeModelParameters(PLL_instance, PLL_partitions, tolerance_);

  // getting the new tree with new branch lengths
  pllTreeToNewick(PLL_instance->tree_string, PLL_instance, PLL_partitions, PLL_instance->start->back, true, true, 0, 0, 0, true, 0,0);
  newickStingForPll.str(PLL_instance->tree_string);
  
  //debug
 // cout << "returned tree from PLL \n" << newickStingForPll.str() << endl;
  
  delete *treeToEvaluate;
  *treeToEvaluate = newickForPll.read(newickStingForPll);
  
  
  //re-rooting if needed
  if(wasRooted)
  {
    // the plyogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (*treeToEvaluate)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
      
      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currSon);

    }
    
    // if not, we will try all the internal branches as potential roots    
    if(!(*treeToEvaluate)->isRooted())
    {
      vector<Node*> outgroupCandidates = (*treeToEvaluate)->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
        
        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currCandidate);
      } 
    }
    
    
    if(!(*treeToEvaluate)->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    } 
    
  }
  
  restoreTreeFromStrict(*treeToEvaluate);
  
  //TODO debug remove
 /* stringstream debugSS2;
  debugTree.write(**treeToEvaluate,debugSS2);
  cout << "Final tree for BPP" << debugSS2.str() << endl;
  */
  return(PLL_instance->likelihood);
}


  double LikelihoodEvaluator::PLL_optimizeBranchLengthsAndParameters(bpp::TreeTemplate<bpp::Node>** treeToEvaluate) {
  
  //TODO debug remove
  Newick debugTree;
  stringstream debugSS;
  debugTree.write(**treeToEvaluate,debugSS);
//  cout << "tree to evaluate:\n" << debugSS.str() << endl;
  
  // preparing the tree
  TreeTemplate<Node>* treeForPLL = (*treeToEvaluate)->clone();
  convertTreeToStrict(treeForPLL);
  
  // getting the root
  bool wasRooted = (treeForPLL->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted)
  {
    Node* root = treeForPLL->getRootNode();
    Node* son1 = root->getSon(0);
    Node* son2 = root->getSon(1);
    vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
    vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
    leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
    leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
    treeForPLL->unroot();
  }
  
  
  Newick newickForPll;
  stringstream newickStingForPll;
  newickForPll.write(*treeForPLL,newickStingForPll);
  PLL_loadNewick_fromString(newickStingForPll.str());
  delete treeForPLL;
  
  // processing by PLL
  PLL_connectTreeAndAlignment();
  pllInitModel(PLL_instance, PLL_partitions, PLL_alignmentData);
//  pllOptimizeBranchLengths (PLL_instance, PLL_partitions, 64);
/*std::cout << "IIIIIHEHEHEHE " << std::endl;
std::cout << "tolerance_ "<< tolerance_ << std::endl;
std::cout << "IIIIIHEHEHEHE 2" << std::endl;*/

  pllOptimizeModelParameters(PLL_instance, PLL_partitions, tolerance_);

  // getting the new tree with new branch lengths
  pllTreeToNewick(PLL_instance->tree_string, PLL_instance, PLL_partitions, PLL_instance->start->back, true, true, 0, 0, 0, true, 0,0);
  newickStingForPll.str(PLL_instance->tree_string);

  //Set the backup of the parameters of the model
  alpha_ = pllGetAlpha ( PLL_partitions, 0);

 // std::cout << "Alpha: "<< alpha_ <<std::endl;
  if (alpha_ < 0.02 ) {
    alpha_ = 0.02;
  } 
  else if ( alpha_ > 1000.0 ) {
    alpha_ = 1000.0;
  }


 /* pllGetBaseFrequencies(PLL_instance,  PLL_partitions, 0, baseFreq_);
  pllGetSubstitutionMatrix(PLL_instance,  PLL_partitions, 0, subsMatrix_);
*/

  //debug
//  cout << "returned tree from PLL \n" << newickStingForPll.str() << endl;
  
  delete *treeToEvaluate;
  *treeToEvaluate = newickForPll.read(newickStingForPll);
  
  
  //re-rooting if needed
  if(wasRooted)
  {
    // the plyogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (*treeToEvaluate)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
      
      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currSon);

    }
    
    // if not, we will try all the internal branches as potential roots    
    if(!(*treeToEvaluate)->isRooted())
    {
      vector<Node*> outgroupCandidates = (*treeToEvaluate)->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
        
        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currCandidate);
      } 
    }
    
    
    if(!(*treeToEvaluate)->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    } 
    
  }
  
  restoreTreeFromStrict(*treeToEvaluate);
  
  //TODO debug remove
 /* stringstream debugSS2;
  debugTree.write(**treeToEvaluate,debugSS2);
  cout << "Final tree for BPP" << debugSS2.str() << endl;*/
  
  return(PLL_instance->likelihood);

}

double LikelihoodEvaluator::BPP_evaluate(TreeTemplate<Node>** treeToEvaluate)
{ 

  // preparing the tree
  TreeTemplate<Node>* treeForBPP = (*treeToEvaluate)->clone();

  // getting the root
  bool wasRooted = (treeForBPP->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted)
  {
    Node* root = treeForBPP->getRootNode();
    Node* son1 = root->getSon(0);
    Node* son2 = root->getSon(1);
    vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
    vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
    leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
    leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
    treeForBPP->unroot();
  }


  NNIHomogeneousTreeLikelihood * drlk = 0;
  drlk  = new NNIHomogeneousTreeLikelihood (**treeToEvaluate, 
                                              *(nniLk->getData()), 
                                              nniLk->getSubstitutionModel(), 
                                              nniLk->getRateDistribution(), 
                                              true, false);

  drlk->initialize();
  auto_ptr<BackupListener> backupListener;
  int tlEvalMax = 100;
  OutputStream* messageHandler = 0 ; 
  OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (drlk), drlk->getParameters(), backupListener.get(), tolerance_, tlEvalMax, messageHandler, messageHandler, 0);

  delete *treeToEvaluate;
  *treeToEvaluate = static_cast< TreeTemplate<Node>* > (drlk->getTree().clone() );

 //re-rooting if needed
  if(wasRooted)
  {
    // the phylogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (*treeToEvaluate)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
      
      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currSon);

    }
    
    // if not, we will try all the internal branches as potential roots    
    if(!(*treeToEvaluate)->isRooted())
    {
      vector<Node*> outgroupCandidates = (*treeToEvaluate)->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
        
        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currCandidate);
      } 
    }
    
    
    if(!(*treeToEvaluate)->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    } 
    
  }

  return drlk->getValue();
 

}





// CLONABLE?
// LikelihoodEvaluator* LikelihoodEvaluator::clone()
// {
//   return(new LikelihoodEvaluator(this));
// }


LikelihoodEvaluator::~LikelihoodEvaluator()
{
//   if(method == PLL){
//     delete tree;
//   }
//   else
//   {
//     delete nniLk;
//     if(nniLkAlternative)
//       delete nniLkAlternative;
//     //TODO delete all the trees, etc
//   }
}

void LikelihoodEvaluator::initialize()
{
  
  //checking the alignment and the tree contain the same number of sequences
  if(sites->getNumberOfSequences() != tree->getNumberOfLeaves()){
    ostringstream errorMessage;
    errorMessage << "\nNumber of sequences (here: "<< sites->getNumberOfSequences() << "sequences) must match to number of leaves in the tree (here: "<< tree->getNumberOfLeaves() << "leaves). I give up.";
    throw Exception(errorMessage.str());
  }
  
  // ### common requirements for initialization
  if(method == PLL)
    initialize_PLL();
  else
    initialize_BPP_nniLk();
  
  //
  initialized = true;
}


void LikelihoodEvaluator::setAlternativeTree(TreeTemplate< Node >* newAlternative)
{
  if(!initialized)
    throw Exception("Set alternative tree on a non initalized evaluator.");
  if(alternativeTree != 00)
    delete alternativeTree;
  alternativeTree = newAlternative->clone();
  
  if(method == PLL){
    alternativeLogLikelihood = PLL_evaluate( &alternativeTree );
  }
  else
  {
  /*  if ( nniLkAlternative)
      delete nniLkAlternative;*/
    alternativeLogLikelihood = BPP_evaluate( &alternativeTree );

   /* nniLkAlternative =  new NNIHomogeneousTreeLikelihood (*alternativeTree, *sites, substitutionModel, rateDistribution, mustUnrootTrees, verbose);
    alternativeLogLikelihood = nniLkAlternative->getLogLikelihood();*/
  }
}

void LikelihoodEvaluator::acceptAlternativeTree()
{
  delete tree;
  tree = alternativeTree->clone();
  logLikelihood = alternativeLogLikelihood;
}





TreeTemplate< Node >* LikelihoodEvaluator::getAlternativeTree()
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
  return alternativeLogLikelihood;
}

double LikelihoodEvaluator::getLogLikelihood()
{
  return logLikelihood;
}

DiscreteDistribution* LikelihoodEvaluator::getRateDistribution()
{
  return rateDistribution;
}

VectorSiteContainer* LikelihoodEvaluator::getSites()
{
  return sites;
}

TreeTemplate< Node >* LikelihoodEvaluator::getTree()
{
  return tree;
}

SubstitutionModel* LikelihoodEvaluator::getSubstitutionModel()
{
  return substitutionModel;
}


void LikelihoodEvaluator::loadStrictNamesFromAlignment_forPLL()
{
  vector<string> seqNames = sites->getSequencesNames();
  string currName;
  ostringstream currStrictName;
  for(unsigned int currIdx = 0; currIdx != seqNames.size(); currIdx++)
  {
    currName = seqNames.at(currIdx);
    if(currName.at(currName.size()-1) == '\r')
      currName = currName.substr(0,(currName.size()-1));
    
//     cout << "Curr sequence name = ..." << currName << "..." << endl; 
    
    
//     // determining current postfix with letters instead of integers.
//     // Thanks PLL to be so strict!...
//     // (it is not a clean base conversion)
//     ostringstream currPrefix;
//     unsigned int currIntegerPrefix = currIdx+1;
//     while(currIntegerPrefix != 0)
//     {
//       unsigned int unit = currIntegerPrefix % 10;
//       currPrefix << (char)('A' + unit);
//       currIntegerPrefix = currIntegerPrefix / 10;
//     }
//     cout << "curr prefix = " << currPrefix.str() << endl;
    
    currStrictName.clear();
    currStrictName.str("");
    currStrictName << "seq" << currIdx;
    realToStrict[currName] = currStrictName.str();
    strictToReal[currStrictName.str()] = currName;
    
  }
  
}

void LikelihoodEvaluator::convertTreeToStrict(TreeTemplate< Node >* targetTree)
{
  
  vector<Node*> leaves = targetTree->getLeaves();
  for(vector<Node*>::iterator currLeaf = leaves.begin(); currLeaf != leaves.end(); currLeaf++)
  {
    map<string,string>::iterator found = realToStrict.find((*currLeaf)->getName());
    if(found == realToStrict.end())
    {
      cout << "Unable to find sequence named ++" << (*currLeaf)->getName() << "++ in the alignment." << endl;
    }
    else
    {
      (*currLeaf)->setName(found->second);
    }
  }
}

void LikelihoodEvaluator::restoreTreeFromStrict(TreeTemplate< Node >* targetTree)
{
  vector<Node*> leaves = targetTree->getLeaves();
  for(vector<Node*>::iterator currLeaf = leaves.begin(); currLeaf != leaves.end(); currLeaf++){
    (*currLeaf)->setName(strictToReal[(*currLeaf)->getName()]);
  }
}


void LikelihoodEvaluator::writeAlignmentFilesForPLL()
{
  fileNamePrefix = "tmpPLL_" + name + "_" ;
  ofstream alignementFile(string(fileNamePrefix + "alignment.fasta").c_str(), ofstream::out);

  //preparing the file for the alignment
  BasicSequence currSequence(sites->getAlphabet());
  //DEBUG
  cout << "Writing an alignment for PLL with " << sites->getNumberOfSequences() << " sequences" << endl;
  
  for(unsigned int currSeqIndex = 0; currSeqIndex != sites->getNumberOfSequences(); currSeqIndex++)
  {
    currSequence = sites->getSequence(currSeqIndex);
    string currSequenceName = currSequence.getName();
    if(currSequenceName.at(currSequenceName.size()-1) == '\r')
      currSequenceName = currSequenceName.substr(0,(currSequenceName.size()-1));
    alignementFile << ">" << realToStrict[currSequenceName] << "\n" << currSequence.toString() << "\n";
  }
  alignementFile.close();
  ofstream partitionFile(string(fileNamePrefix + "partition.txt").c_str(), ofstream::out);
if (alphabet->getSize() == 4) {
    if (substitutionModel->getName()!="GTR") {
     std::cout << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes GTR for nucleotide models." <<std::endl;
     cout.flush();
     std::cerr << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes GTR for nucleotide models." <<std::endl;
     cerr.flush();   
     MPI::COMM_WORLD.Abort(1);
     exit(-1);
    }
  partitionFile << "DNA, p1=1-" << sites->getNumberOfSites() << "\n";
}
else if (alphabet->getSize() == 20) {
    if (substitutionModel->getName()=="LG08") {
        partitionFile << "LG, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else if (substitutionModel->getName()=="WAG01") {
        partitionFile << "WAG, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else if (substitutionModel->getName()=="JTT92") {
        partitionFile << "JTT, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else {
     std::cout << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes LG08, W1G01, JTT92 for protein models." <<std::endl;
     cout.flush();
     std::cerr << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes LG08, W1G01, JTT92 for protein models." <<std::endl;
     cerr.flush();   
     MPI::COMM_WORLD.Abort(1);
     exit(-1);
    }
}
else {
  std::cout << "Error: alphabet incompatible with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only works with DNA/RNA or Protein alphabets." <<std::endl;
  cout.flush();
  std::cerr << "Error: alphabet incompatible with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only works with DNA/RNA or Protein alphabets." <<std::endl;
  cerr.flush();
  MPI::COMM_WORLD.Abort(1);
  exit(-1);
}
  partitionFile.close();
}

LikelihoodEvaluator::LikelihoodEvaluator(LikelihoodEvaluator const &leval):
params(leval.params), initialized(false), PLL_instance(00), PLL_alignmentData(00), PLL_newick(00), PLL_partitions(00), PLL_partitionInfo(00), tree(00), alternativeTree(00), nniLk(00), nniLkAlternative(00), substitutionModel(00), rateDistribution(00), sites(00), alphabet(00)
{
  
  loadDataFromParams();
  tree = leval.tree;
  sites = leval.sites;
  
  if(leval.initialized)
    initialize();
}

LikelihoodEvaluator* LikelihoodEvaluator::clone()
{
  return new LikelihoodEvaluator(*this);
}

LikelihoodEvaluator::LikelihoodEvaluator(const Tree* tree, const SiteContainer* alignment, SubstitutionModel* model, DiscreteDistribution* rateDistribution, std::map<std::string, std::string> par, bool mustUnrootTrees, bool verbose):
initialized(false), PLL_instance(00), PLL_alignmentData(00), PLL_newick(00), PLL_partitions(00), PLL_partitionInfo(00), tree(00), alternativeTree(00), nniLk(00), nniLkAlternative(00), substitutionModel(00), rateDistribution(00), sites(00), alphabet(00), params(par)
{
  this->tree = dynamic_cast<TreeTemplate<Node> *>(tree->clone());
  this->substitutionModel = model->clone();
  this->rateDistribution = rateDistribution->clone();
  this->sites = dynamic_cast<VectorSiteContainer*>(alignment->clone());
  
  string methodString = ApplicationTools::getStringParameter("likelihood.evaluator",params,"PLL");
  this->method = (methodString == "PLL"? PLL:BPP);
  
  initialize();
  
}
