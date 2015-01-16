/*
 * Copyright or © or Copr. Centre National de la Recherche Scientifique
 * contributor : Bastien Boussau (2009-2013)
 * 
 * bastien.boussau@univ-lyon1.fr
 * 
 * This software is a computer program whose purpose is to simultaneously build 
 * gene and species trees when gene families have undergone duplications and 
 * losses. It can analyze thousands of gene families in dozens of genomes 
 * simultaneously, and was presented in an article in Genome Research. Trees and 
 * parameters are estimated in the maximum likelihood framework, by maximizing 
 * theprobability of alignments given the species tree, the gene trees and the 
 * parameters of duplication and loss.
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

#include <iostream>

#include "Constants.h"
#include "GeneTreeLikelihood.h"

using namespace bpp;



class geneTreeLikelihoodException: public exception
{
  virtual const char* what() const throw()
  {
  WHEREAMI( __FILE__ , __LINE__ );
    return "Can't create a GeneTreeLikelihood object, avoiding family.";
  }
} geneTreeLikelihoodEx;




GeneTreeLikelihood::GeneTreeLikelihood():
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(), spId_() {
  WHEREAMI( __FILE__ , __LINE__ );
  totalIterations_ = 0;
  counter_ = 0;
}


GeneTreeLikelihood::GeneTreeLikelihood(
  std::string file,
  map<string, string> params,
  TreeTemplate<Node> & spTree )
throw (exception):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(), spId_(), 
params_(params), considerSequenceLikelihood_(true)
{
  WHEREAMI( __FILE__ , __LINE__ );
  totalIterations_ = 0;
  counter_ = 0;
  spTree_ = spTree.clone();
  spId_ = computeSpeciesNamesToIdsMap(*spTree_);
  
  bool avoidFamily = false;
  
  // instanciating a new likelihood evaluator with params :
  
  // this will load the tree & the alignemnts
  
 // cout << "new LE 86 - params size = " << params.size() << endl;
  
  levaluator_ = new LikelihoodEvaluator(params_);

  bool qualityFilters = ApplicationTools::getBooleanParameter("use.quality.filters",params_, true);
  
  //TODO: dirty cont to eliminate
  bool cont = true;
  //method to optimize the gene tree root.
  bool rootOptimization = false;
  if (!cont)
      throw(Exception("Unable to load this family"));
  
  std::map<std::string, std::deque<std::string> > spSeq;
  getCorrespondanceSequenceSpeciesFromOptions(params_, cont, seqSp_, spSeq );
  
  if (!cont)
    throw(Exception("Unable to load this family"));
  removeUselessSequencesFromAlignment( spTree_, levaluator_->getSites(), cont , spSeq, file) ;
  
  if (cont) {
    /****************************************************************************
     * Then we need to get the file containing the gene tree, 
     * or build the gene tree.
     *****************************************************************************/
    rootedTree_ = getTreeFromOptions(params_, levaluator_->getAlphabet(), levaluator_->getSites(), levaluator_->getSubstitutionModel(), levaluator_->getRateDistribution(), cont);
  }
  
  if (cont && qualityFilters) 
  { //This family is phylogenetically informative
    qualityControlGeneTree ( rootedTree_, levaluator_->getSites(), cont , file) ;
      }
 
  if (cont) {
    // set the levaluator tree to the modified one
    levaluator_->setTree(rootedTree_);    
  }
  if (cont || !qualityFilters) 
  { //This family is phylogenetically informative
    
    /************************************************************************************************************/
    /********************************************COMPUTING LIKELIHOOD********************************************/
    /************************************************************************************************************/
    bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", params_, true, "", false, false);
    if(!computeLikelihood)
    {
      if (levaluator_)
        delete levaluator_;
      std::cout << "PHYLDOG's done. Bye." << std::endl;
      MPI::COMM_WORLD.Abort(1);
      exit(-1);
    }

    levaluator_->initialize();
    
    
    TreeTemplate<Node> * unrootedGeneTree = rootedTree_->clone();
    if (!rootedTree_->isRooted()) {
      std::cout <<"gene tree is not rooted!!! "<< file <<std::endl;
    }
    unrootedGeneTree->unroot();
    
    //printing the gene trees with the species names instead of the sequence names
    //This is useful to build an input for duptree for instance
    geneTreeWithSpNames_ = unrootedGeneTree->clone();
    std::vector <Node*> leaves = geneTreeWithSpNames_->getLeaves();
    for (unsigned int j =0; j<leaves.size() ; j++) 
    {
      leaves[j]->setName(seqSp_[leaves[j]->getName()]);
    }
    //Now we have a gene tree with species names and bootstrap values and branch lengths
    std::vector <Node *> nodes =  geneTreeWithSpNames_->getNodes();
    for (unsigned int j =0; j<nodes.size() ; j++) 
    {
      if (nodes[j]->hasFather()) 
      {
        nodes[j]->deleteDistanceToFather(); 
      }
      if (nodes[j]->hasBootstrapValue()) 
      {
        nodes[j]->deleteBranchProperty(TreeTools::BOOTSTRAP); 
      }
    }
    
    //Outputting the starting tree, with species names, and with sequence names
    Newick newick(true);
    std::string startingGeneTreeFile = ApplicationTools::getStringParameter("output.starting.gene.tree.file",params_,"startingGeneTree.tree");
    try
    {
      Nhx *nhx = new Nhx();
      annotateGeneTreeWithDuplicationEvents (*spTree_, 
                                             *rootedTree_, 
                                             rootedTree_->getRootNode(), 
                                             seqSp_, spId_); 
      nhx->write(*rootedTree_, startingGeneTreeFile, true);
      
      // newick.write(*geneTree_, startingGeneTreeFile, true);
    }
    catch (IOException e)
    {
      cout << "Problem writing tree to file "<< startingGeneTreeFile <<"\n Is the file path correct and do you have the proper authorizations? "  << endl;
    }
    try
    {
      newick.write(*geneTreeWithSpNames_, startingGeneTreeFile+"_sp_names", true);
    }
    catch (IOException e)
    {
      cout << "Problem writing tree to file "<< startingGeneTreeFile+"_sp_names" <<"\n Is the file path correct and do you have the proper authorizations? "  << endl;
    }
    //std::cout << " Rooted tree? : "<<TreeTemplateTools::treeToParenthesis(*geneTree_, true) << std::endl;
    
    GeneTreeLikelihood* tl;
    std::string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", params_, "no", "", true, false);
    
    sprLimitGeneTree_ = ApplicationTools::getIntParameter("SPR.limit.gene.tree", params_, 2, "", false, false);  
    //ApplicationTools::displayResult("Clock", optimizeClock);
    
    if(optimizeClock == "global")
    {
      std::cout<<"Sorry, clocklike trees have not been implemented yet."<<std::endl;
      MPI::COMM_WORLD.Abort(1);
      exit(0);
    }// This has not been implemented!
    else if(optimizeClock == "no")
    {
     // levaluator_->initialize();
      
      
      
    }
    else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
    
  }
  else 
  {
    throw geneTreeLikelihoodEx;
  }
  
  
  
}

// /**
//  * @brief Build a new ReconciliationTreeLikelihood object.
//  *
//  * @param tree The tree to use.
//  * @param model The substitution model to use.
//  * @param rDist The rate across sites distribution to use.
//  * @param spTree The species tree
//  * @param rootedTree rooted version of the gene tree
//  * @param seqSp link between sequence and species names
//  * @param spId link between species name and species ID
//  * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
//  * @param MLindex ML rooting position
//  * @param checkRooted Tell if we have to check for the tree to be unrooted.
//  * If true, any rooted tree will be unrooted before likelihood computation.
//  * @param verbose Should I display some info?
//  * @throw Exception if an error occured.
//  */
// GeneTreeLikelihood::GeneTreeLikelihood(
//   const Tree & tree,
//   SubstitutionModel * model,
//   DiscreteDistribution * rDist,
//   TreeTemplate<Node> & spTree,  
//   TreeTemplate<Node> & rootedTree, 
//   TreeTemplate<Node> & geneTreeWithSpNames,
//   const std::map <std::string, std::string> seqSp,
//   std::map <std::string,int> spId,
//   int speciesIdLimitForRootPosition,
//   int & MLindex, 
//   bool checkRooted,
//   bool verbose,
//   bool rootOptimization, 
//   bool considerSequenceLikelihood, 
//   unsigned int sprLimit)
// throw (Exception):
// levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(seqSp), spId_(spId)
// {
//   cout << "new LE 266 - params size = " << params_.size() << endl;
// 
//   levaluator_ = new LikelihoodEvaluator(tree,); 
//   spTree_ = spTree.clone();
//   rootedTree_ = rootedTree.clone();
//   geneTreeWithSpNames_ = geneTreeWithSpNames.clone();
//   scenarioLikelihood_ = UNLIKELY;
//   // _sequenceLikelihood = UNLIKELY;
//   MLindex_ = MLindex;
//   rootOptimization_ = rootOptimization; 
//   tentativeMLindex_ = MLindex;
//   totalIterations_ = 0;
//   counter_ = 0;
//   _speciesIdLimitForRootPosition_ = speciesIdLimitForRootPosition;
//   optimizeSequenceLikelihood_ = true;
//   optimizeReconciliationLikelihood_ = true;
//   considerSequenceLikelihood_ = considerSequenceLikelihood;
//   sprLimit_ = sprLimit;
//   // listOfPreviousRoots_ = new std::vector <int> ();
// }


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
  int & MLindex, 
	std::map <std::string, std::string > params,
  bool checkRooted,
  bool verbose, 
  bool rootOptimization, 
  bool considerSequenceLikelihood, 
  unsigned int sprLimitGeneTree)
throw (Exception):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_ (seqSp), spId_(spId), params_(params)
{
  WHEREAMI( __FILE__ , __LINE__ );
  levaluator_ = new LikelihoodEvaluator(&tree, &data, model, rDist, params, false, verbose);
  spTree_ = spTree.clone();
  rootedTree_ = rootedTree.clone();
  geneTreeWithSpNames_ = geneTreeWithSpNames.clone();
  scenarioLikelihood_ = UNLIKELY;
  MLindex_ = MLindex;
  rootOptimization_ = rootOptimization; 
  tentativeMLindex_ = MLindex;
  totalIterations_ = 0; 
  counter_ = 0;
  _speciesIdLimitForRootPosition_ = speciesIdLimitForRootPosition;
  optimizeSequenceLikelihood_ = true;
  optimizeReconciliationLikelihood_ = true;
  considerSequenceLikelihood_ = considerSequenceLikelihood;
  sprLimitGeneTree_ = sprLimitGeneTree;
}


/**
 * @brief Copy constructor.
 */ 
GeneTreeLikelihood::GeneTreeLikelihood(const GeneTreeLikelihood & lik):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_ (lik.seqSp_), spId_(lik.spId_)
{
  WHEREAMI( __FILE__ , __LINE__ );
  levaluator_ = lik.levaluator_->clone(); 
  spTree_ = dynamic_cast<TreeTemplate<Node> *> (lik.spTree_->clone()) ;
  rootedTree_ = dynamic_cast<TreeTemplate<Node> *> (lik.rootedTree_->clone()) ;
  geneTreeWithSpNames_ = dynamic_cast<TreeTemplate<Node> *> (lik.geneTreeWithSpNames_->clone()) ;
  scenarioLikelihood_ = lik.scenarioLikelihood_;
  MLindex_ = lik.MLindex_;
  rootOptimization_ = lik.rootOptimization_; 
  tentativeMLindex_ = lik.MLindex_;
  totalIterations_ = lik.totalIterations_;
  counter_ = lik.counter_;
  _speciesIdLimitForRootPosition_ = lik._speciesIdLimitForRootPosition_;
  nodesToTryInNNISearch_ = lik.nodesToTryInNNISearch_;
  tentativeNodesToTryInNNISearch_ = lik.tentativeNodesToTryInNNISearch_;
  optimizeSequenceLikelihood_ = lik.optimizeSequenceLikelihood_;
  optimizeReconciliationLikelihood_ = lik.optimizeReconciliationLikelihood_ ;
  considerSequenceLikelihood_ = lik.considerSequenceLikelihood_;
  sprLimitGeneTree_ = lik.sprLimitGeneTree_;
}

GeneTreeLikelihood & GeneTreeLikelihood::operator=(const GeneTreeLikelihood & lik)
{
  WHEREAMI( __FILE__ , __LINE__ );
  if (levaluator_) delete levaluator_;
  levaluator_ = lik.levaluator_->clone(); 

  if (spTree_) delete spTree_;
  spTree_ = dynamic_cast<TreeTemplate<Node> *> (lik.spTree_->clone());
  if (rootedTree_) delete rootedTree_;
  rootedTree_= dynamic_cast<TreeTemplate<Node> *> (lik.rootedTree_->clone());
  if (geneTreeWithSpNames_) delete geneTreeWithSpNames_;
  geneTreeWithSpNames_ = dynamic_cast<TreeTemplate<Node> *> (lik.geneTreeWithSpNames_->clone()) ;
  spId_ = lik.spId_;
  scenarioLikelihood_ = lik.scenarioLikelihood_;
  MLindex_ = lik.MLindex_;
  rootOptimization_ = lik.rootOptimization_;
  tentativeMLindex_ = lik.MLindex_;
  totalIterations_ = lik.totalIterations_;
  counter_ = lik.counter_;
  _speciesIdLimitForRootPosition_ = lik._speciesIdLimitForRootPosition_;
  nodesToTryInNNISearch_ = lik.nodesToTryInNNISearch_;
  tentativeNodesToTryInNNISearch_ = lik.tentativeNodesToTryInNNISearch_;
  optimizeSequenceLikelihood_ = lik.optimizeSequenceLikelihood_;
  optimizeReconciliationLikelihood_ = lik.optimizeReconciliationLikelihood_ ;
  considerSequenceLikelihood_ = lik.considerSequenceLikelihood_;
  sprLimitGeneTree_ = lik.sprLimitGeneTree_;
  return *this;
}

void GeneTreeLikelihood::setGeneTree(TreeTemplate<Node>* tree, TreeTemplate<Node>* rootedTree) {
  WHEREAMI( __FILE__ , __LINE__ );
  if (rootedTree_) delete rootedTree_;
  rootedTree_= dynamic_cast<TreeTemplate<Node> *> (rootedTree->clone());
  //recreating geneTreeWithSpNames_
  if (geneTreeWithSpNames_) delete geneTreeWithSpNames_;
 	geneTreeWithSpNames_ = tree->clone();
	std::vector <Node*> leaves = geneTreeWithSpNames_->getLeaves();
	for (unsigned int j =0; j<leaves.size() ; j++) 
	  {
	    leaves[j]->setName(seqSp_[leaves[j]->getName()]);
	  }
	levaluator_->setAlternativeTree(rootedTree);
  levaluator_->acceptAlternativeTree();
}


