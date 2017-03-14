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

#include "Constants.h"
#include "DLGeneTreeLikelihood.h"
#include "GenericTreeExplorationAlgorithms.h"

using namespace bpp;

#define FAST 0


  /******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood( std::string file, 
					     															std::map<std::string, std::string> params, 
					     															TreeTemplate<Node> & spTree ) 
  throw (exception) : GeneTreeLikelihood(file, params, spTree)
  {
  WHEREAMI( __FILE__ , __LINE__ );
  size_t speciesTreeNodeNumber = spTree.getNumberOfNodes();
    for (int i=0; i< speciesTreeNodeNumber ; i++) 
    {
        //We fill the vectors with 0.1s until they are the right sizes.
        //DL model:
        lossExpectedNumbers_.push_back(0.3141514799);
        duplicationExpectedNumbers_.push_back(0.0016105081);
        num0Lineages_.push_back(0);
        num1Lineages_.push_back(0);
        num2Lineages_.push_back(0);
    }
  tentativeNum0Lineages_ =num0Lineages_;
  tentativeNum1Lineages_ =num1Lineages_; 
  tentativeNum2Lineages_ =num2Lineages_;
  }

/******************************************************************************/
// DLGeneTreeLikelihood::DLGeneTreeLikelihood(
//   const Tree & tree,
//   SubstitutionModel * model,
//   DiscreteDistribution * rDist,
//   TreeTemplate<Node> & spTree,
//   TreeTemplate<Node> & rootedTree,
//   TreeTemplate<Node> & geneTreeWithSpNames,
//   const std::map <std::string, std::string> seqSp,
//   std::map <std::string,int> spId,
//   //std::vector <int> & lossNumbers, 
//   std::vector <double> & lossProbabilities, 
//   //std::vector <int> & duplicationNumbers, 
//   std::vector <double> & duplicationProbabilities,
//   //std::vector <int> & branchNumbers,
//   std::vector <int> & num0Lineages,
//   std::vector <int> & num1Lineages,
//   std::vector <int> & num2Lineages, 
//   int speciesIdLimitForRootPosition,
//   int & MLindex, 
//   bool checkRooted,
//   bool verbose, 
//   bool rootOptimization, 
//   bool considerSequenceLikelihood, 
//   bool DLStartingGeneTree, 
//   unsigned int sprLimit)
// throw (Exception):
// GeneTreeLikelihood(tree,
// 		   model,
// 		   rDist,
// 		   spTree,  
// 		   rootedTree, 
// 		   geneTreeWithSpNames,
// 		   seqSp,
// 		   spId,
// 		   speciesIdLimitForRootPosition,
// 		   MLindex, 
// 		   checkRooted,
// 		   verbose,
// 		   rootOptimization, 
// 		   considerSequenceLikelihood, 
// 		   sprLimit)
// 
// {
//   lossExpectedNumbers_ = lossProbabilities;
//   duplicationExpectedNumbers_ = duplicationProbabilities;
//   num0Lineages_=num0Lineages;
//   num1Lineages_=num1Lineages;
//   num2Lineages_=num2Lineages;
//   tentativeNum0Lineages_ =num0Lineages;
//   tentativeNum1Lineages_ =num1Lineages; 
//   tentativeNum2Lineages_ =num2Lineages;
//   tentativeMLindex_ = MLindex;
//   DLStartingGeneTree_ = DLStartingGeneTree;
// }

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  TreeTemplate<Node> & spTree,
  TreeTemplate<Node> & rootedTree,
  TreeTemplate<Node> & geneTreeWithSpNames,
  const std::map <std::string, std::string> seqSp,
  std::map <std::string,int> spId,
  std::vector <double> & lossProbabilities,
  std::vector <double> & duplicationProbabilities,
  std::vector <int> & num0Lineages,
  std::vector <int> & num1Lineages,
  std::vector <int> & num2Lineages,
  int speciesIdLimitForRootPosition,
  int & MLindex,
  std::map <std::string, std::string > params,
  bool checkRooted,
  bool verbose,
  bool rootOptimization,
  bool considerSequenceLikelihood,
  bool DLStartingGeneTree,
  unsigned int sprLimitGeneTree)
throw (Exception):
GeneTreeLikelihood(tree,
		   data,
		   model,
		   rDist,
		   spTree,  
		   rootedTree, 
		   geneTreeWithSpNames,
		   seqSp,
		   spId,
		   speciesIdLimitForRootPosition,
		   MLindex, 
		   params,
		   checkRooted,
		   verbose,
		   rootOptimization, 
		   considerSequenceLikelihood, 
		   sprLimitGeneTree)
{
  WHEREAMI( __FILE__ , __LINE__ );
  lossExpectedNumbers_ = lossProbabilities;
  duplicationExpectedNumbers_ = duplicationProbabilities; 
  num0Lineages_=num0Lineages;
  num1Lineages_=num1Lineages;
  num2Lineages_=num2Lineages;
  tentativeNum0Lineages_ =num0Lineages;
  tentativeNum1Lineages_ =num1Lineages; 
  tentativeNum2Lineages_ =num2Lineages;
  DLStartingGeneTree_ = DLStartingGeneTree;
}

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(const DLGeneTreeLikelihood & lik):
GeneTreeLikelihood(lik)
{
  WHEREAMI( __FILE__ , __LINE__ );
  lossExpectedNumbers_ = lik.lossExpectedNumbers_;
  duplicationExpectedNumbers_ = lik.duplicationExpectedNumbers_; 
  num0Lineages_=lik.num0Lineages_;
  num1Lineages_=lik.num1Lineages_;
  num2Lineages_=lik.num2Lineages_;
  tentativeNum0Lineages_ =lik.tentativeNum0Lineages_;
  tentativeNum1Lineages_ =lik.tentativeNum1Lineages_;
  tentativeNum2Lineages_ =lik.tentativeNum2Lineages_;
  DLStartingGeneTree_ = lik.DLStartingGeneTree_;
}

/******************************************************************************/

DLGeneTreeLikelihood & DLGeneTreeLikelihood::operator=(const DLGeneTreeLikelihood & lik)
{
  WHEREAMI( __FILE__ , __LINE__ );
  GeneTreeLikelihood::operator=(lik);
  lossExpectedNumbers_ = lik.lossExpectedNumbers_;
  duplicationExpectedNumbers_ = lik.duplicationExpectedNumbers_;
  num0Lineages_=lik.num0Lineages_;
  num1Lineages_=lik.num1Lineages_;
  num2Lineages_=lik.num2Lineages_;
  tentativeNum0Lineages_ =lik.tentativeNum0Lineages_;
  tentativeNum1Lineages_ =lik.tentativeNum1Lineages_;
  tentativeNum2Lineages_ =lik.tentativeNum2Lineages_;
  DLStartingGeneTree_ = lik.DLStartingGeneTree_;
  return *this;
}




/******************************************************************************/

DLGeneTreeLikelihood::~DLGeneTreeLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  if (levaluator_) delete levaluator_;
  if (spTree_) delete spTree_;
  if (rootedTree_) delete rootedTree_;
  if (geneTreeWithSpNames_) delete geneTreeWithSpNames_;
}



/******************************************************************************/

void DLGeneTreeLikelihood::initParameters()
{
  WHEREAMI( __FILE__ , __LINE__ );
  if (considerSequenceLikelihood_) {
    //TODO What to do here with levaluator?
    //nniLk_->initParameters();
  }
      //Here we check whether there are several trees inside the gene tree file (if init.gene.tree was user).
      //If there are several trees, then we test each of them, and we will choose the tree
      //with the highest total likelihood as a starting point.
      vector<Tree*> trees;
      string initTree = ApplicationTools::getStringParameter ( "init.gene.tree", params_, "user", "", false, false );
      if (initTree == "user") {
          std::string geneTree_File =ApplicationTools::getStringParameter ( "gene.tree.file",params_,"none" );
          IMultiTree* treeReader;
          treeReader = new Newick(true);
          treeReader->read(geneTree_File, trees);
          delete treeReader;
      }
      //std::cout << "FIRST : leval loglk: "<< levaluator_->getLogLikelihood() << " and scenario loglk: "<< scenarioLikelihood_<< std::endl;

      if (trees.size() > 1) {
          //We have several trees, we need to choose which one is the best!
          double bestSequenceLogL;
	  findBestGeneTreeAmongCandidates(trees, rootedTree_, bestSequenceLogL, scenarioLikelihood_);
	  //std::cout << "SECOND: leval loglk: "<< levaluator_->getLogLikelihood() << " and scenario loglk: "<< scenarioLikelihood_<< std::endl;

      }
      else {
      scenarioLikelihood_ = findMLReconciliationDR (spTree_, rootedTree_,
                                                    seqSp_, spId_, lossExpectedNumbers_,
                                                    duplicationExpectedNumbers_, MLindex_,
                                                    num0Lineages_, num1Lineages_,
                                                    num2Lineages_, nodesToTryInNNISearch_);
      }
      //std::cout << "THIRD: leval loglk: "<< levaluator_->getLogLikelihood() << " and scenario loglk: "<< scenarioLikelihood_<< std::endl;

      //Rooting the tree properly:
      vector<Node*> nodes = rootedTree_->getNodes();
      for (unsigned int j = 0 ; j < nodes.size() ; j++) {
          
          if (nodes[j]->hasNodeProperty("outgroupNode")) {
              
              if (rootedTree_->getRootNode() == nodes[j]) {
                  if (j < nodes.size()-1)
                  {
                      rootedTree_->rootAt(nodes[nodes.size()-1]);   
                  }
                  else {
                      rootedTree_->rootAt(nodes[nodes.size()-2]);
                  }
              };
              rootedTree_->newOutGroup( nodes[j] );
              break;
          }
      }
  MLindex_ = -1;

}



/******************************************************************************/

void DLGeneTreeLikelihood::resetMLindex() {
  WHEREAMI( __FILE__ , __LINE__ );
  MLindex_ = -1;
}




/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double DLGeneTreeLikelihood::getLogLikelihood() const
{
  WHEREAMI( __FILE__ , __LINE__ );
  double ll = 0;

  
  if (considerSequenceLikelihood_) {
    //TODO: debug to remove
    //cout << "(DTL) sequence likelihood=" <<  levaluator_->getLogLikelihood() << "/nscenario likelihood=" << scenarioLikelihood_ << endl;
    ll = levaluator_->getLogLikelihood() + scenarioLikelihood_;
  }
  else {
    ll = scenarioLikelihood_;
  }
  return ll;
}
/******************************************************************************/
//returns -loglikelihood
//As a reminder: _sequenceLikelihood < 0, scenarioLikelihood_ < 0, getValue >0
double DLGeneTreeLikelihood::getValue() const  
throw (Exception)
{
  {
  WHEREAMI( __FILE__ , __LINE__ );
    if (considerSequenceLikelihood_) {
      return (- levaluator_->getLogLikelihood() - scenarioLikelihood_);
    }
    else {
      return (- scenarioLikelihood_);
    }
  }
}



/******************************************************************************/
/*
 * was for BPP optimizable objects
 * 

// void DLGeneTreeLikelihood::fireParameterChanged(const ParameterList & params)
// {
//   if (considerSequenceLikelihood_) {
//     levaluator_->getNniLk()->applyParameters();
//     levaluator_->getNniLk()->matchParametersValues(params);
//   }
//   //If we need to update the reconciliation likelihood
//   if (optimizeReconciliationLikelihood_) {
//     computeReconciliationLikelihood();
//   }
//   
//   
// }*/


/******************************************************************************/


void DLGeneTreeLikelihood::computeSequenceLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  if ( considerSequenceLikelihood_ && (optimizeSequenceLikelihood_==true) ) {
//     with likelihood evaluator, likelihood are always computed when the tree is modified
//     nniLk()->computeTreeLikelihood ();
  }
}




/******************************************************************************/

void DLGeneTreeLikelihood::computeReconciliationLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  resetLossesAndDuplications(*spTree_, /*lossNumbers_, */lossExpectedNumbers_, /*duplicationNumbers_, */duplicationExpectedNumbers_);
    scenarioLikelihood_ = findMLReconciliationDR (spTree_, rootedTree_, 
						  seqSp_, spId_, 
						  lossExpectedNumbers_, duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_); 
    MLindex_ = tentativeMLindex_;
    num0Lineages_ = tentativeNum0Lineages_;
    num1Lineages_ = tentativeNum1Lineages_;
    num2Lineages_ = tentativeNum2Lineages_;
    nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
  
}



/******************************************************************************/


void DLGeneTreeLikelihood::computeTreeLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  if (considerSequenceLikelihood_ )
  {
    computeSequenceLikelihood();
  }
  computeReconciliationLikelihood();  
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.
 * 
 */
double DLGeneTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  WHEREAMI( __FILE__ , __LINE__ );
  //If the NNI is around a branch where a duplication was found, 
  //or if we just try all branches because the starting gene trees are parsimonious in
  //numbers of DL.
  if ((nodesToTryInNNISearch_.count(nodeId)==1) /*|| DLStartingGeneTree_*/) {
    TreeTemplate<Node> * treeForNNI = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
    
    tentativeMLindex_ = MLindex_;
    tentativeNum0Lineages_ = num0Lineages_;
    tentativeNum1Lineages_ = num1Lineages_;
    tentativeNum2Lineages_ =num2Lineages_;
    tentativeNodesToTryInNNISearch_.clear();
    
    //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
    //We use the same procedure as in doNNI !
    const Node * son    = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->getNode(nodeId);
    
    
    if(!son->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'son' must not be the root node.", nodeId);
    const Node * parent = son->getFather();
    
    if(!parent->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent->getId());
    const Node * grandFather = parent->getFather();
    
    //From here: Bifurcation assumed.
    //In case of multifurcation, an arbitrary uncle is chosen.
    //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
    unsigned int parentPosition = grandFather->getSonPosition(parent);
    
    Node * sonForNNI    = treeForNNI->getNode(nodeId);
    Node * parentForNNI = sonForNNI->getFather();
    Node * grandFatherForNNI = parentForNNI->getFather();
    Node * uncleForNNI = grandFatherForNNI->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
    parentForNNI->removeSon(sonForNNI);
    grandFatherForNNI->removeSon(uncleForNNI);
    parentForNNI->addSon(uncleForNNI);
    grandFatherForNNI->addSon(sonForNNI);
    
    //Now we root the tree sent to findMLReconciliation as in rootedTree_
    int id = treeForNNI->getRootNode()->getId();
    if(TreeTemplateTools::hasNodeWithId(*(rootedTree_->getRootNode()->getSon(0)),id)) {
      treeForNNI->newOutGroup(rootedTree_->getRootNode()->getSon(1)->getId());
    }
    else {
      treeForNNI->newOutGroup(rootedTree_->getRootNode()->getSon(0)->getId());
    }
    double candidateScenarioLk = 0;
    totalIterations_ = totalIterations_+1;
    

    candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForNNI/*&rootedTree_*/, 
					       seqSp_, spId_, 
					       lossExpectedNumbers_, duplicationExpectedNumbers_, 
					       tentativeMLindex_, 
					       tentativeNum0Lineages_, tentativeNum1Lineages_, 
					       tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_, false); //false so that _tentativeNum*Lineages are not updated
    
    if (considerSequenceLikelihood_ ) 
    {
      if  (candidateScenarioLk >  scenarioLikelihood_) 
      { //If it is worth computing the sequence likelihood
                levaluator_->setAlternativeTree(treeForNNI);             

        double tot = -( candidateScenarioLk + levaluator_->getAlternativeLogLikelihood() ) + ( getSequenceLikelihood() + scenarioLikelihood_ ) ;


        if (tot < 0)
        {
          tentativeScenarioLikelihood_= candidateScenarioLk;
        }
        return tot;
	
      }
      else 
      {
        tentativeMLindex_ = -1;
        return 1;
      }
    }
    else 
    {
      if  (candidateScenarioLk >  scenarioLikelihood_) 
      {
        tentativeScenarioLikelihood_= candidateScenarioLk;
        return (- candidateScenarioLk + scenarioLikelihood_);
      }
      else 
      {
        tentativeMLindex_ = -1;
        return 1;
      }
    }
  }
  else {
    tentativeMLindex_ = -1;
    return 1;
  }
}

/*******************************************************************************/

void DLGeneTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  WHEREAMI( __FILE__ , __LINE__ );
     levaluator_->acceptAlternativeTree();


  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  
  MLindex_ = tentativeMLindex_;
  
  nodesToTryInNNISearch_ = tentativeNodesToTryInNNISearch_;
  
  // std::cout <<"tentativeScenarioLikelihood_: "<<tentativeScenarioLikelihood_<<std::endl;
  scenarioLikelihood_ = tentativeScenarioLikelihood_;// + _brLikFunction->getValue();
  
  //Now we need to update rootedTree_
  TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
  //First we root this temporary tree as in rootedTree_ (same lines as in testNNI)
  int id = tree->getRootNode()->getId();
  if(TreeTemplateTools::hasNodeWithId(*(rootedTree_->getRootNode()->getSon(0)),id)) {
    tree->newOutGroup(rootedTree_->getRootNode()->getSon(1)->getId());
  }
  else {
    tree->newOutGroup(rootedTree_->getRootNode()->getSon(0)->getId());
  }
  //Then we root this tree according to MLindex
  tree->newOutGroup(MLindex_);
  //We update rootedTree_
  if (rootedTree_) delete rootedTree_;
  rootedTree_ = tree->clone();
  delete tree;
  //we need to update the sequence likelihood
  if (considerSequenceLikelihood_ ) 
  {
    OptimizeSequenceLikelihood(true);
  }
  OptimizeReconciliationLikelihood(true);
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getDuplicationNumbers(){
  WHEREAMI( __FILE__ , __LINE__ );
  return duplicationNumbers_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getLossNumbers(){
  WHEREAMI( __FILE__ , __LINE__ );
  return lossNumbers_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getBranchNumbers(){
  WHEREAMI( __FILE__ , __LINE__ );
  return branchNumbers_;
}


/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get0LineagesNumbers() const {
  WHEREAMI( __FILE__ , __LINE__ );
  return num0Lineages_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get1LineagesNumbers() const {
  WHEREAMI( __FILE__ , __LINE__ );
  return num1Lineages_;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get2LineagesNumbers() const {
  WHEREAMI( __FILE__ , __LINE__ );
  return num2Lineages_;
}

/*******************************************************************************/

void DLGeneTreeLikelihood::setExpectedNumbers(std::vector <double> duplicationProbabilities, std::vector <double> lossProbabilities){
  WHEREAMI( __FILE__ , __LINE__ );
  
  lossExpectedNumbers_ = lossProbabilities;
  duplicationExpectedNumbers_ = duplicationProbabilities; 
}

/*******************************************************************************/

int DLGeneTreeLikelihood::getRootNodeindex(){
  WHEREAMI( __FILE__ , __LINE__ );
  return MLindex_;
}

/*******************************************************************************/
/*
 * void DLGeneTreeLikelihood::resetSequenceLikelihood(){
 *    _sequenceLikelihood = UNLIKELY;
 * }
 */
/*******************************************************************************/

double DLGeneTreeLikelihood::getSequenceLikelihood() const {
  WHEREAMI( __FILE__ , __LINE__ );
  return levaluator_->getLogLikelihood(); 
}
/*******************************************************************************/

void DLGeneTreeLikelihood::initialize() {
  WHEREAMI( __FILE__ , __LINE__ );
  levaluator_->initialize();

  // Update the rooted tree so that its branch lengths are as optimized by the levaluator.
    TreeTemplate<Node> * treeTemp = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
    
    
    //Now we root the tree sent to findMLReconciliation as in rootedTree_
    int id = treeTemp->getRootNode()->getId();
    if(TreeTemplateTools::hasNodeWithId(*(rootedTree_->getRootNode()->getSon(0)),id)) {
      treeTemp->newOutGroup(rootedTree_->getRootNode()->getSon(1)->getId());
    }
    else {
      treeTemp->newOutGroup(rootedTree_->getRootNode()->getSon(0)->getId());
    }

    if (rootedTree_) 
      {
	delete rootedTree_;
        rootedTree_ = 0;
      }
    rootedTree_ = treeTemp->clone();
    breadthFirstreNumberAndResetProperties (*rootedTree_);
    delete treeTemp;

  return;
}



/*******************************************************************************/
void DLGeneTreeLikelihood::print () const {
  WHEREAMI( __FILE__ , __LINE__ );
  std::cout << "Species tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getSpTree(), true)<<std::endl;
  std::cout << "Gene family rooted tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
  std::cout << "Gene family tree:"<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (*(levaluator_->getTree()), true)<<std::endl;
  std::cout << "0 lineage numbers"<<std::endl;
  VectorTools::print(get0LineagesNumbers());
  std::cout << "1 lineage numbers"<<std::endl;
  VectorTools::print(get1LineagesNumbers());
  std::cout << "2 lineages numbers"<<std::endl;
  VectorTools::print(get2LineagesNumbers());
  std::cout << "Expected numbers of losses"<<std::endl;
  VectorTools::print(lossExpectedNumbers_);
  std::cout << "Expected numbers of duplications"<<std::endl;
  VectorTools::print(duplicationExpectedNumbers_);
  std::cout << "Root index"<<std::endl;
  std::cout << MLindex_ <<std::endl;
  
  
}





/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * To do all this as fast as possible, we optimize only a few branch lengths on the SPR tree, 
 * and we use a simple recursion for that.
 * WORKS.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRsFast2 (map<string, string> params) {
  WHEREAMI( __FILE__ , __LINE__ );
  
  double startingTime = ApplicationTools::getTime();
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    computeReconciliationLikelihood();
    return;
  }
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  // if (getLogLikelihood()==UNLIKELY) 
  computeReconciliationLikelihood();
  double logL = getLogLikelihood();
  
  double bestlogL = logL;
  double candidateScenarioLk ;
  double bestSequenceLogL = getSequenceLikelihood();
  double bestScenarioLk = getScenarioLikelihood();
  //std::cout << "LOGL: "<<logL << " ScenarioLK: "<< bestScenarioLk <<"; sequenceLK: "<<getSequenceLikelihood() << std::endl;
  unsigned int numIterationsWithoutImprovement = 0;
  breadthFirstreNumber (*rootedTree_);
  
  string parentDup;
  string nodeDup;
  string numLoss = "0";
  
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    	    
  //std::cout << "DEBUG: DLGeneTreeLk BEFORE: "<< TreeTemplateTools::treeToParenthesis(*rootedTree_)<<std::endl;

  
  while (numIterationsWithoutImprovement < rootedTree_->getNumberOfNodes() - 2 && (timeLimit_ == 0 || elapsedTime_ < timeLimit_))
  {
    annotateGeneTreeWithDuplicationEvents (*spTree_, 
					   *rootedTree_, 
					   rootedTree_->getRootNode(), 
					   seqSp_, spId_); 
    
    for (int nodeForSPR=rootedTree_->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
    {
      Node * n = rootedTree_->getNode(nodeForSPR);
     /* if (n->hasBranchProperty("L")) {
        numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
      }*/
      if ( /*numLoss != "0"*/ 1  ) {
	
	buildVectorOfRegraftingNodesGeneTree(*spTree_, *rootedTree_, nodeForSPR, sprLimitGeneTree_, nodeIdsToRegraft);
	betterTree = false;
	for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
	{
	  if (treeForSPR) 
	  {
	    delete treeForSPR;
	    treeForSPR = 0;
	  }
	  treeForSPR = rootedTree_->clone();
	  
	  nodesToUpdate = makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false, true);
	  
	  //Compute the DL likelihood
	  candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
							 seqSp_, spId_, 
						  lossExpectedNumbers_, 
						  duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, 
						  tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, 
						  tentativeNodesToTryInNNISearch_, false); 
	  if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the DL likelihood is not bad
	  {    

	    if (computeSequenceLikelihoodForSPR) {     
	      levaluator_->setAlternativeTree(treeForSPR);
 	      logL = candidateScenarioLk + levaluator_->getAlternativeLogLikelihood();
	    }
	    else {
	      logL = candidateScenarioLk + bestSequenceLogL;
	    }
	  }
	  else { 
	    logL =logL - 10;
	  }
	  
	  //std::cout << "found logLk: "<< logL << "compared to" << bestlogL<< " bestSequenceLogL: "<< bestSequenceLogL << " levaluator_->getAlternativeLogLikelihood(): "<<levaluator_->getAlternativeLogLikelihood() <<" bestScenarioLk: "<< bestScenarioLk << " candidateScenarioLk: "<< candidateScenarioLk <<std::endl;
    
	  //If the candidate tree has a DL + sequence Lk better than the current best
	  if (logL - 0.01 > bestlogL) 
	  {
            levaluator_->acceptAlternativeTree();
      WHEREAMI( __FILE__ , __LINE__ );
	    std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
	    
	    betterTree = true;
	    bestlogL = logL;
	    bestScenarioLk = candidateScenarioLk;
	    if (computeSequenceLikelihoodForSPR) {
	      bestSequenceLogL = levaluator_->getAlternativeLogLikelihood();
	    }
	    if (bestTree) {
	      delete bestTree;
	      bestTree = 0;
	    }
	    
	    bestTree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();

	    //Rooting bestTree as in TreeForSPR:
	    vector<Node*> rlkNodes = bestTree->getNodes();
	    for (unsigned int j = 0 ; j < rlkNodes.size() ; j++) {
	      if (rlkNodes[j]->hasNodeProperty("outgroupNode")) {
		if (bestTree->getRootNode() == rlkNodes[j]) {
		  if (j < rlkNodes.size()-1) 
		  {
		    bestTree->rootAt(rlkNodes[rlkNodes.size()-1]);   
		  }
		  else {
		    bestTree->rootAt(rlkNodes[rlkNodes.size()-2]);
		  }
		};
		bestTree->newOutGroup( rlkNodes[j] );
		//  std::cout << "FOUND"<<std::endl;
		break;
	      }
	    }
	    writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
	    
	  }
	  
	}
	if (betterTree) //If, among all the SPRs tried, a better tree has been found 
	{
	  
	  logL = bestlogL; 
	  numIterationsWithoutImprovement = 0;
	  if (treeForSPR) 
	  {
	    delete treeForSPR;
	    treeForSPR = 0;
	  }
	  if (rootedTree_) 
	  {
	    delete rootedTree_;
	    rootedTree_ = 0;
	  }
	  rootedTree_ = bestTree->clone();
	  breadthFirstreNumberAndResetProperties (*rootedTree_);
	  
	  if (bestTree) {
	    delete bestTree;
	    bestTree = 0;
	  }
	  
	  
	}
	else 
	{
	  numIterationsWithoutImprovement++;
	}
	
	if (treeForSPR) 
	{
	  delete treeForSPR;
	  treeForSPR = 0;
	}
	
	
      }
      else {
	numIterationsWithoutImprovement++;
      }
    }
    elapsedTime_ += (ApplicationTools::getTime() - startingTime);
    startingTime = ApplicationTools::getTime();
  }
  
  rootedTree_->resetNodesId();
  
  
  //One more reconciliation, to update the "_num*Lineages" vectors.
  computeReconciliationLikelihood();
  
  Nhx *nhx = new Nhx();
  annotateGeneTreeWithDuplicationEvents (*spTree_, 
					 *rootedTree_, 
					 rootedTree_->getRootNode(), 
					 seqSp_, spId_); 

  cout << "Reconciled tree: "<<endl;
  nhx->write(*rootedTree_, cout);

  if (bestTree) {
    delete bestTree;
    bestTree = 0;
  }
  
  if (nhx) delete nhx;
  
}


/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * To do all this as fast as possible, we compute the sequence likelihood only for the most promising SPRs
 * given their reconciliation score.
 * We optimize all branch lengths on the SPR tree.
 * We use a simple recursion for that.
 * WORKS.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRsFast3 (map<string, string> params) {
  WHEREAMI( __FILE__ , __LINE__ );
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    //std::cout << "WE DONT DO SPRS"<<std::endl;
    computeReconciliationLikelihood();
    return;
  }
  //   std::cout << "WE DO SPRS; seq lk: "<< getSequenceLikelihood()  <<std::endl;
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
  // if (getLogLikelihood()==UNLIKELY) 
  computeReconciliationLikelihood();
  double logL = getLogLikelihood();

  double bestlogL = logL;
  double candidateScenarioLk ;
  double bestSequenceLogL = getSequenceLikelihood();
  double bestScenarioLk = getScenarioLikelihood();

  unsigned int numIterationsWithoutImprovement = 0;
  breadthFirstreNumber (*rootedTree_);
  
  
  string parentDup;
  string nodeDup;
  string numLoss = "0";
  std::map < double, TreeTemplate<Node> * >  treesToOptimizeSeqLk ;
  std::map < double, TreeTemplate<Node> * >::reverse_iterator it;
  double bestCurrentCandidateScenarioLk;
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
  
  
  while (numIterationsWithoutImprovement < rootedTree_->getNumberOfNodes() - 2)
  {
    
    annotateGeneTreeWithDuplicationEvents (*spTree_, 
					   *rootedTree_, 
					   rootedTree_->getRootNode(), 
					   seqSp_, spId_); 
    
    for (unsigned int nodeForSPR=rootedTree_->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
    {
      Node * n = rootedTree_->getNode(nodeForSPR);
      if (n->hasBranchProperty("L")) {
	numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
      }
      if ( numLoss != "0"  ) {
	treesToOptimizeSeqLk.clear();
	buildVectorOfRegraftingNodesGeneTree(*spTree_, *rootedTree_, nodeForSPR, sprLimitGeneTree_, nodeIdsToRegraft);
	betterTree = false;
	for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
	{
	  if (treeForSPR) 
	  {
	    delete treeForSPR;
	    treeForSPR = 0;
	  }
	  treeForSPR = rootedTree_->clone();
	  nodesToUpdate = makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false, true);
	  
	  //Compute the DL likelihood
	  candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
							 seqSp_, spId_, 
						  lossExpectedNumbers_, 
						  duplicationExpectedNumbers_, 
						  tentativeMLindex_, 
						  tentativeNum0Lineages_, 
						  tentativeNum1Lineages_, 
						  tentativeNum2Lineages_, 
						  tentativeNodesToTryInNNISearch_, false); 
	  while (treesToOptimizeSeqLk.count(candidateScenarioLk) > 0 )
	  {
	    candidateScenarioLk = candidateScenarioLk + NumConstants::SMALL();
	  }
	  treesToOptimizeSeqLk[candidateScenarioLk] = treeForSPR->clone();
	}
	if (treesToOptimizeSeqLk.size() > 0)
	  bestCurrentCandidateScenarioLk = (* (treesToOptimizeSeqLk.rbegin() ) ).first;
	for ( it=treesToOptimizeSeqLk.rbegin() ; it != treesToOptimizeSeqLk.rend(); it++ ) {
	  candidateScenarioLk = (*it).first;
	  if ( ( candidateScenarioLk > bestScenarioLk) && ( candidateScenarioLk < bestCurrentCandidateScenarioLk + 10 * NumConstants::SMALL()) )//If the likelihood is better than the current best, and one of the best ones to test
	  {
	    //We compute the sequence likelihood
	    if (computeSequenceLikelihoodForSPR) {
	      if (treeForSPR) 
	      {
		delete treeForSPR;
		treeForSPR = 0;
	      }
	      treeForSPR = (*it).second->clone();
	      
              levaluator_->setAlternativeTree(treeForSPR);
              
	      logL = candidateScenarioLk + levaluator_->getAlternativeLogLikelihood();
	      
	      
	      //If the candidate tree has a DL + sequence Lk better than the current best
	      if (logL - 0.01 > bestlogL) 
	      {
                // levaluator: since it the best tree, setting it as the current one
                levaluator_->acceptAlternativeTree();
    WHEREAMI( __FILE__ , __LINE__ );
		std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
		betterTree = true;
		bestlogL = logL;
		bestScenarioLk = candidateScenarioLk;
		if (computeSequenceLikelihoodForSPR) {
		  bestSequenceLogL = levaluator_->getLogLikelihood();
		}
		if (bestTree) {
		  delete bestTree;
		  bestTree = 0;
		}
		bestTree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
		//Rooting bestTree as in TreeForSPR:
		vector<Node*> rlkNodes = bestTree->getNodes();
		for (unsigned int j = 0 ; j < rlkNodes.size() ; j++) {
		  if (rlkNodes[j]->hasNodeProperty("outgroupNode")) {
		    if (bestTree->getRootNode() == rlkNodes[j]) {
		      if (j < rlkNodes.size()-1) 
		      {
			bestTree->rootAt(rlkNodes[rlkNodes.size()-1]);   
		      }
		      else {
			bestTree->rootAt(rlkNodes[rlkNodes.size()-2]);
		      }
		    };
		    bestTree->newOutGroup( rlkNodes[j] );
		    //  std::cout << "FOUND"<<std::endl;
		    break;
		  }
		}
		writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
		break; //If we have found one topology better than the current one for seqlk+scenlk
	      }
	    }
	    else {
	      logL = candidateScenarioLk - bestSequenceLogL;
	    }
	    
	  }
	  else {
	    break; //No more likelihoods better than the current best
	  }
	}
	if (betterTree) //If, among all the SPRs tried, a better tree has been found 
	{
	  
	  logL = bestlogL; 
	  numIterationsWithoutImprovement = 0;
	  if (treeForSPR) 
	  {
	    delete treeForSPR;
	    treeForSPR = 0;
	  }
	  if (rootedTree_) 
	  {
	    delete rootedTree_;
	    rootedTree_ = 0;
	  }
	  rootedTree_ = bestTree->clone();
	  scenarioLikelihood_ = bestScenarioLk;
	  //  breadthFirstreNumber (*rootedTree_);
	  breadthFirstreNumberAndResetProperties (*rootedTree_);
	  
	  if (bestTree) {
	    delete bestTree;
	    bestTree = 0;
	  }
	  
	}
	else 
	{
	  numIterationsWithoutImprovement++;
	}
	
	if (treeForSPR) 
	{
	  delete treeForSPR;
	  treeForSPR = 0;
	}
      }
      else
      {
	numIterationsWithoutImprovement++;
      }
    }
  }
  
  rootedTree_->resetNodesId();
    
  //One more reconciliation, to update the "_num*Lineages" vectors.
  computeReconciliationLikelihood();
  
  Nhx *nhx = new Nhx();
  annotateGeneTreeWithDuplicationEvents (*spTree_, 
					 *rootedTree_, 
					 rootedTree_->getRootNode(), 
					 seqSp_, spId_); 
  cout << "Reconciled tree: "<<endl;
  nhx->write(*rootedTree_, cout);
  
  if (bestTree) {
    delete bestTree;
    bestTree = 0;
  }
  
  if (nhx) delete nhx;
  
}


/************************************************************************
 * Implements the moves by Muffato and Crollius. Goes through the tree
 * while computing the number of genes ng and the number of species ns downstream, 
 * and in cases where a duplication is predicted but ng << 2 * ns, 
 * rearranges the tree to resolve the duplication.
 * First classical double-recursive tree traversals are used to find the root of the tree,
 * then another post-order tree traversal is done to compute ng and ns, 
 * and at the same time produce rearranged trees. Then the DL likelihood is computed
 * for all rearranged trees, and the best ones have their sequence likelihood computed.
 * 
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeMuffato (map<string, string> params) {
  WHEREAMI( __FILE__ , __LINE__ );

  computeReconciliationLikelihood();
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do SPRs
    return;
  }
  std::vector<Node*> nodesToUpdate;
  std::vector <int> nodeIdsToRegraft;
  bool betterTree = false;
  TreeTemplate<Node> * treeForSPR = 0;
  TreeTemplate<Node> * bestTree = 0;
    
  double candidateScenarioLk ;
  double bestSequenceLogL = getSequenceLikelihood();
  double bestScenarioLk = getScenarioLikelihood();
  candidateScenarioLk = bestScenarioLk;
  unsigned int numIterationsWithoutImprovement = 0;
  breadthFirstreNumber (*rootedTree_);
  
    
  string parentDup;
  string nodeDup;
  string numLoss = "0";
  std::map < double, TreeTemplate<Node> * >  treesToOptimizeSeqLk ;
  std::map < double, TreeTemplate<Node> * >::reverse_iterator it;
  //double bestCurrentCandidateScenarioLk;
  bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
  
  double logL = 0;
  double bestlogL = 0;
  
  while (1)
  {
    std::cout << "\t\tStarting Muffato optimization loop"<<std::endl;
    
    annotateGeneTreeWithScoredDuplicationEvents (*spTree_, 
						 *rootedTree_, 
						 rootedTree_->getRootNode(), 
						 seqSp_, spId_); 
    
    double editionThreshold = ApplicationTools::getDoubleParameter("muffato.edition.threshold", params, 0.3, "", false, false);
    
    treeForSPR = rootedTree_->clone();
    bool edited = editDuplicationNodesMuffato(*spTree_, *treeForSPR, editionThreshold);

    //editDuplicationNodesMuffato(*spTree_, *treeForSPR, treeForSPR->getRootNode(), editionThreshold);
    //    std::cout <<"Here 13"<<std::endl;
    if (edited) {
    //Compute the DL likelihood
    candidateScenarioLk =  findMLReconciliationDR (spTree_, treeForSPR, 
						   seqSp_, spId_, 
						   lossExpectedNumbers_, 
						   duplicationExpectedNumbers_, 
						   tentativeMLindex_, 
						   tentativeNum0Lineages_, 
						   tentativeNum1Lineages_, 
						   tentativeNum2Lineages_, 
						   tentativeNodesToTryInNNISearch_, false); 
    std::cout <<"Candidate scenario lk for Muffato rearranged tree: "<< candidateScenarioLk <<std::endl;
    
    
    //We compute the sequence likelihood
    if (computeSequenceLikelihoodForSPR) {
      
      unsigned int numEval = 0;
      
      levaluator_->setAlternativeTree(treeForSPR);
      
      logL = candidateScenarioLk + levaluator_->getAlternativeLogLikelihood();
      
      //If the candidate tree has a DL + sequence Lk better than the current best
      if (logL - 0.1 > bestlogL) 
      {
        // then we accept the tree
        levaluator_->acceptAlternativeTree();
  WHEREAMI( __FILE__ , __LINE__ );      
	std::cout << "Better tree overall: "<<logL << " compared to "<<bestlogL<<std::endl;
	betterTree = true;
	bestlogL = logL;
	bestScenarioLk = candidateScenarioLk;
	if (computeSequenceLikelihoodForSPR) {
	  bestSequenceLogL = levaluator_->getAlternativeLogLikelihood();
	}
	
	if (bestTree) {
	  delete bestTree;
	  bestTree = 0;
	}
	bestTree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
	//Rooting bestTree as in TreeForSPR:
	vector<Node*> rlkNodes = bestTree->getNodes();
	
	for (unsigned int j = 0 ; j < rlkNodes.size() ; j++) {
	  if (rlkNodes[j]->hasNodeProperty("outgroupNode")) {
	    if (bestTree->getRootNode() == rlkNodes[j]) {
	      if (j < rlkNodes.size()-1) 
	      {
		bestTree->rootAt(rlkNodes[rlkNodes.size()-1]);   
	      }
	      else {
		bestTree->rootAt(rlkNodes[rlkNodes.size()-2]);
	      }
	    };
	    bestTree->newOutGroup( rlkNodes[j] );
	    break;
	  }
	}
	writeReconciledGeneTree ( params, dynamic_cast<const TreeTemplate<Node> *> ((bestTree))->clone(), spTree_, seqSp_, true ) ;
	
      }
      else {
	break;
	
      }
    }
    else {
      
      logL = candidateScenarioLk - bestSequenceLogL;
    }
  }
    if (betterTree) //If, among all the SPRs tried, a better tree has been found 
    {
      
      logL = bestlogL; 
      numIterationsWithoutImprovement = 0;
      if (treeForSPR) 
      {
	delete treeForSPR;
	treeForSPR = 0;
      }
      
      if (rootedTree_) 
      {
	delete rootedTree_;
	rootedTree_ = 0;
      }
      
      rootedTree_ = bestTree->clone();
      scenarioLikelihood_ = bestScenarioLk;
      breadthFirstreNumberAndResetProperties (*rootedTree_);
      
      if (bestTree) {
	delete bestTree;
	bestTree = 0;
      }
      
    }
    else {
      break;
    }
  }
  
  rootedTree_->resetNodesId();
  
  
  //One more reconciliation, to update the "_num*Lineages" vectors.
  computeReconciliationLikelihood();
  
  Nhx *nhx = new Nhx();
  annotateGeneTreeWithDuplicationEvents (*spTree_, 
					 *rootedTree_, 
					 rootedTree_->getRootNode(), 
					 seqSp_, spId_); 
  cout << "Muffato reconciled tree: "<<endl;
  nhx->write(*rootedTree_, cout);
  
  if (bestTree) {
    delete bestTree;
    bestTree = 0;
  }

  if (nhx) delete nhx;
  
}








/************************************************************************
 * Tries all NNIs, and accepts NNIs that improve the likelihood as soon as
 * they have been tried.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeNNIs(map<string, string> params, unsigned int verbose ) {
  WHEREAMI( __FILE__ , __LINE__ );
  
  double startingTime = 0;
  
  if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
    //We don't do NNIs
    computeReconciliationLikelihood();
    return;
  }
  bool test = true;
  do
  { 
    TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();
    vector<Node *> nodes = tree->getNodes();
    
    vector<Node *> nodesSub = nodes;
    for(unsigned int i = nodesSub.size(); i>0; i--)
    {
      // !!! must not reach i==0 because of unsigned int
      if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.      
      else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node. 
    }
    
    // Test all NNIs:
    test = false;
    for(unsigned int i = 0; !test && i < nodesSub.size(); i++)
    {
      Node* node = nodesSub[i];
      double diff = testNNI(node->getId());
      if(verbose >= 3)
      {
	ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
	+ " at " + TextTools::toString(node->getFather()->getId()),
					TextTools::toString(diff));
      }
      
      if(diff < 0.)
      { //Good NNI found...
	if(verbose >= 2)
	{
	  ApplicationTools::displayResult("   Swapping node " + TextTools::toString(node->getId())
	  + " at " + TextTools::toString(node->getFather()->getId()),
					  TextTools::toString(diff));
	}
	doNNI(node->getId());
	test = true;
//	nniLk_->topologyChangeTested(   TopologyChangeEvent ());
//	nniLk_->topologyChangeSuccessful(   TopologyChangeEvent ());
	if(verbose >= 1)
	  ApplicationTools::displayResult("   Current value", TextTools::toString(- getValue(), 10));
      }
    }
    elapsedTime_ += (ApplicationTools::getTime() - startingTime);
    startingTime = ApplicationTools::getTime();
  }
  while(test && (timeLimit_ == 0 || elapsedTime_ < timeLimit_));
}


/************************************************************************
 * Tells if the gene family is single copy (1 gene per sp)
 ************************************************************************/
bool DLGeneTreeLikelihood::isSingleCopy() {
  WHEREAMI( __FILE__ , __LINE__ );
  vector<string> names = TreeTemplateTools::getLeavesNames(*(geneTreeWithSpNames_->getRootNode() ) );
  vector<string> uniqueNames = VectorTools::unique(names);
  if (uniqueNames.size() == names.size() ) {
    return true;
  }
  return false;
}



/************************************************************************
 * Evaluate vector of gene trees, and returns the index of the best one
 ************************************************************************/

size_t DLGeneTreeLikelihood::findBestGeneTreeAmongCandidates(vector<Tree*> &trees,
                                                             TreeTemplate < Node > *bestTree,
                                                             double bestSequenceLogL,
                                                             double bestScenarioLk) {
    
    double candidateScenarioLk = 0.0;
    double candidateSequenceLk = 0.0;
    size_t bestI = 0;
//TreeTemplate < Node > * bestTree = 0;
    TreeTemplate < Node > * candidateTree = 0;
    
    for (size_t i = 0; i < trees.size() ; ++i) {
        candidateTree = dynamic_cast < TreeTemplate < Node > * > (trees[i]);
        candidateScenarioLk =  findMLReconciliationDR (spTree_, candidateTree,
                                                       seqSp_, spId_,
                                                       lossExpectedNumbers_, duplicationExpectedNumbers_,
                                                       tentativeMLindex_,
                                                       tentativeNum0Lineages_, tentativeNum1Lineages_,
                                                       tentativeNum2Lineages_, tentativeNodesToTryInNNISearch_, false); //false so that _tentativeNum*Lineages are not updated
        
        if (considerSequenceLikelihood_ )
        {
            levaluator_->setAlternativeTree(candidateTree);
            candidateSequenceLk = levaluator_->getAlternativeLogLikelihood();
        }
        std::cout << "Tree "<< i <<" : scenario LogLk: "<< candidateScenarioLk <<" ; sequence logLk : " << candidateSequenceLk <<std::endl;
        
        if (i == 0) {
            bestScenarioLk = candidateScenarioLk;
            bestSequenceLogL = candidateSequenceLk;
            if (bestTree ) {
                delete bestTree;
            }
            bestTree = dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();// candidateTree->clone(); 
        }
        else if (candidateScenarioLk + candidateSequenceLk > bestScenarioLk + bestSequenceLogL) {
            bestScenarioLk = candidateScenarioLk;
            bestSequenceLogL = candidateSequenceLk;
            levaluator_->acceptAlternativeTree();
            if (bestTree ) {
                delete bestTree;
            }
            bestTree =  dynamic_cast<const TreeTemplate<Node> *> (levaluator_->getTree())->clone();// candidateTree->clone();
            bestI = i;

	//Rooting bestTree properly:
	vector<Node*> rlkNodes = bestTree->getNodes();
	
	for (unsigned int j = 0 ; j < rlkNodes.size() ; j++) {
	  if (rlkNodes[j]->hasNodeProperty("outgroupNode")) {
	    if (bestTree->getRootNode() == rlkNodes[j]) {
	      if (j < rlkNodes.size()-1) 
	      {
		bestTree->rootAt(rlkNodes[rlkNodes.size()-1]);   
	      }
	      else {
		bestTree->rootAt(rlkNodes[rlkNodes.size()-2]);
	      }
	    };
	    bestTree->newOutGroup( rlkNodes[j] );
	    break;
	  }
	}
	    
        }
        delete candidateTree;
    }
    std::cout << "Best Tree: "<< bestI <<" : scenario LogLk: "<< bestScenarioLk <<" ; sequence logLk : " << bestSequenceLogL <<std::endl;

    return bestI;
}


