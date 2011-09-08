/*
 *  RearrangeGeneTreeDTL.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 29/06/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

/******************************************************************************/

#include "DTLGeneTreeLikelihood.h"

namespace mpi = boost::mpi;

//#include "mpi.h" 
//using namespace std;
using namespace bpp;



void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "RearrangeGeneTreeDTL parameter1_name=parameter1_value parameter2_name=parameter2_value"   ).endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "output.infos                      | file where to write site infos").endLine();
  (*ApplicationTools::message << "output.estimates                  | file where to write estimated parameter values").endLine();
  (*ApplicationTools::message << "species.tree.file                 | Path to a species tree" ).endLine();
  (*ApplicationTools::message << "init.gene.tree                     | user, bionj or phyml").endLine();
  (*ApplicationTools::message << "gene.tree.file                     | file containing a list of gene option files to analyse").endLine();
  (*ApplicationTools::message << "branchProbabilities.optimization  | average, branchwise, average_then_branchwise or no: how we optimize duplication, transfer and loss probabilities").endLine();
  (*ApplicationTools::message << "spr.limit                         | integer giving the breadth of SPR movements, in number of nodes. 0.1* number of nodes in the species tree might be OK.").endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of supplementary options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}





/**************************************************************************
 * This function optimizes a gene tree based on DTL and sequence likelihoods for 
 * a fixed species tree. It uses SPRs and NNIs.
 **************************************************************************/

double refineGeneTreeDTL (Species_tree * scoringTree, 
                          TreeTemplate<Node> *& geneTree, 
                          std::map<std::string, std::string > seqSp,
                          std::map<std::string, int > spID, 
                          pair <double, string> MLValueAndTree,
                          int & MLindex, 
                          string scoringMethod)
{
  
 /* TreeTemplate<Node> *tree = 0;
  TreeTemplate<Node> *bestTree = 0;
  TreeTemplate<Node> *currentTree = 0;
  currentTree = geneTree->clone();
  breadthFirstreNumber (*currentTree);//, duplicationExpectedNumbers, lossExpectedNumbers);
  int sprLimit = 10; //Arbitrary.
  std::vector <int> nodeIdsToRegraft; 
  double bestlogL;
  double logL;
  bool betterTree;
  int numIterationsWithoutImprovement = 0;
  double startingML = findMLReconciliationDR (scoringTree, 
                                              currentTree, 
                                              seqSp,
                                              spID,
                                              lossExpectedNumbers, 
                                              duplicationExpectedNumbers, 
                                              MLindex, 
                                              num0lineages, 
                                              num1lineages, 
                                              num2lineages, 
                                              nodesToTryInNNISearch, false);
  tree = currentTree->clone();
  bestTree = currentTree->clone();
  bestlogL = startingML;
  string strTree = "";
  
  while (numIterationsWithoutImprovement < geneTree->getNumberOfNodes()-1) {
    for (int nodeForSPR=geneTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) {
      betterTree = false;
      tree = currentTree->clone();
      buildVectorOfRegraftingNodesLimitedDistance(*tree, nodeForSPR, sprLimit, nodeIdsToRegraft);
      for (int i =0 ; i<nodeIdsToRegraft.size() ; i++) {
        if (tree) {
          delete tree;
        }
        tree = currentTree->clone();
        makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i], false);
        strTree = TreeTemplateTools::treeToParenthesis(*tree);
        logL = findMLReconciliationDR (spTree, 
                                       tree, 
                                       seqSp,
                                       spID,
                                       lossExpectedNumbers, 
                                       duplicationExpectedNumbers, 
                                       MLindex, 
                                       num0lineages, 
                                       num1lineages, 
                                       num2lineages, 
                                       nodesToTryInNNISearch, false);
        if (logL-0.01>bestlogL) {
          betterTree = true;
          bestlogL =logL;
          if (bestTree)
            delete bestTree;
          bestTree = tree->clone();  
          //std::cout << "Gene tree SPR: Better candidate tree likelihood : "<<bestlogL<< std::endl;
          //std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
        }
      }
      if (betterTree) {
        logL = bestlogL; 
        if (currentTree)
          delete currentTree;
        currentTree = bestTree->clone();
        breadthFirstreNumber (*currentTree);//, duplicationExpectedNumbers, lossExpectedNumbers);
                                            //std::cout <<"NEW BETTER TREE: \n"<< TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;
        numIterationsWithoutImprovement = 0;
      }
      else {
        logL = bestlogL; 
        if (currentTree)
          delete currentTree;
        currentTree = bestTree->clone(); 
        numIterationsWithoutImprovement++;
      }
    }
  }
  if (geneTree)
    delete geneTree;
  geneTree = bestTree->clone();
  if (tree) delete tree;
  if (bestTree) delete bestTree;
  if (currentTree) delete currentTree;  
  std::cout << "DTL initial likelihood: "<< startingML << "; Optimized DTL log likelihood "<< bestlogL <<" for family: " << file << std::endl; 

  return bestlogL;*/
}










/******************************************************************************/
// Program to rearrange a gene tree taking into account sequence likelihood and 
// DTL likelihood.
// Algorithm: we generate a lot of gene tree topologies, and compute the DTL scores
// for these. We compute the sequence lk for those that have higher DTL lks than 
// the current topology, select the best one, and start again from this new one.
//
/******************************************************************************/



int main(int args, char ** argv)
{
	if(args == 1)
    {
		help();
		exit(0);
    }
  try {
		ApplicationTools::startTimer();
    std::cout << "******************************************************************" << std::endl;
    std::cout << "*         DTL Gene Tree Rearrangement Program, version 1.0       *" << std::endl;
    std::cout << "* Authors: G. Szollosi, B. Boussau            Created 29/06/2011 *" << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;
    
    std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);
    
    DTLGeneTreeLikelihood DTLgtl = DTLGeneTreeLikelihood(params);
    DTLgtl.initialize();

    //Now we optimize the gene tree.
    DTLgtl.MLSearch();
    
    //Outputting the resulting gene tree
    DTLgtl.printGeneTree();
    
    std::cout << "RearrangeGeneTreeDTL's done. Bye." << std::endl;
    ApplicationTools::displayTime("Total execution time:");
	}
	catch(std::exception & e)
	{
  std::cout << e.what() << std::endl;
  exit(-1);
	}
	return (0);
}