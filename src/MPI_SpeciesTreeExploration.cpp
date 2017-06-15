/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a computer program whose purpose is to simultaneously build
gene and species trees when gene families have undergone duplications and
losses. It can analyze thousands of gene families in dozens of genomes
simultaneously, and was presented in an article in Genome Research. Trees and
parameters are estimated in the maximum likelihood framework, by maximizing
theprobability of alignments given the species tree, the gene trees and the
parameters of duplication and loss.

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
/* This file contains various functions useful for the search of the species tree*/

#include "Constants.h"
#include "MPI_SpeciesTreeExploration.h"




/************************************************************************
 * Procedure that makes NNIs and ReRootings by calling makeDeterministicNNIsAndRootChangesOnly.
************************************************************************/


void localOptimizationWithNNIsAndReRootings(const mpi::communicator& world,
                                            TreeTemplate<Node> *& tree,
                                            TreeTemplate<Node> *& bestTree,
                                            unsigned int &index,
                                            unsigned int &bestIndex,
                                            bool &stop,
                                            int timeLimit,
                                            double &logL,
                                            double &bestlogL,
                                            std::vector<int> &num0Lineages,
                                            std::vector<int> &num1Lineages,
                                            std::vector<int> &num2Lineages,
                                            std::vector<int> &bestNum0Lineages,
                                            std::vector<int> &bestNum1Lineages,
                                            std::vector<int> &bestNum2Lineages,
                                            std::vector< std::vector<int> > &allNum0Lineages,
                                            std::vector< std::vector<int> > &allNum1Lineages,
                                            std::vector< std::vector<int> > &allNum2Lineages,
                                            std::vector<double> &lossExpectedNumbers,
                                            std::vector<double> &duplicationExpectedNumbers,
                                            std::vector<unsigned int> &num12Lineages,
                                            std::vector<unsigned int> &num22Lineages,
                                            std::vector<unsigned int> &bestNum12Lineages,
                                            std::vector<unsigned int> &bestNum22Lineages,
                                            std::vector<double> &coalBls,
                                            std::string &reconciliationModel,
                                            bool rearrange,
                                            unsigned int &numIterationsWithoutImprovement,
                                            unsigned int server,
                                            size_t & nodeForNNI,
                                            size_t & nodeForRooting,
                                            map <string, double > treesToLogLk,
                                            std::string & branchExpectedNumbersOptimization,
                                            std::map < std::string, int> genomeMissing,
                                            std::vector < double >  &NNILks,
                                            std::vector<double> &rootLks, unsigned int currentStep,
											 const bool fixedOutgroupSpecies_,
											 const std::vector < std::string > outgroupSpecies_)
{
    std::vector <double> bestDupProba;
    std::vector<double> bestLossProba;
    std::vector <double> bestCoalBls;
    if (reconciliationModel == "DL") {
        breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
        bestDupProba=duplicationExpectedNumbers;
        bestLossProba=lossExpectedNumbers;
    }
    else if (reconciliationModel == "COAL") {
        breadthFirstreNumber (*tree, coalBls);
        bestCoalBls=coalBls;
    }

  std::string currentSpeciesTree;
  while (!stop)
    {
    //If we have already computed the value for the NNI or the rerooting
    //we are about to make, we don't make it!
	/*Seems buggy, does not do what it's supposed to do. Replaced below.
     *	if (fixedOutgroupSpecies_)
			nodeForRooting = tree->getNumberOfNodes();
		stop = checkChangeHasNotBeenDone(*tree, bestTree, nodeForNNI, nodeForRooting, NNILks, rootLks);
    */

    if (!stop)
      {
		  //Here we have implemented functions to not change the root if fixedOutgroupSpecies_==true (untested)
      makeDeterministicNNIsAndRootChangesOnly(*tree, nodeForNNI, nodeForRooting, fixedOutgroupSpecies_);
      double cachedLogLk = checkChangeHasNotBeenDone(*tree, treesToLogLk);
        if (cachedLogLk == 0.0) {
		 /* if (reconciliationModel == "DL") {
			  breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
		  }
		  else if (reconciliationModel == "COAL") {
			  breadthFirstreNumber (*tree, coalBls);
		  }*///TEST170413

		  breadthFirstreNumber (*tree);
		    //     std::cout << "HEHE 1 " << TreeTemplateTools::treeToParenthesis(*tree, true) <<std::endl;

		  //We set preliminary loss and duplication rates, correcting for genome coverage
		  computeDuplicationAndLossRatesForTheSpeciesTreeInitially(branchExpectedNumbersOptimization,
									    num0Lineages,
									    num1Lineages,
									    num2Lineages,
									    lossExpectedNumbers,
									    duplicationExpectedNumbers,
									    genomeMissing,
									    *tree);

      computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop,
                                                                         logL, num0Lineages,
                                                                         num1Lineages,num2Lineages,
                                                                         allNum0Lineages, allNum1Lineages,
                                                                         allNum2Lineages, lossExpectedNumbers,
                                                                         duplicationExpectedNumbers,
                                                                         num12Lineages, num22Lineages, coalBls,
                                                                         reconciliationModel, rearrange,
                                                                         server, branchExpectedNumbersOptimization,
                                                                         genomeMissing, *tree, bestlogL, currentStep);

      if ((nodeForNNI <3) /*&& (nodeForRooting == 4)*/) //A NNI has been done
        {
        if (logL < rootLks[nodeForRooting-1])
          {
          rootLks[nodeForRooting-1] = logL;
          }
        }
      else {
        if (logL < NNILks[nodeForNNI-1])
          {
          NNILks[nodeForNNI-1] = logL;
          }
      }
        //int branchId = bestTree->getNode(nodeForNNI-1)->getFather()->getId();
      /*  if (logL < NNILks[nodeForNNI])
          {
          NNILks[nodeForNNI] = logL;
          }
        }
      else //A rerooting has been done
        {
        if (logL < rootLks[nodeForRooting])
          {
          rootLks[nodeForRooting] = logL;
          }
        }*/

              treesToLogLk[TreeTemplateTools::treeToParenthesis( *tree, false )] = logL;

          if (logL+0.01<bestlogL)
          {
              std::cout << "\t\tNNIs or Root changes: Improvement: new total Likelihood value "<<logL<<" compared to the best log Likelihood : "<<bestlogL<< std::endl;
              numIterationsWithoutImprovement = 0;
              bestlogL =logL;
			 // breadthFirstreNumber (*tree);
              if (reconciliationModel == "DL") {
				 // breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
                  bestDupProba=duplicationExpectedNumbers;
                  bestLossProba=lossExpectedNumbers;
                  bestNum0Lineages = num0Lineages;
                  bestNum1Lineages = num1Lineages;
                  bestNum2Lineages = num2Lineages;
              }
              else if (reconciliationModel == "COAL") {
				 // breadthFirstreNumber (*tree, coalBls);
                  bestCoalBls=coalBls;
                  bestNum12Lineages = num12Lineages;
                  bestNum22Lineages = num22Lineages;
              }
			  if (bestTree)
              {
                  delete bestTree;
                  bestTree = 0;
              }
              bestTree = tree->clone();
              std::cout << "Improved species tree: "<<TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;
              bestIndex = index;
              for (unsigned int i = 0 ; i< NNILks.size() ; i++ )
              {
                  NNILks[i]=NumConstants::VERY_BIG();
                  rootLks[i]=NumConstants::VERY_BIG();
              }
              if (ApplicationTools::getTime() >= timeLimit)
              {
                  stop = true;
                  lastCommunicationsServerClient (world,
                                                  server,
                                                  stop,
                                                  bestIndex);
              }
          }
      else
        {
            numIterationsWithoutImprovement++;
            std::cout <<"\t\tNNIs or Root changes: Number of iterations without improvement: "<<numIterationsWithoutImprovement<< std::endl;
            if (tree) delete tree;
            tree = bestTree->clone();
            if (reconciliationModel == "DL") {
            duplicationExpectedNumbers = bestDupProba;
            lossExpectedNumbers = bestLossProba;
            }
            else if (reconciliationModel == "COAL") {
                coalBls = bestCoalBls;
            }
            if ( ( numIterationsWithoutImprovement > 2*tree->getNumberOfNodes() ) || (ApplicationTools::getTime() >= timeLimit) )
            {
                stop = true;
                if (ApplicationTools::getTime() >= timeLimit) {
                    lastCommunicationsServerClient (world,
                                                    server,
                                                    stop,
                                                    bestIndex);
                }
            }
        }
      }
      else { //The tree has already been found.
            std::cout << "This species tree has already been tried. "<<std::endl;
            numIterationsWithoutImprovement++;
            std::cout <<"\t\tNNIs or Root changes: Number of iterations without improvement: "<<numIterationsWithoutImprovement<< std::endl;
            if (tree) delete tree;
            tree = bestTree->clone();
            if (reconciliationModel == "DL") {
            duplicationExpectedNumbers = bestDupProba;
            lossExpectedNumbers = bestLossProba;
            }
            else if (reconciliationModel == "COAL") {
                coalBls = bestCoalBls;
            }
            if ( ( numIterationsWithoutImprovement > 2*tree->getNumberOfNodes() ) || (ApplicationTools::getTime() >= timeLimit) )
            {
                stop = true;
                if (ApplicationTools::getTime() >= timeLimit) {
                    lastCommunicationsServerClient (world,
                                                    server,
                                                    stop,
                                                    bestIndex);
                }
            }
      }
      }
    else //stop is true
      {
        stop = false;
        computeSpeciesTreeLikelihood(world, index, stop,
                                     logL, num0Lineages,
                                     num1Lineages,num2Lineages,
                                     allNum0Lineages, allNum1Lineages,
                                     allNum2Lineages, lossExpectedNumbers,
                                     duplicationExpectedNumbers,
                                     num12Lineages, num22Lineages,
                                     coalBls, reconciliationModel,
                                     rearrange, server, branchExpectedNumbersOptimization,
                                     genomeMissing, *bestTree, currentStep);
    /*    stop = true;
        broadcast(world, stop, server);
        broadcast(world, bestIndex, server);*/
      }
    }
}


/************************************************************************
 * Only optimizes duplication and loss rates and sends signals to end the program.
 ************************************************************************/
void optimizeOnlyDuplicationAndLossRates(const mpi::communicator& world,
                                         TreeTemplate<Node> *& tree,
                                         TreeTemplate<Node> *& bestTree,
                                         unsigned int &index,
                                         unsigned int &bestIndex,
                                         bool &stop,
                                         int timeLimit,
                                         double &logL,
                                         double &bestlogL,
                                         std::vector<int> &num0Lineages,
                                         std::vector<int> &num1Lineages,
                                         std::vector<int> &num2Lineages,
                                         std::vector<int> &bestNum0Lineages,
                                         std::vector<int> &bestNum1Lineages,
                                         std::vector<int> &bestNum2Lineages,
                                         std::vector< std::vector<int> > &allNum0Lineages,
                                         std::vector< std::vector<int> > &allNum1Lineages,
                                         std::vector< std::vector<int> > &allNum2Lineages,
                                         std::vector<double> &lossExpectedNumbers,
                                         std::vector<double> &duplicationExpectedNumbers,
                                         std::vector<unsigned int> &num12Lineages,
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange,
                                         unsigned int &numIterationsWithoutImprovement,
                                         unsigned int server,
                                         std::string & branchExpectedNumbersOptimization,
                                         std::map < std::string,
                                         int> genomeMissing, unsigned int currentStep) {
  computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop,
                                                                     logL, num0Lineages,
                                                                     num1Lineages,num2Lineages,
                                                                     allNum0Lineages, allNum1Lineages,
                                                                     allNum2Lineages, lossExpectedNumbers,
                                                                     duplicationExpectedNumbers,
                                                                     num12Lineages, num22Lineages,
                                                                     coalBls, reconciliationModel, rearrange,
                                                                     server, branchExpectedNumbersOptimization,
                                                                     genomeMissing, *tree, bestlogL, currentStep);
  bestIndex = index;
    if (bestTree) {
        delete bestTree;
        bestTree=0;
    }
  bestTree = tree->clone();
  bestNum0Lineages = num0Lineages;
  bestNum1Lineages = num1Lineages;
  bestNum2Lineages = num2Lineages;

    stop = true;
    lastCommunicationsServerClient (world,
                                    server,
                                    stop,
                                    bestIndex);
}



/************************************************************************
 * Tries all reRootings of the species tree, and executes the one with the highest likelihood.
 Here we do not use branchwise rates of duplications and losses,
 and we do not optimise these rates often.
 However, there are lots of useless identical std::vector copies...
 ************************************************************************/
void fastTryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world,
                                                TreeTemplate<Node> *& currentTree,
                                                TreeTemplate<Node> *& bestTree,
                                                unsigned int &index, unsigned int &bestIndex,
                                                bool stop, int timeLimit,
                                                double &logL, double &bestlogL,
                                                std::vector<int> &num0Lineages, std::vector<int> &num1Lineages,
                                                std::vector<int> &num2Lineages,
                                                std::vector<int> &bestNum0Lineages,
                                                std::vector<int> &bestNum1Lineages,
                                                std::vector<int> &bestNum2Lineages,
                                                std::vector< std::vector<int> > &allNum0Lineages,
                                                std::vector< std::vector<int> > &allNum1Lineages,
                                                std::vector< std::vector<int> > &allNum2Lineages,
                                                std::vector<double> &lossExpectedNumbers,
                                                std::vector<double> &duplicationExpectedNumbers,
                                                std::vector<unsigned int> &num12Lineages,
                                                std::vector<unsigned int> &num22Lineages,
                                                std::vector<unsigned int> &bestNum12Lineages,
                                                std::vector<unsigned int> &bestNum22Lineages,
                                                std::vector<double> &coalBls,
                                                std::string &reconciliationModel,
                                                bool rearrange, unsigned int &numIterationsWithoutImprovement,
                                                unsigned int server, std::string &branchExpectedNumbersOptimization,
                                                std::map < std::string, int> genomeMissing,
                                                bool optimizeRates, unsigned int currentStep) {
    std::vector <double> backupDupProba;
    std::vector <double> backupLossProba;
    std::vector <double> bestDupProba;
    std::vector <double> bestLossProba;
    std::vector <double> backupCoalBls;
    std::vector <double> bestCoalBls;
    if (reconciliationModel == "DL") {
        breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers);
        backupDupProba=duplicationExpectedNumbers;
        backupLossProba=lossExpectedNumbers;
        bestDupProba=duplicationExpectedNumbers;
        bestLossProba=lossExpectedNumbers;
    }
    else if (reconciliationModel == "COAL") {
        breadthFirstreNumber (*currentTree, coalBls);
        backupCoalBls=coalBls;
        bestCoalBls=coalBls;
    }
	if (bestTree) delete bestTree;
	bestTree = currentTree->clone();
  TreeTemplate<Node> *tree;
  tree = currentTree->clone();
  bool betterTree = false;
  std::vector <int> nodeIds = tree->getNodesId();
  for (unsigned int i =0 ; i<nodeIds.size() ; i++) {
    if ((nodeIds[i]!=0)&&(nodeIds[i]!=1)&&(nodeIds[i]!=2)) { //We do not want to try the root we're already at
      if (i!=0) {
        if (tree) delete tree;
        tree = currentTree->clone();
      }
      changeRoot(*tree, nodeIds[i]);

		breadthFirstreNumber (*tree);
		 //      std::cout << "HEHE 2" << TreeTemplateTools::treeToParenthesis(*tree, true) <<std::endl;

		//We set preliminary loss and duplication rates, correcting for genome coverage
		computeDuplicationAndLossRatesForTheSpeciesTreeInitially(branchExpectedNumbersOptimization,
																 num0Lineages,
																 num1Lineages,
																 num2Lineages,
																 lossExpectedNumbers,
																 duplicationExpectedNumbers,
																 genomeMissing,
																 *tree);

		/*
		if (reconciliationModel == "DL") {
			breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
		}
		else if (reconciliationModel == "COAL") {
			breadthFirstreNumber (*tree, coalBls);
		}*///TEST170413
   /*   if (optimizeRates)
        {*/
            computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL,
                                                                               num0Lineages,
                                                                               num1Lineages,num2Lineages,
                                                                               allNum0Lineages, allNum1Lineages,
                                                                               allNum2Lineages, lossExpectedNumbers,
                                                                               duplicationExpectedNumbers,
                                                                               num12Lineages, num22Lineages,
                                                                               coalBls, reconciliationModel,
                                                                               rearrange, server,
                                                                               branchExpectedNumbersOptimization,
                                                                               genomeMissing, *tree,
                                                                               bestlogL, currentStep);
/*        }
      else
        {
          computeSpeciesTreeLikelihood(world, index, stop, logL, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossExpectedNumbers, duplicationExpectedNumbers, rearrange, server, branchExpectedNumbersOptimization, genomeMissing, *tree, currentStep);
      }*/
      if (logL+0.01<bestlogL) {
        numIterationsWithoutImprovement = 0;
        betterTree = true;
        bestlogL =logL;
          if (bestTree) {
              delete bestTree;
              bestTree=0;
          }
        bestTree = tree->clone();
          if (reconciliationModel == "DL") {
              bestDupProba=duplicationExpectedNumbers;
              bestLossProba=lossExpectedNumbers;
              bestNum0Lineages = num0Lineages;
              bestNum1Lineages = num1Lineages;
              bestNum2Lineages = num2Lineages;
          }
          else if (reconciliationModel == "COAL") {
              bestCoalBls=coalBls;
              bestNum12Lineages = num12Lineages;
              bestNum22Lineages = num22Lineages;
          }

        bestIndex = index;
       std::cout <<"ReRooting: Improvement! : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
       std::cout << "Better candidate tree likelihood : "<<bestlogL<< std::endl;
       std::cout << TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;
		/*  //TEMP PRINTING
		  //For loss rates
		  for (unsigned int i =0; i<num0Lineages.size() ; i++ )
		  {
			  tree->getNode(i)->setBranchProperty("LOSSES", Number<double>(lossExpectedNumbers[i]));
			  if (tree->getNode(i)->hasFather())
			  {
				  tree->getNode(i)->setDistanceToFather(lossExpectedNumbers[i]);
			  }
		  }
		  std::cout <<"\n\n\t\t Better Species Tree found after rerooting, with Losses: "<<std::endl;
		  //			std::cout << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "LOSSES")<<std::endl;
		  std::cout << TreeTemplateTools::treeToParenthesis(*tree, false)<<std::endl;
		 */


      }
      else {
        numIterationsWithoutImprovement++;
		 duplicationExpectedNumbers = backupDupProba ;
		 lossExpectedNumbers = backupLossProba ;
        std::cout <<"ReRooting: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
      }
    }
    if (ApplicationTools::getTime() >= timeLimit)
      {
      stop = true;
          lastCommunicationsServerClient (world,
                                          server,
                                          stop,
                                          bestIndex);
      break;
      }
  }

  if (betterTree) {
    logL = bestlogL;
      if (reconciliationModel == "DL") {
          duplicationExpectedNumbers = bestDupProba;
          lossExpectedNumbers = bestLossProba;
          num0Lineages = bestNum0Lineages;
          num1Lineages = bestNum1Lineages;
          num2Lineages = bestNum2Lineages;
      }
      else if (reconciliationModel == "COAL") {
          coalBls = bestCoalBls;
          num12Lineages = bestNum12Lineages;
          num22Lineages = bestNum22Lineages;
      }
    if (currentTree) delete currentTree;
    currentTree = bestTree->clone();
   std::cout << "\t\tServer: tryAllPossibleReRootingsAndMakeBestOne: new total Likelihood value "<<logL<< std::endl;
   std::cout << TreeTemplateTools::treeToParenthesis(*currentTree, true)<< std::endl;

	  /*
	  //TEMP PRINTING
	  //For loss rates
	  for (unsigned int i =0; i<num0Lineages.size() ; i++ )
	  {
		  currentTree->getNode(i)->setBranchProperty("LOSSES", Number<double>(lossExpectedNumbers[i]));
		  if (currentTree->getNode(i)->hasFather())
		  {
			  currentTree->getNode(i)->setDistanceToFather(lossExpectedNumbers[i]);
		  }
	  }
	  std::cout <<"\n\n\t\t Better Species Tree found after rerooting, with Losses: "<<std::endl;
	  //			std::cout << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "LOSSES")<<std::endl;
	  std::cout << TreeTemplateTools::treeToParenthesis(*currentTree, false)<<std::endl;
	   */

    std::string currentSpeciesTree = TreeTemplateTools::treeToParenthesis(*currentTree, true);
//TEST  WEIRD
	  //breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers); //TEST170413
      //breadthFirstreNumber (*currentTree);
    if (ApplicationTools::getTime() < timeLimit)
      {
   /*TEST 08022012
          broadcastsAllInformation(world, server, stop, rearrange, lossExpectedNumbers, duplicationExpectedNumbers, currentSpeciesTree, currentStep);
      //COMPUTATION IN CLIENTS
      index++;
      bestIndex = index;
      std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
          gathersInformationFromClients (world,
                                         server,
                                         server,
                                         logL,
                                         num0Lineages,
                                         num1Lineages,
                                         num2Lineages,
                                         allNum0Lineages,
                                         allNum1Lineages,
                                         allNum2Lineages);*/

      /*logL = 0.0;
      std::vector<double> logLs;
      resetVector(num0Lineages);
      resetVector(num1Lineages);
      resetVector(num2Lineages);
      gather(world, logL, logLs, server);
      logL =  VectorTools::sum(logLs);
      gather(world, num0Lineages, allNum0Lineages, server);
      gather(world, num1Lineages, allNum1Lineages, server);
      gather(world, num2Lineages, allNum2Lineages, server);*/
      }
  }
  else {
   std::cout<< "No improvement in fastTryAllPossibleReRootingsAndMakeBestOne"<< std::endl;
    logL = bestlogL;
      if (reconciliationModel == "DL") {
          duplicationExpectedNumbers = bestDupProba;
          lossExpectedNumbers = bestLossProba;
          num0Lineages = bestNum0Lineages;
          num1Lineages = bestNum1Lineages;
          num2Lineages = bestNum2Lineages;
      }
      else if (reconciliationModel == "COAL") {
          coalBls = bestCoalBls;
          num12Lineages = bestNum12Lineages;
          num22Lineages = bestNum22Lineages;
      }
      if (currentTree) delete currentTree;
    currentTree = bestTree->clone();
  }

  if (tree) delete tree;
}


/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR,
 * and executes the ones with the highest likelihood.
 ************************************************************************/
void fastTryAllPossibleSPRs(const mpi::communicator& world, TreeTemplate<Node> *& currentTree,
                            TreeTemplate<Node> *& bestTree, unsigned int &index, unsigned int &bestIndex,
                            bool stop, int timeLimit, double &logL, double &bestlogL,
                            std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages,
                            std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages,
                            std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages,
                            std::vector< std::vector<int> > &allNum2Lineages,
                            std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers,
                            std::vector<unsigned int> &num12Lineages,
                            std::vector<unsigned int> &num22Lineages,
                            std::vector<unsigned int> &bestNum12Lineages,
                            std::vector<unsigned int> &bestNum22Lineages,
                            std::vector<double> &coalBls,
                            std::string &reconciliationModel,
                            bool rearrange, unsigned int &numIterationsWithoutImprovement, unsigned int server,
                            std::string &branchExpectedNumbersOptimization, std::map < std::string, int> genomeMissing,
                            int sprLimit, bool optimizeRates, unsigned int currentStep,
							const bool fixedOutgroupSpecies_, const std::vector < std::string > outgroupSpecies_) {

    std::vector <double> bestDupProba;
    std::vector <double> bestLossProba;
    std::vector <double> backupDupProba;
    std::vector <double> backupLossProba;

    std::vector <double> bestCoalBls;
    if (reconciliationModel == "DL") {
        breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers);
        bestDupProba=duplicationExpectedNumbers;
        bestLossProba=lossExpectedNumbers;
		backupDupProba=duplicationExpectedNumbers;
        backupLossProba=lossExpectedNumbers;

    }
    else if (reconciliationModel == "COAL") {
        breadthFirstreNumber (*currentTree, coalBls);
        bestCoalBls=coalBls;
    }
    std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> *tree = 0;
	//Main loop
  for (unsigned int nodeForSPR=currentTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) {
    buildVectorOfRegraftingNodesLimitedDistance(*currentTree, nodeForSPR, sprLimit, nodeIdsToRegraft);
    betterTree = false;
    for (size_t i =0 ; i<nodeIdsToRegraft.size() ; i++) {
      if (tree) {
        delete tree;
        tree = 0;
      }
      tree = currentTree->clone();

      makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i]);
		if (!fixedOutgroupSpecies_ || (fixedOutgroupSpecies_ && isTreeRootedWithOutgroup (*tree, outgroupSpecies_) ) )
		{
       breadthFirstreNumber (*tree); //TEST170413
      // std::cout << "HEHE 0 "<< TreeTemplateTools::treeToParenthesis(*tree, true) <<std::endl;
			//We set preliminary loss and duplication rates, correcting for genome coverage
			computeDuplicationAndLossRatesForTheSpeciesTreeInitially(branchExpectedNumbersOptimization,
										num0Lineages,
										num1Lineages,
										num2Lineages,
										lossExpectedNumbers,
										duplicationExpectedNumbers,
										genomeMissing,
										*tree);

		/*	if (reconciliationModel == "DL") {
				breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
			}
			else if (reconciliationModel == "COAL") {
				breadthFirstreNumber (*tree, coalBls);
			}*/
   /*   if (optimizeRates)
        {*/
          computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index,
                                                                             stop, logL,
                                                                             num0Lineages, num1Lineages,num2Lineages,
                                                                             allNum0Lineages, allNum1Lineages, allNum2Lineages,
                                                                             lossExpectedNumbers, duplicationExpectedNumbers,
                                                                             num12Lineages, num22Lineages,
                                                                             coalBls, reconciliationModel,
                                                                             rearrange, server, branchExpectedNumbersOptimization,
                                                                             genomeMissing, *tree, bestlogL, currentStep);

   /*     }
      else
        {
          computeSpeciesTreeLikelihood(world, index,
                                       stop, logL,
                                       num0Lineages, num1Lineages,num2Lineages,
                                       allNum0Lineages, allNum1Lineages, allNum2Lineages,
                                       lossExpectedNumbers, duplicationExpectedNumbers,
                                       rearrange, server, branchExpectedNumbersOptimization,
                                       genomeMissing, *tree, currentStep);

        }*/
		}
		else {
			//The tree to test is not properly rooted, we do not compute its likelihood
			logL = bestlogL + 2000 ;
		}
      if (logL+0.01<bestlogL) {
        betterTree = true;
        bestlogL =logL;
          if (bestTree) {
              delete bestTree;
              bestTree=0;
          }
        bestTree = tree->clone();
          if (reconciliationModel == "DL" ) {
              bestDupProba = duplicationExpectedNumbers;
              bestLossProba = lossExpectedNumbers;
              bestNum0Lineages = num0Lineages;
              bestNum1Lineages = num1Lineages;
              bestNum2Lineages = num2Lineages;
              /* VectorTools::print(duplicationExpectedNumbers);
               VectorTools::print(lossExpectedNumbers);*/
          }
          else if (reconciliationModel == "COAL" ) {
              bestCoalBls = coalBls;
              bestNum12Lineages = num12Lineages;
              bestNum22Lineages = num22Lineages;
          }
        bestIndex = index;

       std::cout << "SPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
       std::cout << TreeTemplateTools::treeToParenthesis(*bestTree, true)<< std::endl;
		  /*
		  //TEMP PRINTING
		  //For loss rates
		  for (unsigned int i =0; i<num0Lineages.size() ; i++ )
		  {
			  bestTree->getNode(i)->setBranchProperty("LOSSES", Number<double>(lossExpectedNumbers[i]));
			  if (bestTree->getNode(i)->hasFather())
			  {
				  bestTree->getNode(i)->setDistanceToFather(lossExpectedNumbers[i]);
			  }
		  }
		  std::cout <<"\n\n\t\t Better Species Tree found, with Losses: "<<std::endl;
		  //			std::cout << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "LOSSES")<<std::endl;
		  std::cout << TreeTemplateTools::treeToParenthesis(*bestTree, false)<<std::endl;
		  */
      }
	  else { // No improvement with this SPR
		  duplicationExpectedNumbers = backupDupProba ;
		  lossExpectedNumbers = backupLossProba ;
	  }
      if (ApplicationTools::getTime() >= timeLimit)
        {
        stop = true;
            lastCommunicationsServerClient (world,
                                            server,
                                            stop,
                                            bestIndex);
        break;
        }
    }
    if (betterTree) {
      logL = bestlogL;
        numIterationsWithoutImprovement = 0;
        if ( reconciliationModel == "DL" ) {
            duplicationExpectedNumbers = bestDupProba;
            lossExpectedNumbers = bestLossProba;
            num0Lineages = bestNum0Lineages;
            num1Lineages = bestNum1Lineages;
            num2Lineages = bestNum2Lineages;
        }
        else if ( reconciliationModel == "COAL" ) {
            coalBls = bestCoalBls;
            num12Lineages = bestNum12Lineages;
            num22Lineages = bestNum22Lineages;
        }
     // deleteTreeProperties(*currentTree);
      if (currentTree) delete currentTree;
      currentTree = bestTree->clone();

//      breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers); //TEST
     std::cout <<"SPRs: Improvement! : "<<numIterationsWithoutImprovement<< std::endl;
     std::cout << "\t\tServer: SPRs: new total Likelihood value "<<logL<< std::endl;

      if (ApplicationTools::getTime() < timeLimit)
        {


        //Send the new std::vectors to compute the new likelihood of the best tree
        std::string currentSpeciesTree = TreeTemplateTools::treeToParenthesis(*currentTree, true);
        broadcastsAllInformation(world, server, stop,
                                 rearrange,
                                 lossExpectedNumbers, duplicationExpectedNumbers,
                                 coalBls,
                                 currentSpeciesTree, currentStep,
                                 reconciliationModel);
        //COMPUTATION IN CLIENTS
        index++;
        bestIndex = index;
        std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
            gathersInformationFromClients (world,
                                           server,
                                           server,
                                           logL,
                                           num0Lineages,
                                           num1Lineages,
                                           num2Lineages,
                                           allNum0Lineages,
                                           allNum1Lineages,
                                           allNum2Lineages,
                                           num12Lineages,
                                           num22Lineages,
                                           reconciliationModel);
			ApplicationTools::displayTime("Execution time so far:");
        }
      else
        {
        stop = true;
            lastCommunicationsServerClient (world,
                                            server,
                                            stop,
                                            bestIndex);
        break;
        }
    }
    else {
        logL = bestlogL;
        if ( reconciliationModel == "DL" ) {
            duplicationExpectedNumbers = bestDupProba;
            lossExpectedNumbers = bestLossProba;
            num0Lineages = bestNum0Lineages;
            num1Lineages = bestNum1Lineages;
            num2Lineages = bestNum2Lineages;
        }
        else if ( reconciliationModel == "COAL" ) {
            coalBls = bestCoalBls;
            num12Lineages = bestNum12Lineages;
            num22Lineages = bestNum22Lineages;
        }
        if (currentTree) delete currentTree;
        currentTree = bestTree->clone();
        numIterationsWithoutImprovement++;
        std::cout <<"SPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
    }
    if (tree) {
      delete tree;
      tree = 0;
    }
    if (ApplicationTools::getTime() >= timeLimit)
      {
      stop = true;
          lastCommunicationsServerClient (world,
                                          server,
                                          stop,
                                          bestIndex);
      break;
      }
  }
}



/************************************************************************
 * Tries all SPRs and all rerootings with average rates of duplication and loss,
 * not branchwise rates. Only does SPRs at a given distance.
 ************************************************************************/
void fastTryAllPossibleSPRsAndReRootings(const mpi::communicator& world,
                                         TreeTemplate<Node> *& currentTree,
                                         TreeTemplate<Node> *& bestTree,
                                         unsigned int &index, unsigned int &bestIndex,
                                         bool stop,
                                         int timeLimit,
                                         double &logL,
                                         double &bestlogL,
                                         std::vector<int> &num0Lineages,
                                         std::vector<int> &num1Lineages,
                                         std::vector<int> &num2Lineages,
                                         std::vector<int> &bestNum0Lineages,
                                         std::vector<int> &bestNum1Lineages,
                                         std::vector<int> &bestNum2Lineages,
                                         std::vector< std::vector<int> > &allNum0Lineages,
                                         std::vector< std::vector<int> > &allNum1Lineages,
                                         std::vector< std::vector<int> > &allNum2Lineages,
                                         std::vector<double> &lossExpectedNumbers,
                                         std::vector<double> &duplicationExpectedNumbers,
                                         std::vector<unsigned int> &num12Lineages,
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<unsigned int> &bestNum12Lineages,
                                         std::vector<unsigned int> &bestNum22Lineages,
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange,
                                         unsigned int &numIterationsWithoutImprovement,
                                         unsigned int server,
                                         std::string &branchExpectedNumbersOptimization,
                                         std::map < std::string, int> genomeMissing,
                                         int sprLimit,
                                         bool optimizeRates,
                                         unsigned int currentStep,
										 const bool fixedOutgroupSpecies_,
										 const std::vector < std::string > outgroupSpecies_) {
  if (optimizeRates)
    {
    std::cout <<"Making SPRs and Rerootings and optimizing duplication and loss rates."<< std::endl;
    }
  else
    {
    std::cout <<"Making SPRs and Rerootings but NOT optimizing duplication and loss rates."<< std::endl;
    }
  numIterationsWithoutImprovement=0;

  while ( numIterationsWithoutImprovement < 2*currentTree->getNumberOfNodes() ) {
    fastTryAllPossibleSPRs(world, currentTree, bestTree,
                           index, bestIndex,
                           stop, timeLimit,
                           logL, bestlogL,
                           num0Lineages, num1Lineages, num2Lineages,
                           bestNum0Lineages, bestNum1Lineages, bestNum2Lineages,
                           allNum0Lineages, allNum1Lineages, allNum2Lineages,
                           lossExpectedNumbers, duplicationExpectedNumbers,
                           num12Lineages, num22Lineages,
                           bestNum12Lineages, bestNum22Lineages,
                           coalBls, reconciliationModel,
                           rearrange, numIterationsWithoutImprovement,
                           server, branchExpectedNumbersOptimization, genomeMissing,
                           sprLimit, optimizeRates, currentStep,
						   fixedOutgroupSpecies_, outgroupSpecies_);

    if (ApplicationTools::getTime() >= timeLimit)
      {
        stop = true;
          lastCommunicationsServerClient (world,
                                          server,
                                          stop,
                                          bestIndex);
        break;
      }
    if  ( numIterationsWithoutImprovement > 2*currentTree->getNumberOfNodes() )
      {
        break;
      }
	  if ( ! fixedOutgroupSpecies_) {
		  fastTryAllPossibleReRootingsAndMakeBestOne(world, currentTree, bestTree,
													 index, bestIndex,
													 stop, timeLimit,
													 logL, bestlogL,
													 num0Lineages, num1Lineages, num2Lineages,
													 bestNum0Lineages, bestNum1Lineages, bestNum2Lineages,
													 allNum0Lineages, allNum1Lineages, allNum2Lineages,
													 lossExpectedNumbers, duplicationExpectedNumbers,
													 num12Lineages, num22Lineages,
													 bestNum12Lineages, bestNum12Lineages,
													 coalBls, reconciliationModel,
													 rearrange, numIterationsWithoutImprovement,
													 server, branchExpectedNumbersOptimization,
													 genomeMissing, optimizeRates, currentStep);
	  }
	  else {
		  numIterationsWithoutImprovement += 2*bestTree->getNumberOfLeaves() - 2 ;
	  }
   /* std::cout << "after fastTryAllPossibleReRootingsAndMakeBestOne :currentTree: "<< std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis(*currentTree, true)<< std::endl;*/
    if (ApplicationTools::getTime() >= timeLimit)
      {
        stop = true;
          lastCommunicationsServerClient (world,
                                          server,
                                          stop,
                                          bestIndex);
        break;
      }

  }

}


/************************************************************************
 * Computes the likelihood of a species tree, only once.
 ************************************************************************/
std::string computeSpeciesTreeLikelihood(const mpi::communicator& world,
                                         unsigned int &index,
                                         bool stop,
                                         double &logL,
                                         std::vector<int> &num0Lineages,
                                         std::vector<int> &num1Lineages,
                                         std::vector<int> &num2Lineages,
                                         std::vector< std::vector<int> > &allNum0Lineages,
                                         std::vector< std::vector<int> > &allNum1Lineages,
                                         std::vector< std::vector<int> > &allNum2Lineages,
                                         std::vector<double> &lossExpectedNumbers,
                                         std::vector<double> &duplicationExpectedNumbers,
                                         std::vector<unsigned int> &num12Lineages,
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange,
                                         unsigned int server,
                                         std::string &branchExpectedNumbersOptimization,
                                         std::map < std::string, int> genomeMissing,
                                         TreeTemplate<Node> &tree, unsigned int currentStep)
{
 //TEST breadthFirstreNumber (tree, duplicationExpectedNumbers, lossExpectedNumbers);
  std::string currentSpeciesTree = TreeTemplateTools::treeToParenthesis(tree, true);
  computeSpeciesTreeLikelihoodWithGivenStringSpeciesTree(world,index,
                                                         stop, logL,
                                                         num0Lineages, num1Lineages,
                                                         num2Lineages, allNum0Lineages,
                                                         allNum1Lineages, allNum2Lineages,
                                                         lossExpectedNumbers, duplicationExpectedNumbers,
                                                         num12Lineages, num22Lineages,
                                                         coalBls, reconciliationModel,
                                                         rearrange, server,
                                                         branchExpectedNumbersOptimization, genomeMissing,
                                                         tree, currentSpeciesTree, false, currentStep);
  return currentSpeciesTree;
}





/************************************************************************
 * Computes the likelihood of a species tree, only once.
 ************************************************************************/
std::string computeSpeciesTreeLikelihoodWithGivenStringSpeciesTree(const mpi::communicator& world,
                                                                   unsigned int &index,
                                                                   bool stop,
                                                                   double &logL,
                                                                   std::vector<int> &num0Lineages,
                                                                   std::vector<int> &num1Lineages,
                                                                   std::vector<int> &num2Lineages,
                                                                   std::vector< std::vector<int> > &allNum0Lineages,
                                                                   std::vector< std::vector<int> > &allNum1Lineages,
                                                                   std::vector< std::vector<int> > &allNum2Lineages,
                                                                   std::vector<double> &lossExpectedNumbers,
                                                                   std::vector<double> &duplicationExpectedNumbers,
                                                                   std::vector<unsigned int> &num12Lineages,
                                                                   std::vector<unsigned int> &num22Lineages,
                                                                   std::vector<double> &coalBls,
                                                                   std::string &reconciliationModel,
                                                                   bool rearrange,
                                                                   unsigned int server,
                                                                   std::string &branchExpectedNumbersOptimization,
                                                                   std::map < std::string, int> genomeMissing,
                                                                   TreeTemplate<Node> &tree,
                                                                   std::string currentSpeciesTree,
                                                                   bool firstTime,
                                                                   unsigned int currentStep)
{
 /* if (!firstTime)
    {*/
    broadcastsAllInformation(world, server, stop,
                             rearrange, lossExpectedNumbers,
                             duplicationExpectedNumbers,
                             coalBls,
                             currentSpeciesTree, currentStep,
                             reconciliationModel);

    index++;
 //   }

    gathersInformationFromClients (world,
                                   server,
                                   server,
                                   logL,
                                   num0Lineages,
                                   num1Lineages,
                                   num2Lineages,
                                   allNum0Lineages,
                                   allNum1Lineages,
                                   allNum2Lineages,
                                   num12Lineages,
                                   num22Lineages,
                                   reconciliationModel);

  /*
  //COMPUTATION IN CLIENTS
  logL = 0.0;
  std::vector<double> logLs;
  resetVector(num0Lineages);
  resetVector(num1Lineages);
  resetVector(num2Lineages);
  gather(world, logL, logLs, server);
  logL =  VectorTools::sum(logLs);
  gather(world, num0Lineages, allNum0Lineages, server);
  gather(world, num1Lineages, allNum1Lineages, server);
  gather(world, num2Lineages, allNum2Lineages, server);
  int temp = allNum0Lineages.size();
  for (int k =0; k<temp ; k++ ) {
    num0Lineages= num0Lineages+allNum0Lineages[k];
    num1Lineages= num1Lineages+allNum1Lineages[k];
    num2Lineages= num2Lineages+allNum2Lineages[k];
  }  */
  return currentSpeciesTree;
}




/************************************************************************
 * Computes the likelihood of a species tree, optimizes the duplication and loss
 * rates, and updates the likelihood of the species tree. This way the likelihood
 * of this species tree is computed with adequate rates.
 ************************************************************************/
void computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(const mpi::communicator& world,
                                                                        unsigned int &index, bool stop, double &logL,
                                                                        std::vector<int> &num0Lineages,
                                                                        std::vector<int> &num1Lineages,
                                                                        std::vector<int> &num2Lineages,
                                                                        std::vector< std::vector<int> > &allNum0Lineages,
                                                                        std::vector< std::vector<int> > &allNum1Lineages,
                                                                        std::vector< std::vector<int> > &allNum2Lineages,
                                                                        std::vector<double> &lossExpectedNumbers,
                                                                        std::vector<double> &duplicationExpectedNumbers,
                                                                        std::vector<unsigned int> &num12Lineages,
                                                                        std::vector<unsigned int> &num22Lineages,
                                                                        std::vector<double> &coalBls,
                                                                        std::string &reconciliationModel,
                                                                        bool rearrange, unsigned int server,
                                                                        std::string &branchExpectedNumbersOptimization,
                                                                        std::map < std::string, int> genomeMissing,
                                                                        TreeTemplate<Node> &tree, double & bestlogL,
                                                                        unsigned int currentStep)
{
	/*if (branchExpectedNumbersOptimization != "no")
	 {
	 //we set the expected numbers to 0.001 uniformly: we start afresh for each new species tree topology.
	 for (unsigned int i = 0 ; i < lossExpectedNumbers.size() ; i++ )
	 {
	 lossExpectedNumbers[i] = 0.001;
	 duplicationExpectedNumbers[i] = 0.001;
	 }
	 }*/

	std::string currentSpeciesTree = computeSpeciesTreeLikelihood(world, index,
																  stop, logL,
																  num0Lineages,
																  num1Lineages,
																  num2Lineages,
																  allNum0Lineages,
																  allNum1Lineages,
																  allNum2Lineages,
																  lossExpectedNumbers,
																  duplicationExpectedNumbers,
																  num12Lineages, num22Lineages,
																  coalBls, reconciliationModel,
																  rearrange, server,
																  branchExpectedNumbersOptimization,
																  genomeMissing, tree, currentStep);

	double currentlogL = -UNLIKELY;
	int i=1;
	if (branchExpectedNumbersOptimization != "no")
    {
		//   std::cout << "Species tree LogLikelihood after the first round: "<< - logL<<std::endl;
		//Then we update duplication and loss rates based on the results of this first
		//computation, until the likelihood stabilizes (roughly)
        std::cout<<";\t\tLogLk value for the species before optimizing DL parameters: "<< - logL<<std::endl;
        while ((i<=1)&&(currentlogL-logL>logL-bestlogL))
            //      while (logL - bestlogL > 0.1)
        {
            currentlogL = logL;
            i++;
            if (reconciliationModel == "DL") {
				computeDuplicationAndLossRatesForTheSpeciesTree (branchExpectedNumbersOptimization,
																 num0Lineages, num1Lineages,
																 num2Lineages, lossExpectedNumbers,
																 duplicationExpectedNumbers,
																 genomeMissing, tree);
            }
            else if (reconciliationModel == "COAL") {
				if (currentStep == 0) {
					std::string temp = "no";
					computeCoalBls (temp,
									num12Lineages,
									num22Lineages,
									coalBls) ;
				}
				else {
					computeCoalBls (branchExpectedNumbersOptimization,
									num12Lineages,
									num22Lineages,
									coalBls) ;
				}
            }
            computeSpeciesTreeLikelihoodWithGivenStringSpeciesTree(world,index,
                                                                   stop, logL,
                                                                   num0Lineages,
                                                                   num1Lineages,
                                                                   num2Lineages,
                                                                   allNum0Lineages,
                                                                   allNum1Lineages,
                                                                   allNum2Lineages,
                                                                   lossExpectedNumbers,
                                                                   duplicationExpectedNumbers,
                                                                   num12Lineages, num22Lineages,
                                                                   coalBls, reconciliationModel,
                                                                   false, server,
                                                                   branchExpectedNumbersOptimization,
                                                                   genomeMissing, tree,
                                                                   currentSpeciesTree, false, currentStep);
        }
        std::cout<<";\t\tLogLk value for the species after optimizing DL parameters: "<< - logL<<std::endl;
    }
	ApplicationTools::displayTime("Execution time so far:");

	// std::cout << i<< " iterations of likelihood computation for optimizing duplication and loss rates have been done."<< std::endl;
}





/************************************************************************
 * Broadcasts all necessary information.
 ************************************************************************/
void broadcastsAllInformation(const mpi::communicator& world, unsigned int server,
                              bool stop, bool &rearrange,
                              std::vector<double> &lossExpectedNumbers,
                              std::vector<double> &duplicationExpectedNumbers,
                              std::vector<double> &coalBls,
                              std::string & currentSpeciesTree,
                              unsigned int &currentStep,
                              std::string &reconciliationModel) {
  //  MPI_Barrier(world);

    broadcast(world, stop, server);

    broadcastsAllInformationButStop(world,server, rearrange,
                                    lossExpectedNumbers,
                                    duplicationExpectedNumbers,
                                    coalBls, currentSpeciesTree,
                                    currentStep, reconciliationModel);
}

void broadcastsAllInformationButStop(const mpi::communicator& world, unsigned int server,
                                     bool &rearrange,
                                     std::vector<double> &lossExpectedNumbers,
                                     std::vector<double> &duplicationExpectedNumbers,
                                     std::vector<double> &coalBls,
                                     std::string &currentSpeciesTree,
                                     unsigned int &currentStep,
                                     std::string &reconciliationModel) {

    broadcast(world, rearrange, server);
    if (reconciliationModel == "DL")
    {
        broadcast(world, lossExpectedNumbers, server);
        broadcast(world, duplicationExpectedNumbers, server);
    }
    else if (reconciliationModel == "COAL")
    {
        broadcast(world, coalBls, server);
    }

    broadcast(world, currentSpeciesTree, server);
    broadcast(world, currentStep, server);
 /*  double t = MPI_Wtime();
    broadcast(world, t, server);
    std::cout << "broadcastsAllInformationButStop: " << currentStep<<" "<< setprecision(30)<< t <<std::endl;
    if (currentStep == 4 ) {
        std::cout << "CURRENT STEP = 4 "<< std::endl;
    }*/
  //  MPI_Barrier(world);
}


/******************************************************************************/
// This function runs the first communications between the server and the clients.
/******************************************************************************/
void firstCommunicationsServerClient (const mpi::communicator & world,
                                      unsigned int & server,
                                      std::vector <unsigned int>  & numbersOfGenesPerClient,
                                      unsigned int & assignedNumberOfGenes,
                                      std::vector<std::string> & assignedFilenames,
                                      std::vector <std::vector<std::string> > & listOfOptionsPerClient,
                                      bool & optimizeSpeciesTreeTopology,
                                      int & SpeciesNodeNumber,
                                      std::vector <double> & lossExpectedNumbers,
                                      std::vector <double> & duplicationExpectedNumbers,
                                      std::vector <int> & num0Lineages,
                                      std::vector <int> & num1Lineages,
                                      std::vector <int> & num2Lineages,
                                      std::vector <unsigned int> & num12Lineages,
                                      std::vector <unsigned int> & num22Lineages,
                                      std::vector < double >& coalBls,
                                      std::string & currentSpeciesTree)
{
    //   MPI_Barrier(world);

    //The server sends the number of genes assigned to each client
    scatter(world , numbersOfGenesPerClient, assignedNumberOfGenes, server);
    //Then the server sends each client the std::vector of filenames it is in charge of
    scatter(world , listOfOptionsPerClient, assignedFilenames, server);
    //Now the server sends to all clients the same information.
    broadcast(world, optimizeSpeciesTreeTopology, server);
    broadcast(world, SpeciesNodeNumber, server);
    broadcast(world, lossExpectedNumbers, server);
    broadcast(world, duplicationExpectedNumbers, server);
    broadcast(world, num0Lineages, server);
    broadcast(world, num1Lineages, server);
    broadcast(world, num2Lineages, server);
    broadcast(world, num12Lineages, server);
    broadcast(world, num22Lineages, server);
    broadcast(world, coalBls, server);

    broadcast(world, currentSpeciesTree, server);
    MPI_Barrier(world);

    //Then the clients will send back how many families they are going to work with.
    return;
}

/******************************************************************************/
// This function runs the second communication between the server and the clients.
// The clients send back how many families they deal with,
// and the server can then decide whether to stop or not.
/******************************************************************************/
void secondCommunicationsServerClient (const mpi::communicator & world ,
                                       unsigned int & server,
                                       unsigned int & whoami,
                                       unsigned int & finalNumberOfGeneFamilies,
                                       vector<unsigned int> & finalNumbersOfGeneFamilies,
                                       vector <double> & lossExpectedNumbers,
                                       vector <double> & duplicationExpectedNumbers,
                                       std::vector<double> &coalBls,
                                       std::string & currentSpeciesTree) {
//    MPI_Barrier(world);

    if (whoami == server)
        gather(world, finalNumberOfGeneFamilies, finalNumbersOfGeneFamilies, server);
    else {
        gather(world, finalNumberOfGeneFamilies, server);
    }

    broadcast(world, lossExpectedNumbers, server);
    broadcast(world, duplicationExpectedNumbers, server);
    broadcast(world, coalBls, server);

    broadcast(world, currentSpeciesTree, server);
 //   MPI_Barrier(world);
    return;
}

/******************************************************************************/
// This function runs the communication between the server and the clients for building a MRP species tree.
// The clients send back gene trees with one gene per species,
// and the server can then build the MRP species tree.
/******************************************************************************/
void mrpCommunicationsServerClient (const mpi::communicator & world,
                               unsigned int & server,
                               unsigned int & whoami,
                               string & trees1PerSpecies,
                               std::vector<string> & allTrees1PerSpecies) {
//    MPI_Barrier(world);

    if (whoami == server)
        gather(world, trees1PerSpecies, allTrees1PerSpecies, server);
    else {
        gather(world, trees1PerSpecies, server);
    }
  //  MPI_Barrier(world);
    return;
}

/******************************************************************************/
// This function runs the communication between the server and the clients for counting gene families after filtering.
/******************************************************************************/

void numberOfFilteredFamiliesCommunicationsServerClient (const mpi::communicator & world, unsigned int & server,
                                                    unsigned int & whoami, unsigned int & numberOfGeneFamilies) {

    if (whoami == server) {
        unsigned int allNumbersOfRemainingFamilies = 0;
        mpi::reduce(world, allNumbersOfRemainingFamilies, numberOfGeneFamilies, std::plus<unsigned int> (), 0);

        bool stop = false;

        if (numberOfGeneFamilies == 0) {
            stop = true;
            broadcast(world, stop, server);
            //MPI::COMM_WORLD.Abort(0);
            exit(0);
        }
        else {
            broadcast(world, stop, server);
            if (numberOfGeneFamilies == 1) {
                std::cout << numberOfGeneFamilies <<" gene family left after filtering. Continuing" <<std::endl;
            }
            else {
                std::cout << numberOfGeneFamilies <<" gene families left after filtering. Continuing" <<std::endl;
            }
        }
    }
    else {
        mpi::reduce(world, numberOfGeneFamilies, std::plus<unsigned int> (), 0);
        bool stop;
        broadcast(world, stop, server);
        if (stop) {
            exit(0);
        }
    }

    return;


}


/******************************************************************************/
// This function runs the communication between the server and the clients to stop the program.
// The server send stop and the index of the best tree found.
/******************************************************************************/
void lastCommunicationsServerClient (const mpi::communicator & world,
                                    unsigned int & server,
                                    bool & stop,
                                    unsigned int & bestIndex) {
   // MPI_Barrier(world);
    broadcast(world, stop, server);
    broadcast(world, bestIndex, server);
    return;
}


/************************************************************************
 * Gathers information from clients.
 ************************************************************************/


/*std::vector< int > PLUS_VEC (std::vector < int > i, std::vector< int > j)
{
    for (unsigned int k =0; k<i.size() ; k++ ) {
        i[k] = i[k]+j[k];
    }
    return i;
}*/


// unused since Boost > 1.49 (not sure but used to work with 1.49, does not work anymore since 1.55)

// struct plus_vec {
//     std::vector< int > operator()(std::vector< int > i, std::vector< int > j) const {
//         for (unsigned int k =0; k<i.size() ; k++ ) {
//             i[k] = i[k]+j[k];
//         }
//         return i;
//     }
// };
//
//
//
// struct plus_vec_unsigned {
//     std::vector< unsigned int > operator()(std::vector< unsigned int > i, std::vector< unsigned int > j) const {
//         for (unsigned int k =0; k<i.size() ; k++ ) {
//             i[k] = i[k]+j[k];
//         }
//         return i;
//     }
// };



/*
namespace boost { namespace mpi {

    //template<>
    struct is_commutative< plus_vec < std::vector < int>, std::vector< int > >,  std::vector< int > > : mpl::true_ { };

} } // end namespace boost::mpi
*/

void gathersInformationFromClients (const mpi::communicator & world,
                                    unsigned int & server,
                                    unsigned int &whoami,
                                    double &logL,
                                    std::vector<int> &num0Lineages,
                                    std::vector<int> &num1Lineages,
                                    std::vector<int> &num2Lineages,
                                    std::vector< std::vector<int> > &allNum0Lineages,
                                    std::vector< std::vector<int> > &allNum1Lineages,
                                    std::vector< std::vector<int> > &allNum2Lineages,
                                    std::vector<unsigned int> &num12Lineages,
                                    std::vector<unsigned int> &num22Lineages,
                                    std::string &reconciliationModel)
{
  //  MPI_Barrier(world);

    if (whoami == server) {
        //AFTER COMPUTATION IN CLIENTS

        logL = 0.0;
        double logLVal = 0.0;
    //    std::vector<double> logLs;
        resetVector(num0Lineages);
        resetVector(num1Lineages);
        resetVector(num2Lineages);
        resetVector(num12Lineages);
        resetVector(num22Lineages);
  //      gather(world, logL, logLs, server);
//        logL =  VectorTools::sum(logLs);
//	std::cout << "TEST: "<< logL << " and " << logLVal << std::endl;
        mpi::reduce(world, logLVal, logL, std::plus<double> (), 0);
	//	std::cout << "TEST BIS: "<< logL << " and " << logLVal << std::endl;


        if (reconciliationModel == "DL") {
            vector< int> tempNums = num0Lineages;
            mpi::reduce(world, &tempNums.front(), tempNums.size(), &num0Lineages.front(), std::plus<int>(), 0);
            mpi::reduce(world, &tempNums.front(), tempNums.size(), &num1Lineages.front(), std::plus<int>(), 0);
            mpi::reduce(world, &tempNums.front(), tempNums.size(), &num2Lineages.front(), std::plus<int>(), 0);
        }
        else if (reconciliationModel == "COAL") {
            vector< unsigned int> tempNums = num12Lineages;
            mpi::reduce(world, &tempNums.front(), tempNums.size(), &num12Lineages.front(), std::plus<unsigned int>(), 0);
            mpi::reduce(world, &tempNums.front(), tempNums.size(), &num22Lineages.front(), std::plus<unsigned int>(), 0);
        }

       /* std::cout <<"LOOK HERE:"<<std::endl;
        VectorTools::print(num0Lineages);
        VectorTools::print(num1Lineages);
        VectorTools::print(num2Lineages);
        resetVector(num0Lineages);
        resetVector(num1Lineages);
        resetVector(num2Lineages);
        gather(world, num0Lineages, allNum0Lineages, server);
        gather(world, num1Lineages, allNum1Lineages, server);
        gather(world, num2Lineages, allNum2Lineages, server);
        unsigned int temp = allNum0Lineages.size();
        for (unsigned int k =0; k<temp ; k++ ) {
            num0Lineages= num0Lineages+allNum0Lineages[k];
            num1Lineages= num1Lineages+allNum1Lineages[k];
            num2Lineages= num2Lineages+allNum2Lineages[k];
        }
        std::cout <<"COMPARED TO:"<<std::endl;
        VectorTools::print(num0Lineages);
        VectorTools::print(num1Lineages);
        VectorTools::print(num2Lineages);
        */
    }
    else {
        //Clients send back stuff to the server.
//        gather(world, logL, server);

        mpi::reduce(world, logL, std::plus<double> (), 0);


        if (reconciliationModel == "DL") {

            mpi::reduce(world, &num0Lineages.front(), num0Lineages.size(), std::plus<int>(), 0);
            mpi::reduce(world, &num1Lineages.front(), num1Lineages.size(), std::plus<int>(), 0);
            mpi::reduce(world, &num2Lineages.front(), num2Lineages.size(),std::plus<int>(), 0);
        }

        else if (reconciliationModel == "COAL") {
          mpi::reduce(world, &num12Lineages.front(), num12Lineages.size(), std::plus<unsigned int>(), 0);
          mpi::reduce(world, &num22Lineages.front(), num22Lineages.size(), std::plus<unsigned int>(), 0);
        }

/*      gather(world, num0Lineages, allNum0Lineages, server);
        gather(world, num1Lineages, allNum1Lineages, server);
        gather(world, num2Lineages, allNum2Lineages, server);*/
    }
  //  std::cout << "GATHERALL: "<<std::endl;
  /*  double t = MPI_Wtime();
    broadcast(world, t, server);
    std::cout << "GATHERALL: "<< setprecision(30)<< t <<std::endl;*/
 //   MPI_Barrier(world);
    return ;
}


/******************************************************************************/
// These functions input and output alternate topologies likelihoods.
/******************************************************************************/

void inputNNIAndRootLks(std::vector <double> & NNILks, std::vector <double> & rootLks, std::map<std::string, std::string> & params, std::string & suffix)
{
  std::string file = ApplicationTools::getAFilePath("alternate.topology.likelihoods", params, false, false, suffix, true);

  std::ifstream in (file.c_str(), std::ios::in);

  std::string line;

  if (in.is_open())
    {

    size_t i=0;
    while ((! in.eof() )  && (i<NNILks.size() ) )
      {

      getline (in,line);
      StringTokenizer st(line, "\t", true, true);
      if (st.getToken(0) != "NNI Likelihoods")
        {
        std::cout <<"Inputing values for node "<<i<<std::endl;
        NNILks[i] = TextTools::toDouble(st.getToken(0));
        rootLks[i] = TextTools::toDouble(st.getToken(1));
        i = i+1;
        }
      }
    in.close();
    }
}


void outputNNIAndRootLks(std::vector <double> & NNILks, std::vector <double> & rootLks, std::map<std::string, std::string> & params, std::string & suffix)
{
  std::string file = ApplicationTools::getAFilePath("alternate.topology.likelihoods", params, false, false, suffix, true);
  std::ofstream out (file.c_str(), std::ios::out);
  out << "NNI Likelihoods"<<"\t"<<"Root likelihoods"<<std::endl;
  for (unsigned int i = 0 ; i < NNILks.size() ; i++)
    {
    out<<TextTools::toString(NNILks[i], 15)<<"\t"<<TextTools::toString(rootLks[i], 15)<<std::endl;
    }
  out.close();
}








/**********************************CRAP******************************************/
