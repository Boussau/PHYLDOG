
/* This file contains various functions useful for the search of the species tree*/

#include "SpeciesTreeExploration.h"




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
                                            int &numIterationsWithoutImprovement, 
                                            unsigned int server, 
                                            int & nodeForNNI, 
                                            int & nodeForRooting, 
                                            std::string & branchExpectedNumbersOptimization, 
                                            std::map < std::string, int> genomeMissing,
                                            std::vector < double >  &NNILks, 
                                            std::vector<double> &rootLks, unsigned int currentStep) 
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
    stop = checkChangeHasNotBeenDone(*tree, bestTree, nodeForNNI, nodeForRooting, NNILks, rootLks);

    if (!stop) 
      {
      makeDeterministicNNIsAndRootChangesOnly(*tree, nodeForNNI, nodeForRooting);  

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
      
      
          if (logL+0.01<bestlogL) 
          {
              std::cout << "\t\tNNIs or Root changes: Improvement: new total Likelihood value "<<logL<<" compared to the best log Likelihood : "<<bestlogL<< std::endl;
              breadthFirstreNumber (*tree);
              std::cout << "Improved species tree: "<<TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;
              numIterationsWithoutImprovement = 0;
              bestlogL =logL;
              if (bestTree)
              {
                  delete bestTree;
                  bestTree = 0;
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
              for (unsigned int i = 0 ; i< NNILks.size() ; i++ ) 
              {
                  NNILks[i]=NumConstants::VERY_BIG;
                  rootLks[i]=NumConstants::VERY_BIG;
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
            if ( (numIterationsWithoutImprovement>2*tree->getNumberOfNodes()) || (ApplicationTools::getTime() >= timeLimit) ) 
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
                                         int &numIterationsWithoutImprovement, 
                                         unsigned int server, 
                                         int & nodeForNNI, 
                                         int & nodeForRooting, 
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
                                                bool rearrange, int &numIterationsWithoutImprovement, 
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
      }
      else {
        numIterationsWithoutImprovement++;
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
    std::string currentSpeciesTree = TreeTemplateTools::treeToParenthesis(*currentTree, true);
//TEST    breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers);
      breadthFirstreNumber (*currentTree);
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
                            bool rearrange, int &numIterationsWithoutImprovement, unsigned int server, 
                            std::string &branchExpectedNumbersOptimization, std::map < std::string, int> genomeMissing, 
                            int sprLimit, bool optimizeRates, unsigned int currentStep) {

    std::vector <double> bestDupProba;
    std::vector <double> bestLossProba;
    std::vector <double> bestCoalBls;
    if (reconciliationModel == "DL") {
        breadthFirstreNumber (*currentTree, duplicationExpectedNumbers, lossExpectedNumbers);
        bestDupProba=duplicationExpectedNumbers;
        bestLossProba=lossExpectedNumbers;
    }
    else if (reconciliationModel == "COAL") {
        breadthFirstreNumber (*currentTree, coalBls);
        bestCoalBls=coalBls;
    }
    std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> *tree = 0;
  for (unsigned int nodeForSPR=currentTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) {
    buildVectorOfRegraftingNodesLimitedDistance(*currentTree, nodeForSPR, sprLimit, nodeIdsToRegraft);
    betterTree = false;
    for (int i =0 ; i<nodeIdsToRegraft.size() ; i++) {
      if (tree) {
        delete tree;
        tree = 0;
      }
      tree = currentTree->clone();

      makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i]);
        breadthFirstreNumber (*tree);
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
       std::cout << TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;
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
                                         int &numIterationsWithoutImprovement, 
                                         unsigned int server, 
                                         std::string &branchExpectedNumbersOptimization, 
                                         std::map < std::string, int> genomeMissing, 
                                         int sprLimit, 
                                         bool optimizeRates, 
                                         unsigned int currentStep) {
  if (optimizeRates)
    {
    std::cout <<"Making SPRs and NNIs and optimizing duplication and loss rates."<< std::endl;
    }
  else 
    {
    std::cout <<"Making SPRs and NNIs but NOT optimizing duplication and loss rates."<< std::endl;
    }
  numIterationsWithoutImprovement=0;
  
  while (numIterationsWithoutImprovement<2*currentTree->getNumberOfNodes()) {
    fastTryAllPossibleSPRs(world, currentTree, bestTree, 
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
                           server, branchExpectedNumbersOptimization, genomeMissing, 
                           sprLimit, optimizeRates, currentStep);

    if (ApplicationTools::getTime() >= timeLimit) 
      {	
        stop = true;
          lastCommunicationsServerClient (world, 
                                          server, 
                                          stop, 
                                          bestIndex);
        break;      
      }    
    if  (numIterationsWithoutImprovement>2*currentTree->getNumberOfNodes()) 
      {	
        break;
      }
    
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
                computeCoalBls (num12Lineages, 
                                     num22Lineages, 
                                     coalBls) ;
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
    }
  std::cout <<"\t\tNumber of species trees tried: "<< index ;
  std::cout<<";\t\tLogLk value for this species tree: "<< - logL<<std::endl;

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
        reduce(world, allNumbersOfRemainingFamilies, numberOfGeneFamilies, std::plus<unsigned int> (), 0);

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
        reduce(world, numberOfGeneFamilies, std::plus<unsigned int> (), 0);
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


struct plus_vec {
    std::vector< int > operator()(std::vector< int > i, std::vector< int > j) const { 
        for (unsigned int k =0; k<i.size() ; k++ ) {
            i[k] = i[k]+j[k];
        }
        return i; 
    }
};



struct plus_vec_unsigned {
    std::vector< unsigned int > operator()(std::vector< unsigned int > i, std::vector< unsigned int > j) const { 
        for (unsigned int k =0; k<i.size() ; k++ ) {
            i[k] = i[k]+j[k];
        }
        return i; 
    }
};



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
        reduce(world, logLVal, logL, std::plus<double> (), 0);
        if (reconciliationModel == "DL") {
            reduce(world, num0Lineages, num0Lineages, plus_vec(), 0);
            reduce(world, num1Lineages, num1Lineages, plus_vec(), 0);
            reduce(world, num2Lineages, num2Lineages, plus_vec(), 0);
        }
        else if (reconciliationModel == "COAL") {
            reduce(world, num12Lineages, num12Lineages, plus_vec_unsigned(), 0);
            reduce(world, num22Lineages, num22Lineages, plus_vec_unsigned(), 0);
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
        reduce(world, logL, std::plus<double> (), 0);
        if (reconciliationModel == "DL") {
            reduce(world, num0Lineages, plus_vec(), 0);
            reduce(world, num1Lineages, plus_vec(), 0);
            reduce(world, num2Lineages, plus_vec(), 0);
        }
        else if (reconciliationModel == "COAL") {
            reduce(world, num12Lineages, plus_vec_unsigned(), 0);
            reduce(world, num22Lineages, plus_vec_unsigned(), 0);
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

    int i=0;
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





/*  
 broadcastsAllInformation(world, server, stop, rearrange, lossExpectedNumbers, duplicationExpectedNumbers, currentSpeciesTree);
 //COMPUTATION IN CLIENTS
 index++;  
 std::cout <<"\t\tNumber of species trees tried: "<<index<< std::endl;
 logL = 0.0;
 std::vector<double> logLs;
 resetVector(num0Lineages);
 resetVector(num1Lineages);
 resetVector(num2Lineages);
 gather(world, logL, logLs, server);
 
 logL =  VectorTools::sum(logLs);
 std::cout<<"New minus logLk value in computeSpeciesTreeLikelihood: "<<logL<<std::endl;
 gather(world, num0Lineages, allNum0Lineages, server);
 gather(world, num1Lineages, allNum1Lineages, server);
 gather(world, num2Lineages, allNum2Lineages, server);
 int temp = allNum0Lineages.size();
 for (int k =0; k<temp ; k++ ) {
 num0Lineages= num0Lineages+allNum0Lineages[k];
 num1Lineages= num1Lineages+allNum1Lineages[k];
 num2Lineages= num2Lineages+allNum2Lineages[k];
 }
 */





/*  std::string currentSpeciesTree = TreeTemplateTools::treeToParenthesis(*tree, true);
 
 while (!stop) {
 std::vector<double> logLs;
 broadcastsAllInformation(world, server, stop, rearrange, lossExpectedNumbers, duplicationExpectedNumbers, currentSpeciesTree);
 //Computation in clients
 index++;  
 std::cout <<"\tNumber of species trees tried : "<<index<< std::endl;
 logL = 0.0;
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
 } 
 
 if (logL+0.01<bestlogL) {
 std::cout << "\t\tDuplication and loss probabilities optimization: Server : new total Likelihood value "<<logL<<" compared to the best log Likelihood : "<<bestlogL<< std::endl;
 std::cout << TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;
 numIterationsWithoutImprovement = 0;
 bestlogL =logL;
 if (bestTree) 	
 {
 deleteTreeProperties(*bestTree);
 delete bestTree;
 }
 bestTree = tree->clone();
 bestIndex = index;
 bestNum0Lineages = num0Lineages;
 bestNum1Lineages = num1Lineages;
 bestNum2Lineages = num2Lineages;
 
 computeDuplicationAndLossRatesForTheSpeciesTree (branchExpectedNumbersOptimization, num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers, genomeMissing, *tree);
 breadthFirstreNumber (*tree, duplicationExpectedNumbers, lossExpectedNumbers);
 
 broadcastsAllInformation(world, server, stop, rearrange, lossExpectedNumbers, duplicationExpectedNumbers, currentSpeciesTree);
 
 //COMPUTATION IN CLIENTS
 index++;  
 bestIndex = index;
 std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
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
 temp = allNum0Lineages.size();
 for (int i =0; i<temp ; i++ ) {
 num0Lineages= num0Lineages+allNum0Lineages[i];
 num1Lineages= num1Lineages+allNum1Lineages[i];
 num2Lineages= num2Lineages+allNum2Lineages[i];
 }
 if ( (ApplicationTools::getTime() > timeLimit) ) {
 stop = true;
 broadcast(world, stop, server); 
 broadcast(world, bestIndex, server);
 }
 }
 else {
 numIterationsWithoutImprovement++;
 std::cout <<"\t\tNo more increase in likelihood "<< std::endl;
 deleteTreeProperties(*tree);
 delete tree;
 tree = bestTree->clone();
 stop = true;
 broadcast(world, stop, server); 
 broadcast(world, bestIndex, server);
 }
 }*/






/************************************************************************
 * Optimizes the duplication and loss rates with a numerical algorithm.
 * As we should start from overestimated values, we start by decreasing all 
 * values, and keep on trying to minimize these rates until the likelihood stops 
 * increasing. 
 ************************************************************************/
/*void numericalOptimizationOfDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers, bool rearrange, int server, std::string &branchExpectedNumbersOptimization, std::map < std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL) {
  std::cout<<"Numerical optimization of Duplication And Loss expected numbers of events"<<std::endl;
  std::cout<<"Before optimization: Duplications"<<std::endl;
  VectorTools::print(duplicationExpectedNumbers);  
  std::cout<<"Before optimization: Losses"<<std::endl;
  VectorTools::print(lossExpectedNumbers);
  
  std::vector <double> formerDupRates = duplicationExpectedNumbers;
  std::vector <double> formerLossRates = lossExpectedNumbers;
  
  std::vector <double> backupDupProba  = duplicationExpectedNumbers;
  std::vector <double> backupLossProba = lossExpectedNumbers;
  
  std::vector <double> newDupRates = duplicationExpectedNumbers;
  std::vector <double> newLossRates = lossExpectedNumbers;
  double currentlogL = -UNLIKELY;
  std::cout <<"HEHEHEHEHcurrentlogL: "<<currentlogL<<" logL: "<<logL<<std::endl;
  for (int i = 1 ; i<10 ; i++) {
    currentlogL = logL;
    double d = (double)i /10.0;
    newDupRates = formerDupRates * d;
    newLossRates = formerLossRates * d;
    
    computeSpeciesTreeLikelihood(world, index, stop, logL, num0Lineages, num1Lineages, num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, newDupRates, newLossRates, rearrange, server, branchExpectedNumbersOptimization, genomeMissing, tree);
    std::cout <<"i: "<<i<<" currentlogL: "<<currentlogL<<" logL: "<<logL<<std::endl;
    //If we have improved the log-likelihood
    if (currentlogL < bestlogL) {
      duplicationExpectedNumbers = newDupRates;
      lossExpectedNumbers = newLossRates;
      bestlogL = currentlogL;
      std::cout<<"Better logLikelihood! logL: "<< currentlogL <<std::endl;
      std::cout<<"After optimization: Duplications"<<std::endl;
      VectorTools::print(duplicationExpectedNumbers);  
      std::cout<<"After optimization: Losses"<<std::endl;
      VectorTools::print(lossExpectedNumbers);
    }
  }
  stop = true;
  broadcast(world, stop, server); 
  broadcast(world, index, server);
  
  return;  
  
}

*/












/************************************************************************
 * Randomly chooses the modification to be done, either SPR or changeRoot. The probabilities of the two events are a function of expectedFrequencySPR and expectedFrequencyChangeRoot.
************************************************************************/
/*
void makeRandomModifications(TreeTemplate<Node> &tree, int expectedFrequencyNNI, int expectedFrequencySPR, int expectedFrequencyChangeRoot) {
  //breadthFirstreNumber (tree);
  int tot = expectedFrequencyNNI + expectedFrequencySPR + expectedFrequencyChangeRoot;
  int ran=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(tot);
  std::vector <int> allNodeIds;
  //TEMPORARY
  //ran = 0;
  if (ran < expectedFrequencyNNI) {//Make a NNI move
    //we benefit from the fact that the tree has been numbered with the lowest indices as closest to the root : we discard 0,1,2    
    for (int i = 3; i<tree.getNumberOfNodes(); i++) {
      allNodeIds.push_back(i);
    }
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size())+3;
    makeNNI(tree, ran);
  }
  else if (ran < expectedFrequencySPR ) { //Make a SPR move ! we do not want to get the root as one of the two nodes involved
    Node * N = tree.getRootNode();
    std::vector <int> firstHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(0)));
    std::vector <int> secondHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(1)));
    
    if (secondHalfNodeIds.size()<firstHalfNodeIds.size()) {
      allNodeIds = firstHalfNodeIds;
      for (int i = 0 ; i< secondHalfNodeIds.size() ; i++ ) {
	allNodeIds.push_back(secondHalfNodeIds[i]);
      }
    }
    else {
      allNodeIds = secondHalfNodeIds;
      for (int i = 0 ; i< firstHalfNodeIds.size() ; i++ ) {
	allNodeIds.push_back(firstHalfNodeIds[i]);
      }
    }
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    int cutNodeId = allNodeIds[ran];
    //TEMPORARY
    //  cutNodeId = 2;
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(cutNodeId)));
    std::vector <int> toRemove;
    for (int i = 0 ; i< allNodeIds.size() ; i++) {
      if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
	toRemove.push_back(i);
      }
    }
    sort(toRemove.begin(), toRemove.end(), cmp);
    for (int i = 0 ; i< toRemove.size() ; i++) {
      std::vector<int>::iterator vi = allNodeIds.begin();
      allNodeIds.erase(vi+toRemove[i]);
    }
    int ran2 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    //TEMPORARY
    // ran2 = 0;
    int newBrotherId = allNodeIds[ran2];
    makeSPR(tree, cutNodeId, newBrotherId);
    }
  else { //Make a ChangeRoot move !
   std::cout << "CHANGE ROOT MOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< std::endl;
    allNodeIds = tree.getNodesId();
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    changeRoot(tree, allNodeIds[ran]);
  }
}






 //TEST
    /*  std::cout << TreeTemplateTools::treeToParenthesis(*currentTree, true)<< std::endl;

    //Send the new std::vectors to compute the new likelihood of the best tree
    computeAverageDuplicationAndlossExpectedNumbersForAllBranches (num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers);
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    broadcast(world, lossExpectedNumbers, server); 
    broadcast(world, duplicationExpectedNumbers, server); 
    currentSpeciesTree = TreeTemplateTools::treeToParenthesis(*currentTree, false);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    */
    //END TEST













