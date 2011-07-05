/*
 *  SpeciesTreeLikelihood.cpp
 *  
 *
 *  Created by boussau on 18/02/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */

#include "SpeciesTreeLikelihood.h"
#include "SpeciesTreeExploration.h"

using namespace bpp;

using namespace std;

/*******************************************************************************/

void SpeciesTreeLikelihood::updateDuplicationAndLossExpectedNumbers() 
{
  double d = getParameterValue("coefDup");
  duplicationExpectedNumbers_ = backupDuplicationExpectedNumbers_ * d;
  double l = getParameterValue("coefLoss");
  lossExpectedNumbers_ = backupLossExpectedNumbers_ * d;
}

/*******************************************************************************/
void SpeciesTreeLikelihood::initialize() 
{
  rearrange_ = false;

  /****************************************************************************
   * First communications between the server and the clients.
   *****************************************************************************/
  firstCommunicationsServerClient (world_ , server_, numbersOfGenesPerClient_, assignedNumberOfGenes_,
                                   assignedFilenames_, listOfOptionsPerClient_, 
                                   optimizeSpeciesTreeTopology_, speciesTreeNodeNumber_, 
                                   lossExpectedNumbers_, duplicationExpectedNumbers_, num0Lineages_, 
                                   num1Lineages_, num2Lineages_, currentSpeciesTree_);
  
  /****************************************************************************
   * First species tree likelihood computation.
   *****************************************************************************/
  breadthFirstreNumber (*tree_, duplicationExpectedNumbers_, lossExpectedNumbers_);
  bestTree_ = tree_->clone();
  currentTree_ = tree_->clone();

  std::cout << TreeTools::treeToParenthesis(*currentTree_, true)<<std::endl;
  
  computeSpeciesTreeLikelihoodWithGivenStringSpeciesTree(world_, 
                                                         index_, 
                                                         stop_, 
                                                         logL_, 
                                                         num0Lineages_, 
                                                         num1Lineages_, 
                                                         num2Lineages_, 
                                                         allNum0Lineages_, 
                                                         allNum1Lineages_, 
                                                         allNum2Lineages_, 
                                                         lossExpectedNumbers_, 
                                                         duplicationExpectedNumbers_, 
                                                         rearrange_, 
                                                         server_, 
                                                         branchExpectedNumbersOptimization_, 
                                                         genomeMissing_, 
                                                         *tree_, 
                                                         currentSpeciesTree_,
                                                         true); 
  std::cout << "\t\tServer: total initial Likelihood value "<<logL_<<std::endl;
  bestlogL_ = logL_;
  ApplicationTools::displayTime("Execution time so far:");
  for (int i =0; i<num0Lineages_.size() ; i++ ) {
    std::cout <<"branch Number#"<< i<<"Expected numbers:  dup: "<< duplicationExpectedNumbers_[i]<<" loss: "<< lossExpectedNumbers_[i]<<std::endl;
  }
  bestNum0Lineages_ = num0Lineages_;
  bestNum1Lineages_ = num1Lineages_;
  bestNum2Lineages_ = num2Lineages_;

  return;    
}



/*******************************************************************************/
void SpeciesTreeLikelihood::computeLogLikelihood() 
{
  computeSpeciesTreeLikelihood(world_, index_,  stop_, logL_, num0Lineages_, 
                               num1Lineages_, num2Lineages_, 
                               allNum0Lineages_, allNum1Lineages_, 
                               allNum2Lineages_, lossExpectedNumbers_, 
                               duplicationExpectedNumbers_, rearrange_, server_, 
                               branchExpectedNumbersOptimization_, genomeMissing_, *tree_);
}




/*******************************************************************************/
void SpeciesTreeLikelihood::parseOptions()
{
  std::vector <std::string> spNames;
  /****************************************************************************
   //First, we need to get the species tree.
   *****************************************************************************/
  // Get the initial tree
  std::string initTree = ApplicationTools::getStringParameter("init.species.tree", params_, "user", "", false, false);
  ApplicationTools::displayResult("Input species tree", initTree);
  // A given species tree
  if(initTree == "user")
    {
    std::string spTreeFile =ApplicationTools::getStringParameter("species.tree.file",params_,"none");
    if (spTreeFile=="none" )
      {
      std::cout << "\n\nNo Species tree was provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file.\nIf you do not have a species tree file, the program can start from a random tree, if you set init.species.tree at random, and give a list of species names as species.names.file\n\n" << std::endl;
      exit(-1);
      }
    ApplicationTools::displayResult("Species Tree file", spTreeFile);
    Newick newick(true);
    tree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
    if (!tree_->isRooted()) 
      {
      std::cout << "The tree is not rooted, midpoint-rooting it!\n";
      TreeTools::midpointRooting(*tree_);
      }
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree_->getNumberOfLeaves()));
    spNames=tree_->getLeavesNames();
    }
  //We build a random species tree
  else if(initTree == "random")
    {
    //Then we need a file containing species names
    std::string spNamesFile =ApplicationTools::getStringParameter("species.names.file",params_,"none");
    if (spNamesFile=="none" )
      {
      std::cout << "\n\nNo Species names were provided. The option init.species.tree is set to random, which means that the option species.names.file must be filled with the path of a file containing a list of species names. \n\n" << std::endl;
      exit(-1);
      }
    std::ifstream inListNames (spNamesFile.c_str());
    std::string line;
    while(getline(inListNames,line)) 
      {
      spNames.push_back(line);
      }
    cleanVectorOfOptions(spNames, false);	
    if (spNames.size()<2) 
      {
      std::cout <<"No more than two species: no need to make a species tree! Quitting."<<std::endl;
      exit (0);
      }
    tree_ = TreeTemplateTools::getRandomTree(spNames);
    tree_->setBranchLengths(1.0);
    TreeTools::midpointRooting(*tree_);
    std::cout <<"Initial species tree<<std::endl";
    std::cout <<TreeTools::treeToParenthesis(*tree_, false)<<std::endl;
    }
  else throw Exception("Unknown init.species.tree method.");
  assignArbitraryBranchLengths(*tree_);
  //We write the starting species tree to a file
  std::string file = ApplicationTools::getStringParameter("starting.tree.file", params_, "starting.tree");

  Newick newick;

  newick.write(*tree_, file, true);
  // Try to write the current tree to file. This will be overwritten by the optimized tree,
  // but allows to check file existence before running optimization!

  PhylogeneticsApplicationTools::writeTree(*tree_, params_, "", "", true, false, true);

  speciesTreeNodeNumber_ = tree_->getNumberOfNodes();
  

  /****************************************************************************
   * branchExpectedNumbersOptimization_ sets the type of optimization that is applied to duplication and loss rates:
   * "average": all branches have the same average rates. These rates are optimized during the course of the program.
   * "branchwise": each branch has its own set of duplication and loss rates. These rates are optimized during the course of the program.
   * "average_then_branchwise": at the beginning, all branches have the same average rates, and then they are individualised.
   *****************************************************************************/
  optimizeSpeciesTreeTopology_ = ApplicationTools::getBooleanParameter("optimization.topology", params_, false, "", true, false);
  branchExpectedNumbersOptimization_ = ApplicationTools::getStringParameter("branchProbabilities.optimization",params_,"average");
  std::cout << "Branch expected numbers of duplications and losses: "<<branchExpectedNumbersOptimization_ <<std::endl;
  if ((branchExpectedNumbersOptimization_!="average")&&(branchExpectedNumbersOptimization_!="branchwise")&&
      (branchExpectedNumbersOptimization_!="average_then_branchwise")&&(branchExpectedNumbersOptimization_!="no")) 
    {
    std::cout << "branchProbabilities.optimization is not properly set; please set to either 'average', 'branchwise', 'average_then_branchwise', or 'no'. "<<std::endl;
    exit(-1);
    }
  //When doing a spr, how far from the original position can we regraft the pruned subtree?
  //Hordjik and Gascuel (2005) consider 10% of the total number of hedges is already large.
  sprLimit_=ApplicationTools::getIntParameter("spr.limit",params_,4);
  

  /****************************************************************************
   // We get percent coverage of the genomes under study.
   // We produce a genomeMissing std::map that associates species names to percent of missing data.
   *****************************************************************************/
  std::string percentCoverageFile = ApplicationTools::getStringParameter("genome.coverage.file",params_,"none");
  for(std::vector<std::string >::iterator it = spNames.begin(); it != spNames.end(); it++)
    {
    genomeMissing_[*it]=0;
    }
  
  if (percentCoverageFile=="none" )
    {
    std::cout << "No file for genome.coverage.file, we assume that all genomes are covered at 100%."<<std::endl;
    }
  else 
    {
    std::ifstream inCoverage (percentCoverageFile.c_str());
    std::vector <std::string> listCoverages;
    std::string line;
    while(getline(inCoverage,line)) 
      {
      listCoverages.push_back(line);
      }
    cleanVectorOfOptions(listCoverages, false);
    for(std::vector<std::string >::iterator it = listCoverages.begin(); it != listCoverages.end(); it++)
      {
      StringTokenizer st1 = StringTokenizer::StringTokenizer (*it, ":", true);
      std::map<std::string,int>::iterator iter = genomeMissing_.find(st1.getToken(0));
      if (iter != genomeMissing_.end() ) 
        {
        iter->second = 100-TextTools::toInt(st1.getToken(1));
        }
      else 
        {
        std::cout <<"Coverage from species "<<st1.getToken(0)<<" is useless as long as species "<<st1.getToken(0)<<" is not included in the dataset..."<<std::endl;
        }
      }
    }

  /****************************************************************************
   // Then, we get the maximum time allowed, in hours. 1h before the limit, 
   // the program will nicely stop and save the current species tree and 
   // the step we are at in a new option file, so that next time the program
   // runs, we approximately start from the same place.
   *****************************************************************************/
  timeLimit_ = ApplicationTools::getIntParameter("time.limit",params_,1000);
  //we convert hours in seconds. 
  //1h before the end, we will launch the stop signal.
  if (timeLimit_ <= 1 )
    {
    timeLimit_ = 2;
    std::cout<<"Setting time.limit to 2, otherwise the program will exit after 0 (=(1-1) * 3600) second !."<<std::endl;
    }
  timeLimit_ = (timeLimit_-1) * 3600 ; 
  std::cout<<"The program will exit after "<<timeLimit_<<" seconds."<<std::endl;
  //Current step: information to give when the run is stopped for timing reasons, 
  //or to get if given.
  currentStep_ = ApplicationTools::getIntParameter("current.step",params_,0);
  
  /****************************************************************************
   * We can also get duplication and loss expected numbers, 
   * if the respective trees have been provided. Otherwise, 
   * we set the expected numbers of duplications to arbitrary values,
   * we set the numbers of duplications and losses to 0, 
   * we take into account genome coverage.
   *****************************************************************************/
  
  for (int i=0; i<speciesTreeNodeNumber_; i++) 
    {
    //We fill the vectors with 0s until they are the right sizes.
    lossExpectedNumbers_.push_back(0);
    duplicationExpectedNumbers_.push_back(0);
    num0Lineages_.push_back(0);
    num1Lineages_.push_back(0);
    num2Lineages_.push_back(0);
    }
  
  resetLossesAndDuplications(*tree_, lossExpectedNumbers_, duplicationExpectedNumbers_);
  breadthFirstreNumber (*tree_);
  std::string spTreeDupFile =ApplicationTools::getStringParameter("species.duplication.tree.file",params_,"none");
  std::string spTreeLossFile =ApplicationTools::getStringParameter("species.loss.tree.file",params_,"none");
  
  if ( (spTreeDupFile == "none") && (spTreeLossFile == "none") ) {
   //We set preliminary loss and duplication rates, correcting for genome coverage
  computeDuplicationAndLossRatesForTheSpeciesTreeInitially(branchExpectedNumbersOptimization_, 
                                                           num0Lineages_, 
                                                           num1Lineages_, 
                                                           num2Lineages_, 
                                                           lossExpectedNumbers_, 
                                                           duplicationExpectedNumbers_, 
                                                           genomeMissing_, 
                                                           *tree_);
    
  }
  else {
    if ( (spTreeDupFile == "none") || (spTreeLossFile == "none") ) {
      std::cout <<"Sorry, you have to input both loss and duplication rates, or none of them."<<std::endl;
      exit(-1);
    }
    TreeTemplate<Node> * treeD = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeDupFile));
    breadthFirstreNumber (*treeD);
    TreeTemplate<Node> * treeL = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeLossFile));
    breadthFirstreNumber (*treeL);

    for (int i = 0 ; i < treeD->getNumberOfNodes() ; i++)
      {
       if (treeD->getNode(i)->hasFather()) 
         {
         duplicationExpectedNumbers_[i] = treeD->getNode(i)->getDistanceToFather();
         }
      if (treeL->getNode(i)->hasFather()) 
        {
        lossExpectedNumbers_[i] = treeL->getNode(i)->getDistanceToFather();
        }
      }
    backupDuplicationExpectedNumbers_ = duplicationExpectedNumbers_;
    backupLossExpectedNumbers_ = lossExpectedNumbers_;
  }

  currentSpeciesTree_ = TreeTools::treeToParenthesis(*tree_, true);
  bestlogL_ = -UNLIKELY;
  numIterationsWithoutImprovement_ = 0;
  bestIndex_ = 0;
  index_ = 0; 
  //vector to keep NNIs likelihoods, for making aLRTs.
  for (int i = 0 ; i <tree_->getNumberOfNodes() ; i ++)
    {
    NNILks_.push_back(NumConstants::VERY_BIG);
    rootLks_.push_back(NumConstants::VERY_BIG);
    }
  
  
  
  /****************************************************************************
   // Then, we handle all gene families.
   // First, we get and read the file containing the list of all gene family files.
   //In this file, the format is expected to be as follows :
   optionFile1 for gene1
   optionFile2 for gene2
   optionFile3 for gene3
   ...
   //Alternatively, there can be sizes for the gene families, in which case the format is:
   optionFile1:size1
   optionFile2:size2
   optionFile3:size3
   ...
   *****************************************************************************/
  std::string listGeneFile = ApplicationTools::getStringParameter("genelist.file",params_,"none");
  if (listGeneFile=="none" )
    {
    std::cout << "\n\nNo list of genes was provided. Cannot compute a reconciliation between a species tree and gene families if you do not tell me where the options for gene families are ! Use the option genelist.file, which should be a file containing a list of gene option files.\n" << std::endl;
    std::cout << "ReconcileDuplications species.tree.file=bigtree taxaseq.file=taxaseqlist genelist.file=geneList output.tree.file=outputtree\n"<<std::endl;
    exit(-1);
    }

  //Addresses to these option files are expected to be absolute: 
  //may be improved, for instance through the use of global variables, as in other files.
  std::ifstream inListOpt (listGeneFile.c_str());
  std::string line;
  std::vector<std::string> listOptions;
  while(getline(inListOpt,line)) 
    {
    listOptions.push_back(line);
    }
  cleanVectorOfOptions(listOptions, true);
  std::cout <<"Using "<<listOptions.size()<<" gene families."<<std::endl;
  std::cout <<"Using " <<size_<<" nodes"<<std::endl;
  if (listOptions.size()==0) 
    {
    std::cout << "\n\nThere should be at least one option file specified in the list of option files !"<<std::endl;
    exit(-1);
    }
  else if (listOptions.size() < size_ -1) 
    {
    std::cout << "You want to use more nodes than (gene families+1). This is not possible (and useless). Please decrease the number of processors to use."<<std::endl;
    exit(-1);
    }
  //We compute and send the number of genes per client.
  generateListOfOptionsPerClient(listOptions, size_, listOfOptionsPerClient_, numbersOfGenesPerClient_);
  
  suffix_ = ApplicationTools::getStringParameter("output.file.suffix", params_, "", "", false, false);
  return;  
}


/*******************************************************************************/
void SpeciesTreeLikelihood::MLsearch()
{
  //Indices used in the exploration
  int nodeForNNI = 0;
  int nodeForRooting = 4;
  bool noMoreSPR;
  stop_ = false;
  if(optimizeSpeciesTreeTopology_) 
    {
    std::cout <<"Optimizing the species tree topology"<<std::endl;
    noMoreSPR=false; 
    }
  else 
    {
    noMoreSPR=true;
    }  
  
  /****************************************************************************
   ****************************************************************************
   * Species tree likelihood optimization: main loop.
   ****************************************************************************
   *****************************************************************************/
  while (!stop_) 
    {
    //Using deterministic SPRs first, and then NNIs
    //Making SPRs, from leaves to deeper nodes (approximately), plus root changes
    if (currentStep_ >= 3)
      {
      noMoreSPR = true;
      numIterationsWithoutImprovement_=2*speciesTreeNodeNumber_;
      }
    if (!noMoreSPR) 
      {
      if (currentStep_ == 0) 
        {
        /****************************************************************************
         * Doing SPRs and rerootings. This first function does not optimize duplication 
         * and loss (D and L) expected numbers, as proved by the last "false" argument.
         *****************************************************************************/            
        fastTryAllPossibleSPRsAndReRootings(world_, currentTree_, bestTree_, 
                                            index_, bestIndex_, stop_, timeLimit_, 
                                            logL_, bestlogL_, 
                                            num0Lineages_, num1Lineages_, num2Lineages_, 
                                            bestNum0Lineages_, bestNum1Lineages_, 
                                            bestNum2Lineages_, 
                                            allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                            lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                            rearrange_, numIterationsWithoutImprovement_, 
                                            server_, branchExpectedNumbersOptimization_, 
                                            genomeMissing_, sprLimit_, false);
        if (ApplicationTools::getTime() < timeLimit_) 
          {
          std::cout << "\n\n\t\t\tFirst fast step of SPRs and rerootings over.\n\n"<< std::endl;
          ApplicationTools::displayTime("Execution time so far:");
          currentStep_ = 1;
          }  
        else 
          {
          std::cout << "\n\n\t\t\tFirst fast step of SPRs and rerootings is not over yet. The program exits because of the time limit.\n\n"<< std::endl;
          ApplicationTools::displayTime("Execution time so far:");
          if (stop_==false)
            {
            stop_ = true;
            broadcast(world_, stop_, server_); 
            broadcast(world_, bestIndex_, server_);                      
            }              
          }
        }
      
      /****************************************************************************
       * Doing SPRs and rerootings. Now we update expected numbers of D and L.
       *****************************************************************************/            
      std::cout <<"Before updating expected numbers of Duplication and loss; current Likelihood "<<bestlogL_  <<" and logL: "<<logL_<<std::endl;
      backupLossExpectedNumbers_ = lossExpectedNumbers_;
      backupDuplicationExpectedNumbers_ = duplicationExpectedNumbers_;
      if ((currentStep_ == 1) && (ApplicationTools::getTime() < timeLimit_) )
        {
        computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world_, index_, 
                                                                           stop_, logL_, 
                                                                           num0Lineages_, num1Lineages_, num2Lineages_, 
                                                                           allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                                                           lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                                                           rearrange_, server_, 
                                                                           branchExpectedNumbersOptimization_, genomeMissing_, 
                                                                           *currentTree_, bestlogL_);
        std::cout <<"After updating expected numbers of Duplication and loss; current Likelihood "<<logL_<<std::endl;
        if (logL_+0.01<bestlogL_) 
          {
          bestlogL_ =logL_;
          bestNum0Lineages_ = num0Lineages_;
          bestNum1Lineages_ = num1Lineages_;
          bestNum2Lineages_ = num2Lineages_;
          bestIndex_ = index_;
          std::cout << "Updating duplication and loss expected numbers yields a better candidate tree likelihood : "<<bestlogL_<<std::endl;
          std::cout << TreeTools::treeToParenthesis(*currentTree_, true)<<std::endl;
          } 
        else 
          {
          std::cout << "No improvement: we keep former expected numbers. "<<bestlogL_<<std::endl;
          lossExpectedNumbers_ = backupLossExpectedNumbers_;
          duplicationExpectedNumbers_ = backupDuplicationExpectedNumbers_;
          }
        for (int i =0; i<num0Lineages_.size() ; i++ ) 
          {
          std::cout <<"branch Number#"<< i<<"Expected numbers:  dup: "<< duplicationExpectedNumbers_[i]<<" loss: "<< lossExpectedNumbers_[i]<<std::endl;
          }
        if (ApplicationTools::getTime() < timeLimit_) 
          currentStep_ = 2;
        }
      if ( (currentStep_ == 2) && (ApplicationTools::getTime() < timeLimit_) )
        {
        fastTryAllPossibleSPRsAndReRootings(world_, currentTree_, bestTree_, 
                                            index_, bestIndex_, stop_, timeLimit_, 
                                            logL_, bestlogL_, 
                                            num0Lineages_, num1Lineages_, num2Lineages_, 
                                            bestNum0Lineages_, bestNum1Lineages_, bestNum2Lineages_, 
                                            allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                            lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                            rearrange_, numIterationsWithoutImprovement_, 
                                            server_, branchExpectedNumbersOptimization_, 
                                            genomeMissing_, sprLimit_, true);
        if (ApplicationTools::getTime() < timeLimit_) 
          {
          std::cout << "\n\n\t\t\tStep of SPRs and rerootings with optimization of the duplication and loss parameters over.\n\n"<< std::endl;
          ApplicationTools::displayTime("Execution time so far:");
          currentStep_ = 3;
          }
        else 
          {
          std::cout << "\n\n\t\t\tStep of SPRs and rerootings with optimization of the duplication and loss parameters is not over yet. The program exits because of the time limit.\n\n"<< std::endl;
          ApplicationTools::displayTime("Execution time so far:");
          if (stop_==false)
            {
            stop_ = true;
            broadcast(world_, stop_, server_); 
            broadcast(world_, bestIndex_, server_);                      
            }
          }
        
        noMoreSPR=true;
        
        std::cout << "\n\n\t\t\tNow entering the final optimization steps, without SPRs\n\n"<<std::endl;
        std::cout << TreeTools::treeToParenthesis(*currentTree_, true)<<std::endl;
        }
      }        
    else 
      { 
        rearrange_ = true; //Now we rearrange gene trees
        /****************************************************************************
         * noMoreSPR==true, we thus only make NNIs and root changes, 
         * and we rearrange gene trees.
         *****************************************************************************/            
        if ( ( ((!rearrange_)&&(numIterationsWithoutImprovement_>=2*speciesTreeNodeNumber_)) || ((!rearrange_)&&(!optimizeSpeciesTreeTopology_)) || (currentStep_ == 3) ) && (ApplicationTools::getTime() < timeLimit_)  )
          {
          //rearrange_ = true; //Now we rearrange gene trees
          numIterationsWithoutImprovement_ = 0;
          //This first computation is done without rearranging gene trees
          computeSpeciesTreeLikelihood(world_, index_, 
                                       stop_, logL_, 
                                       num0Lineages_, num1Lineages_, num2Lineages_, 
                                       allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                       lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                       rearrange_, server_, 
                                       branchExpectedNumbersOptimization_, genomeMissing_, 
                                       *currentTree_);

          std::cout<<"Species tree likelihood without gene tree rearrangement: "<< - logL_<<std::endl;
          //Now gene trees are really rearranged.
          computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world_, index_, 
                                                                             stop_, logL_, 
                                                                             num0Lineages_, num1Lineages_, num2Lineages_, 
                                                                             allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                                                             lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                                                             rearrange_, server_, 
                                                                             branchExpectedNumbersOptimization_, genomeMissing_, 
                                                                             *currentTree_, bestlogL_);

          std::cout << "\t\tSpecies tree likelihood with gene tree optimization and new branch probabilities: "<< - logL_<<" compared to the former log-likelihood : "<< - bestlogL_<<std::endl;
          numIterationsWithoutImprovement_ = 0;
          bestlogL_ =logL_;
          bestIndex_ = index_;
          bestNum0Lineages_ = num0Lineages_;
          bestNum1Lineages_ = num1Lineages_;
          bestNum2Lineages_ = num2Lineages_;
          }
       
        if ( (currentStep_ == 3) && (ApplicationTools::getTime() < timeLimit_) || (!optimizeSpeciesTreeTopology_))
          {
          if (optimizeSpeciesTreeTopology_) 
            {
            std::cout<<"Reading previous topology likelihoods"<<std::endl;
            inputNNIAndRootLks(NNILks_, rootLks_, params_, suffix_);

            std::cout <<"\tNNIs or Root changes: Number of iterations without improvement : "<<numIterationsWithoutImprovement_<<std::endl;
            localOptimizationWithNNIsAndReRootings(world_, currentTree_, bestTree_, 
                                                   index_, bestIndex_, 
                                                   stop_, timeLimit_, 
                                                   logL_, bestlogL_, 
                                                   num0Lineages_, num1Lineages_, num2Lineages_, 
                                                   bestNum0Lineages_, bestNum1Lineages_, bestNum2Lineages_, 
                                                   allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                                   lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                                   rearrange_, numIterationsWithoutImprovement_, 
                                                   server_, nodeForNNI, nodeForRooting, 
                                                   branchExpectedNumbersOptimization_, genomeMissing_, NNILks_, rootLks_);
            }
          else 
            {
            std::cout <<"\tDuplication and loss rate optimization: Number of iterations without improvement : "<<numIterationsWithoutImprovement_<<std::endl;
            optimizeOnlyDuplicationAndLossRates(world_, currentTree_, bestTree_, 
                                                index_, bestIndex_, 
                                                stop_, timeLimit_, 
                                                logL_, bestlogL_, 
                                                num0Lineages_, num1Lineages_, num2Lineages_, 
                                                bestNum0Lineages_, bestNum1Lineages_, bestNum2Lineages_, 
                                                allNum0Lineages_, allNum1Lineages_, allNum2Lineages_, 
                                                lossExpectedNumbers_, duplicationExpectedNumbers_, 
                                                rearrange_, numIterationsWithoutImprovement_, 
                                                server_, nodeForNNI, nodeForRooting, 
                                                branchExpectedNumbersOptimization_, genomeMissing_);
            }
          if (ApplicationTools::getTime() < timeLimit_) 
            {
            std::cout << "\n\n\t\t\tStep of final optimization over.\n\n"<< std::endl;
            ApplicationTools::displayTime("Execution time so far:");
            currentStep_ = 4;
            }
          else 
            {
            std::cout << "\n\n\t\t\tStep of final optimization is not over yet. The program exits because of the time limit 1.\n\n"<< std::endl;
            ApplicationTools::displayTime("Execution time so far:");
            if (stop_==false)
              {
              stop_ = true;
             broadcast(world_, stop_, server_); 
              broadcast(world_, bestIndex_, server_);                     
              }
            }
          }
        else if (ApplicationTools::getTime() >= timeLimit_)
          {
          std::cout << "\n\n\t\t\tStep of final optimization is not over yet. The program exits because of the time limit 2.\n\n"<< std::endl;
          ApplicationTools::displayTime("Execution time so far:");
          if (stop_==false)
            {
            stop_ = true;
           broadcast(world_, stop_, server_); 
            broadcast(world_, bestIndex_, server_);                    
            }          
          }
      }
    }      
  
  
  
  
  /****************************************************************************
   ****************************************************************************
   * End of the main loop: outputting results
   * If we have to stop the program before the end, we need to:
   * - save the species tree to a given file to start from the good species tree
   * - it would be nice to make client option files for the next run to not recompute starting gene trees. Not for now, though.
   ****************************************************************************
   *****************************************************************************/
  outputNNIAndRootLks(NNILks_, rootLks_, params_, suffix_);

  
  if (ApplicationTools::getTime() >= timeLimit_) 
    {
    
    /****************************************************************************
     * Run not finished yet, outputting information for next run. 
     *****************************************************************************/            
    if (rearrange_)
      broadcast(world_, rearrange_, server_); 
    std::cout <<"\n\n\t\t\tNo time to finish the run. "<<std::endl;
    std::cout <<"\t\tSaving the current species tree for the next run."<<std::endl;
    std::string tempSpTree = ApplicationTools::getStringParameter("output.temporary.tree.file", params_, "CurrentSpeciesTree.tree", "", false, false);
    tempSpTree = tempSpTree + suffix_;
    std::ofstream out (tempSpTree.c_str(), std::ios::out);
    out << TreeTools::treeToParenthesis(*bestTree_, false)<<std::endl;
    out.close();  
    std::cout<<"\n\n\t\t\tAdd:\ninit.species.tree=user\nspecies.tree.file="<<tempSpTree<<"\ncurrent.step="<<currentStep_<<"\nto your options.\n"<<std::endl;
    }
  else 
    {
    
    /****************************************************************************
     * Run finished, outputting end results. 
     *****************************************************************************/            
    //COMPUTING ALRTs
    //In Anisimova and Gascuel, the relevant distribution is a mixture of chi^2_1 and chi^2_0.
    //Here, I am not sure one can do the same. We stick with the classical chi^2_1 distribution, more conservative.
    
    for (int i = 0; i<bestTree_->getNumberOfNodes() ; i++ ) 
      {
      if ((! bestTree_->getNode(i)->isLeaf()) && (bestTree_->getNode(i)->hasFather())) 
        {
        // double proba=(1/2)*(RandomTools::pChisq(2*(NNILks_[i] - bestlogL_), 1)+1); //If one wants to use the mixture.
        double proba=RandomTools::pChisq(2*(NNILks_[i] - bestlogL_), 1);
        //Bonferroni correction
        proba = 1-3*(1-proba);
        if (proba<0) 
          {
          std::cout <<"Negative aLRT!"<<std::endl;
          proba == 0;
          }
        std::cout <<"Branch "<<i<<" second best Lk: "<< NNILks_[i]<< ";Lk difference: "<< NNILks_[i] - bestlogL_ <<"; aLRT: "<<proba<<std::endl;
        bestTree_->getNode(i)->setBranchProperty("ALRT", Number<double>(proba));
        }
      else 
        {
        bestTree_->getNode(i)->setBranchProperty("ALRT", Number<double>(1));
        }
      }
    std::cout <<"\n\n\t\tBest Species Tree found, with node Ids: "<<std::endl;
    std::cout << TreeTools::treeToParenthesis (*bestTree_, true)<<std::endl;
    std::cout <<"\n\n\t\tBest Species Tree found, with aLRTs: "<<std::endl;
    std::cout << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "ALRT")<<std::endl;
    //Here we output the species tree with rates of duplication and loss
    //For duplication rates
    for (int i =0; i<num0Lineages_.size() ; i++ ) 
      {
      bestTree_->getNode(i)->setBranchProperty("DUPLICATIONS", Number<double>( duplicationExpectedNumbers_[i]));
      if (bestTree_->getNode(i)->hasFather()) 
        {
        bestTree_->getNode(i)->setDistanceToFather(duplicationExpectedNumbers_[i]);
        }
      }
    std::string dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", params_, "AllDuplications.tree", "", false, false); 
    dupTree = dupTree + suffix_;
    std::ofstream out (dupTree.c_str(), std::ios::out);

    out << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "DUPLICATIONS")<<std::endl;
    out.close();

    //For loss rates
    for (int i =0; i<num0Lineages_.size() ; i++ ) 
      {
      bestTree_->getNode(i)->setBranchProperty("LOSSES", Number<double>(lossExpectedNumbers_[i]));
      if (bestTree_->getNode(i)->hasFather()) 
        {
        bestTree_->getNode(i)->setDistanceToFather(lossExpectedNumbers_[i]);
        }
      }

    std::string lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", params_, "AllLosses.tree", "", false, false);
    lossTree = lossTree + suffix_;
    out.open (lossTree.c_str(), std::ios::out);
    out << treeToParenthesisWithDoubleNodeValues(*bestTree_, false, "LOSSES")<<std::endl;

    out.close();
    std::string numTree = ApplicationTools::getStringParameter("output.numbered.tree.file", params_, "ServerNumbered.tree", "", false, false);
    numTree = numTree + suffix_;
    out.open (numTree.c_str(), std::ios::out);
    out << TreeTools::treeToParenthesis (*bestTree_, true)<<std::endl;
    out.close();

    //Here we output the species tree with numbers of times 
    //a given number of lineages has been found per branch.
    assignNumLineagesOnSpeciesTree(*bestTree_, 
                                   num0Lineages_, 
                                   num1Lineages_, 
                                   num2Lineages_);

    std::string lineagesTree = ApplicationTools::getStringParameter("output.lineages.tree.file", params_, "lineageNumbers.tree", "", false, false); 
    lineagesTree = lineagesTree + suffix_;
    out.open (lineagesTree.c_str(), std::ios::out);
    out << TreeTools::treeToParenthesis(*bestTree_, false, NUMLINEAGES)<<std::endl;
    out.close();

    std::cout <<"Number of species trees tried : "<<index_<<std::endl;        
    PhylogeneticsApplicationTools::writeTree(*bestTree_, params_, "", "", true, false, false);
    std::cout << "\t\tServer : best found logLikelihood value : "<< - bestlogL_<<std::endl;
    }
}



