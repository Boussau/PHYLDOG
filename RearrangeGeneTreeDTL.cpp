/*
 *  RearrangeGeneTreeDTL.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 29/06/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

/******************************************************************************/

#include "DTL.h"


void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "RearrangeGeneTreeDTL parameter1_name=parameter1_value parameter2_name=parameter2_value"   ).endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  
  SequenceApplicationTools::printInputAlignmentHelp();
  PhylogeneticsApplicationTools::printInputTreeHelp();
  PhylogeneticsApplicationTools::printSubstitutionModelHelp();
  PhylogeneticsApplicationTools::printRateDistributionHelp();
  PhylogeneticsApplicationTools::printCovarionModelHelp();
  PhylogeneticsApplicationTools::printOptimizationHelp(true, false);
  PhylogeneticsApplicationTools::printOutputTreeHelp();
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
 * This function produces a string version of a gene tree, 
 * replaced by species names.
 **************************************************************************/
string geneTreeToParenthesisWithSpeciesNames (TreeTemplate<Node> * geneTree,
                                          std::map<std::string, std::string > seqSp) {
  
}

string parenthesisWithSpeciesNamesToGeneTree (TreeTemplate<Node> * geneTree,
                                          std::map<std::string, std::string > seqSp) {
  
}


/**************************************************************************
 * This function optimizes a gene tree based on the reconciliation score only.
 * It uses SPRs and NNIs, and calls findMLReconciliationDR to compute the likelihood.
 **************************************************************************/

double refineGeneTreeDTL (Species_tree * scoringTree, 
                          TreeTemplate<Node> *& geneTree, 
                          std::map<std::string, std::string > seqSp,
                          std::map<std::string, int > spID, 
                          pair <double, string> MLValueAndTree,
                          int & MLindex)
{
  TreeTemplate<Node> *tree = 0;
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
  double startingML = findMLReconciliationDR (spTree, 
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
          /*std::cout << "Gene tree SPR: Better candidate tree likelihood : "<<bestlogL<< std::endl;
           std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;*/
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

  return bestlogL;
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
		std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);
    std::cout << "******************************************************************" << std::endl;
    std::cout << "*         DTL Gene Tree Rearrangement Program, version 1.0       *" << std::endl;
    std::cout << "* Authors: G. Szollosi, B. Boussau            Created 29/06/2011 *" << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;
    
    
    ////////////////////////////////////////////////
    //Species tree related options:
    ////////////////////////////////////////////////
    // Get the initial tree
    std::string initTree = ApplicationTools::getStringParameter("init.species.tree", params_, "user", "", false, false);
    ApplicationTools::displayResult("Input species tree", initTree);
    // A given species tree
    if(initTree == "user")
      {
      std::string spTreeFile =ApplicationTools::getStringParameter("species.tree.file",params_,"none");
      if (spTreeFile=="none" )
        {
        std::cout << "\n\nNo Species tree was provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file.\n\n" << std::endl;
        exit(-1);
        }
      ApplicationTools::displayResult("Species Tree file", spTreeFile);
      Newick newick(true);
      spTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
      if (!tree_->isRooted()) 
        {
        std::cout << "The tree is not rooted, midpoint-rooting it!\n";
        TreeTools::midpointRooting(*spTree);
        }
      ApplicationTools::displayResult("Number of leaves", TextTools::toString(spTree->getNumberOfLeaves()));
      spNames=spTree->getLeavesNames();
      }
    else {
      std::cout << "\n\nA user-defined Species tree needs to be provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file.\n\n" << std::endl;
      exit(-1);
    }
    
    
    ////////////////////////////////////////////////
    //DTL related options:
    ////////////////////////////////////////////////
    double delta=ApplicationTools::getDoubleParameter("duplication.rate", params, 0.2, "", false, false);
    double tau=ApplicationTools::getDoubleParameter("transfer.rate", params, 0.1, "", false, false);
    double delta=ApplicationTools::getDoubleParameter("loss.rate", params, 0.7, "", false, false);    
    
    
    ////////////////////////////////////////////////
    //Gene family related options:
    ////////////////////////////////////////////////
    //Sequences and model of evolution
    Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);
    std::string seqFile = ApplicationTools::getStringParameter("input.sequence.file",params,"none");
    if(!FileTools::fileExists(seqFile))
      {
      std::cerr << "Error: Sequence file "<< seqFile <<" not found." << std::endl;
      exit(-1);
      }
    VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);       
    VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, params);     
    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
    
    /****************************************************************************
     //Then we need to get the file containing links between sequences and species.
     *****************************************************************************/
    std::string taxaseqFile = ApplicationTools::getStringParameter("taxaseq.file",params,"none");
    if (taxaseqFile=="none" ){
      std::cout << "\n\nNo taxaseqfile was provided. Cannot compute a reconciliation between a species tree and a gene tree using sequences if the relation between the sequences and the species is not explicit !\n" << std::endl;
      std::cout << "ReconcileDuplications species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
      exit(-1);
    }
    if(!FileTools::fileExists(taxaseqFile))
      {
      std::cerr << "Error: taxaseqfile "<< taxaseqFile <<" not found." << std::endl;
      exit(-1);
      }
    
    //Getting the relations between species and sequence names
    //In this file, the format is expected to be as follows :
    /*
     SpeciesA:sequence1;sequence2
     SpeciesB:sequence5
     SpeciesC:sequence3;sequence4;sequence6
     ...
     */
    //We use a std::map to record the links between species names and sequence names
    //For one species name, we can have several sequence names
    std::map<std::string, std::deque<std::string> > spSeq;
    //We use another std::map to store the link between sequence and species.
    std::map<std::string, std::string> seqSp;
    
    std::ifstream inSpSeq (taxaseqFile.c_str());
    std::string line;
    while(getline(inSpSeq,line)) {
      //We divide the line in 2 : first, the species name, second the sequence names
      StringTokenizer st1 = StringTokenizer::StringTokenizer (line, ":", true);
      //Then we divide the sequence names
      if (st1.numberOfRemainingTokens ()>1) {
        StringTokenizer st2 = StringTokenizer::StringTokenizer (st1.getToken(1), ";", true);
        if (spSeq.find(st1.getToken(0)) == spSeq.end())
          spSeq.insert( make_pair(st1.getToken(0),st2.getTokens()));
        else {
          for (int i = 0 ; i < (st2.getTokens()).size() ; i++)
            spSeq.find(st1.getToken(0))->second.push_back(st2.getTokens()[i]);
        }
      }
    }
    //Printing the contents and building seqSp 
    //At the same time, we gather sequences we will have to remove from the 
    //alignment and from the gene tree
    std::vector <std::string> spNamesToTake = tree->getLeavesNames(); 
    std::vector <std::string> seqsToRemove;
    for(std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++){
      spNames.push_back(it->first);
      for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
        seqSp.insert(make_pair(*it2, it->first));
      }
      if (!VectorTools::contains(spNamesToTake,it->first)) {
        for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
          seqsToRemove.push_back(*it2);
        }
      }
    }
    std::map <std::string, std::string> spSelSeq;
    //If we need to remove all sequences or all sequences except one, 
    //better remove the gene family
    if (seqsToRemove.size()>=sites->getNumberOfSequences()-1) {
      numDeletedFamilies = numDeletedFamilies+1;
      avoidFamily=true;
      std::cout <<"All or almost all sequences have been removed: avoiding family "<<assignedFilenames[i-numDeletedFamilies+1]<<std::endl;
    }
    
    if (!avoidFamily) {
      //We need to prune the alignment so that they contain
      //only sequences from the species under study.
      for (int j =0 ; j<seqsToRemove.size(); j++) 
        {
        std::vector <std::string> seqNames = sites->getSequencesNames();
        if ( VectorTools::contains(seqNames, seqsToRemove[j]) ) 
          {
          sites->deleteSequence(seqsToRemove[j]);
          }
        else 
          std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene alignment."<<std::endl;
        }
      
      /****************************************************************************
       //Then we need to get the substitution model.
       *****************************************************************************/
      
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params); 
      if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
        {
        // Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
        }
      else
        {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
        }
      
      /****************************************************************************
       //Then we need to get the file containing the gene tree.
       *****************************************************************************/
      // Get the initial gene tree
      initTree = ApplicationTools::getStringParameter("init.gene.tree", params, "user", "", false, false);
      ApplicationTools::displayResult("Input gene tree", initTree);
      try 
      {
      if(initTree == "user")
        {
        std::string geneTreeFile =ApplicationTools::getStringParameter("gene.tree.file",params,"none");
        if (geneTreeFile=="none" )
          {
          std::cout << "\n\nNo Gene tree was provided. The option init.gene.tree is set to user (by default), which means that the option gene.tree.file must be filled with the path of a valid tree file. \nIf you do not have a gene tree file, the program can start from a random tree, if you set init.gene.tree at random, or can build a gene tree with BioNJ or a PhyML-like algorithm with options bionj or phyml.\n\n" << std::endl;
          exit(-1);
          }
        Newick newick(true);
        if(!FileTools::fileExists(geneTreeFile))
          {
          std::cerr << "Error: geneTreeFile "<< geneTreeFile <<" not found." << std::endl;
          std::cerr << "Building a bionj tree instead for gene " << geneTreeFile << std::endl;
          unrootedGeneTree = buildBioNJTree (params, sites, model, rDist, file, alphabet);
          if (geneTree) 
            {
            delete geneTree;
            }            
          geneTree = unrootedGeneTree->clone();
          geneTree->newOutGroup(0); 
          //exit(-1);
          }
        else {
          geneTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTreeFile));
        }
        if (!geneTree->isRooted()) 
          {
          unrootedGeneTree = geneTree->clone();
          std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
          geneTree->newOutGroup(0);
          }
        else 
          {
          unrootedGeneTree = geneTree->clone();
          unrootedGeneTree->unroot();
          }
        ApplicationTools::displayResult("Gene Tree file", geneTreeFile);
        ApplicationTools::displayResult("Number of leaves", TextTools::toString(geneTree->getNumberOfLeaves()));
        }
      else if ( (initTree == "bionj") || (initTree == "phyml") ) //build a BioNJ starting tree, and possibly refine it using PhyML algorithm
        {
        unrootedGeneTree = buildBioNJTree (params, sites, model, rDist, file, alphabet);
        
        if (initTree == "phyml")//refine the tree using PhyML algorithm (2003)
          { 
            refineGeneTreeUsingSequenceLikelihoodOnly (params, unrootedGeneTree, sites, model, rDist, file, alphabet);
          }
        
        if (geneTree) 
          {
          delete geneTree;
          }        
        
        geneTree = unrootedGeneTree->clone(); 
        geneTree->newOutGroup(0); 
        
        }
      else throw Exception("Unknown init gene tree method. init.gene.tree should be 'user', 'bionj', or 'phyml'.");
      // Try to write the current tree to file. This will be overwritten by the optimized tree,
      // but allows to check file existence before running optimization!
      PhylogeneticsApplicationTools::writeTree(*tree_, params_, "", "", true, false, true);      
      }
      catch (std::exception& e)
      {
      std::cout << e.what() <<"; Unable to get a proper gene tree for family "<<file<<"; avoiding this family."<<std::endl;
      numDeletedFamilies = numDeletedFamilies+1;
      avoidFamily=true;
      }
    }
    else {
      delete sites;
      delete alphabet;
      throw Exception("\nWe do not rearrange this gene tree.\n");
    }

    //This family is phylogenetically informative
    //Pruning sequences from the gene tree
    for (int j =0 ; j<seqsToRemove.size(); j++) 
      {
      std::vector <std::string> leafNames = geneTree->getLeavesNames();
      if ( VectorTools::contains(leafNames, seqsToRemove[j]) )
        {
        removeLeaf(*geneTree, seqsToRemove[j]);
        unrootedGeneTree = geneTree->clone();
        if (!geneTree->isRooted()) {
          std::cout <<"gene tree is not rooted!!! "<< taxaseqFile<<std::endl;
        }
        unrootedGeneTree->unroot();
        }
      else 
        std::cout<<"Sequence "<<seqsToRemove[j] <<" is not present in the gene tree."<<std::endl;
      }
    
    ////////////////////////////////////////////
    // How do we score gene trees?
    // Summing over roots and reconciliations? "integral"
    // Summing over reconciliations only and maximizing over roots? "ML.root"
    // Maximizing over both roots and reconciliations? "ML.root.reconciliation"
    ////////////////////////////////////////////
    
    std::string scoringMethod =ApplicationTools::getStringParameter("gene.tree.scoring.method",params_,"ML.root");

    ////////////////////////////////////////////
    //Real Work Starting Now
    ////////////////////////////////////////////
    // We construct a species tree object that will calculate the probability of G-s in the vector trees according to the scoring method scoringMethod:
    Species_tree * scoringTree = new Species_tree(spTree);
    //init_treewise sets up different pieces of the calculation, and LL_treewise computes likelihoods.
    vector < pair<double,string> > treeDTLlogLks;
    string strTree = TreeTemplateTools::treeToParenthesis(*unrootedGeneTree);
    pair<double, string> MLValueAndTree;
    if (scoringMethod == "integral") {
      std::cout << "DTL likelihood summed over roots and reconciliations."<<std::endl;
      scoringTree->init_treewise(spTree,delta,tau,lambda,"sum");
      treeDTLlogLks=scoringTree->LL_treewise(strTree);
      MLValueAndTree = VectorTools::max(treeDTLlogLks);
      std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
    }
    else if (scoringMethod == "ML.root") {
      std::cout << "DTL likelihood summed over reconciliations, using ML root."<<std::endl;
      scoringTree->init_treewise(spTree,delta,tau,lambda,"root");
      treeDTLlogLks=scoringTree->LL_treewise(strTree);
      MLValueAndTree = VectorTools::max(treeDTLlogLks);
      std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
    }
    else if (scoringMethod == "ML.root.reconciliation") {
      std::cout << "DTL likelihood using ML reconciliations and ML root."<<std::endl;
      scoringTree->init_treewise(spTree,delta,tau,lambda,"tree");
      treeDTLlogLks=scoringTree->LL_treewise(strTree);
      MLValueAndTree = VectorTools::max(treeDTLlogLks);
      std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
    }
    else {
      throw Exception("Unknown DTL scoring method: "<< scoringMethod <<". It should be 'integral', 'ML.root', or 'ML.root.reconciliation'.\n");
    }
    
    //Now we optimize the gene tree.
    
    refineGeneTreeDTL (scoringTree, 
                       unrootedGeneTree, 
                       seqSp,
                       spID, 
                       MLValueAndTree,
                       MLindex);
      
    
    
    
    
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