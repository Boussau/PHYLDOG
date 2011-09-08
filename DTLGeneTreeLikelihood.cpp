/*
 *  DTLGeneTreeLikelihood.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 28/07/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

#include "DTLGeneTreeLikelihood.h"


using namespace bpp;

using namespace std;


void DTLGeneTreeLikelihood::parseOptions() {
  std::vector <std::string> spNames;
  bool avoidFamily=false;
  
  //Some default options, to avoid warnings:
  if(params_.find(std::string("optimization.verbose")) == params_.end()) {
    params_[ std::string("optimization.verbose")] = "0";
  }
  if(params_.find(std::string("optimization.reparametrization")) == params_.end()) {
    params_[ std::string("optimization.reparametrization")] = "0";
  }
  if(params_.find(std::string("optimization.final")) == params_.end()) {
    params_[ std::string("optimization.final")] = "none";
  }
  if(params_.find(std::string("optimization.max_number_f_eval")) == params_.end()) {
    params_[ std::string("optimization.max_number_f_eval")] = "100000";
  }
  if(params_.find(std::string("optimization.tolerance")) == params_.end()) {
    params_[ std::string("optimization.tolerance")] = "0.00001";
  }
  
  
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
    spTree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
    if (!spTree_->isRooted()) 
      {
      std::cout << "The tree is not rooted, midpoint-rooting it!\n";
      TreeTools::midpointRooting(*spTree_);
      }
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(spTree_->getNumberOfLeaves()));
    spNames=spTree_->getLeavesNames();
    }
  else {
    std::cout << "\n\nA user-defined Species tree needs to be provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file.\n\n" << std::endl;
    exit(-1);
  }
  
  
  ////////////////////////////////////////////////
  //DTL related options:
  ////////////////////////////////////////////////
  delta_=ApplicationTools::getDoubleParameter("duplication.rate", params_, 0.2, "", false, false);
  tau_=ApplicationTools::getDoubleParameter("transfer.rate", params_, 0.1, "", false, false);
  lambda_=ApplicationTools::getDoubleParameter("loss.rate", params_, 0.7, "", false, false);    
  
  
  ////////////////////////////////////////////////
  //Gene family related options:
  ////////////////////////////////////////////////
  //Sequences and model of evolution
  alphabet_ = SequenceApplicationTools::getAlphabet(params_, "", false);
  std::string seqFile = ApplicationTools::getStringParameter("input.sequence.file",params_,"none");
  if(!FileTools::fileExists(seqFile))
    {
    std::cerr << "Error: Sequence file "<< seqFile <<" not found." << std::endl;
    exit(-1);
    }
  VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet_, params_);       
  sites_ = SequenceApplicationTools::getSitesToAnalyse(*allSites, params_);     
  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites_->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites_->getNumberOfSites()));
  
  /****************************************************************************
   //Then we need to get the file containing links between sequences and species.
   *****************************************************************************/
  std::string taxaseqFile = ApplicationTools::getStringParameter("taxaseq.file",params_,"none");
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
  //Printing the contents and building seqSp_ 
  //At the same time, we gather sequences we will have to remove from the 
  //alignment and from the gene tree
  std::vector <std::string> spNamesToTake = spTree_->getLeavesNames(); 
  std::vector <std::string> seqsToRemove;
  for(std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++){
    spNames.push_back(it->first);
    for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
      seqSp_.insert(make_pair(*it2, it->first));
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
  if (seqsToRemove.size()>=sites_->getNumberOfSequences()-1) {
    avoidFamily=true;
    std::cout <<"All or almost all sequences have been removed: avoiding family "<<std::endl;
  }
  
  if (!avoidFamily) {
    //We need to prune the alignment so that they contain
    //only sequences from the species under study.
    for (int j =0 ; j<seqsToRemove.size(); j++) 
      {
      std::vector <std::string> seqNames = sites_->getSequencesNames();
      if ( VectorTools::contains(seqNames, seqsToRemove[j]) ) 
        {
        sites_->deleteSequence(seqsToRemove[j]);
        }
      else 
        std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene alignment."<<std::endl;
      }
    
    /****************************************************************************
     //Then we need to get the substitution model.
     *****************************************************************************/
    
    model_ = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet_, sites_, params_); 
    if (model_->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    if (model_->getNumberOfStates() >= 2 * model_->getAlphabet()->getSize())
      {
      // Markov-modulated Markov model!
      rDist_ = new ConstantDistribution(1.);
      }
    else
      {
      rDist_ = PhylogeneticsApplicationTools::getRateDistribution(params_);
      }
    
    /****************************************************************************
     //Then we need to get the file containing the gene tree.
     *****************************************************************************/
    // Get the initial gene tree
    initTree = ApplicationTools::getStringParameter("init.gene.tree", params_, "user", "", false, false);
    ApplicationTools::displayResult("Input gene tree", initTree);
    try 
    {
    if(initTree == "user")
      {
      std::string geneTreeFile =ApplicationTools::getStringParameter("gene.tree.file",params_,"none");
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
        unrootedGeneTree_ = buildBioNJTree (params_, sites_, model_, rDist_, alphabet_);
        if (geneTree_) 
          {
          delete geneTree_;
          }            
        geneTree_ = unrootedGeneTree_->clone();
        geneTree_->newOutGroup(0); 
        //exit(-1);
        }
      else {
        geneTree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTreeFile));
      }
      if (!geneTree_->isRooted()) 
        {
        unrootedGeneTree_ = geneTree_->clone();
        std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
        geneTree_->newOutGroup(0);
        }
      else 
        {
        unrootedGeneTree_ = geneTree_->clone();
        unrootedGeneTree_->unroot();
        }
      ApplicationTools::displayResult("Gene Tree file", geneTreeFile);
      ApplicationTools::displayResult("Number of leaves", TextTools::toString(geneTree_->getNumberOfLeaves()));
      }
    else if ( (initTree == "bionj") || (initTree == "phyml") ) //build a BioNJ starting tree, and possibly refine it using PhyML algorithm
      {
      unrootedGeneTree_ = buildBioNJTree (params_, sites_, model_, rDist_, alphabet_);
      
      if (initTree == "phyml")//refine the tree using PhyML algorithm (2003)
        { 
          refineGeneTreeUsingSequenceLikelihoodOnly (params_, unrootedGeneTree_, sites_, model_, rDist_, seqFile, alphabet_);
        }
      
      if (geneTree_) 
        {
        delete geneTree_;
        }        
      
      geneTree_ = unrootedGeneTree_->clone(); 
      geneTree_->newOutGroup(0); 
      
      }
    else throw Exception("Unknown init gene tree method. init.gene.tree should be 'user', 'bionj', or 'phyml'.");
    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allows to check file existence before running optimization!
    // PhylogeneticsApplicationTools::writeTree(*geneTree_, params_, "", "", true, false, true);      
    }
    catch (std::exception& e)
    {
    std::cout << e.what() <<"; Unable to get a proper gene tree for family; avoiding this family."<<std::endl;
    avoidFamily=true;
    }
  }
  else {
    throw Exception("\nWe do not rearrange this gene tree.\n");
  }
  
  //This family is phylogenetically informative
  //Pruning sequences from the gene tree
  for (int j =0 ; j<seqsToRemove.size(); j++) 
    {
    std::vector <std::string> leafNames = geneTree_->getLeavesNames();
    if ( VectorTools::contains(leafNames, seqsToRemove[j]) )
      {
      removeLeaf(*geneTree_, seqsToRemove[j]);
      unrootedGeneTree_ = geneTree_->clone();
      if (!geneTree_->isRooted()) {
        std::cout <<"gene tree is not rooted!!! "<< taxaseqFile<<std::endl;
      }
      unrootedGeneTree_->unroot();
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
  
  scoringMethod_ =ApplicationTools::getStringParameter("gene.tree.scoring.method",params_,"ML.root");
  sprLimit_=ApplicationTools::getIntParameter("spr.limit",params_,4);
  
  return;
  
}

/************************************************************************
 * Computes the DTL likelihood of the current unrootedGeneTree_.
 *
 ************************************************************************/

pair<double, string> DTLGeneTreeLikelihood::computeDTLLikelihood () {
  string strTree = geneTreeToParenthesisPlusSpeciesNames (unrootedGeneTree_, seqSp_ );
  vector <string> strTrees;
  strTrees.push_back(strTree);
  vector < pair<long double,string> > treeDTLlogLks=DTLScoringTree_->LL_treewise(strTrees);
  pair<double, string> MLValueAndTree = VectorTools::max(treeDTLlogLks);
  DTLlogL_ = MLValueAndTree.first;
  return MLValueAndTree;
}

/************************************************************************
 * Computes the sequence likelihood of the current unrootedGeneTree_.
 *
 ************************************************************************/

double DTLGeneTreeLikelihood::optimizeSequenceLikelihood () {
  //Compute sequence likelihood and improve branch lengths on the gene tree.
  sequencelogL_ = refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (params_, 
                                                                          unrootedGeneTree_, 
                                                                          sites_, 
                                                                          model_, 
                                                                          rDist_, 
                                                                          "", alphabet_, true);
  return sequencelogL_;
}

/************************************************************************
 * Computes sequence and DTL likelihood of the current unrootedGeneTree_.
 *
 ************************************************************************/

/*pair<double, string> DTLGeneTreeLikelihood::computeSequenceAndDTLLikelihood () {
  string strTree = geneTreeToParenthesisPlusSpeciesNames (unrootedGeneTree_, seqSp_ );
  std::cout << strTree <<std::endl;
  vector <string> strTrees;
  strTrees.push_back(strTree);
  vector < pair<long double,string> > treeDTLlogLks=DTLScoringTree_->LL_treewise(strTrees);
  pair<double, string> MLValueAndTree = VectorTools::max(treeDTLlogLks);
  return MLValueAndTree;
}*/

/************************************************************************
 * Initializes the object.
 *
 ************************************************************************/


void DTLGeneTreeLikelihood::initialize() {
  
  numIterationsWithoutImprovement_ = 0;
  //Annotating the leaves of the gene tree with species names
  annotateGeneTreeWithSpeciesNames (geneTree_, seqSp_ );
  if (unrootedGeneTree_)
    delete unrootedGeneTree_;
  unrootedGeneTree_ = geneTree_->clone();
  unrootedGeneTree_->unroot();
  
  //Computation of the sequence likelihood of the gene tree
  /*    vector <Node *> nodes = unrootedGeneTree_->getNodes();
   for (int i = 0 ; i < nodes.size() ; i++) {
   if (nodes[i]->hasFather()) {
   nodes[i]->setDistanceToFather(0.1);
   } 
   }
  ApplicationTools::displayTime ("Time Before optimization using mapping: ") ;
    refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (params_, 
   unrootedGeneTree_, 
   sites_, 
   model, 
   rDist, 
   "", alphabet, true);
  ApplicationTools::displayTime ("Time After optimization using mapping: ") ;*/
  
 // std::cout << TreeTools::treeToParenthesis(*unrootedGeneTree_, false)<< std::endl;
  
  /*
   nodes = unrootedGeneTree_->getNodes();
   for (int i = 0 ; i < nodes.size() ; i++) {
   if (nodes[i]->hasFather()) {
   nodes[i]->setDistanceToFather(0.1);
   } 
   }
   std::cout << TreeTools::treeToParenthesis(*unrootedGeneTree_, false)<< std::endl;
   
   
   ApplicationTools::displayTime ("Time Before: ") ;
   refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (params_, 
   unrootedGeneTree_, 
   sites_, 
   model, 
   rDist, 
   "", alphabet, false);
   ApplicationTools::displayTime ("Time After: ") ;
   std::cout << TreeTools::treeToParenthesis(*unrootedGeneTree_, false)<< std::endl;
   */
  
  ////////////////////////////////////////////
  //Real Work Starting Now
  ////////////////////////////////////////////
  // We construct a species tree object that will calculate the probability of G-s in the vector trees according to the scoring method scoringMethod_:
  DTLScoringTree_ = new Species_tree(spTree_);

  //init_treewise sets up different pieces of the calculation, and LL_treewise computes likelihoods.
  //    string strTree = geneTreeToParenthesisWithSpeciesNames(unrootedGeneTree_, seqSp_);
  // annotateGeneTreeWithSpeciesNames (unrootedGeneTree_, seqSp_ );
 
  
  //Compute sequence likelihood and improve branch lengths on the gene tree.
  optimizeSequenceLikelihood();
  bestSequencelogL_ = sequencelogL_;
  pair<double, string> MLValueAndTree;
  if (scoringMethod_ == "integral") {
    std::cout << "DTL likelihood summed over roots and reconciliations."<<std::endl;
    DTLScoringTree_->init_treewise(spTree_,delta_,tau_,lambda_,"sum");
    MLValueAndTree = computeDTLLikelihood ();
   // std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
  }
  else if (scoringMethod_ == "ML.root") {
    std::cout << "DTL likelihood summed over reconciliations, using ML root."<<std::endl;
    DTLScoringTree_->init_treewise(spTree_,delta_,tau_,lambda_,"root");
    MLValueAndTree = computeDTLLikelihood ();
   // std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
  }
  else if (scoringMethod_ == "ML.root.reconciliation") {
    std::cout << "DTL likelihood using ML reconciliations and ML root."<<std::endl;
    DTLScoringTree_->init_treewise(spTree_,delta_,tau_,lambda_,"tree");
    MLValueAndTree = computeDTLLikelihood ();
   // std::cout << "Initial DTL likelihood: "<< MLValueAndTree.first << ", for tree: "<< MLValueAndTree.second << std::endl;
  }
  else {
    throw Exception("\t\t\tUnknown DTL scoring method: " + scoringMethod_ + ". It should be 'integral', 'ML.root', or 'ML.root.reconciliation'.\n");
  }
  
  bestDTLlogL_ = DTLlogL_;
  
  std::cout << "\t\t\tInitial gene tree with total lolk=" << bestDTLlogL_ + bestSequencelogL_ <<"; DTL loglk: "<< bestDTLlogL_ << "; sequence logLk: " << bestSequencelogL_ << endl;
  
  return;
}


/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 * Uses only average rates of duplication and loss, not branchwise rates.
 ************************************************************************/
void DTLGeneTreeLikelihood::MLSearch() {
  std::cout <<"\t\t\tStarting MLSearch : current tree : "<< std::endl;
  std::cout<< TreeTools::treeToParenthesis(*geneTree_, true)<< std::endl;
  breadthFirstreNumber (*geneTree_);
  std::vector <int> nodeIdsToRegraft;
  bool betterTree;
  TreeTemplate<Node> *tree = 0;
  double logL = DTLlogL_ + sequencelogL_;
  double bestlogL = logL;
  string bestTree;
  pair<double, string> MLValueAndTree;
  while (numIterationsWithoutImprovement_ < geneTree_->getNumberOfNodes())
    {
    for (int nodeForSPR=geneTree_->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
      {
      buildVectorOfRegraftingNodesLimitedDistance(*geneTree_, nodeForSPR, sprLimit_, nodeIdsToRegraft);
      
      betterTree = false;
      for (int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
        {
        if (tree) 
          {
          delete tree;
          tree = 0;
          }
        tree = geneTree_->clone();
        makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i]);
        if (unrootedGeneTree_) 
          {
          delete unrootedGeneTree_;
          unrootedGeneTree_ = 0;
          }      
        unrootedGeneTree_ = tree->clone();
        unrootedGeneTree_->unroot();
        
        //Compute sequence likelihood and improve branch lengths on the gene tree.
        optimizeSequenceLikelihood();
        
        
        //Compute DTL likelihood
        MLValueAndTree = computeDTLLikelihood ();
        
        //Now we should have a well rooted tree and its DTL likelihood in MLValueAndTree.
        DTLlogL_ = MLValueAndTree.first;
        
        logL = DTLlogL_ + sequencelogL_;
        
        
        
        if (logL+0.01 > bestlogL) 
          {
          betterTree = true;
          bestlogL =logL;
          bestDTLlogL_ = DTLlogL_;
          bestSequencelogL_ = sequencelogL_;
          bestTree = MLValueAndTree.second;
          bestIndex_ = index_;
          std::cout << "\t\t\tSPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
          std::cout << "\t\t\tTotal likelihood: "<<logL <<"; DTL likelihood: "<< DTLlogL_ << ", Sequence likelihood: "<< sequencelogL_ <<", for tree: "<< MLValueAndTree.second << std::endl;
          std::cout << bestTree << std::endl;
          }
        else {
          std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
        }
        }
      if (betterTree) 
        {
        logL = bestlogL; 
        numIterationsWithoutImprovement_ = 0;
        delete geneTree_;
        //geneTree_ = TreeTemplateTools::parenthesisToTree(bestTree); 
        geneTree_ = parenthesisPlusSpeciesNamesToGeneTree  (bestTree);
        if (scoringMethod_ == "integral") geneTree_->newOutGroup(0);
        breadthFirstreNumber (*geneTree_);
        std::cout <<"\t\t\tSPRs: Improvement! : "<<numIterationsWithoutImprovement_<< std::endl;
        std::cout << "\t\t\tNew total Likelihood value "<<logL<< std::endl;
        index_++;  
        bestIndex_ = index_;
        std::cout <<"\t\t\tNumber of species trees tried : "<<index<< std::endl;
        }
      else 
        {
        logL = bestlogL;  
        numIterationsWithoutImprovement_++;
        std::cout <<"\t\t\tSPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement_<< std::endl;
        }
      if (tree) 
        {
        delete tree;
        tree = 0;
        }
      }
    }
}

/************************************************************************
 * Outputs gene tree.
 ************************************************************************/

void DTLGeneTreeLikelihood::printGeneTree()
{
  std::cout << "Best gene tree found with total lolk=" << bestDTLlogL_ + bestSequencelogL_ <<"; DTL loglk: "<< bestDTLlogL_ << "; sequence logLk: " << bestSequencelogL_ << endl;
  string outputTreeFile = ApplicationTools::getStringParameter("output.tree.file", params_, "output.tree", "", false, false);
  //std::cout << TreeTools::treeToParenthesis(*geneTree_, false)<< std::endl;
  std::ofstream out;
  out.open (outputTreeFile.c_str(), std::ios::out);
  Nhx *nhx = new Nhx();
  nhx->write(*geneTree_, out);
  out.close();
  // Write resulting tree:
  //PhylogeneticsApplicationTools::writeTree(*geneTree_, params_);

  return;
}

  
  
