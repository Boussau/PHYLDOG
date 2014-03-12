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

#include <iostream>
#include "GeneTreeLikelihood.h"

using namespace bpp;



class geneTreeLikelihoodException: public exception
{
  virtual const char* what() const throw()
  {
    return "Can't create a GeneTreeLikelihood object, avoiding family.";
  }
} geneTreeLikelihoodEx;




GeneTreeLikelihood::GeneTreeLikelihood():
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(), spId_(), heuristicsLevel_(0) {
    totalIterations_ = 0;
    counter_ = 0;
}


GeneTreeLikelihood::GeneTreeLikelihood(std::string file , 
		   map<string, string> params, 
		   TreeTemplate<Node> & spTree ) throw (exception):
    levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(), spId_(), 
    params_(params), heuristicsLevel_(0), considerSequenceLikelihood_(true) {
    totalIterations_ = 0;
    counter_ = 0;
    spTree_ = spTree.clone();
    spId_ = computeSpeciesNamesToIdsMap(*spTree_);
    SubstitutionModel*    model    = 0;
    DiscreteDistribution* rDist    = 0;
    TreeTemplate<Node> * unrootedGeneTree = 0;
    bool avoidFamily = false;
    std::vector <std::string> spNames;

            //Sequences and model of evolution
      Alphabet *alphabet = SequenceApplicationTools::getAlphabet(params_, "", false);
      std::string seqFile = ApplicationTools::getStringParameter("input.sequence.file",params_,"none");
      if(!FileTools::fileExists(seqFile))
      {
          std::cerr << "Error: Sequence file "<< seqFile <<" not found." << std::endl;
          MPI::COMM_WORLD.Abort(1);
          exit(-1);
      }
      VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params_);       
      
      unsigned int numSites = allSites->getNumberOfSites();
      ApplicationTools::displayResult("Number of sequences", TextTools::toString(allSites->getNumberOfSequences()));
      ApplicationTools::displayResult("Number of sites", TextTools::toString(numSites));
    
    std::vector <std::string> seqsToRemove;
    VectorSiteContainer * sites;
    
    if (numSites == 0 ) {
      std::cout<<"WARNING: Discarding a family whose alignment is 0 site long: "<< seqFile <<std::endl;
      avoidFamily = true;
      delete allSites;   
    }
    if (!avoidFamily) {
      unsigned int minPercentSequence = ApplicationTools::getIntParameter("sequence.removal.threshold",params_,0);
      unsigned int threshold = (int) ((double)minPercentSequence * (double)numSites / 100 );
      
      
      if (minPercentSequence > 0) {
        for ( int j = allSites->getNumberOfSequences()-1 ; j >= 0 ; j--) {
          if (SequenceTools::getNumberOfCompleteSites(allSites->getSequence(j) ) < threshold ) {
            ApplicationTools::displayResult("Removing a short sequence:", allSites->getSequence(j).getName()  );
            // allSites->deleteSequence(i);
            seqsToRemove.push_back(allSites->getSequence(j).getName());
          }
        }
      }
      
      for (unsigned int j =0 ; j<seqsToRemove.size(); j++) 
      {
        std::vector <std::string> seqNames = allSites->getSequencesNames();
        if ( VectorTools::contains(seqNames, seqsToRemove[j]) ) 
        {
          allSites->deleteSequence(seqsToRemove[j]);
        }
        else 
          std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene alignment."<<std::endl;
      }
      
      
      ApplicationTools::displayResult("# sequences post size-based removal:", TextTools::toString(allSites->getNumberOfSequences()));
      seqsToRemove.clear();
      
      if (allSites->getNumberOfSequences() <= 1 ) {
        std::cout << "Only one sequence left: discarding gene family "<< file<<std::endl;
        avoidFamily = true;
      }
      sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, params_);     
      delete allSites;   
    }
      
      //method to optimize the gene tree root; only useful if heuristics.level!=0.
      bool rootOptimization = false;
      
      /****************************************************************************
       //Then we need to get the file containing links between sequences and species.
       *****************************************************************************/
      std::string taxaseqFile = ApplicationTools::getStringParameter("taxaseq.file",params_,"none");
      if (!avoidFamily) {
          if (taxaseqFile=="none" ){
              std::cout << "\n\nNo taxaseqfile was provided. Cannot compute a reconciliation between a species tree and a gene tree using sequences if the relation between the sequences and the species is not explicit !\n" << std::endl;
              std::cout << "phyldog species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
              MPI::COMM_WORLD.Abort(1);
              exit(-1);
          }
          if(!FileTools::fileExists(taxaseqFile))
          {
              std::cerr << "Error: taxaseqfile "<< taxaseqFile <<" not found." << std::endl;
              MPI::COMM_WORLD.Abort(1);
              exit(-1);
          }
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
      if (!avoidFamily) {
          while(getline(inSpSeq,line)) {
              //We divide the line in 2 : first, the species name, second the sequence names
              StringTokenizer st1 (line, ":", true);
              //Then we divide the sequence names
              if (st1.numberOfRemainingTokens ()>1) {
                  StringTokenizer st2(st1.getToken(1), ";", true);
                  if (spSeq.find(st1.getToken(0)) == spSeq.end())
                      spSeq.insert( make_pair(st1.getToken(0),st2.getTokens()));
                  else {
                      for (unsigned int j = 0 ; j < (st2.getTokens()).size() ; j++)
                          spSeq.find(st1.getToken(0))->second.push_back(st2.getTokens()[j]);
                  }
              }
          }
      }
      //Printing the contents and building seqSp_ 
      //At the same time, we gather sequences we will have to remove from the 
      //alignment and from the gene tree
      std::vector <std::string> spNamesToTake = spTree_->getLeavesNames(); 
      if (!avoidFamily) {
          for(std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++){
              spNames.push_back(it->first);
              for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
                  seqSp_.insert(make_pair(*it2, it->first));
              }
              if (!VectorTools::contains(spNamesToTake,it->first)) {
                  for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
                      std::cout<<"Removing sequence of species not considered"<<std::endl;
                      seqsToRemove.push_back(*it2);
                  }
              }
          }
      }
      std::map <std::string, std::string> spSelSeq;
      
      if (!avoidFamily) {
          //If we need to remove all sequences or all sequences except one, 
          //better remove the gene family
          if (seqsToRemove.size()>=sites->getNumberOfSequences()-1) {
              avoidFamily=true;
              std::cout <<"All or almost all sequences have been removed: avoiding family "<< file <<std::endl;
          }
      }
      
      if (!avoidFamily) {
          //We need to prune the alignment so that they contain
          //only sequences from the species under study.
          for (unsigned int j =0 ; j<seqsToRemove.size(); j++) 
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

          model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, 00, sites, params_); 

          if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

          if (model->getNumberOfStates() > model->getAlphabet()->getSize())
          {
              // Markov-modulated Markov model!
              rDist = new ConstantDistribution(1.);
          }
          else
          {
              rDist = PhylogeneticsApplicationTools::getRateDistribution(params_);
          }
          
          /****************************************************************************
           * Then we need to get the file containing the gene tree, 
           * or build the gene tree.
           *****************************************************************************/
          // Get the initial gene tree
          string initTree = ApplicationTools::getStringParameter("init.gene.tree", params_, "user", "", false, false);
          ApplicationTools::displayResult("Input gene tree", initTree);
          try 
          {
              if(initTree == "user")
              {
                  std::string geneTree_File =ApplicationTools::getStringParameter("gene.tree.file",params_,"none");
                  if (geneTree_File=="none" )
                  {
                      std::cout << "\n\nNo Gene tree was provided. The option init.gene.tree is set to user (by default), which means that the option gene.tree.file must be filled with the path of a valid tree file. \nIf you do not have a gene tree file, the program can start from a random tree, if you set init.gene.tree at random, or can build a gene tree with BioNJ or a PhyML-like algorithm with options bionj or phyml.\n\n" << std::endl;
                      MPI::COMM_WORLD.Abort(1);
                      exit(-1);
                  }
                  Newick newick(true);
                  if(!FileTools::fileExists(geneTree_File))
                  {
                      std::cerr << "Error: geneTree_File "<< geneTree_File <<" not found." << std::endl;
                      std::cerr << "Building a bionj tree instead for gene " << geneTree_File << std::endl;
                      unrootedGeneTree = buildBioNJTree (params_, sites, model, rDist, alphabet);
                      if (rootedTree_) 
                      {
                          delete rootedTree_;
                          rootedTree_ =0;
                      }            
                      rootedTree_ = unrootedGeneTree->clone();
                      rootedTree_->newOutGroup(0); 
                      //exit(-1);
                  }
                  else {
                      if (rootedTree_) 
                      {
                          delete rootedTree_;
                          rootedTree_ =0;
                      }            
                      rootedTree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTree_File));
                  }
                  if (!rootedTree_->isRooted()) 
                  {
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = rootedTree_->clone();
                      //std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
                      rootedTree_->newOutGroup(0);
                  }
                  else 
                  {
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = rootedTree_->clone();
                      unrootedGeneTree->unroot();
                  }
                  ApplicationTools::displayResult("Gene Tree file", geneTree_File);
                  ApplicationTools::displayResult("Number of leaves", TextTools::toString(rootedTree_->getNumberOfLeaves()));
              }
              else if ( (initTree == "bionj") || (initTree == "phyml") ) //build a BioNJ starting tree, and possibly refine it using PhyML algorithm
              {
                  unrootedGeneTree = buildBioNJTree (params_, sites, model, rDist, alphabet);

                  if (initTree == "phyml")//refine the tree using PhyML algorithm (2003)
                  { 
                      refineGeneTreeUsingSequenceLikelihoodOnly (params_, unrootedGeneTree, sites, model, rDist, file, alphabet);
                  }
                  if (rootedTree_) 
                  {
                      delete rootedTree_;
                      rootedTree_ = 0;
                  }        
                  rootedTree_ = unrootedGeneTree->clone(); 
          breadthFirstreNumber (*rootedTree_);
        //  std::cout << " Problem tree? : "<<TreeTemplateTools::treeToParenthesis(*rootedTree_, true) << std::endl;
          
          rootedTree_->newOutGroup( rootedTree_->getLeavesId()[0] ); 
              }
              else throw Exception("Unknown init gene tree method. init.gene.tree should be 'user', 'bionj', or 'phyml'.");
          }
          catch (std::exception& e)
          {
              std::cout << e.what() <<"; Unable to get a proper gene tree for family "<<file<<"; avoiding this family."<<std::endl;
              avoidFamily=true;
          }
      }

      if (!avoidFamily) 
      { //This family is phylogenetically informative

          //Going through the gene tree to see if leaves have branches that are too long.
          std::vector <Node*> leaves = rootedTree_->getLeaves();
         // std::cout << "leaves.size(): "<<leaves.size() <<std::endl;
          for (unsigned int j = 0 ; j < leaves.size() ; j++) {
              if ( leaves[j] -> hasFather() && leaves[j]->getDistanceToFather() >= 2.0 ) {
                  std::cout << "WARNING: Removing sequence "<< leaves[j]->getName() <<" from family "<<file<< " because its branch is unreasonably long (>=2.0)."<<std::endl;
                  seqsToRemove.push_back( leaves[j]->getName() );
                  //removing the corresponding sequence, if present
                  if (sites->hasSequence(leaves[j]->getName() ) )
                      sites->deleteSequence( leaves[j]->getName() );
              }
          }

          if (sites->getNumberOfSequences() >1) {
              //Pruning sequences from the gene tree
              for (unsigned int j =0 ; j<seqsToRemove.size(); j++) 
              {
                  std::vector <std::string> leafNames = rootedTree_->getLeavesNames();
                  if ( VectorTools::contains(leafNames, seqsToRemove[j]) )
                  {
                      removeLeaf(*rootedTree_, seqsToRemove[j]);
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = rootedTree_->clone();
                      if (!rootedTree_->isRooted()) {
                          std::cout <<"gene tree is not rooted!!! "<< taxaseqFile<<std::endl;
                      }
                      unrootedGeneTree->unroot();
                  }
                  else 
                      std::cout<<"Sequence "<<seqsToRemove[j] <<" is not present in the gene tree."<<std::endl;
              }
          }
          //If we have only one sequence in the end, we do not make a tree
          else {
              avoidFamily=true;
              std::cout <<"All or almost all sequences have been removed: avoiding family "<< file <<std::endl;
          }
      }
      if (!avoidFamily) 
      { //This family is phylogenetically informative

      
          /************************************************************************************************************/
          /********************************************COMPUTING LIKELIHOOD********************************************/
          /************************************************************************************************************/
          bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", params_, true, "", false, false);
          if(!computeLikelihood)
          {
              if (alphabet)
                  delete alphabet;
              if (sites)
                  delete sites;
              if (spTree_)
                  delete spTree_;
              std::cout << "PHYLDOG's done. Bye." << std::endl;
              MPI::COMM_WORLD.Abort(1);
              exit(-1);
          }
          
          
          /****************************************************************************
           //Then we can change the gene tree so that its topology minimizes the number of duplications and losses.
           *****************************************************************************/
	  /*
          std::string alterStartingTopologyWithDL = ApplicationTools::getStringParameter("DL.starting.gene.tree.optimization", params_, "no", "", true, false);
          
          
          bool DLStartingGeneTree;
          if (alterStartingTopologyWithDL == "yes") 
          {
              DLStartingGeneTree = true;
          }
          else 
          {
              DLStartingGeneTree = false;
          }
          
          if (DLStartingGeneTree)
          {
              //we temporarily build a ReconciliationTreeLikelihood object, 
              //but won't consider sequences, to save computational time
              std::cout << "Changing the starting gene tree to minimize the numbers of duplications and losses"<<std::endl;
              std::set <int> nodesToTryInNNISearch;
              refineGeneTreeDLOnly (tree_, 
                                    geneTree_, 
                                    seqSp,
                                    spId_,
                                    lossExpectedNumbers_, 
                                    duplicationExpectedNumbers_, 
                                    MLindex_, 
                                    allNum0Lineages_[i-numDeletedFamilies_], 
                                    allNum1Lineages_[i-numDeletedFamilies_], 
                                    allNum2Lineages_[i-numDeletedFamilies_], 
                                    nodesToTryInNNISearch);
              if (unrootedGeneTree) {
                  delete unrootedGeneTree;
                  unrootedGeneTree = 0;
              }
              unrootedGeneTree = geneTree_->clone();
              unrootedGeneTree->unroot();
              
              //Refining branch lengths for this altered gene tree.
              refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (params_, unrootedGeneTree, sites, model, rDist, file, alphabet);
          }
          */
	  
	  
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
          std::string startingGeneTreeFile =ApplicationTools::getStringParameter("output.starting.gene.tree.file",params_,"none");
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
	    levaluator_ = new NNIHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, true, true); 

              
          }
          else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
 
      }
      else 
      {
          if (sites)
              delete sites;
          if (alphabet)
              delete alphabet;
	  throw geneTreeLikelihoodEx;
      }

      
  
}

/**
 * @brief Build a new ReconciliationTreeLikelihood object.
 *
 * @param tree The tree to use.
 * @param model The substitution model to use.
 * @param rDist The rate across sites distribution to use.
 * @param spTree The species tree
 * @param rootedTree rooted version of the gene tree
 * @param seqSp link between sequence and species names
 * @param spId link between species name and species ID
 * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
 * @param heuristicsLevel type of heuristics used
 * @param MLindex ML rooting position
 * @param checkRooted Tell if we have to check for the tree to be unrooted.
 * If true, any rooted tree will be unrooted before likelihood computation.
 * @param verbose Should I display some info?
 * @throw Exception if an error occured.
 */
GeneTreeLikelihood::GeneTreeLikelihood(
                   const Tree & tree,
                   SubstitutionModel * model,
                   DiscreteDistribution * rDist,
                   TreeTemplate<Node> & spTree,  
                   TreeTemplate<Node> & rootedTree, 
                   TreeTemplate<Node> & geneTreeWithSpNames,
                   const std::map <std::string, std::string> seqSp,
                   std::map <std::string,int> spId,
                   int speciesIdLimitForRootPosition,
                   int heuristicsLevel,
                   int & MLindex, 
                   bool checkRooted,
                   bool verbose,
                   bool rootOptimization, 
                   bool considerSequenceLikelihood, 
                   unsigned int sprLimit)
throw (Exception):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_(seqSp), spId_(spId), heuristicsLevel_(0)
{
    
    levaluator_ = new LikelihoodEvaluator(tree, model, rDist, checkRooted, verbose); 
    spTree_ = spTree.clone();
    rootedTree_ = rootedTree.clone();
    geneTreeWithSpNames_ = geneTreeWithSpNames.clone();
    scenarioLikelihood_ = UNLIKELY;
    // _sequenceLikelihood = UNLIKELY;
    MLindex_ = MLindex;
    rootOptimization_ = rootOptimization; 
    tentativeMLindex_ = MLindex;
    totalIterations_ = 0;
    counter_ = 0;
    _speciesIdLimitForRootPosition_ = speciesIdLimitForRootPosition;
    heuristicsLevel_ = heuristicsLevel;
    optimizeSequenceLikelihood_ = true;
    optimizeReconciliationLikelihood_ = true;
    considerSequenceLikelihood_ = considerSequenceLikelihood;
    sprLimit_ = sprLimit;
    // listOfPreviousRoots_ = new std::vector <int> ();
}


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
 * @param heuristicsLevel type of heuristics used
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
                   int heuristicsLevel,
                   int & MLindex, 
                   bool checkRooted,
                   bool verbose, 
                   bool rootOptimization, 
                   bool considerSequenceLikelihood, 
                   unsigned int sprLimit)
throw (Exception):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_ (seqSp), spId_(spId), heuristicsLevel_(0)
{
    levaluator_ = new LikelihoodEvaluator(tree, data, model, rDist, checkRooted, verbose);
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
    heuristicsLevel_ = heuristicsLevel;
    optimizeSequenceLikelihood_ = true;
    optimizeReconciliationLikelihood_ = true;
    considerSequenceLikelihood_ = considerSequenceLikelihood;
    sprLimit_ = sprLimit;
}


/**
 * @brief Copy constructor.
 */ 
GeneTreeLikelihood::GeneTreeLikelihood(const GeneTreeLikelihood & lik):
levaluator_(00), spTree_(00), rootedTree_(00), geneTreeWithSpNames_(00), seqSp_ (lik.seqSp_), spId_(lik.spId_), heuristicsLevel_(0)
{
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
    heuristicsLevel_ = lik.heuristicsLevel_;
    nodesToTryInNNISearch_ = lik.nodesToTryInNNISearch_;
    tentativeNodesToTryInNNISearch_ = lik.tentativeNodesToTryInNNISearch_;
    optimizeSequenceLikelihood_ = lik.optimizeSequenceLikelihood_;
    optimizeReconciliationLikelihood_ = lik.optimizeReconciliationLikelihood_ ;
    considerSequenceLikelihood_ = lik.considerSequenceLikelihood_;
    sprLimit_ = lik.sprLimit_;
}

GeneTreeLikelihood & GeneTreeLikelihood::operator=(const GeneTreeLikelihood & lik)
{
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
    heuristicsLevel_ = lik.heuristicsLevel_;
    nodesToTryInNNISearch_ = lik.nodesToTryInNNISearch_;
    tentativeNodesToTryInNNISearch_ = lik.tentativeNodesToTryInNNISearch_;
    optimizeSequenceLikelihood_ = lik.optimizeSequenceLikelihood_;
    optimizeReconciliationLikelihood_ = lik.optimizeReconciliationLikelihood_ ;
    considerSequenceLikelihood_ = lik.considerSequenceLikelihood_;
    sprLimit_ = lik.sprLimit_;
    return *this;
}




