//           ClientComputingGeneLikelihoods.cpp
//  Mon March 03 11:28:42 2014
//  Copyright  2014  boussau
//  <user@host>

#include "ClientComputingGeneLikelihoods.h"

namespace mpi = boost::mpi;


/*

const mpi::communicator & world_, 
                               unsigned int & server_, 
                               unsigned int & rank_,
                               std::vector<std::string> & assignedFilenames_, 
                               std::map<std::string, std::string> & params_, 
                               unsigned int & numDeletedFamilies_, TreeTemplate<Node> * geneTree_, 
                               TreeTemplate<Node>* &tree, std::vector <int> & num0Lineages_, 
                               std::vector <std::vector<int> > & allNum0Lineages_, 
                               std::vector <std::vector<int> > & allNum1Lineages_, 
                               std::vector <std::vector<int> > & allNum2Lineages_, 
                               std::vector <double> & lossExpectedNumbers_, 
                               std::vector <double> & duplicationExpectedNumbers_,
                               std::vector < std::vector < std::vector < std::vector <unsigned int> > > > coalCounts_,
                               std::vector <double> & coalBls_,
                               std::vector <std::vector<unsigned int> > & allNum12Lineages_, 
                               std::vector <std::vector<unsigned int> > & allNum22Lineages_, 
                               std::map <std::string, int> & spId_, 
                               int speciesIdLimitForRootPosition_, 
                               int heuristicsLevel_, int MLindex_, 
                               std::vector <double> & allLogLs_, 
                               std::vector <GeneTreeLikelihood *> & treeLikelihoods_,
                               std::vector <std::map<std::string, std::string> > & allParams_, 
                               std::vector <Alphabet *> & allAlphabets_, 
                               std::vector <VectorSiteContainer *> & allDatasets_, 
                               std::vector <SubstitutionModel *> & allModels_, 
                               std::vector <DiscreteDistribution *> & allDistributions_, 
                               std::vector <TreeTemplate<Node> *> & allGeneTrees_, 
                               std::vector <TreeTemplate<Node> *> & allUnrootedGeneTrees_, 
				std::vector < std::map <std::string, std::string> >& allSeqSps_, 
				std::vector<unsigned int>& allSprLimitGeneTree_,
                               string &reconciliationModel_*/

                               
                               
/******************************************************************************/
// Function to initialize various parameters in the client, using the options.
/******************************************************************************/
void ClientComputingGeneLikelihoods::parseOptions()  {
              /****************************************************************************
             * First communications between the server and the clients.
             *****************************************************************************/      
//             bool optimizeSpeciesTreeTopology;
//             std::vector <std::string> assignedFilenames;
//             std::vector <double> lossExpectedNumbers;
//             std::vector <double> duplicationExpectedNumbers;
//             std::vector <double> backupLossExpectedNumbers;
//             std::vector <double> backupDuplicationExpectedNumbers;
//             std::vector <int> num0Lineages;
//             std::vector <int> num1Lineages; 
//             std::vector <int> num2Lineages;
//             std::vector <std::vector<int> > allNum0Lineages;
//             std::vector <std::vector<int> > allNum1Lineages;
//             std::vector <std::vector<int> > allNum2Lineages;
//             std::vector < std::vector < std::vector < std::vector< unsigned int > > > > coalCounts;
//             std::vector < double > coalBls;
//             std::vector <unsigned int> num22Lineages; 
//             std::vector <unsigned int> num12Lineages; 
//             std::vector < std::vector <unsigned int> > allNum22Lineages; 
//             std::vector < std::vector <unsigned int> > allNum12Lineages; 
// 
//             std::string currentSpeciesTree_;
//             std::string reconciliationModel;
//             
             unsigned int assignedNumberOfGenes;
             std::vector <unsigned int> numbersOfGenesPerClient;
             std::vector <std::vector<std::string> > listOfOptionsPerClient;
             int SpeciesNodeNumber;

            firstCommunicationsServerClient (world_, server_, numbersOfGenesPerClient, assignedNumberOfGenes,
                                             assignedFilenames_, listOfOptionsPerClient, 
                                             optimizeSpeciesTreeTopology_, SpeciesNodeNumber, 
                                             lossExpectedNumbers_, duplicationExpectedNumbers_, num0Lineages_, 
                                             num1Lineages_, num2Lineages_, 
                                             num12Lineages_, num22Lineages_, coalBls_, currentSpeciesTree_);
            
            /****************************************************************************
             * Various initializations.
             *****************************************************************************/  
           
            //First we read the species tree from the char[] sent by the server
            tree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
 
            resetLossesAndDuplications(*tree_, lossExpectedNumbers_, duplicationExpectedNumbers_);
            //To make the correspondance between species name and id:
            spId_ = computeSpeciesNamesToIdsMap(*tree_);
            
            /****************************************************************************
             * HISTORY: Meaning of heuristicsLevel:
             * 0: exact double-recursive algorithm. 
             All possible root likelihoods are computed with only 2 tree traversals (default).
             * 1: fastest heuristics : only a few nodes are tried for the roots of the gene trees 
             (the number of these nodes tried depends upon speciesIdLimitForRootPosition), 
             and for each root tried, the events are re-computed only for a subset of the tree.
             * 2: All roots are tried, and for each root tried, 
             the events are re-computed only for a subset of the tree.
             * 3: All roots are tried, and for each root tried, 
             the events are re-computed for all nodes of the tree 
             (which should be useless unless there is a bug in the selection of the subset of the nodes.
             * WARNING: options 1, 2, 3 are probably buggy now, and should not be used.
             *****************************************************************************/
            
            /****************************************************************************
             * Gene family parsing and first likelihood computation.
             *****************************************************************************/            
            string toPrint = "";
            std::cout <<"Client  of rank "<<rank_ <<" with PID "<< TextTools::toString((int)getpid()) <<" is in charge of " << assignedFilenames_.size()<<" gene families:"<<std::endl;
            for (unsigned int i = 0 ; i< assignedFilenames_.size() ; i++ ) {
                toPrint = toPrint + assignedFilenames_[i] + "\t";
            }
            std::cout << toPrint <<std::endl;
            //Gets gene family-specific options, and computes the likelihood
            parseAssignedGeneFamilies();
            std::vector <std::string> t;  
            allParamsBackup_ = allParams_;
            resetGeneTrees_ = ApplicationTools::getBooleanParameter("reset.gene.trees",params_,true);
            currentStep_ = ApplicationTools::getIntParameter("current.step",params_,0);

            for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
            {
                reconciledTrees_.push_back(t);
                duplicationTrees_.push_back(t);
                lossTrees_.push_back(t);
                //This is to avoid optimizing gene tree parameters in the first steps of the program, 
                //if we optimize the species tree topology.
                if (optimizeSpeciesTreeTopology_ && currentStep_ < 3) 
                { 
                    if (ApplicationTools::getBooleanParameter("optimization.topology", allParams_[i], false, "", true, false))
                    {
                        allParams_[i][ std::string("optimization.topology")] = "false";
                    }
                    allParams_[i][ std::string("optimization")] = "None"; //Quite extreme, but the sequence likelihood has no impact on the reconciliation !
                    treeLikelihoods_[i]->OptimizeSequenceLikelihood(false);
                }
            }

          //  bool firstTimeImprovingGeneTrees = false; //When for the first time we optimize gene trees, we set it at true
            if (optimizeSpeciesTreeTopology_ && currentStep_ < 3)
            {//At the beginning, we do not record the gene trees.
                recordGeneTrees_ = false;
            }
            else {
                recordGeneTrees_ = true;
               // firstTimeImprovingGeneTrees=true;
            }
            startRecordingTreesFrom_ = 0; //This int is incremented until the gene trees start to be backed-up, when we start the second phase of the algorithm.
          //  MPI_Barrier(world);
            broadcast(world_, stop_, server_);
            broadcastsAllInformationButStop(world_, server_, rearrange_, 
                                            lossExpectedNumbers_, 
                                            duplicationExpectedNumbers_, 
                                            coalBls_,
                                            currentSpeciesTree_,
                                            currentStep_, 
                                            reconciliationModel_
                                            );
            if (assignedFilenames_.size()-numDeletedFamilies_ > 0) 
            {
                if (tree_) delete tree_;
                tree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
                spId_ = computeSpeciesNamesToIdsMap(*tree_);
            }
            startRecordingTreesFrom_ = 1;
            for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
            {
                treeLikelihoods_[i]->setSpTree(*tree_);
                treeLikelihoods_[i]->setSpId(spId_);
                if (reconciliationModel_ == "DL") {
                    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->setExpectedNumbers(duplicationExpectedNumbers_, lossExpectedNumbers_);
                }
                else if (reconciliationModel_ == "COAL") {
                    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->setCoalBranchLengths(coalBls_);
                }
            }

  
}

                               
                               
/******************************************************************************/
// This function parses a vector of gene family files, 
// discards the ones that do not pass certain criteria, 
// and initializes the others. Used by clients.
/******************************************************************************/
void ClientComputingGeneLikelihoods::parseAssignedGeneFamilies() 
{
  bool avoidFamily;
  std::string initTree;
  std::map<std::string, std::string> famSpecificParams;
  std::vector <std::string> spNames;

  //Here we are going to get all necessary information regarding all gene families the client is in charge of.
  for (unsigned int i = 0 ; i< assignedFilenames_.size() ; i++) 
  { //For each file
      double startingFamilyTime = ApplicationTools::getTime();  
      std::cout <<"Examining family "<<assignedFilenames_[i]<<std::endl;
      avoidFamily = false;
      std::string file =assignedFilenames_[i];
      TreeTemplate<Node> * unrootedGeneTree = 0;
      SubstitutionModel*    model    = 0;
      DiscreteDistribution* rDist    = 0;
      if(!FileTools::fileExists(file))
      {
          std::cerr << "Error: Parameter file "<< file <<" not found." << std::endl;
          fflush(0);
          MPI::COMM_WORLD.Abort(1);
          exit(-1);
      }
      else
      {
          famSpecificParams = AttributesTools::getAttributesMapFromFile(file, "=");
          AttributesTools::resolveVariables(famSpecificParams);
      }
      
      AttributesTools::actualizeAttributesMap(params_, famSpecificParams);
      //COAL or DL?
      reconciliationModel_ = ApplicationTools::getStringParameter("reconciliation.model", params_, "DL", "", true, false);

      //Sequences and model of evolution
      Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params_, "", false);
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
          numDeletedFamilies_ = numDeletedFamilies_+1;
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
        numDeletedFamilies_ = numDeletedFamilies_+1;
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
      //We use another std::map to store the link between sequence and species.
      std::map<std::string, std::string> seqSp;
      
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
      //Printing the contents and building seqSp 
      //At the same time, we gather sequences we will have to remove from the 
      //alignment and from the gene tree
      std::vector <std::string> spNamesToTake = tree_->getLeavesNames(); 
      if (!avoidFamily) {
          for(std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++){
              spNames.push_back(it->first);
              for( std::deque<std::string >::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){
                  seqSp.insert(make_pair(*it2, it->first));
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
              numDeletedFamilies_ = numDeletedFamilies_+1;
              avoidFamily=true;
              std::cout <<"All or almost all sequences have been removed: avoiding family "<<assignedFilenames_[i-numDeletedFamilies_+1]<<std::endl;
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
          initTree = ApplicationTools::getStringParameter("init.gene.tree", params_, "user", "", false, false);
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
                      if (geneTree_) 
                      {
                          delete geneTree_;
                          geneTree_ =0;
                      }            
                      geneTree_ = unrootedGeneTree->clone();
                      geneTree_->newOutGroup(0); 
                      //exit(-1);
                  }
                  else {
                      if (geneTree_) 
                      {
                          delete geneTree_;
                          geneTree_ =0;
                      }            
                      geneTree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTree_File));
                  }
                  if (!geneTree_->isRooted()) 
                  {
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = geneTree_->clone();
                      //std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
                      geneTree_->newOutGroup(0);
                  }
                  else 
                  {
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = geneTree_->clone();
                      unrootedGeneTree->unroot();
                  }
                  ApplicationTools::displayResult("Gene Tree file", geneTree_File);
                  ApplicationTools::displayResult("Number of leaves", TextTools::toString(geneTree_->getNumberOfLeaves()));
              }
              else if ( (initTree == "bionj") || (initTree == "phyml") ) //build a BioNJ starting tree, and possibly refine it using PhyML algorithm
              {
                  unrootedGeneTree = buildBioNJTree (params_, sites, model, rDist, alphabet);

                  if (initTree == "phyml")//refine the tree using PhyML algorithm (2003)
                  { 
                      refineGeneTreeUsingSequenceLikelihoodOnly (params_, unrootedGeneTree, sites, model, rDist, file, alphabet);
                  }
                  if (geneTree_) 
                  {
                      delete geneTree_;
                      geneTree_ = 0;
                  }        
                  geneTree_ = unrootedGeneTree->clone(); 
          breadthFirstreNumber (*geneTree_);
          std::cout << " Problem tree? : "<<TreeTemplateTools::treeToParenthesis(*geneTree_, true) << std::endl;
          
          geneTree_->newOutGroup( geneTree_->getLeavesId()[0] ); 
              }
              else throw Exception("Unknown init gene tree method. init.gene.tree should be 'user', 'bionj', or 'phyml'.");
          }
          catch (std::exception& e)
          {
              std::cout << e.what() <<"; Unable to get a proper gene tree for family "<<file<<"; avoiding this family."<<std::endl;
              numDeletedFamilies_ = numDeletedFamilies_+1;
              avoidFamily=true;
          }
      }

      if (!avoidFamily) 
      { //This family is phylogenetically informative

          //Going through the gene tree to see if leaves have branches that are too long.
          std::vector <Node*> leaves = geneTree_->getLeaves();
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
                  std::vector <std::string> leafNames = geneTree_->getLeavesNames();
                  if ( VectorTools::contains(leafNames, seqsToRemove[j]) )
                  {
                      removeLeaf(*geneTree_, seqsToRemove[j]);
                      if (unrootedGeneTree) {
                          delete unrootedGeneTree;
                          unrootedGeneTree = 0;
                      }
                      unrootedGeneTree = geneTree_->clone();
                      if (!geneTree_->isRooted()) {
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
              numDeletedFamilies_ = numDeletedFamilies_+1;
              avoidFamily=true;
              std::cout <<"All or almost all sequences have been removed: avoiding family "<<assignedFilenames_[i-numDeletedFamilies_+1]<<std::endl;
          }
      }
      if (!avoidFamily) 
      { //This family is phylogenetically informative

          /****************************************************************************
           //Then we initialize the losses and duplication numbers on this tree.
           *****************************************************************************/
          std::vector<int> numbers = num0Lineages_;
          allNum0Lineages_.push_back(numbers);
          allNum1Lineages_.push_back(numbers);
          allNum2Lineages_.push_back(numbers);
          std::vector<unsigned int> uNumbers (tree_->getNumberOfNodes(), 0);

          allNum12Lineages_.push_back(uNumbers);
          allNum22Lineages_.push_back(uNumbers);

          resetLossesAndDuplications(*tree_, lossExpectedNumbers_, duplicationExpectedNumbers_);
          resetVector(allNum0Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum1Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum2Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum12Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum22Lineages_[i-numDeletedFamilies_]);

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
              if (tree_)
                  delete tree_;
              std::cout << "PHYLDOG's done. Bye." << std::endl;
              MPI::COMM_WORLD.Abort(1);
              exit(-1);
          }
          
          
          /****************************************************************************
           //Then we can change the gene tree so that its topology minimizes the number of duplications and losses.
           *****************************************************************************/
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
          
          //printing the gene trees with the species names instead of the sequence names
          //This is useful to build an input for duptree for instance
          TreeTemplate<Node> * treeWithSpNames = unrootedGeneTree->clone();
          std::vector <Node*> leaves = treeWithSpNames->getLeaves();
          for (unsigned int j =0; j<leaves.size() ; j++) 
          {
              leaves[j]->setName(seqSp[leaves[j]->getName()]);
          }
          //Now we have a gene tree with species names and bootstrap values and branch lengths
          std::vector <Node *> nodes =  treeWithSpNames->getNodes();
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
        annotateGeneTreeWithDuplicationEvents (*tree_, 
                           *geneTree_, 
                           geneTree_->getRootNode(), 
                           seqSp, spId_); 
        nhx->write(*geneTree_, startingGeneTreeFile, true);

       // newick.write(*geneTree_, startingGeneTreeFile, true);
      }
      catch (IOException e)
      {
        cout << "Problem writing tree to file "<< startingGeneTreeFile <<"\n Is the file path correct and do you have the proper authorizations? "  << endl;
      }
      try
      {
        newick.write(*treeWithSpNames, startingGeneTreeFile+"_sp_names", true);
      }
      catch (IOException e)
      {
        cout << "Problem writing tree to file "<< startingGeneTreeFile+"_sp_names" <<"\n Is the file path correct and do you have the proper authorizations? "  << endl;
      }
          //std::cout << " Rooted tree? : "<<TreeTemplateTools::treeToParenthesis(*geneTree_, true) << std::endl;
          
          GeneTreeLikelihood* tl;
          std::string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", params_, "no", "", true, false);

          int sprLimitGeneTree = ApplicationTools::getIntParameter("SPR.limit.gene.tree", params_, 2, "", false, false);  
          //ApplicationTools::displayResult("Clock", optimizeClock);
          
          if(optimizeClock == "global")
          {
              std::cout<<"Sorry, clocklike trees have not been implemented yet."<<std::endl;
              MPI::COMM_WORLD.Abort(1);
              exit(0);
          }// This has not been implemented!
          else if(optimizeClock == "no")
          {
              if (reconciliationModel_ == "DL")
              {
                  tl = new DLGeneTreeLikelihood(*unrootedGeneTree, *sites, 
                                                model, rDist, *tree_, 
                                                *geneTree_, *treeWithSpNames, seqSp, spId_, 
                                                lossExpectedNumbers_, 
                                                duplicationExpectedNumbers_, 
                                                allNum0Lineages_[i-numDeletedFamilies_], 
                                                allNum1Lineages_[i-numDeletedFamilies_], 
                                                allNum2Lineages_[i-numDeletedFamilies_], 
                                                speciesIdLimitForRootPosition_, 
                                                heuristicsLevel_, MLindex_, 
                                                true, true, rootOptimization, true, DLStartingGeneTree, sprLimitGeneTree);
                ///  dynamic_cast<DLGeneTreeLikelihood*> (tl)->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
                  
              }
              else if (reconciliationModel_ == "COAL")
              {
                  //coalCounts_: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
                  std::vector< std::vector< std::vector< unsigned int > > > coalCounts_2;
                  std::vector< std::vector<unsigned int> > coalCounts_3;
                  std::vector< unsigned int > coalCounts_4;
                  for (unsigned int j = 0 ; j < 2 ; j++ ) {
                      coalCounts_4.push_back(0);
                  }
                  for (unsigned int j = 0 ; j < tree_->getNumberOfNodes() ; j++ ) {
                      coalCounts_3.push_back(coalCounts_4);
                  }
                  for (unsigned int j = 0 ; j < 3 ; j++ ) {
                      coalCounts_2.push_back(coalCounts_3);
                  }
                  for (unsigned int j = 0 ; j < geneTree_->getNumberOfNodes() ; j++ ) {
                      coalCounts_.push_back(coalCounts_2);
                  }
                  tl = new COALGeneTreeLikelihood(*unrootedGeneTree, *sites, 
                                                  model, rDist, *tree_, 
                                                  *geneTree_, *treeWithSpNames, seqSp, spId_, 
                                                  coalCounts_, coalBls_,  
                                                  speciesIdLimitForRootPosition_, 
                                                  heuristicsLevel_, MLindex_, 
                                                  true, true, rootOptimization, true, sprLimitGeneTree);
              }
              else {
                  std::cerr <<"Unknown reconciliation model: "<< reconciliationModel_ <<std::endl;
                  exit(-1);
              }
              
          }
          else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
          delete treeWithSpNames;
          treeLikelihoods_.push_back(tl);          
          allParams_.push_back(params_); 
          allAlphabets_.push_back(alphabet);
          allDatasets_.push_back(sites);
          allModels_.push_back(model);
          allDistributions_.push_back(rDist);
          allGeneTrees_.push_back(geneTree_->clone());
          allUnrootedGeneTrees_.push_back(unrootedGeneTree->clone());
      allSeqSps_.push_back(seqSp);
      allSprLimitGeneTree_.push_back(sprLimitGeneTree);
    //  delete rDist;
    //  delete model;
      //delete sites;
    
      // delete alphabet;
      
/*          if (tl)
          delete tl;
          if (params_)
          delete params_;
          if (alphabet)
          delete alphabet;
          delete sites;
          delete model;
          delete rDist;
          delete geneTree_;
          delete unrootedGeneTree;*/
      }
      else 
      {
          if (sites)
              delete sites;
          if (alphabet)
              delete alphabet;
      }
      double familyTime = ApplicationTools::getTime() - startingFamilyTime ;  
      std::cout <<"Examined family "<<assignedFilenames_[i] << " in "<<familyTime<<" s."<<std::endl;
      if (unrootedGeneTree) {
          delete unrootedGeneTree;
          unrootedGeneTree = 0;
      }
      if (geneTree_) {
          delete geneTree_;
          geneTree_ = 0;
      }

  }//End for each file
    
    
    if (numDeletedFamilies_ == assignedFilenames_.size()) 
    {
        std::cout<<"WARNING: A client with rank "<< rank_ << " is in charge of 0 gene family after gene family filtering!"<<std::endl;        
        std::cout<<"A processor will be idle most of the time, the load could probably be better distributed."<<std::endl; 
    }

    
    unsigned int numberOfGeneFamilies = assignedFilenames_.size()-numDeletedFamilies_;
    numberOfFilteredFamiliesCommunicationsServerClient (world_, server_, 
                                   rank_, numberOfGeneFamilies);    
    
    //Building a MRP species tree, if the options say so
    initTree = ApplicationTools::getStringParameter("init.species.tree", 
                                                    params_, "user", 
                                                    "", false, false);
    if (initTree == "mrp") {
        //Build a string containing all trees with one gene per species
        string trees1PerSpecies = "";
        vector<string> allTrees1PerSpecies;
        stringstream ss (stringstream::in | stringstream::out);
        for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
        {
            if (reconciliationModel_ == "DL")
            {
            if ( ( dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i]) )->isSingleCopy() ) {
                ss << TreeTemplateTools::treeToParenthesis((dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i]) )->getGeneTreeWithSpNames(), false);
            }
            }
            else {
                if ((dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i]) )->isSingleCopy() ) {
                    ss << TreeTemplateTools::treeToParenthesis((dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i]) )->getGeneTreeWithSpNames(), false);
                }

            }
        }
        trees1PerSpecies = ss.str() ;     
        mrpCommunicationsServerClient (world_, server_, 
                                       rank_, trees1PerSpecies, 
                                       allTrees1PerSpecies);

    }
    vector<unsigned int> numbersOfGeneFamilies;
    
    secondCommunicationsServerClient (world_ , server_, 
                                      rank_, numberOfGeneFamilies, 
                                      numbersOfGeneFamilies, 
                                      lossExpectedNumbers_, duplicationExpectedNumbers_, coalBls_, 
                                      currentSpeciesTree_);

    delete tree_;

    tree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);

    spId_ = computeSpeciesNamesToIdsMap(*tree_);
    
    for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
    {        
        treeLikelihoods_[i]->setSpTree(*tree_);

        treeLikelihoods_[i]->setSpId(spId_);        
        //GeneTreeLikelihood* tl ;
        if (reconciliationModel_ == "DL")
        {
            DLGeneTreeLikelihood* tl = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i]);
            tl->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
            tl->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
            tl->initParameters();
            allLogLs_.push_back(tl->getValue());
            if(std::isinf(allLogLs_[i]))
            {
                // This may be due to null branch lengths, leading to null likelihood!
                ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
                ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
                ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
                std::vector<Node*> nodes = treeLikelihoods_[i]->getRootedTree().getNodes();
                for(unsigned int k = 0; k < nodes.size(); k++)
                {
                    if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
                }
                tl->initParameters();
                // allLogLs_[i-numDeletedFamilies_]= tl->f(tl->getParameters());
                allLogLs_[i]= tl->getValue();
            }
            ApplicationTools::displayResult("Initial sequence and reconciliation likelihood", TextTools::toString(- allLogLs_[i], 15));
            if(std::isinf(allLogLs_[i]))
            {
                ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
                ApplicationTools::displayError("!!! Looking at each site:");
                /*   for(unsigned int k = 0; k < sites->getNumberOfSites(); k++)
                 {
                 (*ApplicationTools::error << "Site " << sites->getSite(k).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(k)).endLine();
                 }*/
                ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularly if datasets are big (>~500 sequences).");
                MPI::COMM_WORLD.Abort(1);
                exit(-1);
            }            
        }
        else if (reconciliationModel_ == "COAL") {
            COALGeneTreeLikelihood* tl = dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i] );
            tl->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
            tl->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
            tl->initParameters();
            allLogLs_.push_back(tl->getValue());
            if(std::isinf(allLogLs_[i]))
            {
                // This may be due to null branch lengths, leading to null likelihood!
                ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
                ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
                ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
                std::vector<Node*> nodes = treeLikelihoods_[i]->getRootedTree().getNodes();
                for(unsigned int k = 0; k < nodes.size(); k++)
                {
                    if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
                }
                tl->initParameters();
                // allLogLs_[i-numDeletedFamilies_]= tl->f(tl->getParameters());
                allLogLs_[i]= tl->getValue();
            }
            ApplicationTools::displayResult("Initial sequence and reconciliation likelihood", TextTools::toString(- allLogLs_[i], 15));
            if(std::isinf(allLogLs_[i]))
            {
                ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
                ApplicationTools::displayError("!!! Looking at each site:");
                /*   for(unsigned int k = 0; k < sites->getNumberOfSites(); k++)
                 {
                 (*ApplicationTools::error << "Site " << sites->getSite(k).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(k)).endLine();
                 }*/
                ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularly if datasets are big (>~500 sequences).");
                MPI::COMM_WORLD.Abort(1);
                exit(-1);
            }

        }
    }

    return;
}


    /****************************************************************************
    ****************************************************************************
    * Main loop: iterative likelihood computations
    ****************************************************************************
    ****************************************************************************/

  void ClientComputingGeneLikelihoods::MLSearch() {
	Nhx *nhx = new Nhx();
	string rearrangementType;
	bool timing = true;
	double startingTime, totalTime;
	while (!stop_)            
	{      
	    logL_=0.0;
	    resetVector(num0Lineages_);
	    resetVector(num1Lineages_);
	    resetVector(num2Lineages_);
	    resetVector(num12Lineages_);
	    resetVector(num22Lineages_);
	    for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
	    {
		if (rearrange_) //(firstTimeImprovingGeneTrees) 
		{
		    treeLikelihoods_[i]->OptimizeSequenceLikelihood(true);
		    allParams_[i][ std::string("optimization.topology")] = "true";
		}
		else {
		    treeLikelihoods_[i]->OptimizeSequenceLikelihood(false);
		    allParams_[i][ std::string("optimization.topology")] = "false";
		}
		// std::cout <<  TreeTemplateTools::treeToParenthesis(*geneTree_, true)<<std::endl;
		
		rearrangementType = ApplicationTools::getStringParameter("rearrangement.gene.tree", allParams_[i], "nni", "", true, false);
		if (rearrangementType == "nni" && currentStep_ !=4 ) {
		    //PhylogeneticsApplicationTools::optimizeParameters(treeLikelihoods_[i], treeLikelihoods_[i]->getParameters(), allParams_[i], "", true, false);
		    if (timing) 
			startingTime = ApplicationTools::getTime();
		    // NNI optimization:  
		    if (reconciliationModel_ == "DL") {
			dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeNNIs(allParams_[i]);
		    }
		    else if (reconciliationModel_ == "COAL") {
			dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeNNIs(allParams_[i]);
		    }
		    if (timing) 
		    {
			totalTime = ApplicationTools::getTime() - startingTime;
			std::cout << "Family "<< assignedFilenames_[i] <<"; Time for NNI exploration: "<<  totalTime << " s." <<std::endl;
		    }
		}
		else {
		    if (timing) 
			startingTime = ApplicationTools::getTime();
		    //SPR optimization:    
		    //std::cout <<"Before optimization: "<<TreeTemplateTools::treeToParenthesis(treeLikelihoods_[i]->getRootedTree(), true)<<std::endl;
		  std::string SPRalgorithm = ApplicationTools::getStringParameter("spr.gene.tree.algorithm", allParams_[i], "normal", "", true, false);
		    if (reconciliationModel_ == "DL") {
	  dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeMuffato(allParams_[i]);
						    if (timing) 
						    {
							    totalTime = ApplicationTools::getTime() - startingTime;
							    if (rearrange_)
							    {
								    std::cout << "Family "<< assignedFilenames_[i] <<"; Time for Muffato exploration: "<<  totalTime << " s." <<std::endl;
							    }
							    startingTime = ApplicationTools::getTime();
						    }
						    if (SPRalgorithm == "normal")
							    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast2(allParams_[i]);
						    else
							    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast3(allParams_[i]);
					    }
					    else if (reconciliationModel_ == "COAL") {
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast(allParams_[i]);
					    }
					    //   treeLikelihoods_[i]->refineGeneTreeSPRs(allParams_[i]);
					    //  treeLikelihoods_[i]->refineGeneTreeSPRs2(allParams_[i]);
					    if (timing) 
					    {
						    totalTime = ApplicationTools::getTime() - startingTime;
						    std::cout << "Family "<< assignedFilenames_[i] <<"; Time for SPR exploration: "<<  totalTime << " s." <<std::endl;
					    }
				    }     
		if (geneTree_) {
		    delete geneTree_;
		    geneTree_ = 0;
		}
		geneTree_ = new TreeTemplate<Node>(treeLikelihoods_[i]->getRootedTree());

		///LIKELIHOOD OPTIMIZED
	      // resetLossesAndDuplications(*tree, lossExpectedNumbers, duplicationExpectedNumbers);
		if ( reconciliationModel_ == "DL" ) {
		    allNum0Lineages_[i] = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->get0LineagesNumbers();
		    allNum1Lineages_[i] = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->get1LineagesNumbers();
		    allNum2Lineages_[i] = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->get2LineagesNumbers();
		    num0Lineages_ = num0Lineages_ + allNum0Lineages_[i];
		    num1Lineages_ = num1Lineages_ + allNum1Lineages_[i];
		    num2Lineages_ = num2Lineages_ + allNum2Lineages_[i];
		    MLindex_ = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->getRootNodeindex();
		    allLogLs_[i] = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->getValue();  
		}
		else if (reconciliationModel_ == "COAL") {
		    allNum12Lineages_[i] = dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->getNum12Lineages();
		    allNum22Lineages_[i] = dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->getNum22Lineages();
		    num12Lineages_ = num12Lineages_ + allNum12Lineages_[i];
		    num22Lineages_ = num22Lineages_ + allNum22Lineages_[i];
		    MLindex_ = dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->getRootNodeindex();
		    allLogLs_[i] = dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->getValue();  
		}
		logL_ = logL_ + allLogLs_[i];
		std::cout<<"Gene Family: " << assignedFilenames_[i] << " total logLk: "<< - allLogLs_[i]<< " ; scenario loglk: "<< treeLikelihoods_[i]->getScenarioLikelihood() <<std::endl;

		if (std::isnan(allLogLs_[i])) 
		{
		    std::cout<<TreeTemplateTools::treeToParenthesis (*geneTree_, false, EVENT)<<std::endl;
		    std::cout<<TreeTemplateTools::treeToParenthesis (*tree_, false, DUPLICATIONS)<<std::endl;
		}
		if (recordGeneTrees_) 
		{
		    reconciledTrees_[i].push_back(nhx->treeToParenthesis (*geneTree_));
		    if (reconciliationModel_ == "DL") {
			duplicationTrees_[i].push_back(nhx->treeToParenthesis (*tree_));
			lossTrees_[i].push_back(nhx->treeToParenthesis (*tree_));
		    }
		}
		if (geneTree_) 
		{
		    delete geneTree_;
		    geneTree_ = 0;
		}
	    }//end for each filename
	    if (!recordGeneTrees_) 
	    {
		startRecordingTreesFrom_++;
	    }
			    if (timing) 
			    {
				    startingTime = ApplicationTools::getTime();
			    }
	    //Clients send back stuff to the server.
	    gathersInformationFromClients (world_, 
					  server_,
					  rank_, 
					  logL_, 
					  num0Lineages_, 
					  num1Lineages_, 
					  num2Lineages_, 
					  allNum0Lineages_, 
					  allNum1Lineages_, 
					  allNum2Lineages_, 
					  num12Lineages_, num22Lineages_, reconciliationModel_);
			    if (timing) 
			    {
				    totalTime = ApplicationTools::getTime() - startingTime;
				    std::cout << "Time for gathering information: "<<  totalTime << " s." <<std::endl;
			    }
			    
	    //Should the computations stop? The server tells us.
	    broadcast(world_, stop_, server_);
			    if (!stop_) 
			    { // we continue the loop
				    //We get the new values from the server
				    if (timing) 
				    {
					    startingTime = ApplicationTools::getTime();
				    }
				    broadcastsAllInformationButStop(world_, server_, rearrange_, 
								    lossExpectedNumbers_, 
								    duplicationExpectedNumbers_, 
								    coalBls_,
								    currentSpeciesTree_,
								    currentStep_, 
								    reconciliationModel_);
				    if (timing) 
				    {
					    totalTime = ApplicationTools::getTime() - startingTime;
					    std::cout << "Time for broadcasting information: "<<  totalTime << " s." <<std::endl;
					    startingTime = ApplicationTools::getTime();
				    }

				    if (assignedFilenames_.size()-numDeletedFamilies_ > 0) 
				    {
					    if (tree_) delete tree_;
					    tree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
					    spId_ = computeSpeciesNamesToIdsMap(*tree_);
				    }
		//if we reset the gene trees by resetting treeLikelihoods_:
		//we always start from ML trees according to sequences only
		//when we optimize dl expected numbers, we may not want to 
		//rearrange the gene trees between the first computation of the species tree
		//lk and the second one; in this case we set rearrange to false.
		//Then there is no need to reset the gene tree!
		if (resetGeneTrees_ && currentStep_ !=4 && rearrange_ == true) {
		    for (unsigned int i =0 ; i<treeLikelihoods_.size() ; i++) {
			if    (treeLikelihoods_[i])
			    delete treeLikelihoods_[i];
		    }
		    treeLikelihoods_.clear();
					    if (reconciliationModel_ == "DL")
					    {
						    for (unsigned int i=0 ; i<allDatasets_.size() ; i++) 
						    {
							    TreeTemplate<Node> * treeWithSpNames = allUnrootedGeneTrees_[i]->clone();
							    std::vector <Node*> leaves = treeWithSpNames->getLeaves();
							    for (unsigned int j =0; j<leaves.size() ; j++) 
							    {
								    leaves[j]->setName(allSeqSps_[i][leaves[j]->getName()]);
							    }
							    treeLikelihoods_.push_back( new DLGeneTreeLikelihood(*(allUnrootedGeneTrees_[i]), *(allDatasets_[i]), 
														allModels_[i], allDistributions_[i], *tree_, 
														*allGeneTrees_[i], *treeWithSpNames, allSeqSps_[i], spId_, 
														lossExpectedNumbers_, 
														duplicationExpectedNumbers_, 
														allNum0Lineages_[i], 
														allNum1Lineages_[i], 
														allNum2Lineages_[i], 
														speciesIdLimitForRootPosition_, 
														heuristicsLevel_, MLindex_, 
														true, true, true, true, false, allSprLimitGeneTree_[i]) );
							    delete treeWithSpNames;
						    }
					    }
					    else if (reconciliationModel_ == "COAL")
					    {
						    for (unsigned int i=0 ; i<allDatasets_.size() ; i++) 
						    {
							    //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
							    std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
							    std::vector< std::vector<unsigned int> > coalCounts3;
							    std::vector< unsigned int > coalCounts4;
							    for (unsigned int j = 0 ; j < 2 ; j++ ) {
								    coalCounts4.push_back(0);
							    }
							    for (unsigned int j = 0 ; j < tree_->getNumberOfNodes() ; j++ ) {
								    coalCounts3.push_back(coalCounts4);
							    }
							    for (unsigned int j = 0 ; j < 3 ; j++ ) {
								    coalCounts2.push_back(coalCounts3);
							    }
							    for (unsigned int j = 0 ; j < geneTree_->getNumberOfNodes() ; j++ ) {
								    coalCounts_.push_back(coalCounts2);
							    }
							    TreeTemplate<Node> * treeWithSpNames = allUnrootedGeneTrees_[i]->clone();
							    std::vector <Node*> leaves = treeWithSpNames->getLeaves();
							    for (unsigned int j =0; j<leaves.size() ; j++) 
							    {
								    leaves[j]->setName(allSeqSps_[i][leaves[j]->getName()]);
							    }

							    treeLikelihoods_.push_back( new COALGeneTreeLikelihood(*allUnrootedGeneTrees_[i], *allDatasets_[i], 
														  allModels_[i], allDistributions_[i], *tree_, 
														  *allGeneTrees_[i], *treeWithSpNames, allSeqSps_[i], spId_, 
														  coalCounts_, coalBls_,  
														  speciesIdLimitForRootPosition_, 
														  heuristicsLevel_, MLindex_, 
														  true, true, true, true, allSprLimitGeneTree_[i]) );
							    delete treeWithSpNames;

						    }
					    }
				    }
				    if (currentStep_>=3)// rearrange_) //?
				    {
					    allParams_ = allParamsBackup_;
					    if (recordGeneTrees_==false) 
					    {
						    recordGeneTrees_=true;
					    }
				    } 
				    for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
				    {
					    treeLikelihoods_[i]->setSpTree(*tree_);
					    treeLikelihoods_[i]->setSpId(spId_);
					    if (reconciliationModel_ == "DL") {
						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->setExpectedNumbers(duplicationExpectedNumbers_, lossExpectedNumbers_);
						    //If not using the backuplks
						    if (! dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->isInitialized() ) {
							    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
						    }
						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->initParameters();
					    }
					    else if (reconciliationModel_ == "COAL") {
						    if (! dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->isInitialized() ) {
							    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->setCoalBranchLengths(coalBls_);
						    }
						    //If not using the backuplks
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->initParameters();
					    }
				    }
				    if (timing) 
				    {
					    totalTime = ApplicationTools::getTime() - startingTime;
					    std::cout << "Time for resetting gene tree likelihoods: "<<  totalTime << " s." <<std::endl;
				    }
			    }
			    else 
			    { 
				    /****************************************************************************
				    * The end, outputting the results.
				    *****************************************************************************/      
				    if (recordGeneTrees_) 
				    {
				      unsigned int bestIndex;
					    broadcast(world_, bestIndex, server_);
					   // std::cout << "bestIndex: "<<bestIndex<<" startRecordingTreesFrom: "<<startRecordingTreesFrom_<<std::endl;
					    outputGeneTrees( bestIndex );
				    }
				    break;
			    }
		    }//End while, END OF MAIN LOOP
		    delete nhx;
  }


  
  
  
  
  
/******************************************************************************/
// This function outputs gene trees from the clients.
/******************************************************************************/
void ClientComputingGeneLikelihoods::outputGeneTrees ( unsigned int & bestIndex )
{
  std::string suffix;
  std::string reconcTree;
  std::ofstream out;
  std::string dupTree;
  std::string lossTree;
  int rightIndex = bestIndex-startRecordingTreesFrom_ ;
  for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
    {
    Nhx *nhx = new Nhx();
    string temp = reconciledTrees_[i][rightIndex];                                                                        
  if (reconciliationModel_ == "DL") {
      TreeTemplate<Node> * geneTree=nhx->parenthesisToTree(temp);
      temp = duplicationTrees_[i][rightIndex]; 
      TreeTemplate<Node> * spTree=nhx->parenthesisToTree( temp);                                                                  
      writeReconciledGeneTree ( allParams_[i], geneTree,  spTree, treeLikelihoods_[i]->getSeqSp(), false ); 

            dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", allParams_[i], "duplications.tree", "", false, false);

            dupTree = dupTree + suffix;
            out.open (dupTree.c_str(), std::ios::out);
            out << duplicationTrees_[i][rightIndex]<<std::endl;
            out.close();
            lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", allParams_[i], "losses.tree", "", false, false);
            lossTree = lossTree + suffix;
            out.open (lossTree.c_str(), std::ios::out);
            out << lossTrees_[i][rightIndex]<<std::endl;
            delete geneTree;
            delete spTree;
        }
        else if (reconciliationModel_ == "COAL") {

    TreeTemplate<Node> * geneTree=TreeTemplateTools::parenthesisToTree(temp);
            out.open (reconcTree.c_str(), std::ios::out);
            nhx->write(*geneTree, out);
            out.close();
            delete geneTree;

        }

        out.close();
        delete nhx;

    }
  return;
}


