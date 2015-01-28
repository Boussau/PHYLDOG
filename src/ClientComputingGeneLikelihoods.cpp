//           ClientComputingGeneLikelihoods.cpp
//  Mon March 03 11:28:42 2014
//  Copyright  2014  boussau
//  <user@host>

#include "Constants.h"
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
                                int MLindex_, 
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
    WHEREAMI( __FILE__ , __LINE__ );

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
            spTree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
 
            resetLossesAndDuplications(*spTree_, lossExpectedNumbers_, duplicationExpectedNumbers_);
            //To make the correspondance between species name and id:
            spId_ = computeSpeciesNamesToIdsMap(*spTree_);
            
            /****************************************************************************
             * Gene family parsing and first likelihood computation.
             *****************************************************************************/            
            string toPrint = "";
            std::cout <<"Client  of rank "<<rank_ <<" with PID "<< TextTools::toString((int)getpid()) <<" is in charge of " << assignedFilenames_.size()<<" gene families:"<<std::endl;
            for (unsigned int i = 0 ; i< assignedFilenames_.size() ; i++ ) {
                toPrint = toPrint + assignedFilenames_[i] + "\t";
            }
            std::cout << toPrint <<std::endl;
            //Gets gene family-specific options, builds the GeneTreeLikelihood objects, and computes the likelihood
            parseAssignedGeneFamilies();
            std::vector <std::string> t;  
            allParamsBackup_ = allParams_;
            resetGeneTrees_ = ApplicationTools::getBooleanParameter("reset.gene.trees",params_,true ); 
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
                if (spTree_) delete spTree_;
                spTree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
                spId_ = computeSpeciesNamesToIdsMap(*spTree_);
            }
            startRecordingTreesFrom_ = 1;
            for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
            {
                treeLikelihoods_[i]->setSpTree(*spTree_);
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
  WHEREAMI( __FILE__ , __LINE__ );
  bool avoidFamily;
  std::string initTree;
  std::map<std::string, std::string> famSpecificParams;

  //Here we are going to get all necessary information regarding all gene families the client is in charge of.
  for (unsigned int i = 0 ; i< assignedFilenames_.size() ; i++) 
  { //For each file
      double startingFamilyTime = ApplicationTools::getTime();  
      std::cout <<"Examining family "<<assignedFilenames_[i]<<std::endl;
      avoidFamily = false;
      std::string familySpecificOptionsFile = assignedFilenames_[i];
      if(!FileTools::fileExists(familySpecificOptionsFile))
      {
          std::cerr << "Error: Parameter file "<< familySpecificOptionsFile <<" not found." << std::endl;
          fflush(0);
          MPI::COMM_WORLD.Abort(1);
          exit(-1);
      }
      else
      {
          famSpecificParams = AttributesTools::getAttributesMapFromFile(familySpecificOptionsFile, "=");
          AttributesTools::resolveVariables(famSpecificParams);
      }
      
      AttributesTools::actualizeAttributesMap( params_, famSpecificParams);
      //COAL or DL?
      reconciliationModel_ = ApplicationTools::getStringParameter("reconciliation.model", params_, "DL", "", true, false);
      GeneTreeLikelihood *tl = 0; 
      //HEHEHEHEHEHEHEHEEHEHEH
      if (reconciliationModel_ == "DL")
              {
								try {
											tl = new DLGeneTreeLikelihood(familySpecificOptionsFile, params_, *spTree_);
								}
								catch (exception& e)
									{
		  							cout << e.what() << '\n';
		  							avoidFamily = true;
		  							numDeletedFamilies_ = numDeletedFamilies_ +1;
									}
									/*
                  tl = new DLGeneTreeLikelihood(*unrootedGeneTree, *sites, 
                                                model, rDist, *spTree_, 
                                                *geneTree_, *treeWithSpNames, seqSp, spId_, 
                                                lossExpectedNumbers_, 
                                                duplicationExpectedNumbers_, 
                                                allNum0Lineages_[i-numDeletedFamilies_], 
                                                allNum1Lineages_[i-numDeletedFamilies_], 
                                                allNum2Lineages_[i-numDeletedFamilies_], 
                                                speciesIdLimitForRootPosition_, 
                                                 MLindex_, 
                                                true, true, rootOptimization, true, DLStartingGeneTree, sprLimitGeneTree);*/
                ///  dynamic_cast<DLGeneTreeLikelihood*> (tl)->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
                  
              }
              else if (reconciliationModel_ == "COAL")
              {
		try {
		 tl = new COALGeneTreeLikelihood( familySpecificOptionsFile, params_, *spTree_ );
		 }
		catch (exception& e)
		{
		  cout << e.what() << '\n';
		  avoidFamily = true;
		  numDeletedFamilies_ = numDeletedFamilies_ +1;
		}
                /*  tl = new COALGeneTreeLikelihood(*unrootedGeneTree, *sites, 
                                                  model, rDist, *spTree_, 
                                                  *geneTree_, *treeWithSpNames, seqSp, spId_, 
                                                  coalCounts_, coalBls_,  
                                                  speciesIdLimitForRootPosition_, 
                                                   MLindex_, 
                                                  true, true, rootOptimization, true, sprLimitGeneTree);*/
              }
              else {
                  std::cerr <<"Unknown reconciliation model: "<< reconciliationModel_ <<std::endl;
                  exit(-1);
              }
      
      
      double familyTime = ApplicationTools::getTime() - startingFamilyTime ;  
      std::cout <<"Examined family "<<assignedFilenames_[i] << " in "<<familyTime<<" s."<<std::endl;
 
      if (!avoidFamily) {
	  // We need to fill these vectors so that we can quickly re-create *GeneTreeLikelihood objects 
	  // in the course of the algorithm.
          treeLikelihoods_.push_back(tl);          
          allParams_.push_back( params_ ); 
         // allAlphabets_.push_back(alphabet);
          allDatasets_.push_back(tl->getSequenceLikelihoodObject()->getSites()->clone());
          allModels_.push_back(tl->getSequenceLikelihoodObject()->getSubstitutionModel ());
          allDistributions_.push_back(tl->getSequenceLikelihoodObject()->getRateDistribution ());
          allGeneTrees_.push_back(tl->getRootedTree().clone());
          allUnrootedGeneTrees_.push_back(new TreeTemplate<Node>(*(tl->getSequenceLikelihoodObject()->getTree())) );
	  allSeqSps_.push_back( tl->getSeqSp() );
	  allSprLimitGeneTree_.push_back(tl->getSprLimitGeneTree() );
	  
	   /****************************************************************************
           //Then we initialize the losses and duplication numbers on this tree for this family.
           *****************************************************************************/
          std::vector<int> numbers = num0Lineages_;
          allNum0Lineages_.push_back(numbers);
          allNum1Lineages_.push_back(numbers);
          allNum2Lineages_.push_back(numbers);
          std::vector<unsigned int> uNumbers (spTree_->getNumberOfNodes(), 0);

          allNum12Lineages_.push_back(uNumbers);
          allNum22Lineages_.push_back(uNumbers);

          resetLossesAndDuplications(*spTree_, lossExpectedNumbers_, duplicationExpectedNumbers_);
          resetVector(allNum0Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum1Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum2Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum12Lineages_[i-numDeletedFamilies_]);
          resetVector(allNum22Lineages_[i-numDeletedFamilies_]);
      }

	  
	  /*
     if (unrootedGeneTree) {
          delete unrootedGeneTree;
          unrootedGeneTree = 0;
      }*/
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

    delete spTree_;

    spTree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);

    spId_ = computeSpeciesNamesToIdsMap(*spTree_);
    
    //std::cout << "THET "<< assignedFilenames_.size() << " and "<< numDeletedFamilies_ << " anana " << reconciliationModel_ <<std::endl;
    
    for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
    {        
        treeLikelihoods_[i]->setSpTree(*spTree_);

        treeLikelihoods_[i]->setSpId(spId_);        
        //GeneTreeLikelihood* tl ;
        if (reconciliationModel_ == "DL")
        {
            DLGeneTreeLikelihood* tl = dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i]);
            tl->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
//             tl->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
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
            std::cout <<"Family "<< assignedFilenames_[i]  << " ; Initial sequence and reconciliation likelihood: " << TextTools::toString(- allLogLs_[i], 15) <<std::endl;
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
//             tl->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
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
  WHEREAMI( __FILE__ , __LINE__ );
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
		
		rearrangementType = ApplicationTools::getStringParameter("rearrangement.gene.tree", allParams_[i], "spr", "", true, false);
		if (rearrangementType == "nni" || currentStep_ !=4 ) {
		    //PhylogeneticsApplicationTools::optimizeParameters(treeLikelihoods_[i], treeLikelihoods_[i]->getParameters(), allParams_[i], "", true, false);
      NNIRearrange(timing, i, startingTime, totalTime);
 		}
		else {
		    if (timing) 
                startingTime = ApplicationTools::getTime();
		    //SPR optimization:    
		    //std::cout <<"Before optimization: "<<TreeTemplateTools::treeToParenthesis(treeLikelihoods_[i]->getRootedTree(), true)<<std::endl;
        std::string SPRalgorithm = ApplicationTools::getStringParameter("spr.gene.tree.algorithm", allParams_[i], "normal", "", true, false);
                  std::cout << "rearrangementType "<< rearrangementType << std::endl;
          std::cout << "SPRalgorithm "<< SPRalgorithm << std::endl;

        if (reconciliationModel_ == "DL") {
  WHEREAMI( __FILE__ , __LINE__ );
          if (rearrangementType == "nni") {
              WHEREAMI( __FILE__ , __LINE__ );

            NNIRearrange(timing, i, startingTime, totalTime);
          }
          else {
                //dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeMuffato(allParams_[i]);
                NNIRearrange(timing, i, startingTime, totalTime);
  WHEREAMI( __FILE__ , __LINE__ );

                if (timing && rearrange_) 
                {
                  totalTime = ApplicationTools::getTime() - startingTime;
                  /*if (rearrange_)
                  {
                    std::cout << "Family "<< assignedFilenames_[i] <<"; Time for Muffato exploration: "<<  totalTime << " s." <<std::endl;
                  }*/
                  startingTime = ApplicationTools::getTime();
                }
                  WHEREAMI( __FILE__ , __LINE__ );

                if (SPRalgorithm == "normal")
                  dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast2(allParams_[i]);
                else
                  dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast3(allParams_[i]);
                }
                  WHEREAMI( __FILE__ , __LINE__ );

          }
          else if (reconciliationModel_ == "COAL") {
              WHEREAMI( __FILE__ , __LINE__ );

            if (rearrangementType == "nni") {
              NNIRearrange(timing, i, startingTime, totalTime);
            }
            else {
              dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeSPRsFast(allParams_[i]);
            }
          }
					    //   treeLikelihoods_[i]->refineGeneTreeSPRs(allParams_[i]);
					    //  treeLikelihoods_[i]->refineGeneTreeSPRs2(allParams_[i]);
          if (timing && rearrange_) 
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
		    std::cout<<TreeTemplateTools::treeToParenthesis (*spTree_, false, DUPLICATIONS)<<std::endl;
		}
		if (recordGeneTrees_) 
		{
		    reconciledTrees_[i].push_back(nhx->treeToParenthesis (*geneTree_));
		    if (reconciliationModel_ == "DL") {
/*			duplicationTrees_[i].push_back(nhx->treeToParenthesis (*spTree_));
			lossTrees_[i].push_back(nhx->treeToParenthesis (*spTree_));*/
      duplicationTrees_[i].push_back(TreeTemplateTools::treeToParenthesis (*spTree_, false, DUPLICATIONS));
      lossTrees_[i].push_back(TreeTemplateTools::treeToParenthesis (*spTree_, false, LOSSES));
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
					    if (spTree_) delete spTree_;
					    spTree_=TreeTemplateTools::parenthesisToTree(currentSpeciesTree_, false, "", true);
					    spId_ = computeSpeciesNamesToIdsMap(*spTree_);
				    }
		//if we reset the gene trees by resetting treeLikelihoods_:
		//we always start from ML trees according to sequences only
		//when we optimize dl expected numbers, we may not want to 
		//rearrange the gene trees between the first computation of the species tree
		//lk and the second one; in this case we set rearrange to false.
		//Then there is no need to reset the gene tree!
		if (resetGeneTrees_ && currentStep_ !=4 && rearrange_ == true) {
					    if (reconciliationModel_ == "DL")
					    {
						    	for (unsigned int i=0 ; i<allDatasets_.size() ; i++) 
						    	{
										string methodString =  treeLikelihoods_[i]->getLikelihoodMethod () ;
										if (methodString == "PLL") {
											treeLikelihoods_[i]->setGeneTree( allUnrootedGeneTrees_[i], allGeneTrees_[i] );
										}
										else {
						    /*	for (unsigned int i=0 ; i<allDatasets_.size() ; i++) 
						    	{*/
										std::map <std::string, std::string > params = treeLikelihoods_[i]->getParams();
										if    (treeLikelihoods_[i])
								    delete treeLikelihoods_[i];
									  TreeTemplate<Node> * treeWithSpNames = allUnrootedGeneTrees_[i]->clone();
							  	  std::vector <Node*> leaves = treeWithSpNames->getLeaves();
							  	  for (unsigned int j =0; j<leaves.size() ; j++) 
							  	  {
									    leaves[j]->setName(allSeqSps_[i][leaves[j]->getName()]);
							  	  }
							  	  treeLikelihoods_[i] =  new DLGeneTreeLikelihood(*(allUnrootedGeneTrees_[i]),                              *(allDatasets_[i]), 
														allModels_[i], allDistributions_[i], *spTree_, 
														*allGeneTrees_[i], *treeWithSpNames, allSeqSps_[i], spId_, 
														lossExpectedNumbers_, 
														duplicationExpectedNumbers_, 
														allNum0Lineages_[i], 
														allNum1Lineages_[i], 
														allNum2Lineages_[i], 
														speciesIdLimitForRootPosition_, 
														 MLindex_, params, 
														true, true, true, true, false, allSprLimitGeneTree_[i]) ;
          /*                          treeLikelihoods_.push_back( new DLGeneTreeLikelihood(*(allUnrootedGeneTrees_[i]), *(allDatasets_[i]), 
                                                        allModels_[i], allDistributions_[i], *spTree_, 
                                                        *allGeneTrees_[i], *treeWithSpNames, allSeqSps_[i], spId_, 
                                                        lossExpectedNumbers_, 
                                                        duplicationExpectedNumbers_, 
                                                        allNum0Lineages_[i], 
                                                        allNum1Lineages_[i], 
                                                        allNum2Lineages_[i], 
                                                        speciesIdLimitForRootPosition_, 
                                                         MLindex_, params, 
                                                        true, true, true, true, false, allSprLimitGeneTree_[i]) );*/
                                  
                                  
                                  
							  	  delete treeWithSpNames;
									}
								}
					    }
					    else if (reconciliationModel_ == "COAL")
					    {
								for (unsigned int i =0 ; i<treeLikelihoods_.size() ; i++) {
									std::map <std::string, std::string > params = treeLikelihoods_[i]->getParams();
									if    (treeLikelihoods_[i])
								    delete treeLikelihoods_[i];
		  					//}
		    			//	treeLikelihoods_.clear();
						 /*   for (unsigned int i=0 ; i<allDatasets_.size() ; i++) 
						    {*/
							    //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
							    std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
							    std::vector< std::vector<unsigned int> > coalCounts3;
							    std::vector< unsigned int > coalCounts4;
							    for (unsigned int j = 0 ; j < 2 ; j++ ) {
								    coalCounts4.push_back(0);
							    }
							    for (unsigned int j = 0 ; j < spTree_->getNumberOfNodes() ; j++ ) {
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
														  allModels_[i], allDistributions_[i], *spTree_, 
														  *allGeneTrees_[i], *treeWithSpNames, allSeqSps_[i], spId_, 
														  coalCounts_, coalBls_,  
														  speciesIdLimitForRootPosition_, 
														   MLindex_, params,
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
					    treeLikelihoods_[i]->setSpTree(*spTree_);
					    treeLikelihoods_[i]->setSpId(spId_);
					    if (reconciliationModel_ == "DL") {
						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->setExpectedNumbers(duplicationExpectedNumbers_, lossExpectedNumbers_);
						    //If not using the backuplks
						    if (! dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->isInitialized() ) {
							    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
						    }
// 						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
						    dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->initParameters();
					    }
					    else if (reconciliationModel_ == "COAL") {
						    if (! dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->isInitialized() ) {
							    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->setCoalBranchLengths(coalBls_);
						    }
						    //If not using the backuplks
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
// 						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->optimizeNumericalParameters(params_); //Initial optimization of all numerical parameters
						    dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->initParameters();
					    }
				    }
				    if (timing && resetGeneTrees_) 
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
  WHEREAMI( __FILE__ , __LINE__ );
  std::string suffix;
  std::string reconcTree;
  std::ofstream out;
  std::string dupTree;
  std::string lossTree;
  int rightIndex = bestIndex-startRecordingTreesFrom_ ;
  for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
    {
    Nhx *nhx = new Nhx();
    if (reconciliationModel_ == "DL") {
        writeReconciledGeneTree ( allParams_[i], treeLikelihoods_[i]->getRootedTree().clone(), spTree_, treeLikelihoods_[i]->getSeqSp(), false ) ;

            /*
     *                                                                          
        TreeTemplate<Node> * geneTree=nhx->parenthesisToTree(temp);
        writeReconciledGeneTree ( allParams_[i], geneTree,  spTree, treeLikelihoods_[i]->getSeqSp(), false ); 
      temp = duplicationTrees_[i][rightIndex]; 

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
            delete spTree;*/
        }
        else if (reconciliationModel_ == "COAL") {
string temp = reconciledTrees_[i][rightIndex];
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


void ClientComputingGeneLikelihoods::NNIRearrange (bool timing, size_t i, double& startingTime, double& totalTime) {
    if (timing && rearrange_) 
      startingTime = ApplicationTools::getTime();
    // NNI optimization:  
    if (reconciliationModel_ == "DL") {
      dynamic_cast<DLGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeNNIs(allParams_[i]);
    }
    else if (reconciliationModel_ == "COAL") {
      dynamic_cast<COALGeneTreeLikelihood*> (treeLikelihoods_[i])->refineGeneTreeNNIs(allParams_[i]);
    }
    if (timing && rearrange_) 
    {
      totalTime = ApplicationTools::getTime() - startingTime;
      std::cout << "Family "<< assignedFilenames_[i] <<"; Time for NNI exploration: "<<  totalTime << " s." <<std::endl;
    }
    return;
}