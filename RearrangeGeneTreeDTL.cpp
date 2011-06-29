/*
 *  RearrangeGeneTreeDTL.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 29/06/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

/******************************************************************************/

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


/******************************************************************************/
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
    std::cout << "*       Bio++ Gene Tree Rearrangement Program, version 1.0       *" << std::endl;
    std::cout << "* Authors: G. Szollosi, B. Boussau            Created 29/06/2011 *" << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;
    
    
    
    ///////////////////////////////
    //Gene family-related options:
    ///////////////////////////////
    
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
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
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
      }
      catch (std::exception& e)
      {
      std::cout << e.what() <<"; Unable to get a proper gene tree for family "<<file<<"; avoiding this family."<<std::endl;
      numDeletedFamilies = numDeletedFamilies+1;
      avoidFamily=true;
      }
    }
      
    
    
    
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