//
// File: ReconcileDuplications.cpp
// Created by: Bastien Boussau
// Created on: Aug 04 23:36 2007
//
/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
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

// From the STL:
#include <iostream>
#include <iomanip>
#include <algorithm>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/HomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
//#include <Phyl/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/ClockTreeLikelihood.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>


// From NumCalc:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From Utils:
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/Text/KeyvalTools.h>


//#include <Phyl/SAHomogeneousTreeLikelihood.h>
#include "ReconciliationTools.h"
#include "ReconciliationTreeLikelihood.h"
#include "SpeciesTreeExploration.h"
#include "SpeciesTreeLikelihood.h"


//From the BOOST library 
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;

//#include "mpi.h" 
//using namespace std;
using namespace bpp;




/*
This program takes a species tree file, sequence alignments and files describing the relations between sequences and species.
To compile it, use the Makefile !
*/

/**************************************************************************
 * This function creates a sequence tree from a species tree and a std::map containing the link between the species and their sequences.
 **************************************************************************/

TreeTemplate<Node> * buildARandomSequenceTreeFromASpeciesTree (std::map <std::string, std::deque<std::string> > & spSeqs, TreeTemplate<Node> & tree, std::map <std::string, std::string> & spSelectedSeq) {

	TreeTemplate<Node> * seqTree ;
	seqTree = new TreeTemplate<Node>(tree);
	//seqTree = tree;
	spSelectedSeq.clear();

	for(std::map<std::string, std::deque<std::string> >::iterator it = spSeqs.begin(); it != spSeqs.end(); it++){
		int numSeq = (it->second).size();
		//    std::cout <<"Numseqs : "<<numSeq<<std::endl;
		if (numSeq > 0) {
		 std::string chosenSequence = (it->second).at(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(numSeq));
			int leafId = seqTree->getLeafId(it->first);
			seqTree->setNodeName(leafId, chosenSequence);
			spSelectedSeq.insert( make_pair(it->first,chosenSequence));
		}
		else {
			int leafId = seqTree->getLeafId(it->first);
			seqTree->getNode(leafId)->getFather()->removeSon(seqTree->getNode(leafId));
		}
	}
	return seqTree;

}


/**************************************************************************
 * This function creates a sequence tree from a species tree and the std::map containing the link between the species and the putative orthologous sequence.
 **************************************************************************/

TreeTemplate<Node> * buildASequenceTreeFromASpeciesTreeAndCorrespondanceMap (TreeTemplate<Node> & tree, std::map <std::string, std::string> & spSelectedSeq) {

	TreeTemplate<Node> * seqTree ;
	seqTree = new TreeTemplate<Node>(tree);

	std::vector <std::string> names = seqTree->getLeavesNames();
	for(std::vector < std::string >::iterator it = names.begin(); it != names.end(); it++){
		if (spSelectedSeq.count(*it)>0) {
		 std::map < std::string, std::string >::iterator iter = spSelectedSeq.find(*it);
		 std::string chosenSequence = iter->second;
			int leafId = seqTree->getLeafId(*it);
			seqTree->setNodeName(leafId, chosenSequence);
		}
		else {
			int leafId = seqTree->getLeafId(*it);
			seqTree->getNode(leafId)->getFather()->removeSon(seqTree->getNode(leafId));
		}
	}

	return seqTree;

}


/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "ReconcileDuplications parameter1_name=parameter1_value parameter2_name=parameter2_value"   ).endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  
  /*SequenceApplicationTools::printInputAlignmentHelp();
    PhylogeneticsApplicationTools::printInputTreeHelp();
    PhylogeneticsApplicationTools::printSubstitutionModelHelp();
    PhylogeneticsApplicationTools::printRateDistributionHelp();
    PhylogeneticsApplicationTools::printCovarionModelHelp();
    PhylogeneticsApplicationTools::printOptimizationHelp(true, false);
    PhylogeneticsApplicationTools::printOutputTreeHelp();*/
  (*ApplicationTools::message << "output.infos                      | file where to write site infos").endLine();
  (*ApplicationTools::message << "output.estimates                  | file where to write estimated parameter values").endLine();
  (*ApplicationTools::message << "starting.tree.file                | file where to write the initial tree").endLine();
  (*ApplicationTools::message << "init.species.tree                 | user or random").endLine();
  (*ApplicationTools::message << "species.tree.file                 | if the former is set to \"user\", path to a species tree" ).endLine();
  (*ApplicationTools::message << "species.names.file                | if instead it was set to \"random\", path to a list of species names ").endLine();
  (*ApplicationTools::message << "heuristics.level                  | 0, 1, 2, or 3; 0: DR exact algorithm (default); 1 fast heuristics; 2: exhaustive and fast; 3: exhaustive and slow").endLine();
  (*ApplicationTools::message << "species.id.limit.for.root.position| Threshold for trying root positions").endLine();
  (*ApplicationTools::message << "genelist.file                     | file containing a list of gene option files to analyse").endLine();
  (*ApplicationTools::message << "branchProbabilities.optimization  | average, branchwise, average_then_branchwise or no: how we optimize duplication and loss probabilities").endLine();
  (*ApplicationTools::message << "genome.coverage.file              | file giving the percent coverage of the genomes used").endLine();
  (*ApplicationTools::message << "spr.limit                         | integer giving the breadth of SPR movements, in number of nodes. 0.1* number of nodes in the species tree might be OK.").endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of supplementary options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}


/******************************************************************************/
/***************************Adding two std::vectors*********************************/


std::vector <int> operator + (std::vector <int> x, std::vector <int> y) {
	int temp=x.size();
	if (temp!=y.size()) {
		std::cout <<"problem : adding two std::vectors of unequal sizes."<<std::endl;
		exit (-1);
	}
 std::vector<int> result;
	for(int i = 0 ; i<temp ; i++){
		result.push_back(x[i]+y[i]);
	}
	return result;
}

std::vector <double> operator + (std::vector <double> x, std::vector <double> y) {
	int temp=x.size();
	if (temp!=y.size()) {
		std::cout <<"problem : adding two std::vectors of unequal sizes."<<std::endl;
		exit (-1);
	}
 std::vector<double> result;
	for(int i = 0 ; i<temp ; i++){
		result.push_back(x[i]+y[i]);
	}
	return result;
}




/******************************************************************************/
// This function parses a vector of gene family files, 
// discards the ones that do not pass certain criteria, 
// and initializes the others. Used by clients.
/******************************************************************************/


void parseAssignedGeneFamilies(std::vector<std::string> & assignedFilenames, std::map<std::string, std::string> & params, int & numDeletedFamilies, TreeTemplate<Node> * geneTree, TreeTemplate<Node> * tree, std::vector <int> & num0Lineages, std::vector <std::vector<int> > & allNum0Lineages, std::vector <std::vector<int> > & allNum1Lineages, std::vector <std::vector<int> > & allNum2Lineages, std::vector <double> & lossExpectedNumbers, std::vector <double> & duplicationExpectedNumbers, std::map <std::string, int> & spId, int speciesIdLimitForRootPosition, int heuristicsLevel, int MLindex, std::vector <double> & allLogLs, std::vector <ReconciliationTreeLikelihood *> & treeLikelihoods, std::vector <ReconciliationTreeLikelihood *> & backupTreeLikelihoods, std::vector <std::map<std::string, std::string> > & allParams, std::vector <Alphabet *> allAlphabets, std::vector <VectorSiteContainer *> & allDatasets, std::vector <SubstitutionModel *> & allModels, std::vector <DiscreteDistribution *> & allDistributions, std::vector <TreeTemplate<Node> *> & allGeneTrees, std::vector <TreeTemplate<Node> *> & allUnrootedGeneTrees) 
{
  bool avoidFamily;
  std::string initTree;
  std::map<std::string, std::string> famSpecificParams;
  std::vector <std::string> spNames;


  //Here we are going to get all necessary information regarding all gene families the client is in charge of.
  for (int i = 0 ; i< assignedFilenames.size() ; i++) 
    { //For each file
      std::cout <<"Examining family "<<assignedFilenames[i]<<std::endl;
      avoidFamily = false;
      std::string file =assignedFilenames[i];
      TreeTemplate<Node> * unrootedGeneTree = NULL;
      SubstitutionModel*    model    = 0;
      SubstitutionModelSet* modelSet = 0;
      DiscreteDistribution* rDist    = 0;
      if(!FileTools::fileExists(file))
        {
        std::cerr << "Parameter file "<< file <<" not found." << std::endl;
        exit(-1);
        }
      else
        {
        famSpecificParams = AttributesTools::getAttributesMapFromFile(file, "=");
        AttributesTools::resolveVariables(famSpecificParams);
        }
      
      AttributesTools::actualizeAttributesMap(params, famSpecificParams);
      
      //Sequences and model of evolution
      Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);
      VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);       
      VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, params);     
      delete allSites;             
      //method to optimize the gene tree root; only useful if heuristics.level!=0.
      bool rootOptimization;
      if (ApplicationTools::getStringParameter("root.optimization",params,"normal")=="intensive") 
        {
        rootOptimization=true;
        }
      else 
        {
        rootOptimization=false;
        }
      
      /****************************************************************************
       //Then we need to get the file containing links between sequences and species.
       *****************************************************************************/
      std::string taxaseqFile = ApplicationTools::getStringParameter("taxaseq.file",params,"none");
      if (taxaseqFile=="none" ){
        std::cout << "\n\nNo taxaseqfile was provided. Cannot compute a reconciliation between a species tree and a gene tree using sequences if the relation between the sequences and the species is not explicit !\n" << std::endl;
        std::cout << "ReconcileDuplications species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
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
          spSeq.insert( make_pair(st1.getToken(0),st2.getTokens()));
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
        std::cout <<std::endl;
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
          geneTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTreeFile));
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
          DistanceEstimation distEstimation(model, rDist, sites, 1, false);
          BioNJ * bionj = new BioNJ();      
          bionj->outputPositiveLengths(true);   
          std::string type = ApplicationTools::getStringParameter("bionj.optimization.method", params, "init");
          if(type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
          else if(type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
          else if(type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
          else throw Exception("Unknown parameter estimation procedure '" + type + "'.");
          // Should I ignore some parameters?
          ParameterList allParameters = model->getParameters();
          allParameters.addParameters(rDist->getParameters());
          ParameterList parametersToIgnore;
          std::string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", "", true, false);
          bool ignoreBrLen = false;
          StringTokenizer st(paramListDesc, ",");
          while(st.hasMoreToken())
            {
            try
              {
              std::string param = st.nextToken();
              if(param == "BrLen")
                ignoreBrLen = true;
              else
                {
                if (allParameters.hasParameter(param))
                  {
                  Parameter* p = &allParameters.getParameter(param);
                  parametersToIgnore.addParameter(*p);
                  }
                else ApplicationTools::displayWarning("Parameter '" + param + "' not found."); 
                }
              } 
            catch(ParameterNotFoundException pnfe)
              {
              ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
              }
            }
          double tolerance = ApplicationTools::getDoubleParameter("bionj.optimization.tolerance", params, .000001);
          unrootedGeneTree = OptimizationTools::buildDistanceTree(distEstimation, *bionj, parametersToIgnore, !ignoreBrLen, false, type, tolerance);
          std::vector<Node*> nodes = unrootedGeneTree->getNodes();
          for(unsigned int k = 0; k < nodes.size(); k++)
            {
            if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
            if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() >= 5) throw Exception("Found a very large branch length in a gene tree; will avoid this family.");
            }
          if (initTree == "phyml")//refine the tree using PhyML algorithm (2003)
            { 
              std::string backupParamOptTopo = params[ std::string("optimization.topology.algorithm_nni.method")];
              params[ std::string("optimization.topology.algorithm_nni.method")] = "phyml";
              NNIHomogeneousTreeLikelihood * tl = new NNIHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, true, true);
              tl->initialize();//Only initializes the parameter list, and computes the likelihood
              double logL = tl->getValue();
              if (std::isinf(logL))
                {
                // This may be due to null branch lengths, leading to null likelihood!
                ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
                ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
                ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
                ParameterList pl = tl->getBranchLengthsParameters();
                for (unsigned int i = 0; i < pl.size(); i++)
                  {
                  if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
                  }
                tl->matchParametersValues(pl);
                logL = tl->getValue();
                }
              ApplicationTools::displayResult("Initial log likelihood for family", file+": "+TextTools::toString(-logL, 15)); 
              if (std::isinf(logL))
                {
                ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
                CodonAlphabet *pca=dynamic_cast<CodonAlphabet*>(alphabet);
                if (pca){
                  bool f=false;
                  unsigned int  s;
                  for (unsigned int i = 0; i < sites->getNumberOfSites(); i++){
                    if (std::isinf(tl->getLogLikelihoodForASite(i))){
                      const Site& site=sites->getSite(i);
                      s=site.size();
                      for (unsigned int j=0;j<s;j++)
                        if (pca->isStop(site.getValue(j))){
                          (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                          f=true;
                        }
                    }
                  }
                  if (f)
                    exit(-1);
                }
                ApplicationTools::displayError("!!! Looking at each site:");
                for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
                  {
                  (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
                  }
                ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
                exit(-1);
                }
                           tl = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
                                                                           PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params));
              //PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params, "", true, false);
              params[ std::string("optimization.topology.algorithm_nni.method")] = backupParamOptTopo ;
              delete unrootedGeneTree;
              unrootedGeneTree = new TreeTemplate<Node> ( tl->getTree() );
              delete tl;
            }
          geneTree = unrootedGeneTree;
          geneTree->newOutGroup(0);
          delete bionj;
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
      
      if (!avoidFamily) 
        { //This family is phylogenetically informative
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
              std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene tree."<<std::endl;
            }
          
          //printing the gene trees with the species names instead of the sequence names
          //This is useful to build an input for duptree for instance
          TreeTemplate<Node> * treeWithSpNames = unrootedGeneTree->clone();
          std::vector <Node *> leaves = treeWithSpNames->getLeaves();
          for (int j =0; j<leaves.size() ; j++) 
            {
            leaves[j]->setName(seqSp[leaves[j]->getName()]);
            }
          std::vector <Node *> nodes =  treeWithSpNames->getNodes();
          for (int j =0; j<nodes.size() ; j++) 
            {
            if (nodes[j]->hasFather()) 
              {
              nodes[j]->deleteDistanceToFather(); 
              }
            if (nodes[j]->hasBootstrapValue()) 
              {
              nodes[j]->removeBranchProperty(TreeTools::BOOTSTRAP); 
              }
            }
          Newick newick(true);
          std::string geneTreeFile =ApplicationTools::getStringParameter("gene.tree.file",params,"none");
          newick.write(*treeWithSpNames, geneTreeFile+"SPNames", true);
          delete treeWithSpNames;
          
          //Outputting the starting tree
          std::string startingGeneTreeFile =ApplicationTools::getStringParameter("output.starting.gene.tree.file",params,"none");
          newick.write(*geneTree, startingGeneTreeFile, true);
          
          /****************************************************************************
           //Then we initialize the losses and duplication numbers on this tree.
           *****************************************************************************/
          std::vector<int> numbers = num0Lineages;
          allNum0Lineages.push_back(numbers);
          allNum1Lineages.push_back(numbers);
          allNum2Lineages.push_back(numbers);
          resetLossesAndDuplications(*tree, lossExpectedNumbers, duplicationExpectedNumbers);
          resetVector(allNum0Lineages[i-numDeletedFamilies]);
          resetVector(allNum1Lineages[i-numDeletedFamilies]);
          resetVector(allNum2Lineages[i-numDeletedFamilies]);
          
          /************************************************************************************************************/
          /********************************************COMPUTING LIKELIHOOD********************************************/
          /************************************************************************************************************/
          bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", params, true, "", false, false);
          if(!computeLikelihood)
            {
            delete alphabet;
            delete sites;
            delete tree;
            std::cout << "ReconcileDuplication's done. Bye." << std::endl;
            exit(-1);
            }
          
          DiscreteRatesAcrossSitesTreeLikelihood* tl;
          
          
          
          std::string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", params, "no", "", true, false);
          ApplicationTools::displayResult("Clock", optimizeClock);
          
          
          if(optimizeClock == "global")
            {
            std::cout<<"Sorry, clocklike trees have not been implemented yet."<<std::endl;
            exit(0);
            }// This has not been implemented!
          else if(optimizeClock == "no")
            {
            tl = new ReconciliationTreeLikelihood(*unrootedGeneTree, *sites, 
                                                  model, rDist, *tree, 
                                                  *geneTree, seqSp, spId, 
            /*allLossNumbers[i-numDeletedFamilies], */
                                                  lossExpectedNumbers, 
            /*allDuplicationNumbers[i-numDeletedFamilies], */
                                                  duplicationExpectedNumbers, 
            /*allBranchNumbers[i-numDeletedFamilies], */
                                                  allNum0Lineages[i-numDeletedFamilies], 
                                                  allNum1Lineages[i-numDeletedFamilies], 
                                                  allNum2Lineages[i-numDeletedFamilies], 
                                                  speciesIdLimitForRootPosition, 
                                                  heuristicsLevel, MLindex, 
                                                  true, true, rootOptimization);
            }
          else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
          
          tl->initialize();//Only initializes the parameter list, and computes the likelihood through fireParameterChanged
          allLogLs.push_back(tl->getValue());
          if(std::isinf(allLogLs[i-numDeletedFamilies]))
            {
            // This may be due to null branch lengths, leading to null likelihood!
            ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
            ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
            ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
            std::vector<Node*> nodes = tree->getNodes();
            for(unsigned int k = 0; k < nodes.size(); k++)
              {
              if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
              }
            dynamic_cast<ReconciliationTreeLikelihood*>(tl)->initParameters();
            
            allLogLs[i-numDeletedFamilies]= tl->f(tl->getParameters());
            }
          ApplicationTools::displayResult("Initial likelihood", TextTools::toString(allLogLs[i-numDeletedFamilies], 15));
          if(std::isinf(allLogLs[i-numDeletedFamilies]))
            {
            ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
            ApplicationTools::displayError("!!! Looking at each site:");
            for(unsigned int k = 0; k < sites->getNumberOfSites(); k++)
              {
              (*ApplicationTools::error << "Site " << sites->getSite(k).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(k)).endLine();
              }
            ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularly if datasets are big (>~500 sequences).");
            exit(-1);
            }
          treeLikelihoods.push_back(dynamic_cast<ReconciliationTreeLikelihood*>(tl));
          allParams.push_back(params); 
          allAlphabets.push_back(alphabet);
          allDatasets.push_back(sites);
          allModels.push_back(model);
          allDistributions.push_back(rDist);
          allGeneTrees.push_back(geneTree);
          allUnrootedGeneTrees.push_back(unrootedGeneTree);
          // allLikelihoodFiles.push_back(likelihoodFile);
        }
      else 
        {
        delete sites;
        delete alphabet;
        }
      std::cout <<"Examined family "<<assignedFilenames[i]<<std::endl;
    }//End for each file
  if (numDeletedFamilies == assignedFilenames.size()) 
    {
    std::cout<<"WARNING: A client is in charge of 0 gene family after gene family filtering!"<<std::endl;        
    std::cout<<"A processor will be idle most of the time, the load could probably be better distributed."<<std::endl; 
    }
  
}


/******************************************************************************/
// This function outputs gene trees from the clients.
/******************************************************************************/
void outputGeneTrees (std::vector<std::string> & assignedFilenames, 
                      std::vector <ReconciliationTreeLikelihood *> treeLikelihoods, 
                      std::vector< std::map<std::string, std::string> > & allParams, 
                      int & numDeletedFamilies, 
                      std::vector <std::vector <std::string> > & reconciledTrees,                       
                      std::vector <std::vector <std::string> > & duplicationTrees, 
                      std::vector <std::vector <std::string> > & lossTrees, 
                      int & bestIndex, int & startRecordingTreesFrom)
{
  std::string suffix;
  std::string reconcTree;
  std::ofstream out;
  std::string dupTree;
  std::string lossTree;
  int rightIndex = bestIndex-startRecordingTreesFrom;
  for (int i = 0 ; i< assignedFilenames.size()-numDeletedFamilies ; i++) 
    {

      suffix = ApplicationTools::getStringParameter("output.file.suffix", allParams[i], "", "", false, false);
      reconcTree = ApplicationTools::getStringParameter("output.reconciled.tree.file", allParams[i], "reconciled.tree", "", false, false);
      reconcTree = reconcTree + suffix;
      Nhx *nhx = new Nhx();
      string temp = reconciledTrees[i][rightIndex];
      TreeTemplate<Node> * geneTree=nhx->parenthesisToTree(temp);
      temp = duplicationTrees[i][rightIndex];
      TreeTemplate<Node> * spTree=nhx->parenthesisToTree( temp);
      breadthFirstreNumber (*spTree);
      std::map <std::string, int> spId = computeSpeciesNamesToIdsMap(*spTree);
      std::map <std::string, std::string> seqSp = treeLikelihoods[i]->getSeqSp();
      annotateGeneTreeWithDuplicationEvents (*spTree, 
                                             *geneTree, 
                                             geneTree->getRootNode(), 
                                             treeLikelihoods[i]->getSeqSp(),
                                             spId); 
    out.open (reconcTree.c_str(), std::ios::out);
    nhx->write(*geneTree, out);
    out.close();
    dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", allParams[i], "duplications.tree", "", false, false);
    dupTree = dupTree + suffix;
    out.open (dupTree.c_str(), std::ios::out);
    out << duplicationTrees[i][rightIndex]<<std::endl;
    out.close();
    lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", allParams[i], "losses.tree", "", false, false);
    lossTree = lossTree + suffix;
    out.open (lossTree.c_str(), std::ios::out);
    out << lossTrees[i][rightIndex]<<std::endl;
    out.close();
    }
  return;
}






/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
//////////////////////////////////////////////////////MAIN/////////////////////////////////////////////////
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


int main(int args, char ** argv)
{
	if(args == 1)
	{
		help();
		exit(0);
	}

	int rank, size;
	int server = 0;

	//Using BOOST :
	mpi::environment env(args, argv);
	mpi::communicator world;
	rank = world.rank();
	size = world.size();
	std::string line;

  if (size==1) {
    std::cout <<"\n\n\n\t\tError: this program can only run if 2 or more processes are used."<<std::endl;
    std::cout <<"\t\tUse 'mpirun -np k ReconcileDuplications ...', where k>=2"<<std::endl;
    exit(-1);
  }
  
  
	try {
		ApplicationTools::startTimer();
		std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);

    
    //##################################################################################################################
    //##################################################################################################################
    //############################################# IF AT THE SERVER NODE ##############################################
    //##################################################################################################################
    //##################################################################################################################
        
		if (rank == server) { 
      
      std::cout << "******************************************************************" << std::endl;
      std::cout << "*       Bio++ Tree Reconciliation Program, version 1.1.0      *" << std::endl;
      std::cout << "* Author: B. Boussau                        Created 16/07/07 *" << std::endl;
      std::cout << "******************************************************************" << std::endl;
      std::cout << std::endl;
      
      
      
      SpeciesTreeLikelihood spTL = SpeciesTreeLikelihood(world, server, size, params);
      spTL.initialize();
      
      spTL.MLsearch();
          
      
			std::cout << "ReconcileDuplication's done. Bye." << std::endl;
			ApplicationTools::displayTime("Total execution time:");
		}//End if at the server node
    
    
    
    
    
    
		//##################################################################################################################
    //##################################################################################################################
    //############################################# IF AT A CLIENT NODE ################################################
    //##################################################################################################################
    //##################################################################################################################
    if (rank >server)
      {

      /****************************************************************************
       * First communications between the server and the clients.
       *****************************************************************************/      
      bool optimizeSpeciesTreeTopology;
      int assignedNumberOfGenes;
      std::vector <std::string> assignedFilenames;
      std::vector <int> numbersOfGenesPerClient;
      std::vector <std::vector<std::string> > listOfOptionsPerClient;
      int SpeciesNodeNumber;
      std::vector <double> lossExpectedNumbers;
      std::vector <double> duplicationExpectedNumbers;
      std::vector <double> backupLossExpectedNumbers;
      std::vector <double> backupDuplicationExpectedNumbers;
      std::vector <int> num0Lineages;
      std::vector <int> num1Lineages; 
      std::vector <int> num2Lineages;
      std::vector <std::vector<int> > allNum0Lineages;
      std::vector <std::vector<int> > allNum1Lineages;
      std::vector <std::vector<int> > allNum2Lineages;
      std::string currentSpeciesTree;

      firstCommunicationsServerClient (world , server, numbersOfGenesPerClient, assignedNumberOfGenes,
                                       assignedFilenames, listOfOptionsPerClient, 
                                       optimizeSpeciesTreeTopology, SpeciesNodeNumber, 
                                       lossExpectedNumbers, duplicationExpectedNumbers, num0Lineages, 
                                       num1Lineages, num2Lineages, currentSpeciesTree);
      
      /****************************************************************************
       * Various initializations.
       *****************************************************************************/  
      TreeTemplate<Node> * tree = 0;
      
      
      /*std::vector <int> lossNumbers;
       std::vector <int> duplicationNumbers;
       std::vector <int> branchNumbers;*/
      
      
      //  char currentSpeciesTree[MAXSPECIESTREESIZE];
      std::string initTree;
      std::string allFileNames;
      
			//First we read the species tree from the char[] sent by the server
      tree=TreeTemplateTools::parenthesisToTree(currentSpeciesTree, false, "", true);
			resetLossesAndDuplications(*tree, lossExpectedNumbers, duplicationExpectedNumbers);
      //To make the correspondance between species name and id:
      std::map <std::string, int> spId = computeSpeciesNamesToIdsMap(*tree);
      //These std::vectors contain all relevant numbers for all gene families the client is in charge of.
			std::vector <double> allLogLs;
      double logL;
			std::vector <std::map<std::string, std::string> > allParams;
			TreeTemplate<Node> * geneTree = NULL;
			int MLindex = 0;
			std::vector <ReconciliationTreeLikelihood *> treeLikelihoods;
			std::vector <ReconciliationTreeLikelihood *> backupTreeLikelihoods;
			std::cout <<"Client  of rank "<<rank <<" is in charge of " << assignedFilenames.size()<<" gene families."<<std::endl;
			std::vector <Alphabet *> allAlphabets;
			std::vector <VectorSiteContainer *>   allDatasets;
			std::vector <SubstitutionModel *> allModels;
			std::vector <DiscreteDistribution *> allDistributions;
      std::vector <TreeTemplate<Node> *> allGeneTrees;
			std::vector <TreeTemplate<Node> *> allUnrootedGeneTrees;
      int numDeletedFamilies=0;
      int bestIndex = 0;
      bool stop = false; 
      bool rearrange = false;

      
      /****************************************************************************
			 // Then, we get 3 options related to the algorithms used.
			 * Meaning of heuristicsLevel:
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
       * WARNING: options 1, 2, 3 are probably bugged now.
			 *****************************************************************************/
      int heuristicsLevel= ApplicationTools::getIntParameter("heuristics.level", params, 0, "", false, false);
			int speciesIdLimitForRootPosition= ApplicationTools::getIntParameter("species.id.limit.for.root.position", params, 2, "", false, false);
			double optimizationTolerance=ApplicationTools::getDoubleParameter("optimization.tolerance", params, 0.000001, "", false, false);
      
      
      /****************************************************************************
       * Gene family parsing and first likelihood computation.
       *****************************************************************************/            
      parseAssignedGeneFamilies(assignedFilenames, params, numDeletedFamilies, 
                                geneTree, tree, num0Lineages, allNum0Lineages, 
                                allNum1Lineages, allNum2Lineages, 
                                lossExpectedNumbers, duplicationExpectedNumbers, 
                                spId, speciesIdLimitForRootPosition, heuristicsLevel, 
                                MLindex, allLogLs, treeLikelihoods, 
                                backupTreeLikelihoods, allParams, 
                                allAlphabets, allDatasets, allModels, 
                                allDistributions, allGeneTrees, allUnrootedGeneTrees );
			std::vector <std::vector <std::string> > reconciledTrees;
			std::vector <std::vector <std::string> > duplicationTrees;
			std::vector <std::vector <std::string> > lossTrees;
			std::vector <std::string> t;  
			std::vector <std::map<std::string, std::string> > allParamsBackup = allParams;
			for (int i = 0 ; i< assignedFilenames.size()-numDeletedFamilies ; i++) 
        {
				reconciledTrees.push_back(t);
				duplicationTrees.push_back(t);
				lossTrees.push_back(t);
        //This is to avoid optimizing gene tree parameters in the first steps of the program, 
        //if we optimize the species tree topology.
        if (optimizeSpeciesTreeTopology) 
          { 
            if (ApplicationTools::getBooleanParameter("optimization.topology", allParams[i], false, "", true, false))
              {
              allParams[i][ std::string("optimization.topology")] = "false";
              }
            allParams[i][ std::string("optimization")] = "None"; //Quite extreme, but the sequence likelihood has no impact on the reconciliation !
            treeLikelihoods[i]->OptimizeSequenceLikelihood(false);
          }
          std::cout <<"Continued family "<<assignedFilenames.size()-numDeletedFamilies<<std::endl;
        }
			bool recordGeneTrees = false; //At the beginning, we do not record the gene trees.
			int startRecordingTreesFrom = 0; //This int is incremented until the gene trees start to be backed-up, when we start the second phase of the algorithm.
      bool firstTimeImprovingGeneTrees = false; //When for the first time we optimize gene trees, we set it at true
      //We make a backup of the gene tree likelihoods.
      for (int i =0 ; i<treeLikelihoods.size() ; i++) 
        {
        backupTreeLikelihoods.push_back(treeLikelihoods[i]->clone());
        }      
 
      /****************************************************************************
       ****************************************************************************
       * Main loop: iterative likelihood computations
       ****************************************************************************
       *****************************************************************************/
        Nhx *nhx = new Nhx();
      while (!stop)            
        {      
          logL=0.0;
          resetVector(num0Lineages);
          resetVector(num1Lineages);
          resetVector(num2Lineages);
          for (int i = 0 ; i< assignedFilenames.size()-numDeletedFamilies ; i++) 
            {
            if (firstTimeImprovingGeneTrees) 
              {
              treeLikelihoods[i]->OptimizeSequenceLikelihood(true);
              backupTreeLikelihoods[i]->OptimizeSequenceLikelihood(true);
              }
           // std::cout <<  TreeTools::treeToParenthesis(*geneTree, true)<<std::endl;
            PhylogeneticsApplicationTools::optimizeParameters(treeLikelihoods[i], treeLikelihoods[i]->getParameters(), allParams[i], "", true, false); 
            geneTree = new TreeTemplate<Node>(treeLikelihoods[i]->getRootedTree()); //MODIFICATION

            ///LIKELIHOOD OPTIMIZED
                        
            resetLossesAndDuplications(*tree, lossExpectedNumbers, duplicationExpectedNumbers);
            allNum0Lineages[i] = treeLikelihoods[i]->get0LineagesNumbers();
            allNum1Lineages[i] = treeLikelihoods[i]->get1LineagesNumbers();
            allNum2Lineages[i] = treeLikelihoods[i]->get2LineagesNumbers();
            MLindex = treeLikelihoods[i]->getRootNodeindex();
            allLogLs[i] = treeLikelihoods[i]->getValue();  
            logL = logL + allLogLs[i];
            num0Lineages = num0Lineages + allNum0Lineages[i];
            num1Lineages = num1Lineages + allNum1Lineages[i];
            num2Lineages = num2Lineages + allNum2Lineages[i];
            std::cout<<"Gene Family: " <<allParams[i][ std::string("taxaseq.file")] << " logLk: "<< allLogLs[i]<< " scenario lk: "<< treeLikelihoods[i]->getScenarioLikelihood() <<std::endl;
            if (std::isnan(allLogLs[i])) 
              {
              std::cout<<TreeTools::treeToParenthesis (*geneTree, false, EVENT)<<std::endl;
              std::cout<<TreeTools::treeToParenthesis (*tree, false, DUPLICATIONS)<<std::endl;
              }
            if (recordGeneTrees) 
              {
                reconciledTrees[i].push_back(nhx->treeToParenthesis (*geneTree));
                duplicationTrees[i].push_back(nhx->treeToParenthesis (*tree));
                lossTrees[i].push_back(nhx->treeToParenthesis (*tree));
                /*
              reconciledTrees[i].push_back(TreeTools::treeToParenthesis (*geneTree, false, EVENT));
              duplicationTrees[i].push_back(treeToParenthesisWithIntNodeValues (*tree, false, DUPLICATIONS));
              lossTrees[i].push_back(treeToParenthesisWithIntNodeValues (*tree, false, LOSSES));*/
              }
            if (geneTree) 
              {
              delete geneTree;
              }
            }//end for each filename
          if (firstTimeImprovingGeneTrees) 
            {
            firstTimeImprovingGeneTrees = false;
            }
          if (!recordGeneTrees) 
            {
            startRecordingTreesFrom++;
            }
          //Clients send back stuff to the server.
          gather(world, logL, server); 
          gather(world, num0Lineages, allNum0Lineages, server); 
          gather(world, num1Lineages, allNum1Lineages, server);   
          gather(world, num2Lineages, allNum2Lineages, server);	 
          //Should the computations stop? The server tells us.
          broadcast(world, stop, server);
          if (!stop) 
            {	// we continue the loop
              //Reset the gene trees by resetting treeLikelihoods:
              //we always start from ML trees according to sequences only
              for (int i =0 ; i<treeLikelihoods.size() ; i++) 
                delete treeLikelihoods[i];
              treeLikelihoods.clear();
              for (int i=0 ; i<backupTreeLikelihoods.size() ; i++) 
                {
                treeLikelihoods.push_back(backupTreeLikelihoods[i]->clone());
                }          
              if (rearrange) 
                {
                allParams = allParamsBackup;
                if (recordGeneTrees==false) 
                  {
                  firstTimeImprovingGeneTrees = true;
                  recordGeneTrees=true;
                  bestIndex=startRecordingTreesFrom;
                  }
                }	
              broadcastsAllInformationButStop(world, server, rearrange, 
                                              lossExpectedNumbers, 
                                              duplicationExpectedNumbers, 
                                              currentSpeciesTree);
              TreeTemplate<Node> * tree=TreeTemplateTools::parenthesisToTree(currentSpeciesTree, false, "", true);
              spId = computeSpeciesNamesToIdsMap(*tree);
              for (int i = 0 ; i< assignedFilenames.size()-numDeletedFamilies ; i++) 
                {
                treeLikelihoods[i]->setSpTree(*tree);
                treeLikelihoods[i]->setSpId(spId);
                treeLikelihoods[i]->setProbabilities(duplicationExpectedNumbers, lossExpectedNumbers);
                treeLikelihoods[i]->computeTreeLikelihood();
                }
              if (tree) delete tree;
            }
          else 
            { 
              /****************************************************************************
               * The end, outputting the results.
               *****************************************************************************/      
              if (recordGeneTrees) 
                {
                broadcast(world, bestIndex, server);
                std::cout << "bestIndex: "<<bestIndex<<" startRecordingTreesFrom: "<<startRecordingTreesFrom<<std::endl;
                outputGeneTrees(assignedFilenames, treeLikelihoods, allParams, numDeletedFamilies, 
                                reconciledTrees, duplicationTrees, lossTrees, 
                                bestIndex, startRecordingTreesFrom);
                }
              break;
            }
        }//End while, END OF MAIN LOOP
	/*		for (int i = 0 ; i< assignedFilenames.size()-numDeletedFamilies ; i++) 
        {  
          delete allAlphabets[i];
          delete allDatasets[i];
          delete allModels[i];
          delete allDistributions[i];
        }*/
    }//end if a client node
	}
	catch(std::exception & e)
	{
		std::cout << e.what() << std::endl;
		exit(-1);
	}
	return (0);
}










/*
 // This bit of code is useful to use GDB on clients, when put into the client's code:
 //launch the application, which will output the client pid
 //then launch gdb, attach to the given pid ("attach pid" or "gdb ReconcileDuplications pid"), 
 //use "up" to go up the stacks, and set the variable z to !=0 to get out of the loop with "set var z = 8".
 int z = 0;
 //   char hostname[256];
 //gethostname(hostname, sizeof(hostname));
 std::cout <<"PID: "<<getpid()<<std::endl;
 std::cout <<"z: "<<z<<std::endl;
 //printf("PID %d on %s ready for attach\n", getpid(), hostname);
 // fflush(stdout);
 while (0 == z){
 std::cout <<z<<std::endl;
 sleep(5);
 }
 */





/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////RUBBISH////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//ATTEMPT 18 02 2010
//OPTIMIZE RATES
/*  SpeciesTreeLikelihood *slk = new SpeciesTreeLikelihood (world, server, tree, index, bestIndex, stop, logL, bestlogL, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossExpectedNumbers, duplicationExpectedNumbers, rearrange, numIterationsWithoutImprovement, branchProbaOptimization, genomeMissing);
 
 slk->getParameters().printParameters (std::cout);
 
 
 //Optimization of rates:
 PowellMultiDimensions * optimizer = new PowellMultiDimensions(slk);
 optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
 optimizer->setProfiler(NULL);
 optimizer->setMessageHandler(NULL);
 optimizer->setVerbose(0);
 
 optimizer->getStopCondition()->setTolerance(0.1);
 //optimizer_->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
 optimizer->init(slk->getParameters()); 
 optimizer->optimize(); 
 
 std::cout <<"optimization Done: "<<logL<<std::endl;
 
 lossExpectedNumbers = lossExpectedNumbers * optimizer->getParameters().getParameter("coefLoss").getValue();
 duplicationExpectedNumbers = duplicationExpectedNumbers * optimizer->getParameters().getParameter("coefDup").getValue();
 
 std::cout <<"Dup coef: "<<optimizer->getParameters().getParameter("coefDup").getValue() <<" Loss coef: "<<optimizer->getParameters().getParameter("coefLoss").getValue()<<std::endl;
 
 
 delete optimizer;
 delete slk; */


/*
 //ATTEMPT 18 02 2010
 //OPTIMIZE RATES
 SpeciesTreeLikelihood *slk = new SpeciesTreeLikelihood (world, server, tree, index, bestIndex, stop, logL, bestlogL, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossExpectedNumbers, duplicationExpectedNumbers, rearrange, numIterationsWithoutImprovement, branchProbaOptimization, genomeMissing);
 
 slk->getParameters().printParameters (std::cout);
 
 //Optimization:
 PowellMultiDimensions * optimizer = new PowellMultiDimensions(slk);
 optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
 optimizer->setProfiler(NULL);
 optimizer->setMessageHandler(NULL);
 optimizer->setVerbose(0);
 
 optimizer->getStopCondition()->setTolerance(0.1);
 //optimizer_->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
 optimizer->init(slk->getParameters()); 
 optimizer->optimize(); 
 
 std::cout <<"optimization Done: "<<logL<<std::endl;
 
 lossExpectedNumbers = lossExpectedNumbers * optimizer->getParameters().getParameter("coefLoss").getValue();
 duplicationExpectedNumbers = duplicationExpectedNumbers * optimizer->getParameters().getParameter("coefDup").getValue();
 
 std::cout <<"Dup coef: "<<optimizer->getParameters().getParameter("coefDup").getValue() <<" Loss coef: "<<optimizer->getParameters().getParameter("coefLoss").getValue()<<std::endl;
 
 
 delete optimizer;
 delete slk;
 */



/*This code permits outputting trees with numbers of duplication or loss events at branches
 //set total numbers of loss and duplications on branches
 setLossesAndDuplications(*bestTree, bestLossNumbers, bestDuplicationNumbers);
 std::string dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", params, "AllDuplications.tree", "", false, false); 
 std::ofstream out (dupTree.c_str(), std::ios::out);
 out << treeToParenthesisWithIntNodeValues (*bestTree, false, DUPLICATIONS)<<std::endl;
 out.close();
 std::string lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", params, "AllLosses.tree", "", false, false);
 out.open (lossTree.c_str(), std::ios::out);
 out << treeToParenthesisWithIntNodeValues (*bestTree, false, LOSSES)<<std::endl;
 out.close();
 */








