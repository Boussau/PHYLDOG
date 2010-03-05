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
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Phyl/HomogeneousTreeLikelihood.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
//#include <Phyl/NNIHomogeneousTreeLikelihood.h>
#include <Phyl/ClockTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/Newick.h>
#include <Phyl/TreeTools.h>
#include <Phyl/BioNJ.h>
#include <Phyl/OptimizationTools.h>


// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>
#include <NumCalc/RandomTools.h>
#include <NumCalc/NumConstants.h>
#include <NumCalc/PowellMultiDimensions.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/Clonable.h>
#include <Utils/Number.h>
#include <Utils/BppString.h>
#include <Utils/KeyvalTools.h>


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
 * Utilitary fonction
 *************************************************************************/

std::string removeComments(
                      const std::string & s,
                      const std::string & begin,
                      const std::string & end)
{
 std::string r = s;
 std::string::size_type last = 0;
	do
	{
	 std::string::size_type first = r.find(begin, last);
		if(first == std::string::npos) return r; //No shell comment.
		//else:  
		last = r.find(end, first);
		if(last == std::string::npos)
		{
			r.erase(r.begin() + first, r.end());
		}
		else
		{
			r.erase(r.begin() + first, r.begin() + last);
		}
	} while(last != std::string::npos);
	return r;
}


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


/**************************************************************************
 * This function distributes gene families so as to homogenize computational load between clients.
 **************************************************************************/

bool sortMaxFunction (std::pair <std::string, double> i, std::pair <std::string, double> j) 
  { 
    if (i.second > j.second ) 
      {
        return (true);
      }
    else 
      {
        return (false);
      }
  }

bool sortMinFunction (std::pair <std::vector<std::string>, double> i, std::pair <std::vector<std::string>, double> j) 
{ 
  if (i.second < j.second ) 
    {
      return (true);
    }
  else 
    {
      return (false);
    }
}




void generateListOfOptionsPerClient(std::vector <std::string> listOptions, int size, std::vector <std::vector<std::string> > &listOfOptionsPerClient, std::vector <int> &numberOfGenesPerClient) {
  std::vector <std::pair <std::string, double> > elements;
  for (int i = 0; i<listOptions.size() ; i++) {
    StringTokenizer st1 = StringTokenizer::StringTokenizer (listOptions[i], ":", true);
    elements.push_back(std::pair <std::string, double>(st1.getToken(0), TextTools::toDouble(st1.getToken(1))) );
  }

  //Now sort the gene families by their size, in descending order
  sort(elements.begin(), elements.end(), sortMaxFunction);
  //Now we assign gene families to nodes.
  //We start with big families, and then go in decreasing order.
  //First we assign the first families
  std::vector <std::pair <std::vector<std::string>, double> > listOfOptionsAndTotSizePerClient;

  std::vector<std::string> temp2;
  std::pair <std::vector<std::string>, double> temp3;
  int j = 0;

  for (int i = 0; i<size-1 ; i++) {
    listOfOptionsAndTotSizePerClient.push_back(temp3);
    listOfOptionsAndTotSizePerClient[i].first.push_back(elements[j].first);
    listOfOptionsAndTotSizePerClient[i].second = elements[j].second;  
    j = j+1;
  }

  while (j<listOptions.size()) {
    //We sort listOfOptionsAndTotSizePerClient in increasing function
    sort(listOfOptionsAndTotSizePerClient.begin(), listOfOptionsAndTotSizePerClient.end(), sortMinFunction);
    listOfOptionsAndTotSizePerClient[0].first.push_back(elements[j].first);
    listOfOptionsAndTotSizePerClient[0].second = listOfOptionsAndTotSizePerClient[0].second + elements[j].second;
    j= j+1;
  }

  //Now all gene families must have been assigned to nodes.
 // std::vector <std::vector<std::string> > listOfOptionsPerClient;
  listOfOptionsPerClient.push_back(temp2);
  listOfOptionsPerClient[0].push_back(std::string("####"));//For the root node
  numberOfGenesPerClient.push_back(0);
  
  //We print the result of the assignment and fill listOfOptionsPerClient:
  for (int i = 0; i<size-1 ; i++) {
    std::cout <<"Client "<<i<<" is in charge of "<< listOfOptionsAndTotSizePerClient[i].first.size()<<" gene families; Total Weight : "<<listOfOptionsAndTotSizePerClient[i].second<<std::endl;
    listOfOptionsPerClient.push_back(listOfOptionsAndTotSizePerClient[i].first);
    numberOfGenesPerClient.push_back(listOfOptionsAndTotSizePerClient[i].first.size());
  }
  return ;
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

/*********************************************************************************************************/
/*********************************************************************************************************/
//////////////////////////////////////////////////////MAIN/////////////////////////////////////////////////
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
	int MAXFILENAMESIZE = 500;
	int MAXSPECIESTREESIZE = 10000; //size of the species tree, in number of CHARs, as it is written in Newick format

	TreeTemplate<Node> * tree = NULL;
	std::vector <double> lossProbabilities;
  std::vector <double> duplicationProbabilities;
  std::vector <double> backupLossProbabilities;
  std::vector <double> backupDuplicationProbabilities;

  
  /*std::vector <int> lossNumbers;
	std::vector <int> duplicationNumbers;
	std::vector <int> branchNumbers;*/
	std::vector <int> num0Lineages;
	std::vector <int> num1Lineages; 
	std::vector <int> num2Lineages;

	//  char currentSpeciesTree[MAXSPECIESTREESIZE];
	std::string currentSpeciesTree;
	int SpeciesNodeNumber;
	std::string initTree;
	std::string allFileNames;
	std::vector <int> numbersOfGenesPerClient;
	std::vector <std::vector<std::string> > listOfOptionsPerClient;
	double optimizationTolerance;
	//Using BOOST :
	mpi::environment env(args, argv);
	mpi::communicator world;
	rank = world.rank();
	size = world.size();
	std::vector <std::string> spNames;
	std::string line;

	try {

		ApplicationTools::startTimer();

		//We use a std::vector to record the list of gene option files
		std::vector<std::string> listOptions;
		double logL = 0.0;  int z=0;
		std::map<std::string, std::string> params;
		//Matrices to store numbers of duplications and losses from the clients
		/*std::vector <std::vector<int> > AllDuplications;
		std::vector <std::vector<int> > AllLosses;
		std::vector <std::vector<int> > AllBranches;*/
		std::vector <std::vector<int> > allNum0Lineages;
		std::vector <std::vector<int> > allNum1Lineages;
		std::vector <std::vector<int> > allNum2Lineages;
		std::string affectedFilename;
		std::vector <std::string> affectedFilenames;
		int affectedNumberOfGenes;
		int heuristicsLevel;
		int speciesIdLimitForRootPosition;
		bool optimizeSpeciesTreeTopology;
		bool stop = false; 
		bool rearrange = false;
		int bestIndex = 0;
		std::string branchProbaOptimization;
		int sprLimit;

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
      
      
      
		 std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);
			/****************************************************************************
			 //First, we need to get the species tree.
			 *****************************************************************************/
			// Get the initial tree
			// This tree should be a rooted species tree

			initTree = ApplicationTools::getStringParameter("init.species.tree", params, "user", "", false, false);
			ApplicationTools::displayResult("Input species tree", initTree);
      // A given species tree
			if(initTree == "user")
			{
			 std::string spTreeFile =ApplicationTools::getStringParameter("species.tree.file",params,"none");
				if (spTreeFile=="none" ){
					std::cout << "\n\nNo Species tree was provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file. \nIf you do not have a species tree file, the program can start from a random tree, if you set init.species.tree at random, and give a list of species names as species.names.file\n\n" << std::endl;
					exit(-1);
				}
				ApplicationTools::displayResult("Species Tree file", spTreeFile);
				Newick newick(true);
				tree = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
				if (!tree->isRooted()) {
					std::cout << "The tree is not rooted, midpoint-rooting it!\n";
					TreeTools::midpointRooting(*tree);
				}
				ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
				spNames=tree->getLeavesNames();
			}
      //We build a random species tree
			else if(initTree == "random")
			{
        //Then we need a file containing species names
			 std::string spNamesFile =ApplicationTools::getStringParameter("species.names.file",params,"none");
				if (spNamesFile=="none" ){
					std::cout << "\n\nNo Species names were provided. The option init.species.tree is set to random, which means that the option species.names.file must be filled with the path of a file containing a list of species names. \n\n" << std::endl;
					exit(-1);
				}
        std::ifstream inListNames (spNamesFile.c_str());
			 std::string line;
				while(getline(inListNames,line)) {
					spNames.push_back(line);
				}
				int maxStrSize=0;
				std::vector <int> toRemove;
				int i=0;
				for(std::vector<std::string >::iterator it = spNames.begin(); it != spNames.end(); it++){
				 std::string arg = removeComments(*it, std::string("#"), std::string("\n"));//remove shell comments.
					arg = removeComments(arg, std::string("//"), std::string("\n"));//remove C simple comments.
					arg = removeComments(arg, std::string("/*"), std::string("*/"));//remove C multiple comments.
					arg = TextTools::removeWhiteSpaces(arg);
					*it = arg;
					int temp=it->length();
					if (temp==0) {
						toRemove.push_back(i);
					}
					else {
						if (temp>maxStrSize){maxStrSize=temp;};
					}
					i++;
				}
				for(std::vector<int>::reverse_iterator it = toRemove.rbegin(); it != toRemove.rend(); it++){
					spNames.erase(spNames.begin()+*it);
				}
        if (spNames.size()<2) {
          std::cout <<"No more than two species: no need to make a species tree! Quitting."<<std::endl;
          exit (0);
        }
				tree = TreeTemplateTools::getRandomTree(spNames);
        tree->setBranchLengths(1.0);
				TreeTools::midpointRooting(*tree);
        std::cout <<"Initial species tree<<std::endl";
        std::cout <<TreeTools::treeToParenthesis(*tree, false)<<std::endl;

			}
			else throw Exception("Unknown init tree method.");

      //We set arbitrary branch lengths
			assignArbitraryBranchLengths(*tree);

      //We write the starting species tree to a file
		 std::string file = ApplicationTools::getStringParameter("starting.tree.file", params, "starting.tree");
			Newick newick;
			newick.write(*tree, file, true);

			// Try to write the current tree to file. This will be overwritten by the optimized tree,
			// but allows to check file existence before running optimization!
			PhylogeneticsApplicationTools::writeTree(*tree, params);

			/****************************************************************************
			 // Then, we get 5 options.
			 * Meaning of heuristicsLevel:
       * 0: exact double-recursive algorithm. All possible root likelihoods are computed with only 2 tree traversals (default).
			 * 1: fastest heuristics : only a few nodes are tried for the roots of the gene trees (the number of these nodes tried depends upon speciesIdLimitForRootPosition), and for each root tried, the events are re-computed only for a subset of the tree.
			 * 2: All roots are tried, and for each root tried, the events are re-computed only for a subset of the tree.
			 * 3: All roots are tried, and for each root tried, the events are re-computed for all nodes of the tree (which should be useless unless there is a bug in the selection of the subset of the nodes.
       * WARNING: options 1, 2, 3 are probably bugged now.
			 * branchProbaOptimization sets the type of optimization that is applied to duplication and loss rates:
			 * "average": all branches have the same average rates. These rates are optimized during the course of the program.
			 * "branchwise": each branch has its own set of duplication and loss rates. These rates are optimized during the course of the program.
			 * "average_then_branchwise": at the beginning, all branches have the same average rates, and then they are individualised.
			 *****************************************************************************/
			heuristicsLevel= ApplicationTools::getIntParameter("heuristics.level", params, 0, "", false, false);
			speciesIdLimitForRootPosition= ApplicationTools::getIntParameter("species.id.limit.for.root.position", params, 2, "", false, false);
			optimizationTolerance=ApplicationTools::getDoubleParameter("optimization.tolerance", params, 0.000001, "", false, false);
			branchProbaOptimization = ApplicationTools::getStringParameter("branchProbabilities.optimization",params,"average");
      std::cout << "Branch expected numbers of duplications and losses: "<<branchProbaOptimization <<std::endl;
			if ((branchProbaOptimization!="average")&&(branchProbaOptimization!="branchwise")&&(branchProbaOptimization!="average_then_branchwise")&&(branchProbaOptimization!="no")) {
				std::cout << "branchProbabilities.optimization is not properly set; please set to either 'average', 'branchwise', 'average_then_branchwise', or 'no'. "<<std::endl;
				exit(-1);
			}
			sprLimit=ApplicationTools::getIntParameter("spr.limit",params,4);


			/****************************************************************************
			 // We get percent coverage of the genomes under study.
			 // We produce a genomeMissing std::map that associates species names to percent of missing data.
			 *****************************************************************************/
		 std::string percentCoverageFile = ApplicationTools::getStringParameter("genome.coverage.file",params,"none");
			std::vector <int> toRemove;
		 std::map <std::string, int> genomeMissing;
			for(std::vector<std::string >::iterator it = spNames.begin(); it != spNames.end(); it++){
				genomeMissing[*it]=0;
			}

			if (percentCoverageFile=="none" ){
				std::cout << "No file for genome.coverage.file, we assume that all genomes are covered at 100%."<<std::endl;
			}
			else {
			 std::ifstream inCoverage (percentCoverageFile.c_str());
				std::vector <std::string> listCoverages;
				while(getline(inCoverage,line)) {
					listCoverages.push_back(line);
				}
				int i =0;
        for(std::vector<std::string >::iterator it = listCoverages.begin(); it != listCoverages.end(); it++){
				 std::string arg = removeComments(*it, std::string("#"), std::string("\n"));//remove shell comments.
					arg = removeComments(arg, std::string("//"), std::string("\n"));//remove C simple comments.
					                                 arg = removeComments(arg, std::string("/*"), std::string("*/"));//remove C multiple comments.
					                                 arg = TextTools::removeWhiteSpaces(arg);
					                                 *it = arg;
					                                 int temp=it->length();
					                                 if (temp==0) {
														 toRemove.push_back(i);
													 }
					                                 i++;
				                                 }
				for(std::vector<int>::reverse_iterator it = toRemove.rbegin(); it != toRemove.rend(); it++){
					listCoverages.erase(listCoverages.begin()+*it);
				}
				for(std::vector<std::string >::iterator it = listCoverages.begin(); it != listCoverages.end(); it++){
					StringTokenizer st1 = StringTokenizer::StringTokenizer (*it, ":", true);
				 std::map<std::string,int>::iterator iter = genomeMissing.find(st1.getToken(0));
					if (iter != genomeMissing.end() ) {
						iter->second = 100-TextTools::toInt(st1.getToken(1));
					}
					else {
						std::cout <<"Coverage from species "<<st1.getToken(0)<<" is useless as long as species "<<st1.getToken(0)<<" is not included in the dataset..."<<std::endl;
					}
				}
			}
			/*for(std::map<std::string, int >::iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++){
				std::cout <<it->first<<" : "<<it->second<<std::endl;
			}*/

			/****************************************************************************
			 // Then, we handle all gene families.
			 *****************************************************************************/

			//We get the file containing the list of all gene family files
		 std::string listGeneFile = ApplicationTools::getStringParameter("genelist.file",params,"none");
			if (listGeneFile=="none" ){
				std::cout << "\n\nNo list of genes was provided. Cannot compute a reconciliation between a species tree and gene trees if you do not tell me where the gene trees are ! Use the option genelist.file, which should be a file containing a list of gene option files.\n" << std::endl;
				std::cout << "ReconcileDuplications species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
				exit(-1);
			}
			//Getting the list of gene alignments, trees, correspondance files and options.
			//In this file, the format is expected to be as follows :
			/*
			 optionFile1 for gene1
			 optionFile2 for gene2
			 optionFile3 for gene3
			 ...
			 */
      //Alternatively, there can be sizes for the gene families, in which case the format is:
      /*
			 optionFile1:size1
			 optionFile2:size2
			 optionFile3:size3
			 ...
			 */      
      
			//Addresses to these option files are expected to be absolute: may need improvement, for instance through the use of global variables, as in other files.
		 std::ifstream inListOpt (listGeneFile.c_str());
			while(getline(inListOpt,line)) {
				listOptions.push_back(line);
			}
      
      
			int maxStrSize=0;
			toRemove.clear();
			int i=0;
			for(std::vector<std::string >::iterator it = listOptions.begin(); it != listOptions.end(); it++){
			 std::string arg = removeComments(*it, std::string("#"), std::string("\n"));//remove shell comments.
				arg = removeComments(arg, std::string("//"), std::string("\n"));//remove C simple comments.
        arg = removeComments(arg, std::string("/*"), std::string("*/"));//remove C multiple comments.
        arg = TextTools::removeWhiteSpaces(arg);
        *it = arg;
        int temp=it->length();
        if (temp==0) {
          toRemove.push_back(i);
        }
        else {
          if (temp>maxStrSize){maxStrSize=temp;};
        }
        i++;
      }
			for(std::vector<int>::reverse_iterator it = toRemove.rbegin(); it != toRemove.rend(); it++){
				listOptions.erase(listOptions.begin()+*it);
			}
			if(maxStrSize>MAXFILENAMESIZE) {
				std::cout << "\nFile names are too long, please abbreviate ! File names (including the path) need to be < "<< MAXFILENAMESIZE <<" letters long.\n"<<std::endl;
				exit(-1);
			}
      std::cout <<"Using "<<listOptions.size()<<" gene families."<<std::endl;
      std::cout <<"Using " <<size<<" nodes"<<std::endl;

			if (listOptions.size()==0) {
				std::cout << "\n\nThere should be at least one option file specified in the list of option files !"<<std::endl;
				exit(-1);
			}
      else if (listOptions.size() < size -1) {
        std::cout << "You want to use more nodes than (gene families+1). This is not possible. Please decrease the number of processors to use."<<std::endl;
        exit(-1);
      }
			//We compute and send the number of genes per client.
      if (size==1) {
        std::cout <<"\n\n\n\t\tError: this program can only run if 2 or more processes are used."<<std::endl;
        std::cout <<"\t\tUse 'mpirun -np k ReconcileDuplications ...', where k>=2"<<std::endl;
        exit(-1);
      }
      //Here, two alternatives: either we do have information regarding the gene family sizes, or we don't.
      //Is there size information?==Is there a ":" in the first line?
      if (TextTools::hasSubstring(listOptions[0],":")) {
          generateListOfOptionsPerClient(listOptions, size, listOfOptionsPerClient, numbersOfGenesPerClient);
      }
      else {
        int numberOfGenesPerClient = (int)(listOptions.size()) / (size -1); 
        int reste = (listOptions.size()) % (size -1);
        std::cout <<"Number of genes per client : "<<numberOfGenesPerClient<< " and extra "<<reste<<std::endl;
        for (int i = 0 ; i< size-1 ; i++ ) {
          if (i<reste) {
            numbersOfGenesPerClient.push_back(numberOfGenesPerClient+1);
          }
          else {
            numbersOfGenesPerClient.push_back(numberOfGenesPerClient);
          }
        } 
        std::vector<int>::iterator it2 = numbersOfGenesPerClient.begin();
        numbersOfGenesPerClient.insert( it2, int(0) ); //For the server, we insert a "dumb" option file at the beginning of the std::vector, so only clients compute the reconciliation
        int currentFile = 0;
        std::vector<std::string> temp2;
        for (int i = 0 ; i< size ; i++ ) {
          listOfOptionsPerClient.push_back(temp2);
          if (i == 0) {
            listOfOptionsPerClient[i].push_back(std::string("####"));
          }
          else {
            for (int j = 0 ; j < numbersOfGenesPerClient[i] ; j++) {
              listOfOptionsPerClient[i].push_back(listOptions[currentFile]);
              currentFile++;
            }
          }
        }
    }
      
      /*std::cout <<"HEHEH"<<std::endl;
      VectorTools::print(listOfOptionsPerClient[1]);
      std::cout <<"HEHEH2"<<std::endl;
      VectorTools::print(numbersOfGenesPerClient);
      std::cout <<"HEHEH3"<<std::endl;
       */
			/***********************************************************
			 * First we set the numbers of duplications and losses at 0.
			 ***********************************************************/

			std::vector <int> nodes = tree->getNodesId();
			SpeciesNodeNumber = nodes.size();
			for (int i=0; i<SpeciesNodeNumber; i++) {
				//We use values for rates of duplications and losses obtained from a dataset containing 150 genes.
				lossProbabilities.push_back(0.313);
        duplicationProbabilities.push_back(0.0159);
        /*lossNumbers.push_back(0);
				duplicationNumbers.push_back(0);
				branchNumbers.push_back(0);*/
				num0Lineages.push_back(0);
				num1Lineages.push_back(0);
				num2Lineages.push_back(0);
			}
			optimizeSpeciesTreeTopology = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
			resetLossesAndDuplications(*tree, /*lossNumbers, */lossProbabilities, /*duplicationNumbers, */duplicationProbabilities);
			breadthFirstreNumber (*tree);
			//We set preliminary loss and duplication rates, correcting for genome coverage
      computeDuplicationAndLossRatesForTheSpeciesTreeInitially(branchProbaOptimization, 
                                                               num0Lineages, 
                                                               num1Lineages, 
                                                               num2Lineages, 
                                                               lossProbabilities, 
                                                               duplicationProbabilities, 
                                                               genomeMissing, 
                                                               *tree);

    /*  for (int i =0; i<num0Lineages.size() ; i++ ) {
        std::cout <<"branch Number#"<< i<<"  du Rate: "<< duplicationProbabilities[i]<<"  loss Rate: "<< lossProbabilities[i]<<std::endl;
      }
      duplicationProbabilities = duplicationProbabilities * 0.1;
      lossProbabilities = lossProbabilities * 0.1;
      std::cout << "After correction:"<<std::endl;
      for (int i =0; i<num0Lineages.size() ; i++ ) {
        std::cout <<"branch Number#"<< i<<"  du Rate: "<< duplicationProbabilities[i]<<"  loss Rate: "<< lossProbabilities[i]<<std::endl;
      }
      */
      
      //We also need to send the species tree to all clients
			currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
			
      //The server sends the number of genes affected to each client
			scatter(world , numbersOfGenesPerClient, affectedNumberOfGenes, server);

			//Then the server sends each client the std::vector of filenames it is in charge of
			scatter(world , listOfOptionsPerClient, affectedFilenames, server);

      //Now the server sends to all clients the same information.
			broadcast(world, optimizeSpeciesTreeTopology, server);
			broadcast(world, SpeciesNodeNumber, server);
      //broadcast(world, lossNumbers, server);
			broadcast(world, lossProbabilities, server); 
			//broadcast(world, duplicationNumbers, server);
			broadcast(world, duplicationProbabilities, server); 
			//broadcast(world, branchNumbers, server);
			broadcast(world, num0Lineages, server);
			broadcast(world, num1Lineages, server);
			broadcast(world, num2Lineages, server);

			//And we also send the species tree to all nodes
			broadcast(world, currentSpeciesTree, server);
			//And we send the two options read by the server
			broadcast(world, heuristicsLevel, server);
			broadcast(world, speciesIdLimitForRootPosition, server);


			/***********************************The server gathers and sums client likelihoods *************************/

			//It might be cool to be able to use the reduce functions, but I have not been able to do so...
			/*  reduce(world, duplicationNumbers, std::plus<std::vector<int> >(), server);
			reduce(world, lossNumbers, std::plus<std::vector<int> >(), server);
			reduce(world, logL, logL, std::plus<double>(),server);*/

      bool noMoreSPR;
			if(optimizeSpeciesTreeTopology) {
				std::cout <<"Optimizing the species tree topology"<<std::endl;
        noMoreSPR=false; 
			}
			else {
				noMoreSPR=true;
			}
      
			breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
			double bestlogL = -UNLIKELY;
			int numIterationsWithoutImprovement = 0;
			bestIndex = 0;
			int index = 0;
			int nodeForNNI = 0;
			int nodeForSPR = 1;
			int nodeForRooting = 4;
			int limit = 3 * tree->getNumberOfNodes();
			int endlimit = limit+1;
			TreeTemplate<Node> * bestTree = tree->clone();
			TreeTemplate<Node> * currentTree = tree->clone();
		/*	std::vector <int> bestLossNumbers;
			std::vector <int> bestDuplicationNumbers;
			std::vector <int> bestBranchNumbers;*/
			std::vector <int> bestNum0Lineages;
			std::vector <int> bestNum1Lineages;
			std::vector <int> bestNum2Lineages;
      
      //vector to keep NNIs likelihoods, for making aLRTs.
      std::vector <double > NNILks(2*tree->getNumberOfLeaves()-2, NumConstants::VERY_BIG); 
      
			double averageDuplicationProbability;
			double averageLossProbability;
      logL = 0.0;
		 std::vector<double> logLs;
			gather(world, logL, logLs, server);
			logL = - VectorTools::sum(logLs);
			//resetVector(duplicationNumbers);
			//resetVector(lossNumbers);
			//resetVector(branchNumbers);
			resetVector(num0Lineages);
			resetVector(num1Lineages);
			resetVector(num2Lineages);
			//gather(world, duplicationNumbers, AllDuplications, server); 
			//gather(world, lossNumbers, AllLosses, server);
			//gather(world, branchNumbers, AllBranches, server);
			gather(world, num0Lineages, allNum0Lineages, server);
			gather(world, num1Lineages, allNum1Lineages, server);
			gather(world, num2Lineages, allNum2Lineages, server);
			int temp = allNum0Lineages.size();
			for (int i =0; i<temp ; i++ ) {
				//duplicationNumbers= duplicationNumbers+AllDuplications[i];
        //lossNumbers= lossNumbers+AllLosses[i];
        //branchNumbers= branchNumbers+AllBranches[i];
        num0Lineages = num0Lineages+allNum0Lineages[i];
				num1Lineages = num1Lineages+allNum1Lineages[i];
				num2Lineages = num2Lineages+allNum2Lineages[i];        
			} 
      /*bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;
      bestBranchNumbers = branchNumbers;*/
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      
      //TEST 16 02 2010
      /*
      if(branchProbaOptimization=="average_then_branchwise") {
        std::string temp = "average";
        computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
      }
      else {
        computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
      }
      */
			std::cout << "\t\tServer: total initial Likelihood value "<<logL<<std::endl;
			bestlogL = logL;
      
      VectorTools::print(num0Lineages);
      VectorTools::print(num1Lineages);
      VectorTools::print( num2Lineages);
      
      
      
      for (int i =0; i<num0Lineages.size() ; i++ ) {
        std::cout <<"branch Number#"<< i<<" du Rate: "<< duplicationProbabilities[i]<<" loss Rate: "<< lossProbabilities[i]<<std::endl;
      }
      
			while (!stop) {
				//Using deterministic SPRs first, and then NNIs
				//Making SPRs, from leaves to deeper nodes (approximately), plus root changes
				if (!noMoreSPR) {
          //This first function does not optimize duplication and loss rates, 
          //as proved by the last "false" argument.
          fastTryAllPossibleSPRsAndReRootings(world, currentTree, bestTree, index, bestIndex, stop, logL, bestlogL, /*lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, */num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, branchProbaOptimization, genomeMissing, sprLimit, false);
					std::cout <<"Before updating rates; current Likelihood "<<bestlogL  <<" and logL: "<<logL<<std::endl;
          
          backupLossProbabilities = lossProbabilities;
          backupDuplicationProbabilities = duplicationProbabilities;
          
          //Now we update the rates before starting a new round of tree exploration
          computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL, /*lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches,*/ num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *currentTree, bestlogL);
          std::cout <<"After updating rates; current Likelihood "<<logL<<std::endl;
          
          if (logL+0.01<bestlogL) {
            bestlogL =logL;
            /*bestLossNumbers = lossNumbers;
            bestDuplicationNumbers = duplicationNumbers;
            bestBranchNumbers = branchNumbers;*/
            bestNum0Lineages = num0Lineages;
            bestNum1Lineages = num1Lineages;
            bestNum2Lineages = num2Lineages;
            bestIndex = index;
            std::cout << "Updating duplication and loss rates yields a better candidate tree likelihood : "<<bestlogL<<std::endl;
            std::cout << TreeTools::treeToParenthesis(*currentTree, true)<<std::endl;
          } 
          else {
            std::cout << "No improvement: we keep former rates. "<<bestlogL<<std::endl;
            lossProbabilities = backupLossProbabilities;
            duplicationProbabilities = backupDuplicationProbabilities;
          }
          
          
          
          for (int i =0; i<num0Lineages.size() ; i++ ) {
            std::cout <<"branch Number#"<< i<<" dup Rate: "<< duplicationProbabilities[i]<<" loss Rate: "<< lossProbabilities[i]<<std::endl;
          }
          

          //Now for each species tree tried, we optimize the duplication and loss rates.
          //Therefore, for each species tree, the likelihood is computed twice, 
          //once before, and once after parameter tuning.
          fastTryAllPossibleSPRsAndReRootings(world, currentTree, bestTree, index, bestIndex, stop, logL, bestlogL, /*lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, */num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, branchProbaOptimization, genomeMissing, sprLimit, true);
          
          noMoreSPR=true;
          
         // numericalOptimizationOfDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *currentTree, bestlogL);
          
          
          std::cout << "\n\n\t\t\tNow entering the final optimization steps, without SPRs\n\n"<<std::endl;
          std::cout << TreeTools::treeToParenthesis(*currentTree, true)<<std::endl;
				}

				else { //noMoreSPR==true, we thus only make NNIs and root changes
					if (((!rearrange)&&(numIterationsWithoutImprovement>=2*tree->getNumberOfNodes()))||((!rearrange)&&(!optimizeSpeciesTreeTopology))) {
            
            rearrange = true; //Now we rearrange gene trees

            numIterationsWithoutImprovement = 0;
            //We compute a new likelihood with rearrangement and possibly branch specific probabilities of duplications and losses.
            //Communications from the server to clients
            currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
           /*broadcast(world, stop, server);
            broadcast(world, rearrange, server); 
            broadcast(world, lossProbabilities, server); 
            broadcast(world, duplicationProbabilities, server); 
            broadcast(world, currentSpeciesTree, server);*/
            broadcastsAllInformation(world, server, stop, rearrange, lossProbabilities, duplicationProbabilities, currentSpeciesTree);

            //Computation in clients
            index++;  
            std::cout <<"\tNumber of species trees tried : "<<index<<std::endl;
            /*resetVector(duplicationNumbers);
            resetVector(lossNumbers);
            resetVector(branchNumbers);*/  
            resetVector(num0Lineages);
            resetVector(num1Lineages);
            resetVector(num2Lineages);
            logL = 0.0;
            resetVector(logLs);
            gather(world, logL, logLs, server); 
            logL = VectorTools::sum(logLs);
            //Communications from the clients to the server
            /*gather(world, duplicationNumbers, AllDuplications, server); 
            gather(world, lossNumbers, AllLosses, server);
            gather(world, branchNumbers, AllBranches, server);*/
            gather(world, num0Lineages, allNum0Lineages, server);
            gather(world, num1Lineages, allNum1Lineages, server);
            gather(world, num2Lineages, allNum2Lineages, server);
            int temp = allNum0Lineages.size();
            for (int i =0; i<temp ; i++ ) {
              /*duplicationNumbers= duplicationNumbers+AllDuplications[i];
              lossNumbers= lossNumbers+AllLosses[i];
              branchNumbers= branchNumbers+AllBranches[i];*/
              num0Lineages = num0Lineages+allNum0Lineages[i];
              num1Lineages = num1Lineages+allNum1Lineages[i];
              num2Lineages = num2Lineages+allNum2Lineages[i];              
            }
                       
            bestlogL =logL;
            bestIndex = index;
            /*bestLossNumbers = lossNumbers;
            bestDuplicationNumbers = duplicationNumbers;
            bestBranchNumbers = branchNumbers;*/
            bestNum0Lineages = num0Lineages;
            bestNum1Lineages = num1Lineages;
            bestNum2Lineages = num2Lineages;
            
            //The first computation with gene tree rearrangement is done, 
            //now we can optimize duplication and loss rates.
            
            computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    
            for (int i =0; i<num0Lineages.size() ; i++ ) {
              std::cout <<"branch Number#"<< i<<" du Rate: "<< duplicationProbabilities[i]<<" loss Rate: "<< lossProbabilities[i]<<std::endl;
            }
            currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
            /*broadcast(world, stop, server);
            broadcast(world, rearrange, server); 
            broadcast(world, lossProbabilities, server); 
            broadcast(world, duplicationProbabilities, server); 
            broadcast(world, currentSpeciesTree, server);*/
            broadcastsAllInformation(world, server, stop, rearrange, lossProbabilities, duplicationProbabilities, currentSpeciesTree);

            //Computation in clients
            index++;  
            std::cout <<"\tNumber of species trees tried : "<<index<<std::endl;
            /*resetVector(duplicationNumbers);
            resetVector(lossNumbers);
            resetVector(branchNumbers);*/  
            resetVector(num0Lineages);
            resetVector(num1Lineages);
            resetVector(num2Lineages);
            logL = 0.0;
            resetVector(logLs);
            gather(world, logL, logLs, server); 
            logL = VectorTools::sum(logLs);
            /*gather(world, duplicationNumbers, AllDuplications, server); 
            gather(world, lossNumbers, AllLosses, server);
            gather(world, branchNumbers, AllBranches, server);*/
            gather(world, num0Lineages, allNum0Lineages, server);
            gather(world, num1Lineages, allNum1Lineages, server);
            gather(world, num2Lineages, allNum2Lineages, server);
            temp = allNum0Lineages.size();
            for (int i =0; i<temp ; i++ ) {
              /*duplicationNumbers= duplicationNumbers+AllDuplications[i];
              lossNumbers= lossNumbers+AllLosses[i];
              branchNumbers= branchNumbers+AllBranches[i];*/
              num0Lineages = num0Lineages+allNum0Lineages[i];
              num1Lineages = num1Lineages+allNum1Lineages[i];
              num2Lineages = num2Lineages+allNum2Lineages[i];              
            }
              
            std::cout << "\t\tServer : Likelihood value with gene tree optimization and new branch probabilities: "<<logL<<" compared to the former log-likelihood : "<<bestlogL<<std::endl;
            numIterationsWithoutImprovement = 0;
            
            bestlogL =logL;
            bestIndex = index;
            /*bestLossNumbers = lossNumbers;
            bestDuplicationNumbers = duplicationNumbers;
            bestBranchNumbers = branchNumbers;*/
            bestNum0Lineages = num0Lineages;
            bestNum1Lineages = num1Lineages;
            bestNum2Lineages = num2Lineages;
          }
          
          //Final steps in the optimization of the species tree topology
					if (optimizeSpeciesTreeTopology) {
						std::cout <<"\tNNIs or Root changes: Number of iterations without improvement : "<<numIterationsWithoutImprovement<<std::endl;
						localOptimizationWithNNIsAndReRootings(world, currentTree, bestTree, index, bestIndex, stop, logL, bestlogL, /*lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, */num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, nodeForNNI, nodeForRooting, branchProbaOptimization, genomeMissing, NNILks);
            //ATTEMPT 18 02 2010
            //OPTIMIZE RATES
          /*  SpeciesTreeLikelihood *slk = new SpeciesTreeLikelihood (world, server, tree, index, bestIndex, stop, logL, bestlogL, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, numIterationsWithoutImprovement, branchProbaOptimization, genomeMissing);
            
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
            
            lossProbabilities = lossProbabilities * optimizer->getParameters().getParameter("coefLoss").getValue();
            duplicationProbabilities = duplicationProbabilities * optimizer->getParameters().getParameter("coefDup").getValue();
            
            std::cout <<"Dup coef: "<<optimizer->getParameters().getParameter("coefDup").getValue() <<" Loss coef: "<<optimizer->getParameters().getParameter("coefLoss").getValue()<<std::endl;
            
            
            delete optimizer;
            delete slk; */
            
            
					}
					else {
						std::cout <<"\tDuplication and loss rate optimization: Number of iterations without improvement : "<<numIterationsWithoutImprovement<<std::endl;
            /*
            //ATTEMPT 18 02 2010
            //OPTIMIZE RATES
            SpeciesTreeLikelihood *slk = new SpeciesTreeLikelihood (world, server, tree, index, bestIndex, stop, logL, bestlogL, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, numIterationsWithoutImprovement, branchProbaOptimization, genomeMissing);
            
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
            
            lossProbabilities = lossProbabilities * optimizer->getParameters().getParameter("coefLoss").getValue();
            duplicationProbabilities = duplicationProbabilities * optimizer->getParameters().getParameter("coefDup").getValue();
            
            std::cout <<"Dup coef: "<<optimizer->getParameters().getParameter("coefDup").getValue() <<" Loss coef: "<<optimizer->getParameters().getParameter("coefLoss").getValue()<<std::endl;
                        
            
            delete optimizer;
            delete slk;
            */
          optimizeOnlyDuplicationAndLossRates(world, currentTree, bestTree, index, bestIndex, stop, logL, bestlogL, /*lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, */num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, nodeForNNI, nodeForRooting, branchProbaOptimization, genomeMissing);
           /* std::cout<< "STOP VALUE: ";
            if (stop) {
              std::cout <<"STOP ==TRUE!"<<std::endl; 
            }
            else {
              std::cout <<"STOP ==FALSE!"<<std::endl;
            }
            numericalOptimizationOfDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *currentTree, bestlogL);*/
					}
				}
			}

      //In Anisimova and Gascuel, the relevant distribution is a mixture of chi^2_1 and chi^2_0.
      //Here, I am not sure one can do the same. We stick with the classical chi^2_1 distribution, more conservative.
      for (int i = 0; i<num0Lineages.size() ; i++ ) {
        if ((! bestTree->getNode(i)->isLeaf()) && (bestTree->getNode(i)->hasFather())) {
         // double proba=(1/2)*(RandomTools::pChisq(2*(NNILks[i] - bestlogL), 1)+1); //If one wants to use the mixture.
          double proba=RandomTools::pChisq(2*(NNILks[i] - bestlogL), 1);
          //Bonferroni correction
          proba = 1-3*(1-proba);
          if (proba<0) {
            std::cout <<"Negative aLRT!"<<std::endl;
            proba ==0;
          }
          std::cout <<"Branch "<<i<<" second best Lk: "<< NNILks[i]<< "; aLRT: "<<proba<<std::endl;
          bestTree->getNode(i)->setBranchProperty("ALRT", Number<double>(proba));
        }
        else {
          bestTree->getNode(i)->setBranchProperty("ALRT", Number<double>(1));
        }
      }
      
      std::cout <<"\n\n\t\tBest Species Tree found, with node Ids: "<<std::endl;
			std::cout << TreeTools::treeToParenthesis (*bestTree, true)<<std::endl;
      std::cout <<"\n\n\t\tBest Species Tree found, with aLRTs: "<<std::endl;
      std::cout << treeToParenthesisWithDoubleNodeValues(*bestTree, false, "ALRT")<<std::endl;
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
      
      //Here we output the species tree with rates of duplication and loss
      //For duplication rates
      for (int i =0; i<num0Lineages.size() ; i++ ) {
        bestTree->getNode(i)->setBranchProperty("DUPLICATIONS", Number<double>( duplicationProbabilities[i]));
        if (bestTree->getNode(i)->hasFather()) {
          bestTree->getNode(i)->setDistanceToFather(duplicationProbabilities[i]);
        }
      }
      
      //Outputting the results of the algorithm, server side.
      std::string dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", params, "AllDuplications.tree", "", false, false); 
      std::ofstream out (dupTree.c_str(), std::ios::out);
			out << treeToParenthesisWithDoubleNodeValues(*bestTree, false, "DUPLICATIONS")<<std::endl;
			out.close();
      //For loss rates
      for (int i =0; i<num0Lineages.size() ; i++ ) {
        bestTree->getNode(i)->setBranchProperty("LOSSES", Number<double>(lossProbabilities[i]));
        if (bestTree->getNode(i)->hasFather()) {
          bestTree->getNode(i)->setDistanceToFather(lossProbabilities[i]);
        }
      }
      std::string lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", params, "AllLosses.tree", "", false, false);
			out.open (lossTree.c_str(), std::ios::out);
			out << treeToParenthesisWithDoubleNodeValues(*bestTree, false, "LOSSES")<<std::endl;
			out.close();
      
      std::string numTree = ApplicationTools::getStringParameter("output.numbered.tree.file", params, "ServerNumbered.tree", "", false, false);
			out.open (numTree.c_str(), std::ios::out);
			out << TreeTools::treeToParenthesis (*bestTree, true)<<std::endl;
			out.close();
      
      //Here we output the species tree with numbers of times 
      //a given number of lineages has been found per branch.
      
      assignNumLineagesOnSpeciesTree(*bestTree, 
                                     num0Lineages, 
                                     num1Lineages, 
                                     num2Lineages);
      std::string lineagesTree = ApplicationTools::getStringParameter("output.lineages.tree.file", params, "lineageNumbers.tree", "", false, false); 
      out.open (lineagesTree.c_str(), std::ios::out);
			out << TreeTools::treeToParenthesis(*bestTree, false, NUMLINEAGES)<<std::endl;
			out.close();
      
          
      
      
			std::cout <<"Number of species trees tried : "<<index<<std::endl;

			PhylogeneticsApplicationTools::writeTree(*bestTree, params);

			std::cout << "\t\tServer : best found logLikelihood value : "<<bestlogL<<std::endl;

			deleteTreeProperties(*currentTree);
			deleteTreeProperties(*bestTree);
			delete currentTree;
			delete bestTree;
			std::cout << "ReconcileDuplication's done. Bye." << std::endl;
			ApplicationTools::displayTime("Total execution time:");
		}//End if at the server node
    
		//##################################################################################################################
    //##################################################################################################################
    //############################################# IF AT A CLIENT NODE ################################################
    //##################################################################################################################
    //##################################################################################################################
	
    if (rank >server){
			//The server sends the number of genes affected to each client
			scatter(world , numbersOfGenesPerClient, affectedNumberOfGenes, server);

			//Then the server sends each client its std::vector of filenames it is in charge of
			scatter(world , listOfOptionsPerClient, affectedFilenames, server);

			broadcast(world, optimizeSpeciesTreeTopology, server);
			broadcast(world, SpeciesNodeNumber, server);
			//Here we need to send the std::vectors of losses and duplications and branch numbers to all clients.
			//broadcast(world, lossNumbers, server);
			broadcast(world, lossProbabilities, server); 
			//broadcast(world, duplicationNumbers, server);
			broadcast(world, duplicationProbabilities, server); 
			//broadcast(world, branchNumbers, server);
			broadcast(world, num0Lineages, server);
			broadcast(world, num1Lineages, server);
			broadcast(world, num2Lineages, server);
			//And we also send the species tree to all nodes
			broadcast(world, currentSpeciesTree, server);
			//And we send the two options read by the server
			broadcast(world, heuristicsLevel, server);
			broadcast(world, speciesIdLimitForRootPosition, server);

			//First we read the species tree from the char[] sent by the server
      tree=TreeTemplateTools::parenthesisToTree(currentSpeciesTree, false, "", true);
			resetLossesAndDuplications(*tree, /*lossNumbers, */lossProbabilities, /*duplicationNumbers, */duplicationProbabilities);
      //To make the correspondance between species name and id:
      std::map <std::string, int> spId = computeSpeciesNamesToIdsMap(*tree);
      //These std::vectors contain all relevant numbers for all gene families the client is in charge of.
			/*std::vector <std::vector <int> > allDuplicationNumbers;
			std::vector <std::vector <int> > allLossNumbers;
			std::vector <std::vector <int> > allBranchNumbers;*/
			std::vector <std::vector <int> > allNum0Lineages;
			std::vector <std::vector <int> > allNum1Lineages;
			std::vector <std::vector <int> > allNum2Lineages;
			std::vector <double> allLogLs;
			std::vector <std::map<std::string, std::string> > allParams;
			TreeTemplate<Node> * geneTree = NULL;
			int MLindex = 0;

			std::vector <ReconciliationTreeLikelihood *> treeLikelihoods;
			std::vector <ReconciliationTreeLikelihood *> backupTreeLikelihoods;

			std::cout <<"Client  of rank "<<rank <<" is in charge of " << affectedFilenames.size()<<" gene families."<<std::endl;
			std::vector <Alphabet *> allAlphabets;
			std::vector <VectorSiteContainer *>   allDatasets;
			std::vector <SubstitutionModel *> allModels;
			std::vector <DiscreteDistribution *> allDistributions;
      std::vector <TreeTemplate<Node> *> allGeneTrees;
			std::vector <TreeTemplate<Node> *> allUnrootedGeneTrees;

      int numDeletedFamilies=0;
      bool avoidFamily;
      
     /* 
     // This bit of code is useful to use GDB on clients:
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
       
      
      
      
      
      
      
      //Here we are going to get all necessary information regarding all gene families the client is in charge of.
			for (int i = 0 ; i< affectedFilenames.size() ; i++) { //For each file
        avoidFamily = false;
			 std::string file =affectedFilenames[i];
				// std::cout << "affectedfilename # " <<i <<" : "<<file<<std::endl;

				if(!FileTools::fileExists(file))
				{
				 std::cerr << "Parameter file not found." << std::endl;
					exit(-1);
				}
				else
				{
					params = AttributesTools::getAttributesMapFromFile(file, "=");
					AttributesTools::resolveVariables(params);
				}

        //Sequences and model of evolution
				Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);
        VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);       
        VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, params);     
				delete allSites;             
				//method to optimize the gene tree root; only useful if heuristics.level!=0.
        bool rootOptimization;
				if (ApplicationTools::getStringParameter("root.optimization",params,"normal")=="intensive") {
					rootOptimization=true;
				}
				else {
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
     
        
        /****************************************************************************
				 //Then we need to get the substitution model.
				 *****************************************************************************/
        
        
        SubstitutionModel*    model    = 0;
        SubstitutionModelSet* modelSet = 0;
        DiscreteDistribution* rDist    = 0;
        
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
				TreeTemplate<Node> * unrootedGeneTree = NULL;
				initTree = ApplicationTools::getStringParameter("init.gene.tree", params, "user", "", false, false);
				ApplicationTools::displayResult("Input gene tree", initTree);
				if(initTree == "user")
				{
				 std::string geneTreeFile =ApplicationTools::getStringParameter("gene.tree.file",params,"none");
					if (geneTreeFile=="none" ){
						std::cout << "\n\nNo Gene tree was provided. The option init.gene.tree is set to user (by default), which means that the option gene.tree.file must be filled with the path of a valid tree file. \nIf you do not have a gene tree file, the program can start from a random tree, if you set init.gene.tree at random\n\n" << std::endl;
						exit(-1);
					}
					Newick newick(true);
					geneTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTreeFile));
					if (!geneTree->isRooted()) {
						unrootedGeneTree = geneTree->clone();
						std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
						geneTree->newOutGroup(0);
					}
					else {
						unrootedGeneTree = geneTree->clone();
						unrootedGeneTree->unroot();
					}
					ApplicationTools::displayResult("Gene Tree file", geneTreeFile);
					ApplicationTools::displayResult("Number of leaves", TextTools::toString(geneTree->getNumberOfLeaves()));
				}
				else if(initTree == "random")
				{
					std::cout << "\n\nCurrently this option is not available.\n"<<std::endl;
					exit(-1);
				}
        else if (initTree == "bionj") //building a BioNJ starting tree
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

            geneTree = unrootedGeneTree;
            geneTree->newOutGroup(0);

          }
				else throw Exception("Unknown init gene tree method.");
     
        /****************************************************************************
				 //Then we need to prune the gene tree and the alignment so that they contain
         //only sequences from the species under study.
				 *****************************************************************************/
        //If we need to remove all sequences or all sequences except one, 
        //better remove the gene family
        if (seqsToRemove.size()>=sites->getNumberOfSequences()-1) {
          numDeletedFamilies = numDeletedFamilies+1;
          avoidFamily=true;
          std::cout <<"Avoiding family "<<affectedFilenames[i-numDeletedFamilies+1]<<std::endl;
        }
        
        if (!avoidFamily) { //This family is phylogenetically informative
          for (int j =0 ; j<seqsToRemove.size(); j++) {
            removeLeaf(*geneTree, seqsToRemove[j]);
            unrootedGeneTree = geneTree->clone();
            if (!geneTree->isRooted()) {
              std::cout <<"gene tree is not rooted!!! "<< taxaseqFile<<std::endl;
            }
            unrootedGeneTree->unroot();
            sites->deleteSequence(seqsToRemove[j]);
          }
          //printing the gene trees with the species names instead of the sequence names
          //This is useful to build an input for duptree for instance
          TreeTemplate<Node> * treeWithSpNames = unrootedGeneTree->clone();
          std::vector <Node *> leaves = treeWithSpNames->getLeaves();
          for (int j =0; j<leaves.size() ; j++) {
            leaves[j]->setName(seqSp[leaves[j]->getName()]);
          }
          std::vector <Node *> nodes =  treeWithSpNames->getNodes();
          for (int j =0; j<nodes.size() ; j++) {
            if (nodes[j]->hasFather()) {
              nodes[j]->deleteDistanceToFather(); 
            }
            if (nodes[j]->hasBootstrapValue()) {
              nodes[j]->removeBranchProperty(TreeTools::BOOTSTRAP); 
            }
          }
          Newick newick(true);
          std::string geneTreeFile =ApplicationTools::getStringParameter("gene.tree.file",params,"none");
          newick.write(*treeWithSpNames, geneTreeFile+"SPNames", true);
          delete treeWithSpNames;
          
          /****************************************************************************
           //Then we initialize the losses and duplication numbers on this tree.
           *****************************************************************************/
          std::vector<int> numbers = num0Lineages;
          /*std::vector<int> numbers = lossNumbers;
          allDuplicationNumbers.push_back(numbers);
          allLossNumbers.push_back(numbers);
          allBranchNumbers.push_back(numbers);*/
          allNum0Lineages.push_back(numbers);
          allNum1Lineages.push_back(numbers);
          allNum2Lineages.push_back(numbers);
          resetLossesAndDuplications(*tree, /*allLossNumbers[i-numDeletedFamilies], */lossProbabilities, /*allDuplicationNumbers[i-numDeletedFamilies], */duplicationProbabilities);
          //resetVector(allBranchNumbers[i-numDeletedFamilies]);
          resetVector(allNum0Lineages[i-numDeletedFamilies]);
          resetVector(allNum1Lineages[i-numDeletedFamilies]);
          resetVector(allNum2Lineages[i-numDeletedFamilies]);
          
          /********************************************COMPUTING LIKELIHOOD********************************************/
     
          bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", params, true, "", false, false);
          if(!computeLikelihood)
            {
              delete alphabet;
              delete sites;
              delete tree;
              std::cout << "ReconcileDuplication's done. Bye." << std::endl;
              return(0);
            }
          
          // Setting branch lengths?
          std::string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", params, "Input", "", true, false);
          std::string cmdName;
          std::map<std::string, std::string> cmdArgs;
          KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
          if (cmdName == "Input")
            {
              // Do nothing!
            }
          else if (cmdName == "Equal")
            {
              double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, false);
              if (value <= 0)
                throw Exception("Value for branch length must be superior to 0");
              ApplicationTools::displayResult("Branch lengths set to", value);
              tree->setBranchLengths(value);
            }
          else if (cmdName == "Clock")
            {
              TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
            }
          else if (cmdName == "Grafen")
            {
              std::string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "Input", "", true, false);
              double h;
              if (grafenHeight == "input")
                {
                  h = TreeTools::getHeight(*tree, tree->getRootId());
                }
              else
                {
                  h = TextTools::toDouble(grafenHeight);
                  if (h <= 0) throw Exception("Height must be positive in Grafen's method.");
                }
              ApplicationTools::displayResult("Total height", TextTools::toString(h));
              
              double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, false);
              ApplicationTools::displayResult("Grafen's rho", rho);
              TreeTools::computeBranchLengthsGrafen(*tree, rho);
              double nh = TreeTools::getHeight(*tree, tree->getRootId());
              tree->scaleTree(h / nh);
            }
          else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
          ApplicationTools::displayResult("Branch lengths", cmdName);

          //     ReconciliationTreeLikelihood *tl;
          //   NNIHomogeneousTreeLikelihood *tl2;
          DiscreteRatesAcrossSitesTreeLikelihood* tl;
        //  DiscreteRatesAcrossSitesTreeLikelihood* tl2;
          
          
          std::string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", params, "no", "", true, false);
          ApplicationTools::displayResult("Clock", optimizeClock);
          bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
          
           
          if(optimizeClock == "global")
            {
              std::cout<<"Sorry, clocklike trees have not been implemented yet."<<std::endl;
              exit(0);
            }// This has not been implemented!
          else if(optimizeClock == "no")
            {
              
             // tl2 = new NNIHomogeneousTreeLikelihood(*geneTree, *sites, model, rDist, true, true);

              tl = new ReconciliationTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, *tree, *geneTree, seqSp, spId, /*allLossNumbers[i-numDeletedFamilies], */lossProbabilities, /*allDuplicationNumbers[i-numDeletedFamilies], */duplicationProbabilities, /*allBranchNumbers[i-numDeletedFamilies], */allNum0Lineages[i-numDeletedFamilies], allNum1Lineages[i-numDeletedFamilies], allNum2Lineages[i-numDeletedFamilies], speciesIdLimitForRootPosition, heuristicsLevel, MLindex, true, true, rootOptimization);

            }
          else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
         // tl2->initialize();
        //  std::cout<<"Value tl2: "<<tl2->getValue()<<std::endl;

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
        }
			}//End for each file

			
			std::vector <std::vector <std::string> > reconciledTrees;
			std::vector <std::vector <std::string> > duplicationTrees;
			std::vector <std::vector <std::string> > lossTrees;
			std::vector <std::string> t;  
			std::vector <std::map<std::string, std::string> > allParamsBackup = allParams;
			for (int i = 0 ; i< affectedFilenames.size()-numDeletedFamilies ; i++) {

				reconciledTrees.push_back(t);
				duplicationTrees.push_back(t);
				lossTrees.push_back(t);
        //This is to avoid optimizing gene tree parameters in the first steps of the program
				if (ApplicationTools::getBooleanParameter("optimization.topology", allParams[i], false, "", true, false)){
					allParams[i][ std::string("optimization.topology")] = "false";
				}
        allParams[i][ std::string("optimization")] = "false"; //Quite extreme, but the sequence likelihood has no impact on the reconciliation !
        treeLikelihoods[i]->OptimizeSequenceLikelihood(false);
			}

			bool recordGeneTrees = false; //At the beginning, we do not record the gene trees.
			int startRecordingTreesFrom = 0; //This int is incremented until the gene trees start to be backed-up, when we start the second phase of the algorithm.
      bool firstTimeImprovingGeneTrees = false; //When for the first time we optimize gene trees, we set it at true

      //We make a backup of the gene tree likelihoods.
      for (int i =0 ; i<treeLikelihoods.size() ; i++) {
        backupTreeLikelihoods.push_back(treeLikelihoods[i]->clone());
      }      

      while (!stop) {      //MAIN LOOP STARTS HERE

				logL=0.0;
				/*resetVector(duplicationNumbers);
				resetVector(lossNumbers);
				resetVector(branchNumbers);*/
				resetVector(num0Lineages);
				resetVector(num1Lineages);
				resetVector(num2Lineages);
				for (int i = 0 ; i< affectedFilenames.size()-numDeletedFamilies ; i++) {
         // std::cout<< "Here 0 "<<std::endl;

          if (firstTimeImprovingGeneTrees) {
            treeLikelihoods[i]->OptimizeSequenceLikelihood(true);
            backupTreeLikelihoods[i]->OptimizeSequenceLikelihood(true);
            PhylogeneticsApplicationTools::optimizeParameters(treeLikelihoods[i], treeLikelihoods[i]->getParameters(), allParams[i], "", true, false); 
          }
          else {
            PhylogeneticsApplicationTools::optimizeParameters(treeLikelihoods[i], treeLikelihoods[i]->getParameters(), allParams[i], "", true, false); 
          }
         // std::cout<< "Here 1 "<<std::endl;

					/********************************************LIKELIHOOD OPTIMIZED********************************************/
					geneTree = new TreeTemplate<Node>(treeLikelihoods[i]->getTree());
					resetLossesAndDuplications(*tree, /*allLossNumbers[i], */lossProbabilities, /*allDuplicationNumbers[i], */duplicationProbabilities);
					/*allDuplicationNumbers[i] = treeLikelihoods[i]->getDuplicationNumbers();
					allLossNumbers[i] = treeLikelihoods[i]->getLossNumbers();	
					allBranchNumbers[i] = treeLikelihoods[i]->getBranchNumbers();*/
					allNum0Lineages[i] = treeLikelihoods[i]->get0LineagesNumbers();
					allNum1Lineages[i] = treeLikelihoods[i]->get1LineagesNumbers();
					allNum2Lineages[i] = treeLikelihoods[i]->get2LineagesNumbers();
					MLindex = treeLikelihoods[i]->getRootNodeindex();
          allLogLs[i] = treeLikelihoods[i]->getValue();  
					/*duplicationNumbers = duplicationNumbers + allDuplicationNumbers[i];
					lossNumbers = lossNumbers + allLossNumbers[i];
					branchNumbers = branchNumbers + allBranchNumbers[i];*/
					logL = logL + allLogLs[i];
					num0Lineages = num0Lineages + allNum0Lineages[i];
					num1Lineages = num1Lineages + allNum1Lineages[i];
					num2Lineages = num2Lineages + allNum2Lineages[i];
          //setLossesAndDuplications(*tree, allLossNumbers[i], allDuplicationNumbers[i]);
				
					if (recordGeneTrees) {
						reconciledTrees[i].push_back(TreeTools::treeToParenthesis (*geneTree, false, EVENT));
						duplicationTrees[i].push_back(treeToParenthesisWithIntNodeValues (*tree, false, DUPLICATIONS));
						lossTrees[i].push_back(treeToParenthesisWithIntNodeValues (*tree, false, LOSSES));
					}

 					if (geneTree) {
						delete geneTree;
					}
				}//end for each filename
           
        if (firstTimeImprovingGeneTrees) {
          firstTimeImprovingGeneTrees = false;
        }
				if (!recordGeneTrees) {
					startRecordingTreesFrom++;
				}
       // std::cout<< "Here 2"<<std::endl;
        //Clients send back stuff to the server.
				gather(world, logL, server);
				/*gather(world, duplicationNumbers, server);
				gather(world, lossNumbers, server);
				gather(world, branchNumbers, server);*/  
				gather(world, num0Lineages, allNum0Lineages, server);
				gather(world, num1Lineages, allNum1Lineages, server);
				gather(world, num2Lineages, allNum2Lineages, server);	
       // std::cout<< "Here 3"<<std::endl;

        //Should the computations stop? The server tells us.
				broadcast(world, stop, server);
       // std::cout<< "Here 4"<<std::endl;

				if (!stop) {	// we continue the loop
         // std::cout<< "Here 5"<<std::endl;

          //Reset the gene trees by resetting treeLikelihoods
          //We always start from ML trees according to sequences only
          for (int i=0 ; i<treeLikelihoods.size() ; i++) {
            ReconciliationTreeLikelihood * tempL= treeLikelihoods[i];
            treeLikelihoods[i] = backupTreeLikelihoods[i]->clone();
            delete tempL;
          }          
					if (!rearrange) {
						broadcast(world, rearrange, server);
					} 
					else {
						allParams = allParamsBackup;
						if (recordGeneTrees==false) {
              firstTimeImprovingGeneTrees = true;
							recordGeneTrees=true;
							bestIndex=startRecordingTreesFrom;
             // std::cout <<"We start recording gene trees : "<<startRecordingTreesFrom<<" and bestIndex :" <<bestIndex<<" and index "<< index <<std::endl;
						}
						broadcast(world, rearrange, server);
					}	
					broadcast(world, lossProbabilities, server); 
					broadcast(world, duplicationProbabilities, server); 
					broadcast(world, currentSpeciesTree, server);
          tree=TreeTemplateTools::parenthesisToTree(currentSpeciesTree, false, "", true);
          spId = computeSpeciesNamesToIdsMap(*tree);
         // std::cout<< "Here 7"<<std::endl;

          for (int i = 0 ; i< affectedFilenames.size()-numDeletedFamilies ; i++) {
            treeLikelihoods[i]->setSpTree(*tree);
            treeLikelihoods[i]->setSpId(spId);
            treeLikelihoods[i]->setProbabilities(duplicationProbabilities, lossProbabilities);
            treeLikelihoods[i]->computeTreeLikelihood();
          }
          //std::cout<< "Here 8"<<std::endl;

				}
				else { //The end, outputting the results
         // std::cout<< "Here 6"<<std::endl;

					broadcast(world, bestIndex, server);
					for (int i = 0 ; i< affectedFilenames.size() ; i++) {
					 std::string reconcTree = ApplicationTools::getStringParameter("output.reconciled.tree.file", allParams[i], "reconciled.tree", "", false, false);
					 std::ofstream out (reconcTree.c_str(), std::ios::out);
						out << reconciledTrees[i][bestIndex-startRecordingTreesFrom]<<std::endl;
						out.close();
					 std::string dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", allParams[i], "duplications.tree", "", false, false);
						out.open (dupTree.c_str(), std::ios::out);
						out << duplicationTrees[i][bestIndex-startRecordingTreesFrom]<<std::endl;
						out.close();
					 std::string lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", allParams[i], "losses.tree", "", false, false);
						out.open (lossTree.c_str(), std::ios::out);
						out << lossTrees[i][bestIndex-startRecordingTreesFrom]<<std::endl;
						out.close();
					}
					break;
				}
			}//End while, END OF MAIN LOOP

			for (int i = 0 ; i< affectedFilenames.size()-numDeletedFamilies ; i++) {  
				delete allAlphabets[i];
				delete allDatasets[i];
				delete allModels[i];
				delete allDistributions[i];
			        //delete allGeneTrees[i];
				//delete allUnrootedGeneTrees[i];
				//delete treeLikelihoods[i];
			}

    }//end if a client node
		if (tree) {
			//TEST (2)
			//deleteTreeProperties(*tree);
			//delete tree;
		}
	}
	catch(std::exception & e)
	{
		std::cout << e.what() << std::endl;
		exit(-1);
	}
	return (0);
}











