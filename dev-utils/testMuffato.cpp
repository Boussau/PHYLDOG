//
//  testMuffato.cpp
//  phyldog
//
//  Created by Bastien Boussau on 05/10/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "ReconciliationTools.h"
#include "GeneTreeAlgorithms.h"
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>

//Compilation: //mpic++ -g -lbpp-core -lbpp-seq -lbpp-phyl -lboost_serialization -lboost_mpi  ReconciliationTools.cpp testMuffato.cpp  GeneTreeAlgorithms.cpp GenericTreeExplorationAlgorithms.cpp -o testMuffato
//On pbil-deb: mpic++ -I /panhome/boussau/libs/include/ -L/usr/local/lib -L/panhome/boussau/libs/lib/ -g -lbpp-core -lbpp-seq -lbpp-phyl -lboost_serialization -lboost_mpi ReconciliationTools.cpp testMuffato.cpp GeneTreeAlgorithms.cpp GenericTreeExplorationAlgorithms.cpp -o testMuffato



int main(int args, char ** argv)
{
	if(args == 1)
	{
		std::cout << "not enough arguments"<<std::endl;	
		exit(0);
	}
	ApplicationTools::startTimer();
	//All processors parse the main options
	std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);

	Newick *treeReader = new Newick(true);

	//Gene tree:
	std::string geneTreesFile = ApplicationTools::getStringParameter("input.gene.tree", params, "", "", false, false);
	ApplicationTools::displayResult("Gene Tree file", geneTreesFile);
	Newick newick(true);
	TreeTemplate <Node>* gTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(geneTreesFile));

	//Species tree:
	std::string speciesTreeFile = ApplicationTools::getStringParameter("input.species.tree", params, "", "", false, false);
	ApplicationTools::displayResult("Species Tree file", speciesTreeFile);
	TreeTemplate <Node>* spTree  = dynamic_cast < TreeTemplate < Node > * > (newick.read(speciesTreeFile));

	
	
	std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
	std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*gTree, true) << std::endl;

	delete treeReader;
	
	std::map<std::string, std::string > seqSp;
	vector<string> leaves = spTree->getLeavesNames();
	for (unsigned int i = 0 ; i < leaves.size() ; i++) {
		seqSp.insert( pair<std::string, std::string >(leaves[i],leaves[i]) );
	}
	
	breadthFirstreNumber (*spTree);
	std::map<std::string, int > spId = computeSpeciesNamesToIdsMap(*spTree);

	annotateGeneTreeWithScoredDuplicationEvents (*spTree, 
												 *gTree, 
												 gTree->getRootNode(), 
												 seqSp, spId); 
	//Now, we generate a tree topology in which we rearrange the poorly supported duplications.
	//By that, we mean duplications with scores < 0.3 (default value)
	double editionThreshold = ApplicationTools::getDoubleParameter("muffato.edition.threshold", params, 0.3, "", false, false);
	std::cout << "SPTREE: " << TreeTemplateTools::treeToParenthesis(*spTree, false) <<std::endl;
	
	std::cout << "TREE: " << TreeTemplateTools::treeToParenthesis(*gTree, false) <<std::endl;

	
	editDuplicationNodesMuffato(*spTree, *gTree, gTree->getRootNode(), editionThreshold);
	
	std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
	std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*gTree, true) << std::endl;
	
	return 0;

}
