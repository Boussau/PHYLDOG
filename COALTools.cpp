//
//  COALTools.cpp
//  phyldog
//
//  Created by Bastien Boussau on 07/03/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "COALTools.h"





/*****************************************************************************
 * This function performs a postorder tree traversal in order to find 
 * fill up vectors of counts of coalescence events for rootings. 
 * When followed by the preorder tree traversal function, 
 * vectors of counts for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalPostorder(TreeTemplate<Node> & spTree, 
                                 TreeTemplate<Node> & geneTree, 
                                 Node * node, 
                                 const std::map<std::string, std::string > & seqSp, 
                                 const std::map<std::string, int > & spID, 
                                 std::vector < std::vector< std::vector<int> > > & coalCounts,
                                 std::vector <std::vector<int> > & speciesIDs) {
	int id=node->getId();
 	if (node->isLeaf()) {
        //Fill all leaf count vectors with 0s
        initializeCountVectors(coalCounts[id]);
        return();
    }
    else {
        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            computeSubtreeCoalPostorder(spTree, geneTree, 
                                              sons[i], seqSp, 
                                              spID, coalCounts, 
                                              speciesIDs);
        }
        
        int idSon0 = sons[0]->getId();
        int idSon1 = sons[1]->getId();
        unsigned int directionSon0, directionSon1;
        std::vector <Node *> neighbors = sons[0]->getNeighbors();
        for (unsigned int i=0; i<neighbors.size(); i++) {
            if (neighbors[i]==node) {
                directionSon0 = i;
            }
        }
        neighbors = sons[1]->getNeighbors();
        for (unsigned int i=0; i<neighbors.size(); i++) {
            if (neighbors[i]==node) {
                directionSon1 = i;
            }
        }
       
        computeCoalCountsFromSons (spTree, sons, 
                                   speciesIDs[id][0], 
                                   coalCounts[idSon0][directionSon0], 
                                   coalCounts[idSon1][directionSon1]
                                   speciesIDs[idSon0][directionSon0], 
                                   speciesIDs[idSon1][directionSon1], );//To implement...
        
        computeConditionalLikelihoodAndAssignSpId(spTree, sons, 
                                                  likelihoodData[id][0], 
                                                  likelihoodData[idSon0][directionSon0], 
                                                  likelihoodData[idSon1][directionSon1], 
                                                  lossRates, duplicationRates, 
                                                  speciesIDs[id][0], 
                                                  speciesIDs[idSon0][directionSon0], 
                                                  speciesIDs[idSon1][directionSon1], 
                                                  dupData[id][0], 
                                                  dupData[idSon0][directionSon0], 
                                                  dupData[idSon1][directionSon1], 
                                                  TreeTemplateTools::isRoot(*node));
                
        return();
	}
	
	
}


/*****************************************************************************
 * Utilitary functions to initialize vectors of counts at leaves.
 ****************************************************************************/
void initializeCountVectors(std::vector< std::vector<int> > vec) {
    int size = vec[0].size();
    for (unsigned int i = 0 ; i < 3 ; i++ ) {
        initializeCountVector(vec[i]);
    }
    return;
}


void initializeCountVector(std::vector<int>  vec) {
    for (unsigned int j = 0 ; j < size ; j++ ) {
        vec[j] = 0;
    }
    return;
}

