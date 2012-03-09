//
//  COALTools.h
//  phyldog
//
//  Created by Bastien Boussau on 07/03/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#ifndef COALTools_h
#define COALTools_h
#include "ReconciliationTools.h"


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
                                 std::vector <std::vector<int> > & speciesIDs);


/*****************************************************************************
 * Utilitary functions to initialize vectors of counts at leaves.
 ****************************************************************************/
void initializeCountVectors(std::vector< std::vector<int> > vec);
void initializeCountVector(std::vector<int>  vec);



#endif
