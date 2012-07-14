//
//  CoalTools.cpp
//  phyldog
//
//  Created by Bastien Boussau on 07/03/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "COALTools.h"



/*****************************************************************************
 * This function performs a postorder tree traversal in order to  
 * fill up vectors of counts of coalescence events for rootings. 
 * For each branch of the species tree, we need to record how many lineages got in,
 * and how many lineages got out.
 * Thus, for each branch of the species tree, we have 2 ints: 
 * vec[0]: number of incoming lineages
 * vec[1]: number of outgoing lineages
 * When followed by the preorder tree traversal function, 
 * vectors of counts for all rootings are computed.
 * coalCounts contains all lower counts for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalCountsPostorder(TreeTemplate<Node> & spTree, 
									   TreeTemplate<Node> & geneTree, 
									   Node * node, 
									   std::map<std::string, std::string > & seqSp, 
									   std::map<std::string, int > & spID, 
									   std::vector < std::vector< std::vector< std::vector< unsigned int > > > > & coalCounts,
									   std::vector <std::vector<unsigned int> > & speciesIDs) {
/*	if (node ->hasFather() == false ){
		std::cout <<"XXXX BITE DEBUT XXXXX" << std::endl;
		std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(geneTree, true) << std::endl;
		std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(spTree, true) << std::endl;
		std::cout <<"SIZES: "<< coalCounts.size() << " "<< coalCounts[0].size() << " "<< coalCounts[0][0].size()<< " "<< coalCounts[0][0][0].size() <<std::endl;
		std::cout << "geneTree->getRootNode()->getId() " << geneTree.getRootNode()->getId() << std::endl;
		map<std::string, std::string>::iterator it;
		for ( it=seqSp.begin() ; it != seqSp.end(); it++ ) {
			std::cout << (*it).first << " => " << (*it).second << std::endl;
		}
		map<std::string, int >::iterator it2;
		for ( it2=spID.begin() ; it2 != spID.end(); it2++ ) {
			std::cout << (*it2).first << " => " << (*it2).second << std::endl;
		}
		//printCoalCounts(coalCounts);

		for (unsigned int i = 0 ; i < speciesIDs.size() ; i++) {
			for (unsigned int j = 0 ; j < speciesIDs[i].size() ; j++) {
				std::cout << "speciesIDs[i][j]: "<< i << " " << j <<" "<< speciesIDs[i][j]<<std::endl;			
			}
		}
		//printCoalCounts(coalCounts);
		std::cout << node->getId() << std::endl;
		std::cout <<"XXXX BITE FIN XXXXX" << std::endl;

	}*/

	int id=node->getId();

 	if (node->isLeaf()) { 
		//In principle this should be done only once at the beginning of the algorithm, if node ids are kept
        //Fill all leaf count vectors with 0s
        initializeCountVectors(coalCounts[id]);

        speciesIDs[id][0] = assignSpeciesIdToLeaf(node, seqSp, spID);
		speciesIDs[id][1] = speciesIDs[id][0];
		speciesIDs[id][2] = speciesIDs[id][0];

//        int i = spID[ seqSp[ node->getName()] ];
        incrementOutCount(coalCounts[id][0], speciesIDs[id][0]);

        return;
    }
    else {

        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            computeSubtreeCoalCountsPostorder(spTree, geneTree, 
											  sons[i], seqSp, 
											  spID, coalCounts, 
											  speciesIDs);

        }
        
        int idSon0 = sons[0]->getId();
        int idSon1 = sons[1]->getId();
		//Index to define the subtree for son0 and son1
        unsigned int directionSon0, directionSon1;
		//Here, because of our convention of linking the father node to direction 0
		directionSon0 = directionSon1 = 0;
		/*
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
        }*/

		/*
		//The speciesIDs[idSon0][directionSon0] and 
		//speciesIDs[idSon1][directionSon1] should have been filled
		//in by computeSubtreeCoalCountsPostorder.
		//We can use them to fill in speciesIDs[id][directionFather]
		unsigned int directionFather;
		Node * father  = 0;
		if (node ->hasFather() ){
			father = node->getFather();
		}
		neighbors = node->getNeighbors();
		for (unsigned int i=0; i<neighbors.size(); i++) {
            if (neighbors[i]==father) {
                directionFather = i;
            }
        }
		*/

        computeCoalCountsFromSons (spTree, sons, 
                                   speciesIDs[id][0], 
                                   speciesIDs[idSon0][directionSon0], 
                                   speciesIDs[idSon1][directionSon1],
                                   coalCounts[id][0],
                                   coalCounts[idSon0][directionSon0], 
                                   coalCounts[idSon1][directionSon1]
                                    ); 

        return;
	}	
}


/*****************************************************************************
 * Utilitary functions to initialize vectors of counts at leaves.
 ****************************************************************************/
void initializeCountVectors(std::vector< std::vector< std::vector<unsigned int> > > &vec) {
    for (unsigned int i = 0 ; i < 3 ; i++ ) { //It should not be necessary to initialize all vectors, i==0 should be enough
        initializeCountVector(vec[i]);
    }
    return;
}


void initializeCountVector(std::vector<std::vector<unsigned int> >  &vec) {
    for (unsigned int j = 0 ; j < vec.size() ; j++ ) {
         for (unsigned int i = 0 ; i < 2 ; i++ ) {
             vec[j][i] = 0;
         }
    }
    return;
}

/*****************************************************************************
 * Utilitary functions to increment vectors of counts.
 ****************************************************************************/

void incrementOutCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos) {
    vec[pos][1] = vec[pos][1] +1;
    return;
}

void incrementInCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos) {
    vec[pos][0] = vec[pos][0] +1;
    return;
}


/*****************************************************************************
 * Computes a vector of counts from two son vectors, and assigns species ID to the
 * father node.
 ****************************************************************************/

void computeCoalCountsFromSons (TreeTemplate<Node> & tree, std::vector <Node *> sons, 
                                unsigned int & rootSpId, 
                                const unsigned int & son0SpId,
                                const unsigned int & son1SpId,
                                std::vector< std::vector<unsigned int> > & coalCountsFather,
                                std::vector< std::vector<unsigned int> > & coalCountsSon0,
                                std::vector< std::vector<unsigned int> > & coalCountsSon1)
{
//std::cout << "SHOULD be not 0: "<<    son0SpId << " and "<< son1SpId <<"; Son 0 id: "<< sons[0]->getId() << " Son 1 id: "<< sons[1]->getId() <<std::endl;
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
    for (unsigned int i= 0 ; i < coalCountsFather.size() ; i++) { //for each branch of the sp tree
        for (unsigned int j= 0 ; j < 2 ; j++) { //for incoming and outgoing counts
            coalCountsFather[i][j] = coalCountsSon0[i][j] + coalCountsSon1[i][j];
        }
        //std::cout << "coalCountsFather[i][j] For Sp branch "<<i<<" Num coal in: "<< coalCountsFather[i][0] << " Num coal out: "<< coalCountsFather[i][1]<<std::endl;
    }
    
    
    
    Node temp0 = *(tree.getNode(son0SpId));
    Node temp1 = *(tree.getNode(son1SpId));
    
    while (a!=b) { //There is ILS!
        if (a>b) {
            recoverILS(temp0, a, olda, coalCountsFather);
        }
        else {
            recoverILS(temp1, b, oldb, coalCountsFather);
        }
    }
    rootSpId = a;
   // std::cout << "SHOULD be not 0: "<<    son0SpId << " and "<< son1SpId <<"; Son 0 id: "<< sons[0]->getId() << " Son 1 id: "<< sons[1]->getId() << " Father Spid: "<< a<< std::endl;
    return;    
}

/*****************************************************************************
 * Computes a vector of counts from two son vectors, and assigns species ID to the
 * father node.
 ****************************************************************************/

void computeCoalCountsFromSonsAndFillTables (TreeTemplate<Node> & tree, std::vector <Node *> sons, 
                                             unsigned int & rootSpId, 
                                             const unsigned int & son0SpId,
                                             const unsigned int & son1SpId,
                                             std::vector< std::vector<unsigned int> > & coalCountsFather,
                                             std::vector< std::vector<unsigned int> > & coalCountsSon0,
                                             std::vector< std::vector<unsigned int> > & coalCountsSon1, 
                                             std::set<int> & nodesToTryInNNISearch)
{
    
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
    for (unsigned int i= 0 ; i < coalCountsFather.size() ; i++) { //for each branch of the sp tree
        for (unsigned int j= 0 ; j < 2 ; j++) { //for incoming and outgoing counts
            coalCountsFather[i][j] = coalCountsSon0[i][j] + coalCountsSon1[i][j];
        }
	/*	if (i == 46) {
			std::cout << "computeCoalCountsFromSonsAndFillTables 46 Adding: "<< coalCountsSon0[i][1] <<" and "<< coalCountsSon1[i][1] <<std::endl;
		}*/
        //std::cout << "coalCountsFather[i][j] For Sp branch "<<i<<" Num coal in: "<< coalCountsFather[i][0] << " Num coal out: "<< coalCountsFather[i][1]<<std::endl;
    }
    
    
    
    Node temp0 = *(tree.getNode(son0SpId));
    Node temp1 = *(tree.getNode(son1SpId));
    
    while (a!=b) { //There is ILS!
		//DO WE FILL CORRECTLY THE nodesToTryInNNISearch SET ? NEED TO TEST THAT
        if (a>b) {
            nodesToTryInNNISearch.insert(sons[0]->getId());
            recoverILS(temp0, a, olda, coalCountsFather);
        }
        else {
            nodesToTryInNNISearch.insert(sons[1]->getId());
            recoverILS(temp1, b, oldb, coalCountsFather);
        }
    }
    rootSpId = a;
/*std::cout << "rootSpId: "<< rootSpId <<std::endl;
	for (unsigned int i= 0 ; i < coalCountsFather.size() ; i++) { //for each branch of the sp tree
		std::cout << coalCountsFather[i][1] << " ";
	}
	std::cout << std::endl;*/
    return;    
}


/*****************************************************************************
 * This function recovers ILS by comparing a subtree in a gene tree to
 * a species tree.
 * WARNING: MAY NEED TO CHANGE: NEED TO TAKE THE LOWEST POSSIBLE BRANCH FOR DOING THE COALESCENCE
 ****************************************************************************/
void recoverILS(Node & node, int & a, int & olda, 
                std::vector <std::vector < unsigned int > > &vec) {
    olda=a;
    Node* nodeA;
    if (node.hasFather()) {
        nodeA = node.getFather();
    }
    else {
        std::cout <<"Problem in recoverILS , nodeA has no father"<<std::endl; 
    }
    a = nodeA->getId();
    node = *nodeA;
	//if (a == 46) 		std::cout << "recoverILS incrementInCount: "<< olda <<" incrementOutCount : "<< a <<std::endl;	
    incrementInCount(vec, olda);
    incrementOutCount(vec, a);
    
    return;
}

/*****************************************************************************
 * Computes the likelihood using our coalescence model, 
 * given a vector of vector giving, for each branch of the species tree,
 * the number of incoming lineages, and the number of outgoing lineages.
 * Formula from Degnan and Salter (2005), Evolution 59(1), pp. 24-37.
 * 3 versions of the function:
 * - working for lots of gene families at once
 * - working for one gene family
 * - working for one branch of one gene family
 ****************************************************************************/


double computeCoalLikelihood (std::vector < std::vector<std::vector<unsigned int> > > vec, std::vector < double > CoalBl ) 
{
    double logLk = 0;
    for (unsigned int i = 0 ; i < vec.size() ; i++ ) 
    {
        logLk += computeCoalLikelihood (vec[i],  CoalBl );
    }
    return logLk;
}


double computeCoalLikelihood (std::vector < std::vector<unsigned int> > vec, std::vector < double > CoalBl ) 
{
    double logLk = 0;
    for (unsigned int i = 0 ; i < vec.size() ; i++ ) 
    {
        logLk += computeCoalLikelihood (vec[i],  CoalBl[i] );
    }
    return logLk;
}

double computeCoalLikelihood ( std::vector<unsigned int>  vec, double CoalBl ) 
{
    
    double logLk = 0;
    double prod;
    int kminus1;

    double v = (double) vec[0];
    double u = (double) vec[1];

//	std::cout << " v and u:"<< v <<" "<< u <<std::endl;
	
    /*
    for (unsigned int k = vec[0] ; k <= vec[1] ; k++ ) 
    {
        prod = 1.0;
        kminus1 = k-1;
        double dy;
        for (unsigned int y = 0 ; y <= kminus1 ; y++ ) {
            dy = double(y);
            prod *= (v + dy) * (u - dy) / ( u + dy ) ;
        }
        if (prod == 0) {
            prod = 1;
        }
        double dk = (double)k;
        double dkminus1 = (double)kminus1;
        double dkminusv = dk - v ;
       std::cout << "CoalBl "<<CoalBl<<std::endl;
        std::cout << "prod "<< prod <<std::endl;
        std::cout <<"fact "<<NumTools::fact(vec[0])<<std::endl;
        std::cout <<"fact2 "<<NumTools::fact(k - vec[0])<<std::endl;
        std::cout <<"exp (-k * (k-1) * CoalBl/2 ) "<<exp (-(double)k * ((double)k-1.0) * CoalBl/2.0  )<<std::endl;
        std::cout <<"(2*k -1) "<<(2.0*(double)k -1)<<std::endl;
        std::cout <<"( (-1)^(k-vec[0]) ) "<<( pow(-1.0, ( (double)k-(double)vec[0]) ) )<<std::endl;
        std::cout <<"( NumTools::fact(vec[0]) * NumTools::fact(k - vec[0]) * (vec[0] + kminus1 ) ) "<<( (double)NumTools::fact(vec[0]) * (double)NumTools::fact(k - vec[0]) * ((double)vec[0] + kminus1 ) )<<std::endl;
        logLk += ( exp (-dk * dkminus1 * CoalBl/2.0 ) ) * (2.0*dk -1.0) * ( pow(-1.0, dkminusv ) ) / ( NumTools::fact(v) * NumTools::fact(dkminusv) * (v + dkminus1 ) ) * prod;
        //std::cout << "logLk: " << logLk<<std::endl;
    }
	if (logLk <= 0) {
		logLk = NumConstants::VERY_TINY * NumConstants::VERY_TINY;
	}
	 return log(logLk);
	 */
	
	double first = CoalBl/2.0  * (-v*v + v);
	
	for (unsigned int k = (int)v ; k <= (int)u ; k++ ) 
	{
		prod = 1.0;
		kminus1 = k-1;
		double dy;
		for (unsigned int y = 0 ; y <= kminus1 ; y++ ) {
			dy = double(y);
			prod *= (v + dy) * (u - dy) / ( u + dy ) ;
		}
		if (prod == 0) {
			prod = 1;
		}		
		double dk = (double)k;
		double dkminus1 = (double)kminus1;
		double dkminusv = dk - v ;
		logLk += ( exp (dkminusv * (1- 2*v - dkminusv) * CoalBl/2.0 ) ) * (2.0*dk -1.0) * ( pow(-1.0, dkminusv ) ) / ( NumTools::fact(v) * NumTools::fact(dkminusv) * (v + dkminus1 ) ) * prod;
	}
	
	if (logLk < 0) {
		std::cout <<"WARNING: correction"<<std::endl;
		logLk = NumConstants::VERY_TINY;
	}
	
	logLk = log(logLk) + first;

    return logLk;
}







/*****************************************************************************
 * This function performs a preorder tree traversal in order to fill vectors of counts. 
 * When used after the postorder tree traversal function, counts for all rootings are computed.
 * coalCounts contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * node: node of the gene tree used as son of the root.
 ****************************************************************************/

void computeSubtreeCoalCountsPreorder(TreeTemplate<Node> & spTree, 
                                      TreeTemplate<Node> & geneTree, 
                                      Node * node, 
                                      const std::map<std::string, std::string > & seqSp, 
                                      const std::map<std::string, int > & spID, 
                                      std::vector < std::vector< std::vector< std::vector<unsigned  int > > > > & coalCounts, 
                                      std::vector<double> & bls, 
                                      std::vector <std::vector<unsigned int> > & speciesIDs, 
                                      int sonNumber, 
                                      std::map <double, Node*> & LksToNodes) {
//	std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(geneTree, true) << std::endl;
    computeRootingCoalCounts(spTree, node, 
                             coalCounts, bls, speciesIDs, 
                             sonNumber, LksToNodes);
   if (node->isLeaf()) {
        return; 
    }
    Node * son;
    if (sonNumber==1) {
        son= node->getSon(1);
	//	std::cout <<"node son getid: "<< node->getSon(1)->getId() <<std::endl;
    }
    else {
        son= node->getSon(0);
	//	std::cout <<"node son getid: "<< node->getSon(0)->getId() <<std::endl;
    }
    //  for (int i = 0; i< sons.size(); i++){
    for (unsigned int j =0; j<son->getNumberOfSons(); j++) {
        computeSubtreeCoalCountsPreorder(spTree, geneTree, 
                                         son, seqSp, spID, 
                                         coalCounts, 
                                         bls, 
                                         speciesIDs, j, LksToNodes);
    }
    //  }
	return;
	
}




/*****************************************************************************
 * This function computes the Coalescent counts of a rooting. 
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 * node is a node of the gene tree. 
 ****************************************************************************/
void computeRootingCoalCounts(TreeTemplate<Node> & spTree, 
                              Node * node, 
                              std::vector < std::vector< std::vector< std::vector<unsigned  int > > > > & coalCounts, 
                              const std::vector< double> & bls, 
                              std::vector <std::vector<unsigned int> > & speciesIDs, 
                              int sonNumber, 
                              std::map <double, Node*> & LksToNodes) {
    int geneNodeId = node->getId();
	//VectorTools::print(TreeTemplateTools::getLeavesNames(*node));
    int directionForFather;
    std::vector <Node*> nodes;

    nodes.push_back(node->getFather());

    if (sonNumber==0) { //If sonNumber==0, the subtree we're interested in is composed of son 1 and father of node.
        nodes.push_back(node->getSon(1));
    }
    else { //If sonNumber==1, the subtree we're interested in is composed of son 0 and father of node.
        nodes.push_back(node->getSon(0));
    }
    if (node->getFather()->getSon(0)==node) {
        directionForFather = 1; //node #1 is son 0: from node, the useful subtree is in 1
    }
    else if (node->getFather()->getSon(1)==node) {
        directionForFather = 2; //node #2 is son 1: from node, the useful subtree is in 2
    }

    int idNode0, idNode1;
    idNode0 = nodes[0]->getId(); //the father
    idNode1 = nodes[1]->getId(); //the son we picked
    unsigned int directionNode0, directionNode1;
    directionNode0 = directionForFather;
    directionNode1 = 0; //we use the inferior index


    unsigned int directionSon0, directionSon1;
    directionSon0 = sonNumber+1;
    directionSon1 = 0;
	//We will need that for the likelihood computation for this root
	std::vector< std::vector<unsigned  int > > rootCounts = coalCounts[geneNodeId][directionSon0];

    computeCoalCountsFromSons (spTree, nodes, 
                               speciesIDs[geneNodeId][directionSon0], 
                               speciesIDs[idNode0][directionNode0], 
                               speciesIDs[idNode1][directionNode1],
                               coalCounts[geneNodeId][directionSon0],
                               coalCounts[idNode0][directionNode0], 
                               coalCounts[idNode1][directionNode1]
                               );  
    
    //Now we have the counts for the upper subtree, 
    //as well as the counts of the lower subtree (which we already had)
    //We can thus compute the total likelihood of the rooting.
    
    std::vector <Node*> sons;
    sons.push_back(node);
    
    sons.push_back(node->getSon(sonNumber));
    int idSon1 = sons[1]->getId();
    
    double rootLikelihood = 0.0;
    unsigned int rootSpId;
    //int rootDupData = 0;
  
    
/*
	if (node->getSon(sonNumber)->getId() == 4) {
		std::cout << "sonNumber: "<<sonNumber <<std::endl;
    for (unsigned int i = 0 ; i < coalCounts[geneNodeId][directionSon0].size() ; i++) {
        std::cout << "BEFORE PREORDER: Superior subtree For outgroup node  "<< node->getSon(sonNumber)->getId() <<"; Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneNodeId][directionSon0][i][0] << " Num coal out: "<< coalCounts[geneNodeId][directionSon0][i][1]<<std::endl;
    }
		
		for (unsigned int i = 0 ; i < coalCounts[geneNodeId][directionSon0].size() ; i++) {
			std::cout << "BEFORE PREORDER: INferior subtree For outgroup node  "<< node->getSon(sonNumber)->getId() <<"; Sp Branch "<<i<<" Num coal in: "<< coalCounts[idSon1][0][i][0] << " Num coal out: "<< coalCounts[idSon1][0][i][1]<<std::endl;
		}

		
	}*/
	
	


	
    computeCoalCountsFromSons (spTree, sons, 
                               rootSpId, 
                               speciesIDs[geneNodeId][directionSon0], 
                               speciesIDs[idSon1][0],
                               rootCounts,
                               coalCounts[geneNodeId][directionSon0], 
                               coalCounts[idSon1][0]
                               );  

    
    
    //Adding the two tables, from the upper subtree, and from the lower subtree
 /*   for (unsigned int i = 0 ; i < spTree.getNumberOfNodes() ; i++ ) {
        coalCounts[geneNodeId][directionSon0][i][0] = coalCounts[geneNodeId][directionSon0][i][0] + coalCounts[idSon1][0][i][0];
        coalCounts[geneNodeId][directionSon0][i][1] = coalCounts[geneNodeId][directionSon0][i][1] + coalCounts[idSon1][0][i][1];
    }*/
    
    //Add the starting lineage at the root
    for (unsigned int i = 0 ; i < spTree.getNumberOfNodes() ; i++ ) {
        if (rootCounts[i][0]==0 && rootCounts[i][1] !=0)
        {
            rootCounts[i][0]=1; 
            break;
        }
    }
    /*
	if (node->getSon(sonNumber)->getId() == 3) {
    for (unsigned int i = 0 ; i < coalCounts[geneNodeId][directionSon0].size() ; i++) {
        std::cout << "PREORDER: For outgroup node "<< node->getSon(sonNumber)->getId()<<"; Sp Branch "<<i<<" Num coal in: "<< rootCounts[i][0] << " Num coal out: "<< rootCounts[i][1]<<std::endl;
    }
	}*/

    //What to put?
    rootLikelihood = computeCoalLikelihood ( rootCounts, bls ) ;
	
 /*   
    computeConditionalLikelihoodAndAssignSpId(spTree, sons, 
                                              rootLikelihood, 
                                              likelihoodData[idSon0][directionSon0], 
                                              likelihoodData[idSon1][directionSon1], 
                                              lossRates, duplicationRates, 
                                              rootSpId, speciesIDs[idSon0][directionSon0], 
                                              speciesIDs[idSon1][directionSon1], 
                                              rootDupData, dupData[idSon0][directionSon0], 
                                              dupData[idSon1][directionSon1], true);*/
//     std::cout <<"LK FOUND "<<rootLikelihood<<std::endl;
    while (LksToNodes.find(rootLikelihood)!=LksToNodes.end()) {
//         std::cout <<"changing rootLikelihood !!!!!!!!!!!!!!!!!!!"<<std::endl;
        rootLikelihood+=SMALLPROBA;
    }
    LksToNodes[rootLikelihood]=node->getSon(sonNumber);
    return;
}





/*****************************************************************************
 * This function performs a postorder tree traversal in order to  
 * fill up vectors of counts of coalescence events for rootings. 
 * For each branch of the species tree, we need to record how many lineages got in,
 * and how many lineages got out.
 * Thus, for each branch of the species tree, we have 2 ints: 
 * vec[0]: number of incoming lineages
 * vec[1]: number of outgoing lineages
 * When followed by the preorder tree traversal function, 
 * vectors of counts for all rootings are computed.
 * coalCounts contains all lower counts for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalCountsPostorderAndFillTables(TreeTemplate<Node> & spTree, 
                                                    TreeTemplate<Node> & geneTree, 
                                                    Node * node, 
                                                    std::map<std::string, std::string > & seqSp, 
                                                    std::map<std::string, int > & spID, 
                                                    std::vector< std::vector< std::vector< std::vector< unsigned int > > > > & coalCounts,
                                                    std::vector <std::vector<unsigned int> > & speciesIDs, 
                                                    std::set<int> &      nodesToTryInNNISearch      ) {
	int id=node->getId();
 	if (node->isLeaf()) { //In principle this should be done only once at the beginning of the algorithm, if node ids are kept
        //Fill all leaf count vectors with 0s
        initializeCountVectors(coalCounts[id]);
        speciesIDs[id][0] = assignSpeciesIdToLeaf(node, seqSp, spID);
		speciesIDs[id][1] = speciesIDs[id][0] ;
		speciesIDs[id][2] = speciesIDs[id][0] ;
        //        int i = spID[ seqSp[ node->getName()] ];
        incrementOutCount(coalCounts[id][0], speciesIDs[id][0]);
		/*if (speciesIDs[id][0] == 46) {
			std::cout << "computeSubtreeCoalCountsPostorderAndFillTables incrementOutCount46 "<<std::endl;
		}*/
        return;
    }
    else {
        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            computeSubtreeCoalCountsPostorderAndFillTables(spTree, geneTree, 
														   sons[i], seqSp, 
														   spID, coalCounts, 
														   speciesIDs, nodesToTryInNNISearch);
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
		
		/*if (sons[0]->hasName()) {
			std::cout << "HASNMAE: "<< sons[0]->getName() << "ISLEAF?: "<< sons[0]->isLeaf() <<std::endl;
		}
		else {
			std::cout << "NO NMAE: "<< idSon0 <<" Sp: "<< speciesIDs[idSon0][directionSon0] <<std::endl;
			
		}
		if (sons[1]->hasName()) {
			std::cout << "HASNMAE: "<< sons[1]->getName() << "ISLEAF?: "<< sons[1]->isLeaf() <<std::endl;
		}
		else {
			std::cout << "NO NMAE: "<< idSon1<<" Sp: "<< speciesIDs[idSon1][directionSon1] <<std::endl;
			
		}*/

		
        computeCoalCountsFromSonsAndFillTables (spTree, sons, 
                                                speciesIDs[id][0], 
                                                speciesIDs[idSon0][directionSon0], 
                                                speciesIDs[idSon1][directionSon1],
                                                coalCounts[id][0],
                                                coalCounts[idSon0][directionSon0], 
                                                coalCounts[idSon1][directionSon1], 
                                                nodesToTryInNNISearch
                                                );                
        return;
	}	
}


/*****************************************************************************
 * Useful for putting all elements of coalCounts to 0.
 ****************************************************************************/

void resetCoalCounts (std::vector < std::vector < std::vector < std::vector<unsigned int> > > > &coalCounts) {
	for (unsigned int i = 0 ; i < coalCounts.size() ; i++) {
		for (unsigned int j = 0 ; j < 3 ; j++) {
			for (unsigned int k = 0 ; k < coalCounts[i][j].size() ; k++) {
				for (unsigned int l = 0 ; l < 2 ; l++) {
					coalCounts[i][j][k][l] = 0;
				}
			}
		}
	}	
}

void printCoalCounts (std::vector < std::vector < std::vector < std::vector<unsigned int> > > > &coalCounts) {
	for (unsigned int i = 0 ; i < coalCounts.size() ; i++) {
		for (unsigned int j = 0 ; j < 3 ; j++) {
			for (unsigned int k = 0 ; k < coalCounts[i][j].size() ; k++) {
				for (unsigned int l = 0 ; l < 2 ; l++) {
					std::cout << "i: "<<i<<"; j:"<<j<<"; k:"<<k<<"; l:"<< coalCounts[i][j][k][l] << std::endl;
				}
			}
		}
	}	
}




/*****************************************************************************
 * This function aims at finding the most likely coalescent reconciliation, 
 * using a double recursive tree traversal. 
 * The first traversal is post-order, and then the second traversal is pre-order.
 * This is a modification of an algorithm quickly explained in 
 * Chen, Durand, Farach-Colton, J. Comp. Biol. pp429-447, 2000.
 * Conditional likelihoods are recorded in a table. 
 * This table has (number of nodes) elements, and for each node, 
 * contains three conditional likelihoods. 
 * The table is thus (number of nodes)*3 cells. For each node i, 
 * likelihoodData[i][j] contains the conditional likelihood of the subtree 
 * having its root in subtree opposite neighbour j of node i.
 * Node species IDs are also recorded in a (number of nodes)*3 cells table.
 * The boolean "fillTables" is here to tell whether we want to update the vectors num*lineages.
 ****************************************************************************/

double findMLCoalReconciliationDR (TreeTemplate<Node> * spTree, 
                               TreeTemplate<Node> * geneTree, 
                               std::map<std::string, std::string > seqSp,
                               std::map<std::string, int > spID,
                               std::vector< double> coalBl, 
                               int & MLindex, 
                               std::vector < std::vector < std::vector<std::vector<unsigned int> > > > &coalCounts,
                               std::set <int> &nodesToTryInNNISearch, 
                               bool fillTables)
{
	if (!geneTree->isRooted()) {
		std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
		std::cout <<"!!!!!!gene tree is not rooted in findMLCoalReconciliationDR !!!!!!"<<std::endl;
        MPI::COMM_WORLD.Abort(1);
		exit(-1);
    }	
	if (spTree->getRootNode()->getId() != 0) {
		std::cout << TreeTemplateTools::treeToParenthesis (*spTree, true)<<std::endl;
		std::cout <<"!!!!!!Species tree is not properly annotated in findMLCoalReconciliationDR !!!!!!"<<std::endl;
        MPI::COMM_WORLD.Abort(1);
		exit(-1);
    }	

    
    std::vector <unsigned int> nodeSpId(3, 0);

	//This speciesIDs vector contains 3 unsigned ints per node of the gene tree.
	//speciesIDs[i][0]: conditional species index, seen from the father
	//speciesIDs[i][1]: conditional species index, seen from the son 0
	//speciesIDs[i][2]: conditional species index, seen from the son 1
    std::vector <std::vector< unsigned int> > speciesIDs(geneTree->getNumberOfNodes(), nodeSpId);

	//Also, reinitialize the coalCounts
	resetCoalCounts(coalCounts);

	/*
	std::string gStr = "(((((taxon26:0.08316,(taxon30:0.01004,taxon3:0.01004):0.073115):0.090045,taxon2:0.173205):0.172685,taxon15:0.34589):0.15941,(((taxon29:0.23656,(taxon7:0.064225,taxon9:0.064225):0.17233):0.044465,(taxon13:0.09682,taxon33:0.09682):0.1842):0.140795,(((((taxon1:0.183495,((taxon22:0.12459,(taxon21:0.025215,taxon20:0.02522):0.09937):0.01579,taxon24:0.140375):0.04311):0.104935,(taxon16:0.21073,taxon12:0.21073):0.077695):0.01399,((taxon34:0.06071,taxon8:0.06071):0.206685,(taxon36:0.252725,((taxon11:0.0872,(taxon27:0.008555,taxon5:0.008555):0.078645):0.10741,taxon17:0.19461):0.058115):0.014665):0.03502):0.003255,((taxon35:0.18384,(taxon18:0.009025,taxon6:0.009025):0.174815):0.049945,(taxon28:0.174385,taxon23:0.174395):0.0594):0.07188):0.075495,(((taxon19:0.20129,taxon14:0.20129):0.093865,(((taxon25:0.076435,taxon4:0.076435):0.06435,(taxon38:0.118835,taxon32:0.118835):0.02195):0.08733,taxon10:0.22812):0.067035):0.04378,(taxon40:0.30888,(taxon39:0.075125,taxon37:0.075125):0.233765):0.03005):0.042235):0.04065):0.08348):0.00623,taxon31:0.511535);";
	geneTree = TreeTemplateTools::parenthesisToTree(gStr);
*/
	
	
	
	//breadthFirstreNumber (*spTree);
    computeSubtreeCoalCountsPostorder(*spTree, 
                                      *geneTree, 
                                      geneTree->getRootNode(), 
                                      seqSp, 
                                      spID, 
                                      coalCounts,
                                      speciesIDs);

    //Add the starting lineage at the root
    for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
        if (coalCounts[geneTree->getRootNode()->getId()][0][i][0]==0 && coalCounts[geneTree->getRootNode()->getId()][0][i][1] !=0)
        {
            coalCounts[geneTree->getRootNode()->getId()][0][i][0]=1; 
            break;
            
        }
    }
    
/*	for (unsigned int i = 0 ; i < coalCounts[0][0].size() ; i++) {
	//	std::cout << "i: "<<i<<" coalCounts[0][0][i][0]: "<< coalCounts[0][0][i][0] <<" coalCounts[0][0][i][1] " << coalCounts[0][0][i][1] << std::endl;
		std::cout << "i: "<<i<<" coalCounts[geneTree->getRootNode()->getId()][0][i][0]: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] <<" coalCounts[geneTree->getRootNode()->getId()][0][i][1] " << coalCounts[geneTree->getRootNode()->getId()][0][i][1] << std::endl;
	}*/
	
	/*
	 std::cout <<"After First postorder traversal: "<< std::endl;
	for (unsigned int i = 0 ; i < speciesIDs.size() ; i++) {
		for (unsigned int j = 0 ; j < speciesIDs[i].size() ; j++) {
			std::cout << "speciesIDs[i][j]: "<< i << " " << j <<" "<< speciesIDs[i][j]<<std::endl;			
		}
	}
	for (unsigned int i = 0 ; i < coalCounts[geneTree->getRootNode()->getId()][0].size() ; i++) {
		std::cout << "Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] << " Num coal out: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][1]<<std::endl;
	}*/

	
 //   std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
 //   std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*geneTree, true) << std::endl;
    /*    
	std::cout << "Printing result of postorder traversal"<<std::endl;
	std::cout << "In"<<std::endl;
	VectorTools::print (coalCounts[geneTree->getRootNode()->getId()][0][0]);
	std::cout << "Out"<<std::endl;
	VectorTools::print (coalCounts[geneTree->getRootNode()->getId()][0][1]);*/
    double initialLikelihood = computeCoalLikelihood ( coalCounts[geneTree->getRootNode()->getId()][0], coalBl ) ;
    
 //   std::cout << "Initial Likelihood: "<< initialLikelihood <<std::endl;
   
    //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
    std::map <double, Node*> LksToNodes;
    //Now doing the preorder tree traversal
    Node * geneRoot = geneTree->getRootNode();
    std::vector <Node *> sons = geneRoot->getSons();
    if (sons.size()!=2) {
        std::cerr <<"Error: "<< sons.size() << "sons at the root!" <<std::endl; 
    }
    
    LksToNodes[initialLikelihood]=sons[0];
    //We fill the likelihood and species ID data for the root node.
    //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
    //with other nodes.
    coalCounts[geneRoot->getId()][1] = coalCounts[geneRoot->getSon(1)->getId()][0];
    coalCounts[geneRoot->getId()][2] = coalCounts[geneRoot->getSon(0)->getId()][0];
    speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon(1)->getId()][0];
    speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon(0)->getId()][0];
    
    for (unsigned int i = 0; i< sons.size(); i++){
        for (unsigned int j =0; j<sons[i]->getNumberOfSons(); j++) {
            computeSubtreeCoalCountsPreorder(*spTree, 
                                             *geneTree, 
                                             sons[i], 
                                             seqSp, 
                                             spID, 
                                             coalCounts,
                                             coalBl, 
                                             speciesIDs, j, LksToNodes);
        }
    }

    
    
    vector<Node*> nodes = geneTree->getNodes();
    for (unsigned int i = 0 ; i < nodes.size() ; i++ ) {
        if (nodes[i]->hasNodeProperty("outgroupNode") ) {
            nodes[i]->deleteNodeProperty("outgroupNode");
            break;
        }
    }
    
    LksToNodes.rbegin()->second->setNodeProperty("outgroupNode", BppString("here") );
    
  //CHANGE02052012  
	geneTree->newOutGroup(LksToNodes.rbegin()->second);     

/*	TEMP TEST 06072012 geneTree->resetNodesId();
	breadthFirstreNumber (*geneTree);*/

    
    if (fillTables) {
        
        //Now the best root has been found. I can thus run a function with this best root to fill all the needed tables. This additional tree traversal could be avoided.
        //To this end, the needed tables should be filled by the postfix and prefix traversals. This has not been done yet.
        //resetting
        speciesIDs = std::vector<std::vector<unsigned int> > (geneTree->getNumberOfNodes(), nodeSpId);
        
        // Getting a well-rooted tree
        TreeTemplate<Node > * tree = geneTree->clone();
        
        //tree->newOutGroup(LksToNodes.rbegin()->second->getId());
        
        
        //std::cout << TreeTemplateTools::treeToParenthesis (*tree, true)<<std::endl;
        
        nodesToTryInNNISearch.clear();
		//Also, reinitialize the coalCounts
	//	std::cout << "RESETTING: "<< std::endl;
        resetCoalCounts(coalCounts);
	/*	for (unsigned int i = 0 ; i < coalCounts[0][0].size() ; i++) {
			if (spTree->getNode(i)->isLeaf() ) {
				std::cout << "Leaf i: "<<i<<" _coalCounts[0][0][i][0]: "<< coalCounts[0][0][i][0] <<" _coalCounts[0][0][i][1] " << coalCounts[0][0][i][1] << std::endl;
			}
		}*/

        computeSubtreeCoalCountsPostorderAndFillTables (*spTree, *tree, 
                                                        tree->getRootNode(), 
                                                        seqSp, spID, 
                                                        coalCounts, speciesIDs, 
                                                        nodesToTryInNNISearch);
	/*	std::cout << "JUST AFTER computeSubtreeCoalCountsPostorderAndFillTables STILL INSIDE findMLCoalReconciliationDR: "<<std::endl;
		for (unsigned int i = 0 ; i < coalCounts[0][0].size() ; i++) {
			if (spTree->getNode(i)->isLeaf() ) {
				std::cout << "Leaf i: "<<i<<" _coalCounts[0][0][i][0]: "<< coalCounts[0][0][i][0] <<" _coalCounts[0][0][i][1] " << coalCounts[0][0][i][1] << std::endl;
			}
		}*/

/*		std::cout << "Printing result of postorder traversal 2"<<std::endl;
		std::cout << "In"<<std::endl;
		VectorTools::print (coalCounts[geneTree->getRootNode()->getId()][0][0]);
		std::cout << "Out"<<std::endl;
		VectorTools::print (coalCounts[geneTree->getRootNode()->getId()][0][1]);*/

       /* computeNumbersOfLineagesFromRoot(spTree, tree, 
                                         tree->getRootNode(), 
                                         seqSp, spID, 
                                         coalCounts, speciesIDs, 
                                         nodesToTryInNNISearch);*/
        delete tree;
        
    }

    //We return the best likelihood
    MLindex = LksToNodes.rbegin()->second->getId();
	//std::cout <<"Best reconciliation likelihood: "<< LksToNodes.rbegin()->first<<std::endl;
	for(std::map<double, Node* >::iterator it = LksToNodes.begin(); it != LksToNodes.end(); it++){
		std::cout << it->second->getId() <<" : "<<it->first <<std::endl;		
		
	}
	
	
	return LksToNodes.rbegin()->first;
    
    
    
}





void computeCoalBls (std::vector< unsigned int > &  num12Lineages, 
                     std::vector< unsigned int > &  num22Lineages, 
                     std::vector<double> & coalBls) 
{
    double n12;
    double n22;
    for (unsigned int i = 0 ; i < num12Lineages.size() ; i++) {
        n12 = num12Lineages[i];
        n22 = num22Lineages[i];
        if (n22 ==0) {
            n22 = 1;
        }
        if (n12 ==0) {
            n12 = 1;
        }
		//std::cout << "branch " << i << "; n12: "<< n12 <<std::endl;
		//std::cout << "branch " << i << "; n22: "<< n22 <<std::endl;
		
       // estimate = log((n12+n22)/n22);        
        //std::cout <<"\t\tAnalytical estimate: "<< estimate<<std::endl;
        coalBls[i] = log((n12+n22)/n22);
    }
    return;
    
}




void computeCoalBls (std::string& branchExpectedNumberOptimization, 
					 std::vector< unsigned int > &  num12Lineages, 
                     std::vector< unsigned int > &  num22Lineages, 
                     std::vector<double> & coalBls) 
{
    double n12 = 0;
    double n22 = 0;
	
	if (branchExpectedNumberOptimization == "average") {
		double coalBl = 0;
		for (unsigned int i = 0 ; i < num12Lineages.size() ; i++) {
			n12 += num12Lineages[i];
			n22 += num22Lineages[i];
		}
		if (n22 ==0) {
            n22 = 1;
        }
        if (n12 ==0) {
            n12 = 1;
        }
		coalBl = log((n12+n22)/n22);
		for (unsigned int i = 0 ; i < num12Lineages.size() ; i++) {
			coalBls[i] = coalBl;
		}
	}
	else if (branchExpectedNumberOptimization == "no") {
		for (unsigned int i = 0 ; i < num12Lineages.size() ; i++) {
			coalBls[i] = 10.0;
		}
	}
	else {
		for (unsigned int i = 0 ; i < num12Lineages.size() ; i++) {
			n12 = num12Lineages[i];
			n22 = num22Lineages[i];
			if (n22 ==0) {
				n22 = 1;
			}
			if (n12 ==0) {
				n12 = 1;
			}
			
			
			//std::cout << "branch " << i << "; n12: "<< n12 <<std::endl;
			//std::cout << "branch " << i << "; n22: "<< n22 <<std::endl;
			
			// estimate = log((n12+n22)/n22);        
			//std::cout <<"\t\tAnalytical estimate: "<< estimate<<std::endl;
			coalBls[i] = log((n12+n22)/n22);
		}
	}
    return;
    
}






void computeCoalBls (std::vector < std::vector < std::vector< unsigned int > > >&  allGeneCounts , std::vector<double> &coalBls) 
{
    for (unsigned int i = 0 ; i < allGeneCounts.size() ; i++) {
        CoalBranchLikelihood *brLikFunction = new CoalBranchLikelihood(allGeneCounts[i]);
        //Initialize BranchLikelihood:
        coalBls[i] = brLikFunction->estimateBl();
    }
    return;
}







/*******************************************************************************/

//Creating the compressed vector, and the map giving the weights for the "site" patterns.
//The compressed vector is like the compressed alignment in sequence likelihood:
//it gives the weights (=number of occurences) for different patterns.
//Finally it givs an estimate of the branch length in coalescent units, and does so by using the counts
//n12 and n22. One could simply use these counts instead of counting all patterns.
double CoalBranchLikelihood::initModel()
{
    string pattern = "";
    compressedVec_.clear();
    std::vector <unsigned int> v (2, 1);
    //We fill the map that associates patterns to weights.
    //If we want to count everything:
    for (unsigned int i = 0 ; i < vec_.size() ; i++) {
        pattern = TextTools::toString(vec_[i][0]) + TextTools::toString(vec_[i][1]);
        if ( patternToWeights_.find(pattern)!= patternToWeights_.end() ) {
            patternToWeights_[pattern] += 1; 
        }
        else {
            patternToWeights_.insert( pair<std::string, unsigned int >(pattern, 1) );
            v[0] = vec_[i][0];
            v[1] = vec_[i][1];
            compressedVec_.push_back(v);
        }
    }
    lks_.resize(compressedVec_.size(), 0.0);
    for(std::map<std::string, unsigned int >::iterator it = patternToWeights_.begin(); it != patternToWeights_.end(); it++){
        std::cout <<"\t\tit->first: "<< it->first << " : "<< it->second <<std::endl;
    }
    double n12 = (double) patternToWeights_["12"];
    double n22 = (double) patternToWeights_["22"];
    double estimate;
    std::cout <<"n12 "<< n12 <<" n22 "<< n22 <<std::endl;
    if (n22 ==0) {
        n22 = 1;
    }
    if (n12 ==0) {
        n12 = 1;
    }
    estimate = log((n12+n22)/n22);
	return estimate;
}

//Just returning the analytical estimate of the branch length, 
//in coalescent units, using the input vector.
double CoalBranchLikelihood::estimateBl()
{
    //We just count n12 and n22
    double n12;
    double n22;
    for (unsigned int i = 0 ; i < vec_.size() ; i++) {
        if ( vec_[i][0]  == 1 && vec_[i][1] == 2) {
            n12 += 1; 
        }
        else if (vec_[i][0]  == 2 && vec_[i][1] == 2) {
            n22 += 1; 
        }
    }    
    double estimate;
    std::cout <<"n12 "<< n12 <<" n22 "<< n22 <<std::endl;
    if (n22 ==0) {
        n22 = 1;
    }
    if (n12 ==0) {
        n12 = 1;
    }
    estimate = log((n12+n22)/n22);
    
    std::cout <<"\t\tAnalytical estimate: "<< estimate<<std::endl;
    
    return estimate;
}

/*******************************************************************************/

void CoalBranchLikelihood::computeLogLikelihood()
{    
    double l = getParameterValue("BrLen"); 
    lnL_ = 0;   
    string pattern = "";
    for (unsigned int i = 0 ; i < compressedVec_.size() ; i++) {
        pattern = TextTools::toString(vec_[i][0]) + TextTools::toString(vec_[i][1]);
        lks_[i] = computeCoalLikelihood(compressedVec_[i], l);
        lnL_ -= patternToWeights_[pattern] * lks_[i];
    }

    return;
    
}


/*******************************************************************************/





//mpic++ -g -lbpp-core -lbpp-seq -lbpp-phyl -lboost_serialization -lboost_mpi CoalTools.cpp ReconciliationTools.cpp -o CoalTools

/*
 //Test main function
 int main(int args, char ** argv)
 {
 std::string spStr = "(((A,B),C),(D,E));";
 std::string gStr = "((A,(B,C)),(D,E));";
 //    std::string gStr = "(((A,D), (B,C)),E);";
 TreeTemplate <Node>* spTree = TreeTemplateTools::parenthesisToTree(spStr);
 TreeTemplate <Node>* geneTree = TreeTemplateTools::parenthesisToTree(gStr);
 
 std::map<std::string, std::string > seqSp;
 
 //    seqSp.insert( pair<std::string, std::string >("taxon1","taxon1") );
 //    seqSp.insert( pair<std::string, std::string >("taxon2","taxon2") );
 //    seqSp.insert( pair<std::string, std::string >("taxon3","taxon3") );
 //    seqSp.insert( pair<std::string, std::string >("taxon4","taxon4") );
 //    seqSp.insert( pair<std::string, std::string >("taxon5","taxon5") );
 //    seqSp.insert( pair<std::string, std::string >("taxon6","taxon6") );
 //    seqSp.insert( pair<std::string, std::string >("taxon7","taxon7") );
 //    seqSp.insert( pair<std::string, std::string >("taxon8","taxon8") );
 
 
 
 seqSp.insert( pair<std::string, std::string >("A","A") );
 seqSp.insert( pair<std::string, std::string >("B","B") );
 seqSp.insert( pair<std::string, std::string >("C","C") );
 seqSp.insert( pair<std::string, std::string >("D","D") );
 seqSp.insert( pair<std::string, std::string >("E","E") );
 
 //    seqSp.insert( pair<std::string, std::string >("F","F") );
 //    seqSp.insert( pair<std::string, std::string >("G","G") );
 //    seqSp.insert( pair<std::string, std::string >("H","H") );
 
 
 
 
 
 //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
 std::vector < std::vector< std::vector< std::vector< unsigned int > > > > coalCounts;
 std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
 std::vector< std::vector<unsigned int> > coalCounts3;
 std::vector< unsigned int > coalCounts4;
 //speciesIDs: vector of genetreenbnodes vectors of 3 (3 directions) ints
 std::vector <std::vector<unsigned int> > speciesIDs;
 std::vector < unsigned int > speciesIDs2;
 for (unsigned int i = 0 ; i < 2 ; i++ ) {
 coalCounts4.push_back(0);
 }
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 coalCounts3.push_back(coalCounts4);
 }
 for (unsigned int i = 0 ; i < 3 ; i++ ) {
 coalCounts2.push_back(coalCounts3);
 speciesIDs2.push_back(0);
 }
 for (unsigned int i = 0 ; i < geneTree->getNumberOfNodes() ; i++ ) {
 coalCounts.push_back(coalCounts2);
 speciesIDs.push_back(speciesIDs2);
 }
 
 breadthFirstreNumber (*spTree);
 std::map<std::string, int > spID = computeSpeciesNamesToIdsMap(*spTree);
 
 for(std::map<std::string, int >::iterator it = spID.begin(); it != spID.end(); it++){
 std::cout <<"it->first: "<< it->first << " : "<< it->second <<std::endl;
 }
 
 computeSubtreeCoalCountsPostorder(*spTree, 
 *geneTree, 
 geneTree->getRootNode(), 
 seqSp, 
 spID, 
 coalCounts,
 speciesIDs);
 //Add the starting lineage at the root
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 if (coalCounts[geneTree->getRootNode()->getId()][0][i][0]==0 && coalCounts[geneTree->getRootNode()->getId()][0][i][1] !=0)
 {
 //coalCounts[geneTree->getRootNode()->getId()][0][0][0] = 1;
 coalCounts[geneTree->getRootNode()->getId()][0][i][0]=1; 
 break;
 
 }
 }
 
 std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
 std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*geneTree, true) << std::endl;
 
 for (unsigned int i = 0 ; i < coalCounts[geneTree->getRootNode()->getId()][0].size() ; i++) {
 std::cout << "Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] << " Num coal out: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][1]<<std::endl;
 }
 
 std::vector<double> bls (spTree->getNumberOfNodes(), 2.0);
 
 double initialLikelihood = computeCoalLikelihood ( coalCounts[geneTree->getRootNode()->getId()][0], bls ) ;
 
 std::cout << "Initial Likelihood: "<< initialLikelihood <<std::endl;
 
 std::map <double, Node*> LksToNodes;
 //Now doing the preorder tree traversal
 Node * geneRoot = geneTree->getRootNode();
 std::vector <Node *> sons = geneRoot->getSons();
 if (sons.size()!=2) {
 std::cerr <<"Error: "<<sons.size()<< "sons at the root!"<<std::endl; 
 }
 
 LksToNodes[initialLikelihood]=sons[0];
 //We fill the likelihood and species ID data for the root node.
 //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
 //with other nodes.
 coalCounts[geneRoot->getId()][1] = coalCounts[geneRoot->getSon(1)->getId()][0];
 coalCounts[geneRoot->getId()][2] = coalCounts[geneRoot->getSon(0)->getId()][0];
 speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon(1)->getId()][0];
 speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon(0)->getId()][0];
 
 
 
 
 for (unsigned int i = 0; i< sons.size(); i++){
 for (unsigned int j =0; j<sons[i]->getNumberOfSons(); j++) {
 computeSubtreeCoalCountsPreorder(*spTree, 
 *geneTree, 
 sons[i], 
 seqSp, 
 spID, 
 coalCounts,
 bls, 
 speciesIDs, j, LksToNodes);
 }
 }
 for (unsigned int j = 0 ; j < geneTree->getNumberOfNodes() ; j++) {
 std::cout << "Node j: " <<std::endl;
 for (unsigned int i = 0 ; i < coalCounts[j][0].size() ; i++) {
 std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][0][i][0] << " Num coal out: "<< coalCounts[j][0][i][1]<<std::endl;
 std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][1][i][0] << " Num coal out: "<< coalCounts[j][1][i][1]<<std::endl;
 std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][2][i][0] << " Num coal out: "<< coalCounts[j][2][i][1]<<std::endl;
 
 }
 }
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 Newick *treeReader = new Newick(true);
 vector<Tree *> gTrees;
 //treeReader->read("testGeneTrees", gTrees);
 //   treeReader->read("geneTrees7", gTrees);
 //treeReader->read("5000SimulatedTrees8Taxa.trees", gTrees);
 // treeReader->read("mesquite500SimCoalTreesNe10000_Depth100000.proj.phy", gTrees);
 //  treeReader->read("mesquite500SimCoalTreesNe10000_Depth100000.proj.phy", gTrees);
 treeReader->read("40Sp_500GeneTrees.phy", gTrees);
 
 
 vector<Tree *> spTrees;
 //    treeReader->read("spTree", spTrees);
 //    treeReader->read("spTree2.0", spTrees);
 //    treeReader->read("spTree7", spTrees);
 //  treeReader->read("SpTree8.tree", spTrees);
 //    treeReader->read("1SimulatedSpTree8Taxa_Depth100000.tree", spTrees);
 treeReader->read("40Sp_First.phy", spTrees);
 
 delete treeReader;
 //TEMP    TreeTemplate <Node>* spTree = dynamic_cast<TreeTemplate<Node>*>(spTrees[0]);
 vector<string> leaves = spTree->getLeavesNames();
 for (unsigned int i = 0 ; i < leaves.size() ; i++) {
 seqSp.insert( pair<std::string, std::string >(leaves[i],leaves[i]) );
 }
 
 breadthFirstreNumber (*spTree);
 //TEMP   std::map<std::string, int > spID = computeSpeciesNamesToIdsMap(*spTree);
 
 for(std::map<std::string, int >::iterator it = spID.begin(); it != spID.end(); it++){
 std::cout <<"it->first: "<< it->first << " : "<< it->second <<std::endl;
 }
 
 
 
 //     Checking the likelihood computation against results from Rosenberg 2002: it works.
 //    //11
 //    vector<unsigned int > vec(2, 1);
 //    double loglk = computeCoalLikelihood ( vec, 1.0 ) ;
 //    std::cout << "11: "<< loglk<< std::endl;
 //    
 //    vec[1] = 2;
 //     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
 //    std::cout << "12: "<<loglk<< std::endl;
 //
 //    vec[0] = 2;
 //     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
 //    std::cout << "22: "<<loglk<< std::endl;
 //
 //    vec[1] = 3;
 //     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
 //    std::cout << "23: "<<loglk<< std::endl;
 //
 //    vec[0] = 3;
 //     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
 //    std::cout << "33: "<<loglk<< std::endl;
 
 
 
 
 
 //TEMP  std::vector<double> bls (spTree->getNumberOfNodes(), 2.0);
 //    for (unsigned int i = 0 ; i < bls.size() ; i++) {
 //        if ( spTree->getNode(i)->hasFather() )
 //            bls[i] = spTree->getNode(i)->getDistanceToFather();
 //    }
 
 std::cout <<"Species Tree: \n" <<    TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
 //TEMP   TreeTemplate <Node>* geneTree = 0;
 
 //allGeneCounts: nbSpeciesBranches vectors of gTrees.size() vectors of 2 ints
 std::vector < std::vector < std::vector<unsigned int> > >  allGeneCounts;
 std::vector< std::vector< unsigned int > > allGeneCounts2;
 std::vector< unsigned int > allGeneCounts3;
 for (unsigned int i = 0 ; i < 2 ; i++ ) {
 allGeneCounts3.push_back(0);
 }
 for (unsigned int i = 0 ; i < gTrees.size() ; i++ ) {
 allGeneCounts2.push_back(allGeneCounts3);
 }
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 allGeneCounts.push_back(allGeneCounts2);
 }
 
 
 
 for (unsigned int j = 0 ; j < gTrees.size() ; j++) {
 geneTree= dynamic_cast<TreeTemplate<Node>*>(gTrees[j]);
 //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
 std::vector < std::vector< std::vector< std::vector< unsigned int > > > > coalCounts;
 std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
 std::vector< std::vector<unsigned int> > coalCounts3;
 std::vector< unsigned int > coalCounts4;
 //speciesIDs: vector of genetreenbnodes vectors of 3 (3 directions) ints
 std::vector <std::vector<unsigned int> > speciesIDs;
 std::vector < unsigned int > speciesIDs2;
 for (unsigned int i = 0 ; i < 2 ; i++ ) {
 coalCounts4.push_back(0);
 }
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 coalCounts3.push_back(coalCounts4);
 }
 for (unsigned int i = 0 ; i < 3 ; i++ ) {
 coalCounts2.push_back(coalCounts3);
 speciesIDs2.push_back(0);
 }
 for (unsigned int i = 0 ; i < geneTree->getNumberOfNodes() ; i++ ) {
 coalCounts.push_back(coalCounts2);
 speciesIDs.push_back(speciesIDs2);
 }
 
 computeSubtreeCoalCountsPostorder(*spTree, 
 *geneTree, 
 geneTree->getRootNode(), 
 seqSp, 
 spID, 
 coalCounts,
 speciesIDs);
 //Add the starting lineage at the root
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 if (coalCounts[geneTree->getRootNode()->getId()][0][i][0]==0 && coalCounts[geneTree->getRootNode()->getId()][0][i][1] !=0)
 {
 //coalCounts[geneTree->getRootNode()->getId()][0][0][0] = 1;
 coalCounts[geneTree->getRootNode()->getId()][0][i][0]=1; 
 break;
 
 }
 }
 
 //            for (unsigned int i = 0 ; i < coalCounts[geneTree->getRootNode()->getId()][0].size() ; i++) {
 //         std::cout << "Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] << " Num coal out: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][1]<<std::endl;
 //         }
 //Now we compute the likelihood of the rooted gene tree:
 double loglk = computeCoalLikelihood ( coalCounts[geneTree->getRootNode()->getId()][0], bls ) ;
 //  std::cout << "lk: "<<exp(loglk)<<std::endl;
 for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
 allGeneCounts[i][j] = coalCounts[geneTree->getRootNode()->getId()][0][i];
 }
 }
 
 //If we want to do numerical optimization
 // BrentOneDimension *brentOptimizer = new BrentOneDimension();
 vector<double> blAnalytical = bls;
 for (unsigned int i = 0 ; i < allGeneCounts.size() ; i++) {
 CoalBranchLikelihood *brLikFunction = new CoalBranchLikelihood(allGeneCounts[i]);
 //Initialize BranchLikelihood:
 blAnalytical[i] = brLikFunction->initModel();
 
 //If we want to do numerical optimization
 
 //        ParameterList parameters;
 //        Parameter brLen = Parameter("BrLen", blAnalytical[i]);//, ExcludingPositiveReal (0) );
 //        parameters.addParameter(brLen);
 //        brLikFunction->setParameters(parameters);
 //        
 //        //Re-estimate branch length:
 //        brentOptimizer->setVerbose(0);
 //        brentOptimizer->setFunction(brLikFunction);
 //        brentOptimizer->getStopCondition()->setTolerance(0.000001);
 //        brentOptimizer->setInitialInterval(0, brLen.getValue()+2);
 //        brentOptimizer->init(parameters);
 //        brentOptimizer->optimize();
 //        bls[i] =brentOptimizer->getParameters().getParameter("BrLen").getValue();
 //        std::cout <<"Value for Sp branch "<<i<<": "<< bls[i] << std::endl;
 
 }
 for (unsigned int i = 0 ; i < bls.size() ; i++) {
 if ( spTree->getNode(i)->hasFather() )
 spTree->getNode(i)->setDistanceToFather(bls[i]);
 }
 std::cout <<"Species Tree optimized: \n" << TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
 std::cout <<"Species Tree optimized: \n" << TreeTemplateTools::treeToParenthesis(*spTree, false) << std::endl;
 
 for (unsigned int i = 0 ; i < bls.size() ; i++) {
 if ( spTree->getNode(i)->hasFather() )
 spTree->getNode(i)->setDistanceToFather(blAnalytical[i]);
 }
 std::cout <<"Species Tree analytical: \n" << TreeTemplateTools::treeToParenthesis(*spTree, true) << std::endl;
 std::cout <<"Species Tree analytical: \n" << TreeTemplateTools::treeToParenthesis(*spTree, false) << std::endl;
 
 return 0;
 
 }
 */





