/*
 *  GenericTreeExplorationAlgorithms.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 17/12/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */

#include "GenericTreeExplorationAlgorithms.h"


/************************************************************************
 * Change the root of the tree by changing the outgroup.
 ************************************************************************/

void changeRoot(TreeTemplate<Node> &tree, int newOutGroup) 
{
  std::cout <<"Changing the outgroup to node :"<<newOutGroup<< std::endl;
  tree.newOutGroup(newOutGroup); 
}




std::string nodeToParenthesisBIS(const Tree & tree, int nodeId, bool bootstrap) throw (NodeNotFoundException)
{
  if (tree.hasNodeName(nodeId)) { }
  if(!tree.hasNode(nodeId)) throw NodeNotFoundException("nodeToParenthesisBIS", nodeId);
  std::ostringstream s;
  if(tree.isLeaf(nodeId))
    { }
  else
    {
    s << "(";
    std::vector<int> sonsId = tree.getSonsId(nodeId);
    printVector(sonsId);
    std::cout << std::endl;
    s << nodeToParenthesisBIS(tree, sonsId[0], bootstrap);
    for(unsigned int i = 1; i < sonsId.size(); i++)
      {
      s << "," << nodeToParenthesisBIS(tree, sonsId[i], bootstrap);
      }
    s << ")";
    
    if(bootstrap)
      {
      if(tree.hasBranchProperty(nodeId, "BOOTSTRAP"))
        s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(nodeId, "BOOTSTRAP"))->getValue());
      }
    
    }
  if(tree.hasDistanceToFather(nodeId)) s << ":" << tree.getDistanceToFather(nodeId);
  return s.str();  
}


/************************************************************************
 * Make a SPR between two nodes. A subtree is cut at node with Id cutNodeId, 
 * and pasted beneath node with Id newFatherId.
 ************************************************************************/
std::vector<Node*> makeSPR(TreeTemplate<Node> &tree, 
						   int cutNodeId, 
						   int newBrotherId, 
						   bool verbose, 
						   bool returnNodesToUpdate) {
	
	Node *cutNode, *newBrother, *oldFather, *oldGrandFather, *brother, *newBrothersFather, *N;
    std::vector<Node*> nodesToUpdate;
	double dist = 0.1;  
	
	if (verbose)
		std::cout <<"\t\t\tMaking a SPR, moving node "<<cutNodeId<< " as brother of node "<< newBrotherId<< std::endl;
	newBrother = tree.getNode(newBrotherId);
	cutNode = tree.getNode(cutNodeId); 
	
    std::vector <int> nodeIds =tree.getNodesId();
	
	if ((!(cutNode->hasFather()))||(!(newBrother->hasFather()))) {
		std::cout <<"Error in makeSPR"<< std::endl;
		if (!(cutNode->hasFather())) {
			std::cout << " Node "<<cutNodeId<<" has no father"<< std::endl;
		}
		else {
			std::cout << " Node "<<newBrotherId<<" has no father"<< std::endl;
		}
		exit (-1);
	}
	
	oldFather = cutNode->getFather();
	//Get all old brothers ; a binary tree is assumed here (because of the "break")
	for(unsigned int i=0;i<oldFather->getNumberOfSons();i++)
		if(oldFather->getSon(i)!=cutNode){brother=oldFather->getSon(i); break;}
	
	newBrothersFather = newBrother->getFather();
	if (newBrothersFather == oldFather) {
		return nodesToUpdate;
	}
	
	if (!(oldFather->hasFather())) {//we displace the outgroup, need to reroot the tree
		//NB : brother is the other son of the root
		int id0 = oldFather->getId();
		int idBrother = brother->getId();
		N=new Node();
		
		N->addSon(newBrother);
		newBrother->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
				
		//we remove cutNode from its old neighborhood
		for(unsigned int i=0;i<oldFather->getNumberOfSons();i++) {
			if(oldFather->getSon(i)==cutNode){oldFather->removeSon(i); break;} 
		}
		// we move node cutNode
		N->addSon(cutNode);
		cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
		
		// update N neighbours 
		for(unsigned int i=0;i<newBrothersFather->getNumberOfSons();i++)
			if(newBrothersFather->getSon(i)==newBrother){newBrothersFather->setSon(i, N); break;}
		N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
								
		tree.rootAt(brother->getId());
		for(unsigned int i=0;i<brother->getNumberOfSons();i++) {
			if(brother->getSon(i)==oldFather){brother->removeSon(i);break;}
		}
		delete oldFather;
		//We renumber the nodes
		brother->setId(id0);
		N->setId(idBrother);
	}
	else  {  
		int id0 = oldFather->getId();
		//we create a new node N which will be the father of cutNode and newBrother
		N=new Node();
		
		N->addSon(newBrother);
		newBrother->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
		// we move node cutNode
		N->addSon(cutNode);
		
		cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
		
		// update N neighbours    
		for(unsigned int i=0;i<newBrothersFather->getNumberOfSons();i++)
			if(newBrothersFather->getSon(i)==newBrother){newBrothersFather->setSon(i, N); break;}
		N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
		oldGrandFather = oldFather->getFather();
		for(unsigned int i=0;i<oldGrandFather->getNumberOfSons();i++)
			if(oldGrandFather->getSon(i)==oldFather){oldGrandFather->setSon(i, brother); break;}
		brother->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
		delete oldFather;
		N->setId(id0);
	}
    if (returnNodesToUpdate) {
		//nodesToUpdate.push_back(newBrother);
		nodesToUpdate.push_back(cutNode);
		nodesToUpdate.push_back(N);
		//nodesToUpdate.push_back(brother);

		if (cutNode->getNumberOfSons() > 0)
			nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, cutNode->getSons() );
		if (brother->getNumberOfSons() > 0)
			nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, brother->getSons() );
		/*if (newBrothersFather->getNumberOfSons() > 0)
		 nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, newBrothersFather->getSons() );*/
		/*	if (newBrother->getNumberOfSons() > 0)
		 nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, newBrother->getSons() );
		 */
		/*		if (cutNode->hasSons() ) {
		 for(unsigned int i = 0 ; i < cutNode->getNumberOfSons(); i++) {
		 nodesToUpdate.push_back( cutNode->getSon(i) );
		 }
		 }*/
        std::vector<Node*> nodesToUpdate2 = TreeTemplateTools::getPathBetweenAnyTwoNodes (*brother, *cutNode, true);
		nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, nodesToUpdate2);

		// nodesToUpdate.push_back(newBrother);
        //We also need to add the path to the root, if the root is not included already
        if ( ! VectorTools::contains(nodesToUpdate, tree.getRootNode() ) ) {
            nodesToUpdate2 = TreeTemplateTools::getPathBetweenAnyTwoNodes (*brother, *(tree.getRootNode()), true);
            nodesToUpdate = VectorTools::vectorUnion (nodesToUpdate, nodesToUpdate2);
        }
		unsigned int num = nodesToUpdate.size();
	//	std::cout <<"BEFORE: "<< num <<std::endl;
		for (unsigned int i = 0 ; i < num ; i++) {
			for (unsigned int j = 0 ; j < nodesToUpdate[i]->getNumberOfSons(); j++) {
 				if (! VectorTools::contains(nodesToUpdate, nodesToUpdate[i]->getSon(j)) ) {
					nodesToUpdate.push_back(nodesToUpdate[i]->getSon(j));
				}
			}
			if ( nodesToUpdate[i]->hasFather() && !VectorTools::contains(nodesToUpdate, nodesToUpdate[i]->getFather() ) ) {
				nodesToUpdate.push_back(nodesToUpdate[i]->getFather() );
			}
		}
	}
//	std::cout <<"AFTER: "<< nodesToUpdate.size() <<std::endl;
    
    return nodesToUpdate;
}


/************************************************************************
 * Make a SPR between two nodes. A subtree is cut at node with Id cutNodeId, 
 * and pasted beneath node with Id newFatherId. Beware this function does not necessarily return a proper binary tree.
 ************************************************************************/
void makeMuffatoSPR(TreeTemplate<Node> &tree, 
						   Node* cutNode, 
						   Node* newFather, 
                           bool verbose) {
	
	Node *oldFather, *N;
	double dist = 0.1;  
	
	if (verbose)
		std::cout <<"\t\t\tMaking a Muffato SPR, moving node "<< cutNode->getId() << " as son of node "<< newFather->getId() << std::endl;

    //std::vector <int> nodeIds =tree.getNodesId();
	
	if ( ! ( cutNode->hasFather() ) ) {
		std::cout <<"Error in makeMuffatoSPR"<< std::endl;
			std::cout << " Node "<< cutNode->getId() <<" has no father"<< std::endl;
		exit (-1);
	}
	
	oldFather = cutNode->getFather();
/*    std::vector<Node*> oldBrothers;
	//Get all old brothers 
	for(unsigned int i=0;i<oldFather->getNumberOfSons();i++)
		if(oldFather->getSon(i)!=cutNode){oldBrothers.push_back(oldFather->getSon(i));}
	
	newBrothersFather = newBrother->getFather();
	if (newBrothersFather == oldFather) {
		return nodesToUpdate;
	}*/
	
    //If the newFather already has two sons, then we need to create a new node (to keep a binary tree)    
    if ( newFather->getNumberOfSons() > 1 ) {
		//we create a new node N which will be the father of cutNode and newBrother
		N=new Node();
      //  std::cout << "newFather id "<< newFather->getId() << std::endl;

        std::vector< Node* > newBrothers = newFather->getSons();
        size_t siz = newBrothers.size();
    //    std::cout << "number of brothers "<< siz << std::endl;

        for (size_t i = 0 ; i < siz; i++) {
            newFather->removeSon( newBrothers[i] );
       //     std::cout << "newBrothers id "<< newBrothers[i]->getId() << std::endl;
            N->addSon( newBrothers[i] );//TEST PRINT id OF newBrothers[i]
            newBrothers[i]->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        }
        newFather->addSon(N);
        N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    }
    // we move node cutNode
  //  std::cout << "cutNode id "<< cutNode->getId() << std::endl;

    oldFather->removeSon( cutNode );
 //   std::cout << "cutNode id again "<< cutNode->getId() << std::endl;

    newFather->addSon(cutNode);
    cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE

    
    //In principle we may not need to clean up behind us, because we only use this function from editDuplicationNodesMuffato
    
    return ;
}





/************************************************************************
 * Make a NNI around a particular branch. The node with Id nodeId is exchanged with its uncle. 
 ************************************************************************/
void makeNNI(TreeTemplate<Node> &tree, int nodeId) {
  std::cout <<"Making a NNI involving node :"<<nodeId<< std::endl;
  double dist = 0.1;  
  Node * son    = tree.getNode(nodeId);
  
  if(!son->hasFather()) throw NodeException("makeNNI(). Node 'son' must not be the root node.", nodeId);
  Node * parent = son->getFather();
  
  if(!parent->hasFather()) throw NodeException("makeNNI(). Node 'parent' must not be the root node.", parent->getId());
  Node * grandFather = parent->getFather();
  //From here: Bifurcation assumed.
  //In case of multifurcation, an arbitrary uncle is chosen.
  //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  unsigned int parentPosition = grandFather->getSonPosition(parent);
  Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
  
  parent->removeSon(son);
  grandFather->removeSon(uncle);
  parent->addSon(uncle);
  uncle->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
  grandFather->addSon(son);
  son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
  return;
}



/************************************************************************
 * Defines the nodes where a subtree can be regrafted.
 ************************************************************************/
void buildVectorOfRegraftingNodes(TreeTemplate<Node> &tree, int nodeForSPR, std::vector <int> & nodeIdsToRegraft) {
  
  std::vector <int> allNodeIds = tree.getNodesId();
  std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
  forbiddenIds.push_back(tree.getRootNode()->getId());
  int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
  int brotherId;
  //Get one brother ; a binary tree is supposed here (because of the "break")
  for(unsigned int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
    if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
  if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
    forbiddenIds.push_back(brotherId);
  }
  
  std::vector <int> toRemove;
  for (unsigned int i = 0 ; i< allNodeIds.size() ; i++) {
    if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
      toRemove.push_back(i);
    }
  }
  sort(toRemove.begin(), toRemove.end(), cmp);
  for (unsigned int i = 0 ; i< toRemove.size() ; i++) {
    std::vector<int>::iterator vi = allNodeIds.begin();
    allNodeIds.erase(vi+toRemove[i]);
  }
  
  //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
  nodeIdsToRegraft = allNodeIds;
  
}


/************************************************************************
 * Defines the nodes where a subtree can be regrafted, at a given distance (in number of branches) from the node to be regrafted.
 ************************************************************************/
//Utilitary functions

std::vector<int> getRemainingNeighbors(const Node * node1, const Node * node2)
{
  std::vector<const Node *> neighbors = node1->getNeighbors();
  std::vector<int> neighbors2;
  for(unsigned int k = 0; k < neighbors.size(); k++)
    {
    const Node * n = neighbors[k];
    if(n != node2 ) neighbors2.push_back(n->getId());
    }
  return neighbors2;
}


/***************************************************************************************/


void getRemainingNeighborsUntilDistance(TreeTemplate<Node> &tree, const Node * node1, const Node * node2, int distance, int d, std::vector <int> & neighbors)
{
  std::vector<const Node *> neighbors1 = node1->getNeighbors();
  std::vector<int> neighbors2;
  for(unsigned int k = 0; k < neighbors1.size(); k++)
    {
    const Node * n = neighbors1[k];
    if(n != node2 ) neighbors2.push_back(n->getId());
    }
  VectorTools::append(neighbors,neighbors2);
  d=d+1;
  if (d<distance) {
    for(unsigned int k = 0; k < neighbors2.size(); k++)
      {
      getRemainingNeighborsUntilDistance(tree, tree.getNode(neighbors2[k]), node1, distance, d, neighbors);
      }
  }
}

/***************************************************************************************/


void getNeighboringNodesIdLimitedDistance (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds) {
  int d=1;
  std::vector <Node *> neighbors = tree.getNode(nodeId)->getNeighbors();
  for (unsigned int i =0 ; i < neighbors.size(); i++) {
    neighboringNodeIds.push_back(neighbors[i]->getId());
    getRemainingNeighborsUntilDistance(tree, neighbors[i], tree.getNode(nodeId), distance, d, neighboringNodeIds);
  }
}
/***************************************************************************************/


void getNeighboringNodesIdLimitedDistanceLowerNodes (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds) {
    int d=1;
    std::vector <Node *> neighbors = tree.getNode(nodeId)->getNeighbors();
    for (unsigned int i =0 ; i < neighbors.size(); i++) {
        if (neighbors[i] != tree.getNode(nodeId)->getFather()) {
            neighboringNodeIds.push_back(neighbors[i]->getId());
            getRemainingNeighborsUntilDistance(tree, neighbors[i], tree.getNode(nodeId), distance, d, neighboringNodeIds);
        }
        else {
            neighboringNodeIds.push_back(neighbors[i]->getId());
        }
    }
}

/***************************************************************************************/

//This whole function efficiency may well be improved
void buildVectorOfRegraftingNodesLimitedDistance(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft) {
  
//  Node * N = tree.getRootNode();

  // std::vector <int> allNodeIds = tree.getNodesId();
  std::vector <int> allNodeIds;
  getNeighboringNodesIdLimitedDistance(tree, nodeForSPR, distance, allNodeIds);
    
  //std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)->getFather()));
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
    int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
    forbiddenIds.push_back(oldFatherId);

/*
  std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));

  

  int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();

  forbiddenIds.push_back(oldFatherId);
  std::cout <<"FatherID: "<< oldFatherId <<std::endl; 

  int brotherId;
  //Get one brother ; a binary tree is supposed here (because of the "break")
  for(int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
    if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
 // if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
    forbiddenIds.push_back(brotherId);
 // }
  std::cout <<"BrotherID: "<< brotherId <<std::endl; 
*/
  forbiddenIds.push_back(tree.getRootNode()->getId());
  
  //We remove the nodes that are not appropriate for regrafting
  std::vector <int> toRemove;
  for (unsigned int i = 0 ; i< allNodeIds.size() ; i++) {
    if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
      toRemove.push_back(i);
    }
  }

/*  std:: cout <<"nodeForSPR: "<< nodeForSPR <<"; FORBIDDEN IDS: "<<std::endl;
  VectorTools::print(forbiddenIds);*/
  sort(toRemove.begin(), toRemove.end(), cmp);
  /*VectorTools::print(forbiddenIds);
  sort(allNodeIds.begin(), allNodeIds.end(), anticmp);*/
  for (unsigned int i = 0 ; i< toRemove.size() ; i++) {
    std::vector<int>::iterator vi = allNodeIds.begin();
    allNodeIds.erase(vi+toRemove[i]);
  }
  
  //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
  nodeIdsToRegraft = allNodeIds;

}




/***************************************************************************************/

//This whole function efficiency may well be improved
void buildVectorOfRegraftingNodesLimitedDistanceLowerNodes(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft) {
    
    //  Node * N = tree.getRootNode();
    
    // std::vector <int> allNodeIds = tree.getNodesId();
    std::vector <int> allNodeIds;
    getNeighboringNodesIdLimitedDistanceLowerNodes(tree, nodeForSPR, distance, allNodeIds);
    
  /*  std::cout << "Before forbidden: "<< std::endl;
    VectorTools::print(allNodeIds);
    */
    
    //std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)->getFather()));
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
    
    /*
     std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
     
     
     
     int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
     
     forbiddenIds.push_back(oldFatherId);
     std::cout <<"FatherID: "<< oldFatherId <<std::endl; 
     
     int brotherId;
     //Get one brother ; a binary tree is supposed here (because of the "break")
     for(int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
     if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
     // if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
     forbiddenIds.push_back(brotherId);
     // }
     std::cout <<"BrotherID: "<< brotherId <<std::endl; 
     */
    forbiddenIds.push_back(tree.getRootNode()->getId());
    
    //We remove the nodes that are not appropriate for regrafting
    std::vector <int> toRemove;
    for (unsigned int i = 0 ; i< allNodeIds.size() ; i++) {
        if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
            toRemove.push_back(i);
        }
    }
    
    /*  std:: cout <<"nodeForSPR: "<< nodeForSPR <<"; FORBIDDEN IDS: "<<std::endl;
     VectorTools::print(forbiddenIds);*/
    sort(toRemove.begin(), toRemove.end(), cmp);
    /*VectorTools::print(forbiddenIds);
     sort(allNodeIds.begin(), allNodeIds.end(), anticmp);*/
    for (unsigned int i = 0 ; i< toRemove.size() ; i++) {
        std::vector<int>::iterator vi = allNodeIds.begin();
        allNodeIds.erase(vi+toRemove[i]);
    }
    
    //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
    nodeIdsToRegraft = allNodeIds;
    
}




/************************************************************************
 * Makes a modification, knowing what previous modifications have been done.
 ************************************************************************/
void makeDeterministicModifications(TreeTemplate<Node> &tree, size_t & nodeForNNI, size_t & nodeForSPR, size_t & nodeForRooting) {
  if (nodeForNNI < tree.getNumberOfNodes()) {//Make a NNI move
    if (nodeForNNI <3) {
      if (nodeForRooting<tree.getNumberOfNodes()) {
        changeRoot(tree, nodeForRooting);
        nodeForRooting++;
      }
      else {
        nodeForNNI=3;
        nodeForRooting = 1;
        makeNNI(tree, nodeForNNI);
        nodeForNNI++;
      }
    }
    else {
      makeNNI(tree, nodeForNNI);
      nodeForNNI++;
    }
  }
  else {//we make a SPR
    if (nodeForSPR == tree.getNumberOfNodes()) {
      nodeForSPR = 1;
      nodeForNNI = 0;
      changeRoot(tree, nodeForRooting);
      nodeForRooting++;
    }
    std::vector <int> allNodeIds;
    Node * N = tree.getRootNode();
    std::vector <int> firstHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(0)));
    std::vector <int> secondHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(1)));
    
    if (secondHalfNodeIds.size()<firstHalfNodeIds.size()) {
      allNodeIds = firstHalfNodeIds;
      for (unsigned int i = 0 ; i< secondHalfNodeIds.size() ; i++ ) {
        allNodeIds.push_back(secondHalfNodeIds[i]);
      }
    }
    else {
      allNodeIds = secondHalfNodeIds;
      for (unsigned int i = 0 ; i< firstHalfNodeIds.size() ; i++ ) {
        allNodeIds.push_back(firstHalfNodeIds[i]);
      }
    }
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
    
    int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
    int brotherId;
    //Get one brother ; a binary tree is supposed here (because of the "break")
    for(size_t i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
      if( (size_t) tree.getNode(oldFatherId)->getSon(i)->getId() != nodeForSPR ){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
    if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
      forbiddenIds.push_back(brotherId);
    }
    std::vector <int> toRemove;
    for (size_t i = 0 ; i< allNodeIds.size() ; i++) {
      if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
        toRemove.push_back(i);
      }
    }
    sort(toRemove.begin(), toRemove.end(), cmp);
    for (size_t i = 0 ; i< toRemove.size() ; i++) {
      std::vector<int>::iterator vi = allNodeIds.begin();
      allNodeIds.erase(vi+toRemove[i]);
    }
    int ran2 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    int newBrotherId = allNodeIds[ran2];
    makeSPR(tree, nodeForSPR, newBrotherId);
    nodeForSPR++;
  }
  
}






/************************************************************************
 * Makes NNIs and root changings only.
 ************************************************************************/
void makeDeterministicNNIsAndRootChangesOnly(TreeTemplate<Node> &tree, size_t & nodeForNNI, size_t & nodeForRooting, const bool fixedOutgroupSpecies_) {
  if (nodeForNNI < tree.getNumberOfNodes() ) {//Make a NNI or rerooting move
    if (nodeForNNI < 3) {
      if ( !fixedOutgroupSpecies_ && nodeForRooting < tree.getNumberOfNodes() ) {//Make a rerooting move
        changeRoot(tree, nodeForRooting);
        nodeForRooting++;
      }
      else { //Make a NNI move
		  nodeForNNI = 3;
		  if ( !fixedOutgroupSpecies_ ) {
			  makeNNI(tree, nodeForNNI);
			  nodeForNNI++;
		  }
		  else {
			  //All grandsons of the root should not be considered as appropriate for a NNI
			  //because we want to consider the root.
			  Node *n = tree.getNode(nodeForNNI);
			  Node *root = tree.getRootNode() ;
			  while ( n->getFather()->getFather() == root) {
				  nodeForNNI++;
				  n = tree.getNode(nodeForNNI);
			  }
			  makeNNI(tree, nodeForNNI);
			  nodeForNNI++;
		  } 
	  }
    }
    else {
		if ( fixedOutgroupSpecies_ ) {
			//All grandsons of the root should not be considered as appropriate for a NNI
			//because we want to consider the root.
			Node *n = tree.getNode(nodeForNNI);
			Node *root = tree.getRootNode() ;
			while ( n->getFather()->getFather() == root) {
				nodeForNNI++;
				n = tree.getNode(nodeForNNI);
			}
		}
		makeNNI(tree, nodeForNNI);
		nodeForNNI++;
    }
  }
  else {//we reset the loop rooting-NNIs
      nodeForNNI = 0;
	  if ( !fixedOutgroupSpecies_ ) {
		  nodeForRooting = 4; //We don't want to root on nodes 1 or 2, 
		  //the two sons of the root.
		  //We do not want to root on node 3 either, 
		  //as a NNI already provides this tree.
		  changeRoot(tree, nodeForRooting);
		  nodeForRooting++;
	  }
	  else {
		  nodeForNNI = 3;
		  //All grandsons of the root should not be considered as appropriate for a NNI
		  //because we want to consider the root.
		  Node *n = tree.getNode(nodeForNNI);
		  Node *root = tree.getRootNode() ;
		  while ( n->getFather()->getFather() == root) {
			  nodeForNNI++;
			  n = tree.getNode(nodeForNNI);
		  }
	  }
  }
}







/************************************************************************
 * Procedure that makes sure that an NNI or rerooting one is about to make has
 * not been computed already.
 * If all NNIs or rerootings have been made, return true.
 ************************************************************************/

bool checkChangeHasNotBeenDone(TreeTemplate<Node> &tree, TreeTemplate<Node> *bestTree, size_t & nodeForNNI, 
                               size_t & nodeForRooting, std::vector < double >  &NNILks, 
                               std::vector < double >  &rootLks)
{
  if (nodeForNNI < tree.getNumberOfNodes()) 
    {//Make a NNI or rerooting move
      if (nodeForNNI <3) 
        {
        if (nodeForRooting<tree.getNumberOfNodes()) 
          {//Make a rerooting move
            while ((rootLks[nodeForRooting] < NumConstants::VERY_BIG()) && (nodeForRooting < tree.getNumberOfNodes())) {
              //std::cout<<rootLks[nodeForRooting]<<" 1_NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
              
              nodeForRooting++;
            }
          }
        if (nodeForRooting >= tree.getNumberOfNodes()) 
          { //Make a NNI move
            nodeForNNI=3;
            //   while ((NNILks[bestTree->getNode(nodeForNNI-1)->getFather()->getId()] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))
            while ((NNILks[nodeForNNI-1] < NumConstants::VERY_BIG()) && (nodeForNNI < tree.getNumberOfNodes()))
              {
              //std::cout<<NNILks[nodeForNNI-1]<<" NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
              nodeForNNI = nodeForNNI+2;
              }
          }
        }
      else 
        {
        //  while ((NNILks[bestTree->getNode(nodeForNNI-1)->getFather()->getId()] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))       
        while ((NNILks[nodeForNNI-1] < NumConstants::VERY_BIG()) && (nodeForNNI < tree.getNumberOfNodes()))
          {
          //std::cout<<NNILks[nodeForNNI-1]<<" NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
          nodeForNNI = nodeForNNI+2;
          }
        }
    }
  else 
    {//we reset the loop rooting-NNIs
      nodeForNNI = 2;
      while ((rootLks[nodeForRooting] < NumConstants::VERY_BIG()) && (nodeForRooting < tree.getNumberOfNodes())) {
        //std::cout<<rootLks[nodeForRooting]<<" 2_NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
        nodeForRooting++;
      }
    }
  //If we have tried all possible rearrangements:
  if (nodeForRooting == tree.getNumberOfNodes() && nodeForNNI == tree.getNumberOfNodes())
    {
    return true;
    }
  else {
    return false;
  }
}



/************************************************************************
 * Function to removes leaves from a tree.
 ************************************************************************/

void dropLeaves(TreeTemplate<Node> & tree, const std::vector<string> &spToDrop) {
    for (unsigned int i = 0 ; i < spToDrop.size() ; i++) {
        TreeTemplateTools::dropLeaf(tree, spToDrop[i]);
    }
    return;
}

/************************************************************************
 * Function to build a MRP bionj tree from a collection of trees.
 * Faster than TreeTools::MRP as we don't do the NNI exploration part.
 ************************************************************************/
Tree* MRP(const vector<Tree*>& vecTr)
{
    //matrix representation
    VectorSiteContainer* sites = TreeTools::MRPEncode(vecTr);
    
    //starting bioNJ tree
    const DNA* alphabet= dynamic_cast<const DNA*>(sites->getAlphabet());
	auto_ptr<SubstitutionModel>   jc (new JCnuc( alphabet ) );
    //ConstantDistribution constRate1(1.);
	auto_ptr<DiscreteDistribution>   constRate (new ConstantDistribution(1.) );
    DistanceEstimation distFunc(jc.release(), constRate.release(), sites, 0, true);
    BioNJ bionjTreeBuilder;
    bionjTreeBuilder.setDistanceMatrix(*(distFunc.getMatrix()));
    bionjTreeBuilder.computeTree();
    if (ApplicationTools::message) ApplicationTools::message->endLine();
    TreeTemplate<Node>* tree = new TreeTemplate<Node>(*bionjTreeBuilder.getTree());
    return tree;
}


/************************************************************************
 * Function to root a tree based on a list of outgroup taxa.
 ************************************************************************/
void rootTreeWithOutgroup (TreeTemplate<Node> &tree, const std::vector<std::string> outgroupTaxa) throw ( TreeException )
{
	if (outgroupTaxa.size() == 0 ) {
		std::cerr<< "Error: trying to root the species tree with an empty outgroup list!"<<std::endl;
		exit(-1);
	}
	std::vector<std::string> allLeafNames = tree.getLeavesNames();
	std::vector<std::string> ingroupTaxa ; 
	VectorTools::diff <std::string> (allLeafNames, *(const_cast< std::vector<std::string>* > (&outgroupTaxa)), ingroupTaxa) ;
	if ( ingroupTaxa.size() == 0 ) {
		std::cerr<< "Error: trying to root the species tree with an outgroup list that contains all taxa!"<<std::endl;
		exit(-1);
	}
	Node *n = tree.getNode( outgroupTaxa[0] );
	std::vector<std::string> descendantNames = TreeTemplateTools::getLeavesNames( *n );
	while ( descendantNames.size() < outgroupTaxa.size() ) {
		n = n->getFather() ;
		descendantNames = TreeTemplateTools::getLeavesNames( *n );
	}
	VectorTools::diff <std::string> (descendantNames, *(const_cast< std::vector<std::string>* > (&outgroupTaxa)), ingroupTaxa) ;
	if ( ingroupTaxa.size() == 0 ) {
		tree.newOutGroup( n ) ;
	}
	else {
		throw ( TreeException("Error: trying to root the species tree with an outgroup but the outgroup is not monophyletic!", &tree ) );
		exit(-1);
	}
	return;	
}



/************************************************************************
 * Function to test whether a tree is rooted according to a list of outgroup taxa.
 ************************************************************************/
bool isTreeRootedWithOutgroup (const TreeTemplate<Node> &tree, const std::vector<std::string> outgroupTaxa) {
	if (outgroupTaxa.size() == 0 ) {
		std::cerr<< "Error: trying to test whether the species tree is rooted with an empty outgroup list!"<<std::endl;
		exit(-1);
	}
	const Node *n = tree.getNode( outgroupTaxa[0] );
	std::vector<std::string> descendantNames = TreeTemplateTools::getLeavesNames( *n );
	while ( descendantNames.size() < outgroupTaxa.size() || ! n->hasFather() ) {
		n = n->getFather() ;
		descendantNames = TreeTemplateTools::getLeavesNames( *n );
	}
	std::vector<std::string> testTaxa ; 

	VectorTools::diff < std::string >(descendantNames, *(const_cast< std::vector<std::string>* > (&outgroupTaxa)), testTaxa) ;
	if ( testTaxa.size() == 0 ) {
		return true;
	}
	else {
		return false;
	}
}


