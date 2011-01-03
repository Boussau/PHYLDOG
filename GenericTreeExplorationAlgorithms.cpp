/*
 *  GenericTreeExplorationAlgorithms.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 17/12/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */

#include "GenericTreeExplorationAlgorithms.h"





//To sort in ascending order
bool cmp( int a, int b ) {
  return a > b;
}  


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
void makeSPR(TreeTemplate<Node> &tree, int cutNodeId, int newBrotherId) {
  
  Node *cutNode, *newBrother, *oldFather, *oldGrandFather, *brother, *newBrothersFather, *N;
  double dist = 0.1;  
  
  std::cout <<"Making a SPR, moving node "<<cutNodeId<< " as brother of node "<< newBrotherId<< std::endl;
  newBrother = tree.getNode(newBrotherId);
  cutNode = tree.getNode(cutNodeId); 
  std::vector <int> nodeIds =tree.getNodesId();
  
  if ((!(cutNode->hasFather()))||(!(newBrother->hasFather()))) {
    std::cout <<"Error in makeSPR"<< std::endl;
    if (!(cutNode->hasFather())) {
      std::cout << " Node "<<cutNodeId<<"has no father"<< std::endl;
    }
    else {
      std::cout << " Node "<<newBrotherId<<"has no father"<< std::endl;
    }
    exit (-1);
  }
  
  oldFather = cutNode->getFather();
  //Get all old brothers ; a binary tree is supposed here (because of the "break")
  for(int i=0;i<oldFather->getNumberOfSons();i++)
    if(oldFather->getSon(i)!=cutNode){brother=oldFather->getSon(i); break;}
  newBrothersFather = newBrother->getFather();
  
  if (newBrothersFather == oldFather) {
    return;
  }
  
  if (!(oldFather->hasFather())) {//we displace the outgroup, need to reroot the tree
                                  //NB : brother is the other son of the root
    int id0 = oldFather->getId();
    int idBrother = brother->getId();
    N=new Node();
    
    N->addSon(newBrother);
    newBrother->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    
    
    //we remove cutNode from its old neighborhood
    for(int i=0;i<oldFather->getNumberOfSons();i++) {
      if(oldFather->getSon(i)==cutNode){oldFather->removeSon(i); break;} 
    }
    // we move node cutNode
    N->addSon(cutNode);
    cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    
    // update N neighbours 
    for(int i=0;i<newBrothersFather->getNumberOfSons();i++)
      if(newBrothersFather->getSon(i)==newBrother){newBrothersFather->setSon(i, N); break;}
    N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    
    tree.rootAt(brother->getId());
    for(int i=0;i<brother->getNumberOfSons();i++) {
      if(brother->getSon(i)==oldFather){brother->removeSon(i);break;}
    }
    delete oldFather;
    //We renumber the nodes
    brother->setId(id0);
    N->setId(idBrother);
    return;   
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
    for(int i=0;i<newBrothersFather->getNumberOfSons();i++)
      if(newBrothersFather->getSon(i)==newBrother){newBrothersFather->setSon(i, N); break;}
    N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    oldGrandFather = oldFather->getFather();
    for(int i=0;i<oldGrandFather->getNumberOfSons();i++)
      if(oldGrandFather->getSon(i)==oldFather){oldGrandFather->setSon(i, brother); break;}
    brother->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
    delete oldFather;
    N->setId(id0);
    return;
  }
  
  
  
  
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
  
  
  Node * N = tree.getRootNode();
  
  std::vector <int> allNodeIds = tree.getNodesId();
  std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
  forbiddenIds.push_back(tree.getRootNode()->getId());
  int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
  int brotherId;
  //Get one brother ; a binary tree is supposed here (because of the "break")
  for(int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
    if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
  if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
    forbiddenIds.push_back(brotherId);
  }
  
  std::vector <int> toRemove;
  for (int i = 0 ; i< allNodeIds.size() ; i++) {
    if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
      toRemove.push_back(i);
    }
  }
  sort(toRemove.begin(), toRemove.end(), cmp);
  for (int i = 0 ; i< toRemove.size() ; i++) {
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

//This whole function efficiency may well be improved
void buildVectorOfRegraftingNodesLimitedDistance(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft) {
  
  
  Node * N = tree.getRootNode();
  
  // std::vector <int> allNodeIds = tree.getNodesId();
  std::vector <int> allNodeIds;
  getNeighboringNodesIdLimitedDistance(tree, nodeForSPR, distance, allNodeIds);
  std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
  forbiddenIds.push_back(tree.getRootNode()->getId());
  int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
  forbiddenIds.push_back(oldFatherId);
  
	
  int brotherId;
  //Get one brother ; a binary tree is supposed here (because of the "break")
  for(int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
    if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
  if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
    forbiddenIds.push_back(brotherId);
  }
  
  //We remove the nodes that are not appropriate for regrafting
  
  std::vector <int> toRemove;
  for (int i = 0 ; i< allNodeIds.size() ; i++) {
    if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
      toRemove.push_back(i);
    }
  }
  sort(toRemove.begin(), toRemove.end(), cmp);
  for (int i = 0 ; i< toRemove.size() ; i++) {
    std::vector<int>::iterator vi = allNodeIds.begin();
    allNodeIds.erase(vi+toRemove[i]);
  }
  
  //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
  nodeIdsToRegraft = allNodeIds;
  
}







/************************************************************************
 * Makes a modification, knowing what previous modifications have been done.
 ************************************************************************/
void makeDeterministicModifications(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForSPR, int & nodeForRooting) {
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
      for (int i = 0 ; i< secondHalfNodeIds.size() ; i++ ) {
        allNodeIds.push_back(secondHalfNodeIds[i]);
      }
    }
    else {
      allNodeIds = secondHalfNodeIds;
      for (int i = 0 ; i< firstHalfNodeIds.size() ; i++ ) {
        allNodeIds.push_back(firstHalfNodeIds[i]);
      }
    }
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)));
    
    int oldFatherId = tree.getNode(nodeForSPR)->getFather()->getId();
    int brotherId;
    //Get one brother ; a binary tree is supposed here (because of the "break")
    for(int i=0;i<tree.getNode(oldFatherId)->getNumberOfSons();i++)
      if(tree.getNode(oldFatherId)->getSon(i)->getId()!=nodeForSPR){brotherId=tree.getNode(oldFatherId)->getSon(i)->getId(); break;}
    if ((tree.getNode(oldFatherId)->hasFather())||(!tree.getNode(brotherId)->isLeaf())) {
      forbiddenIds.push_back(brotherId);
    }
    std::vector <int> toRemove;
    for (int i = 0 ; i< allNodeIds.size() ; i++) {
      if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
        toRemove.push_back(i);
      }
    }
    sort(toRemove.begin(), toRemove.end(), cmp);
    for (int i = 0 ; i< toRemove.size() ; i++) {
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
void makeDeterministicNNIsAndRootChangesOnly(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForRooting) {
  if (nodeForNNI < tree.getNumberOfNodes()) {//Make a NNI or rerooting move
    if (nodeForNNI <3) {
      if (nodeForRooting<tree.getNumberOfNodes()) {//Make a rerooting move
        changeRoot(tree, nodeForRooting);
        nodeForRooting++;
      }
      else { //Make a NNI move
        nodeForNNI=3;
        nodeForRooting = 4; //We don't want to root on nodes 1 or 2, 
                            //the two sons of the root.
                            //We do not want to root on node 3 either, 
                            //as a NNI already provides this tree.
        makeNNI(tree, nodeForNNI);
        nodeForNNI++;
      }
    }
    else {
      makeNNI(tree, nodeForNNI);
      nodeForNNI++;
    }
  }
  else {//we reset the loop rooting-NNIs
    nodeForNNI = 0;
    changeRoot(tree, nodeForRooting);
    nodeForRooting++;
  }
}







/************************************************************************
 * Procedure that makes sure that an NNI or rerooting one is about to make has
 * not been computed already.
 * If all NNIs or rerootings have been made, return true.
 ************************************************************************/

bool checkChangeHasNotBeenDone(TreeTemplate<Node> &tree, TreeTemplate<Node> *bestTree, int & nodeForNNI, 
                               int & nodeForRooting, std::vector < double >  &NNILks, 
                               std::vector < double >  &rootLks)
{
  if (nodeForNNI < tree.getNumberOfNodes()) 
    {//Make a NNI or rerooting move
      if (nodeForNNI <3) 
        {
        if (nodeForRooting<tree.getNumberOfNodes()) 
          {//Make a rerooting move
            while ((rootLks[nodeForRooting] < NumConstants::VERY_BIG) && (nodeForRooting < tree.getNumberOfNodes())) {
              //std::cout<<rootLks[nodeForRooting]<<" 1_NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
              
              nodeForRooting++;
            }
          }
        if (nodeForRooting >= tree.getNumberOfNodes()) 
          { //Make a NNI move
            nodeForNNI=3;
            //   while ((NNILks[bestTree->getNode(nodeForNNI-1)->getFather()->getId()] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))
            while ((NNILks[nodeForNNI-1] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))
              {
              std::cout<<NNILks[nodeForNNI-1]<<" NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
              nodeForNNI++;
              }
          }
        }
      else 
        {
        //  while ((NNILks[bestTree->getNode(nodeForNNI-1)->getFather()->getId()] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))       
        while ((NNILks[nodeForNNI-1] < NumConstants::VERY_BIG) && (nodeForNNI < tree.getNumberOfNodes()))
          {
          std::cout<<NNILks[nodeForNNI-1]<<" NumConstants::VERY_BIG: "<<NumConstants::VERY_BIG <<std::endl;
          nodeForNNI++;
          }
        }
    }
  else 
    {//we reset the loop rooting-NNIs
      nodeForNNI = 2;
      while ((rootLks[nodeForRooting] < NumConstants::VERY_BIG) && (nodeForRooting < tree.getNumberOfNodes())) {
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


