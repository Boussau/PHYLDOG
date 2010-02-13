
/* This file contains various functions useful for the search of the species tree*/

#include "SpeciesTreeExploration.h"



//To sort in ascending order
bool cmp( int a, int b ) {
  return a > b;
}  









/************************************************************************
 * Change the root of the tree by changing the outgroup.
************************************************************************/

void changeRoot(TreeTemplate<Node> &tree, int newOutGroup) {
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
 
  if(!son->hasFather()) throw NodeException("makeNNI(). Node 'son' must not be the root node.", son);
  Node * parent = son->getFather();

  if(!parent->hasFather()) throw NodeException("makeNNI(). Node 'parent' must not be the root node.", parent);
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
 * Tries all SPRs of the subtree starting in node nodeForSPR, and executes the one with the highest likelihood.
 ************************************************************************/
void tryAllPossibleSPRsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree,int nodeForSPR, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing) {

  std::vector <double> backupDupProba=duplicationProbabilities;
  std::vector<double> backupLossProba=lossProbabilities;
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  TreeTemplate<Node> *tree;
  std::vector <int> nodeIdsToRegraft; 
  tree = currentTree->clone();
  buildVectorOfRegraftingNodes(*tree, nodeForSPR, nodeIdsToRegraft);
  bool betterTree = false;
  for (int i =0 ; i<nodeIdsToRegraft.size() ; i++) {
    if (i!=0) {
      deleteTreeProperties(*tree);
      delete tree;
      tree = currentTree->clone();
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i]);
    duplicationProbabilities = backupDupProba;
    lossProbabilities = backupLossProba;
    breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
   
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    //SHOULD WE USE ACCUMULATE ?
    int temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
  
    //SECOND COMPUTATION, TO COMPUTE WITH GOOD RATES OF LOSSES AND DUPLICATION
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(logLs);
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    
    //SHOULD WE USE ACCUMULATE ?
    temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
    
    //END OF SECOND COMPUTATION

   std::cout<<"second computation ended"<< std::endl;

    if (logL+0.01<bestlogL) {
      betterTree = true;
      bestlogL =logL;
      deleteTreeProperties(*bestTree);
      delete bestTree;
      bestTree = tree->clone(); 
      bestDupProba=duplicationProbabilities;
      bestLossProba=lossProbabilities;
      bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;
      bestBranchNumbers = branchNumbers;
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      bestIndex = index;
     std::cout << "SPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
     std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
    }
  }
  if (betterTree) {
    logL = bestlogL;
    duplicationProbabilities = bestDupProba;
    lossProbabilities = bestLossProba;
    numIterationsWithoutImprovement = 0;
    lossNumbers = bestLossNumbers;
    duplicationNumbers = bestDuplicationNumbers;
    branchNumbers = bestBranchNumbers;
    num0Lineages = bestNum0Lineages;
    num1Lineages = bestNum1Lineages;
    num2Lineages = bestNum2Lineages;
    deleteTreeProperties(*currentTree);
    delete currentTree;
    currentTree = bestTree->clone(); 
   std::cout <<"SPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
   std::cout << "\t\tServer: SPRs: new total Likelihood value "<<logL<< std::endl;
   std::cout <<"Total Losses: "<<VectorTools::sum(lossNumbers) <<"; Total Duplications: "<<VectorTools::sum(duplicationNumbers)<< std::endl; 
    //Send the new std::vectors to compute the new likelihood of the best tree
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    breadthFirstreNumber (*currentTree);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"bestIndex again "<<bestIndex<< std::endl;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
  }
  else {
    numIterationsWithoutImprovement++;
   std::cout <<"SPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
  }
  deleteTreeProperties(*tree);
  delete tree;
}
 
/************************************************************************
 * Tries all reRootings of the species tree, and executes the one with the highest likelihood.
 ************************************************************************/
void tryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing) {

  std::vector <double> backupDupProba=duplicationProbabilities;
  std::vector<double> backupLossProba=lossProbabilities;
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  TreeTemplate<Node> *tree;
  tree = currentTree->clone();
  bool betterTree = false;
  std::vector <int> nodeIds = tree->getNodesId();
  for (int i =0 ; i<nodeIds.size() ; i++) {
    if (i!=0) {
      deleteTreeProperties(*tree);
      delete tree;
      tree = currentTree->clone();
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    changeRoot(*tree, nodeIds[i]);
    duplicationProbabilities = backupDupProba;
    lossProbabilities = backupLossProba;
    breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
   
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    //SHOULD WE USE ACCUMULATE ?
    int temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
    //SECOND COMPUTATION, TO COMPUTE WITH GOOD RATES OF LOSSES AND DUPLICATION
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server);
    breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(logLs);
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    
    //SHOULD WE USE ACCUMULATE ?
    temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
    

    //END OF SECOND COMPUTATION
  
    if (logL+0.01<bestlogL) { 
      numIterationsWithoutImprovement = 0;
      betterTree = true;
      bestlogL =logL;
      if (bestTree) {
        delete bestTree;
      }
      bestTree = tree->clone(); 
      bestDupProba=duplicationProbabilities;
      bestLossProba=lossProbabilities;
      bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;
      bestBranchNumbers = branchNumbers; 
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      bestIndex = index;
     std::cout <<"ReRooting: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
     std::cout << "Better candidate tree likelihood : "<<bestlogL<< std::endl; 
     std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
    }
    else {
      numIterationsWithoutImprovement++;
     std::cout <<"ReRooting: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
    }
  }
  if (betterTree) {
    logL = bestlogL; 
    duplicationProbabilities = bestDupProba;
    lossProbabilities = bestLossProba;
    lossNumbers = bestLossNumbers;
    duplicationNumbers = bestDuplicationNumbers;
    branchNumbers = bestBranchNumbers;
    num0Lineages = bestNum0Lineages;
    num1Lineages = bestNum1Lineages;
    num2Lineages = bestNum2Lineages;
    delete currentTree;
    currentTree = bestTree->clone(); 
   std::cout << "\t\tServer: tryAllPossibleReRootingsAndMakeBestOne: new total Likelihood value "<<logL<< std::endl;
   std::cout << TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;

   
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    breadthFirstreNumber (*currentTree);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
  

  }
  deleteTreeProperties(*tree);
  delete tree;
 std::cout<<"tree deleted !!!"<< std::endl;
  
}
 

/************************************************************************
 * Tries all NNIs, and executes NNIs as soon as they increase the likelihood.
 ************************************************************************/
void tryAndMakeNNIs(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages,  std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing) {

  std::vector <double> backupDupProba=duplicationProbabilities;
  std::vector<double> backupLossProba=lossProbabilities;
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  TreeTemplate<Node> *tree;
  std::vector <int> nodeIds; 
  tree = currentTree->clone();
  nodeIds = currentTree->getNodesId();
  bool betterTree = false;
  for (int i =3 ; i<nodeIds.size() ; i++) {
    if (i!=0) {
      deleteTreeProperties(*tree);
      delete tree;
      tree = currentTree->clone();
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    makeNNI(*tree, i);
    duplicationProbabilities = backupDupProba;
    lossProbabilities = backupLossProba;
    breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
  
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    //SHOULD WE USE ACCUMULATE ?
    int temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
    
    //SECOND COMPUTATION, TO COMPUTE WITH GOOD RATES OF LOSSES AND DUPLICATION
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(logLs);
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    
    //SHOULD WE USE ACCUMULATE ?
    temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[i];
      num1Lineages= num1Lineages+allNum1Lineages[i];
      num2Lineages= num2Lineages+allNum2Lineages[i];
    } 
    

    //END OF SECOND COMPUTATION
    
    if (logL+0.01<bestlogL) {
      betterTree = true;    
      numIterationsWithoutImprovement = 0;
      bestlogL =logL;
      deleteTreeProperties(*bestTree);
      delete bestTree;
      bestTree = tree->clone();  
      bestDupProba=duplicationProbabilities;
      bestLossProba=lossProbabilities;
      bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;
      bestBranchNumbers = branchNumbers;
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      bestIndex = index;
     std::cout <<"tryAndMakeNNIs: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
     std::cout << "tryAndMakeNNIs: Better tree likelihood : "<<bestlogL<< std::endl; 
     std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
    }
    else {
      numIterationsWithoutImprovement++;
     std::cout <<"NNI: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
    }
  }
  if (betterTree) {
    logL = bestlogL; 
    duplicationProbabilities = bestDupProba;
    lossProbabilities = bestLossProba;
    lossNumbers = bestLossNumbers;
    duplicationNumbers = bestDuplicationNumbers;
    branchNumbers = bestBranchNumbers;
    num0Lineages = bestNum0Lineages;
    num1Lineages = bestNum1Lineages;
    num2Lineages = bestNum2Lineages;
    deleteTreeProperties(*currentTree);
    delete currentTree;
    currentTree = bestTree->clone(); 
   std::cout << "\t\tServer: tryAndMakeNNIs: new total Likelihood value "<<logL<< std::endl;
   std::cout << TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;
    
    if ((branchProbaOptimization=="average")||(branchProbaOptimization=="average_then_branchwise")) {
      std::string temp = "average";
      computeDuplicationAndLossRatesForTheSpeciesTree (temp, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    else {
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *currentTree);
    }
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
    breadthFirstreNumber (*currentTree);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
  }
  deleteTreeProperties(*tree);
  delete tree;
  
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
 * Makes a modification, knowing what previous modifications have been done. Everything is deterministic, even SPRs.
 * For the moment we allow SPRs at a distance of 2 or 3 branches. There are 4 possible branching points at a distance of 2 branches, and 8 at a distance of 3 branches, which amounts to 12 branching points to try.

************************************************************************/
void makeTotallyDeterministicModifications(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForSPR, int & nodeForRooting) {
  int numberOfBranchingPoints = 12;
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

     int reste = nodeForSPR % numberOfBranchingPoints;
     int prunedSubtree = nodeForSPR / numberOfBranchingPoints;



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
 * Procedure that makes NNIs and ReRootings by calling makeDeterministicNNIsAndRootChangesOnly.

************************************************************************/


void localOptimizationWithNNIsAndReRootings(const mpi::communicator& world, 
                                            TreeTemplate<Node> *tree, 
                                            TreeTemplate<Node> *bestTree, 
                                            int &index, 
                                            int &bestIndex,  
                                            bool &stop, 
                                            double &logL, 
                                            double &bestlogL, 
                                            std::vector<int> &lossNumbers, 
                                            std::vector<int> &duplicationNumbers, 
                                            std::vector<int> &branchNumbers, 
                                            std::vector<int> &bestLossNumbers, 
                                            std::vector<int> &bestDuplicationNumbers, 
                                            std::vector<int> &bestBranchNumbers,
                                            std::vector< std::vector<int> > AllLosses, 
                                            std::vector< std::vector<int> > AllDuplications, 
                                            std::vector< std::vector<int> > AllBranches,
                                            std::vector<int> &num0Lineages, 
                                            std::vector<int> &num1Lineages,
                                            std::vector<int> &num2Lineages, 
                                            std::vector<int> &bestNum0Lineages, 
                                            std::vector<int> &bestNum1Lineages,
                                            std::vector<int> &bestNum2Lineages, 
                                            std::vector< std::vector<int> > &allNum0Lineages, 
                                            std::vector< std::vector<int> > &allNum1Lineages,
                                            std::vector< std::vector<int> > &allNum2Lineages, 
                                            std::vector<double> &lossProbabilities, 
                                            std::vector<double> &duplicationProbabilities, 
                                            double averageDuplicationProbability, 
                                            double averageLossProbability, 
                                            bool rearrange, 
                                            int &numIterationsWithoutImprovement, 
                                            int server, 
                                            int & nodeForNNI, 
                                            int & nodeForRooting, 
                                            std::string & branchProbaOptimization, 
                                            std::map < std::string, int> genomeMissing,
                                            std::vector < double >  &NNILks) {
 
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  std::string currentSpeciesTree;
  
  while (!stop) {
    makeDeterministicNNIsAndRootChangesOnly(*tree, nodeForNNI, nodeForRooting);
   std::cout << "Proposed species tree: "<<TreeTools::treeToParenthesis(*tree, true)<< std::endl;
    computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *tree, bestlogL);
    if ((nodeForNNI >=3) && (nodeForRooting == 4)) {
      int branchId = bestTree->getNode(nodeForNNI-1)->getFather()->getId();
      if (logL < NNILks[branchId]) {
        NNILks[branchId] = logL;
      }
     
    }
    if (logL+0.01<bestlogL) {
     std::cout << "\t\tNNIs or Root changes: Improvement: new total Likelihood value "<<logL<<" compared to the best log Likelihood : "<<bestlogL<< std::endl;
      numIterationsWithoutImprovement = 0;
      bestlogL =logL;
      if (bestTree){
        deleteTreeProperties(*bestTree);
        delete bestTree;
      }
      bestTree = tree->clone();  
      bestDupProba=duplicationProbabilities;
      bestLossProba=lossProbabilities;
      bestIndex = index;
      bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;	 
      bestBranchNumbers = branchNumbers;
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      bestIndex = index;
      for (int i = 0 ; i< NNILks.size() ; i++ ) {
        NNILks[i]=NumConstants::VERY_BIG;
      }
    }
    else {
      numIterationsWithoutImprovement++;
     std::cout <<"\t\tNNIs or Root changes: Number of iterations without improvement: "<<numIterationsWithoutImprovement<< std::endl;
      deleteTreeProperties(*tree);
      delete tree;
      tree = bestTree->clone();
      duplicationProbabilities = bestDupProba;
      lossProbabilities = bestLossProba;
      if (numIterationsWithoutImprovement>2*tree->getNumberOfNodes()) {
        stop = true;
        broadcast(world, stop, server); 
        broadcast(world, bestIndex, server);
      }
    }
  }
}







/************************************************************************
 * Only optimizes duplication and loss rates
************************************************************************/
void optimizeOnlyDuplicationAndLossRates(const mpi::communicator& world, TreeTemplate<Node> *tree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool &stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, int & nodeForNNI, int & nodeForRooting, std::string & branchProbaOptimization, std::map < std::string, int> genomeMissing) {
  std::string currentSpeciesTree;
  while (!stop) {
    std::vector<double> logLs;
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
    broadcast(world, currentSpeciesTree, server);
    //Computation in clients
    index++;  
   std::cout <<"\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    int temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers+AllDuplications[k];
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[k];
      num1Lineages= num1Lineages+allNum1Lineages[k];
      num2Lineages= num2Lineages+allNum2Lineages[k];
    } 
    
    if (logL+0.01<bestlogL) {
     std::cout << "\t\tDuplication and loss probabilities optimization: Server : new total Likelihood value "<<logL<<" compared to the best log Likelihood : "<<bestlogL<< std::endl;
     std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
      numIterationsWithoutImprovement = 0;
      bestlogL =logL;
      if (bestTree) 	
	{
	  deleteTreeProperties(*bestTree);
	  delete bestTree;
	}
      bestTree = tree->clone();
      bestIndex = index;
      bestLossNumbers = lossNumbers;
      bestDuplicationNumbers = duplicationNumbers;	 
      bestBranchNumbers = branchNumbers;
      bestNum0Lineages = num0Lineages;
      bestNum1Lineages = num1Lineages;
      bestNum2Lineages = num2Lineages;
      
      computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, *tree);
  
      broadcast(world, stop, server);
      broadcast(world, rearrange, server); 
      std::string currentSpeciesTree = TreeTools::treeToParenthesis(*tree, true);
      breadthFirstreNumber (*tree, duplicationProbabilities, lossProbabilities);
      broadcast(world, lossProbabilities, server); 
      broadcast(world, duplicationProbabilities, server); 
      broadcast(world, currentSpeciesTree, server);
      //COMPUTATION IN CLIENTS
      index++;  
      bestIndex = index;
     std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
      logL = 0.0;
      std::vector<double> logLs;
      resetVector(duplicationNumbers);
      resetVector(lossNumbers);
      resetVector(branchNumbers);
      resetVector(num0Lineages);
      resetVector(num1Lineages);
      resetVector(num2Lineages);
      gather(world, logL, logLs, server);
      logL = VectorTools::sum(logLs);
      gather(world, duplicationNumbers, AllDuplications, server); 
      gather(world, lossNumbers, AllLosses, server);
      gather(world, branchNumbers, AllBranches, server);
      gather(world, num0Lineages, allNum0Lineages, server);
      gather(world, num1Lineages, allNum1Lineages, server);
      gather(world, num2Lineages, allNum2Lineages, server);
      temp = allNum0Lineages.size();
      for (int i =0; i<temp ; i++ ) {
	num0Lineages= num0Lineages+allNum0Lineages[i];
	num1Lineages= num1Lineages+allNum1Lineages[i];
	num2Lineages= num2Lineages+allNum2Lineages[i];
      }
    }
    else {
      numIterationsWithoutImprovement++;
     std::cout <<"\t\tNo more increase in likelihood "<< std::endl;
      deleteTreeProperties(*tree);
      delete tree;
      tree = bestTree->clone();
      stop = true;
      broadcast(world, stop, server); 
      broadcast(world, bestIndex, server);
    }
  }
}



/************************************************************************
 * Tries all reRootings of the species tree, and executes the one with the highest likelihood. Here we do not use branchwise rates of duplications and losses, and we do not optimise these rates often. However, there are lots of useless identical std::vector copies...
 ************************************************************************/
void fastTryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, bool optimizeRates) {  
  breadthFirstreNumber (*currentTree, duplicationProbabilities, lossProbabilities);
  std::vector <double> backupDupProba=duplicationProbabilities;
  std::vector<double> backupLossProba=lossProbabilities;
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  TreeTemplate<Node> *tree;
  tree = currentTree->clone();
  bool betterTree = false;
  std::vector <int> nodeIds = tree->getNodesId();
  for (int i =0 ; i<nodeIds.size() ; i++) {
    if ((nodeIds[i]!=0)&&(nodeIds[i]!=1)&&(nodeIds[i]!=2)) { //We do not want to try the root we're alredy at
      if (i!=0) {
        deleteTreeProperties(*tree);
        delete tree;
        tree = currentTree->clone();
      }
      changeRoot(*tree, nodeIds[i]);
      if (optimizeRates) 
        {
          computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *tree, bestlogL);
        }
      else 
        {
          computeSpeciesTreeLikelihood(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, backupDupProba, backupLossProba, *tree);
      }
      if (logL+0.01<bestlogL) { 
        numIterationsWithoutImprovement = 0;
        betterTree = true;
        bestlogL =logL;
        if (bestTree) {
          delete bestTree;
        }
        bestTree = tree->clone();   
        bestDupProba=duplicationProbabilities;
        bestLossProba=lossProbabilities;
        bestLossNumbers = lossNumbers;
        bestDuplicationNumbers = duplicationNumbers;
        bestBranchNumbers = branchNumbers; 
        bestNum0Lineages = num0Lineages;
        bestNum1Lineages = num1Lineages;
        bestNum2Lineages = num2Lineages;
        bestIndex = index;
       std::cout <<"ReRooting: Improvement! : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
       std::cout << "Better candidate tree likelihood : "<<bestlogL<< std::endl; 
       std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
      }
      else {
        numIterationsWithoutImprovement++;
       std::cout <<"ReRooting: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< " logLk: "<<logL<< std::endl;
      }
    }
  }
  
   
  
  
  if (betterTree) {
    logL = bestlogL;  
    duplicationProbabilities = bestDupProba;
    lossProbabilities = bestLossProba;
    lossNumbers = bestLossNumbers;
    duplicationNumbers = bestDuplicationNumbers;
    branchNumbers = bestBranchNumbers;
    num0Lineages = bestNum0Lineages;
    num1Lineages = bestNum1Lineages;
    num2Lineages = bestNum2Lineages;
    delete currentTree;
    currentTree = bestTree->clone(); 
   std::cout << "\t\tServer: tryAllPossibleReRootingsAndMakeBestOne: new total Likelihood value "<<logL<< std::endl;
   std::cout << TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    std::string currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
    breadthFirstreNumber (*currentTree, duplicationProbabilities, lossProbabilities);
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;    
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
        
  }
  else {
   std::cout<< "No improvement in fastTryAllPossibleReRootingsAndMakeBestOne"<< std::endl;
    logL = bestlogL;  
    duplicationProbabilities = bestDupProba;
    lossProbabilities = bestLossProba;
    lossNumbers = bestLossNumbers;
    duplicationNumbers = bestDuplicationNumbers;
    branchNumbers = bestBranchNumbers;
    num0Lineages = bestNum0Lineages;
    num1Lineages = bestNum1Lineages;
    num2Lineages = bestNum2Lineages;
    delete currentTree;
    currentTree = bestTree->clone(); 
  }
  
  deleteTreeProperties(*tree);
  delete tree;  
}


/************************************************************************
 * Tries all SPRs at a distance dist for all possible subtrees of the subtree starting in node nodeForSPR, and executes the ones with the highest likelihood. Uses only average rates of duplication and loss, not branchwise rates.
 ************************************************************************/
void fastTryAllPossibleSPRs(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, int sprLimit, bool optimizeRates) {
 std::cout <<"in fastTryAllPossibleSPRs : currentTree : "<< std::endl;
 std::cout<< TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;
  breadthFirstreNumber (*currentTree, duplicationProbabilities, lossProbabilities);
  std::vector <double> backupDupProba=duplicationProbabilities;
  std::vector<double> backupLossProba=lossProbabilities;
  std::vector <double> bestDupProba=duplicationProbabilities;
  std::vector<double> bestLossProba=lossProbabilities;
  for (int nodeForSPR=currentTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) {
    TreeTemplate<Node> *tree;
    std::vector <int> nodeIdsToRegraft; 
    tree = currentTree->clone();
    buildVectorOfRegraftingNodesLimitedDistance(*tree, nodeForSPR, sprLimit, nodeIdsToRegraft);
    bool betterTree = false;
    for (int i =0 ; i<nodeIdsToRegraft.size() ; i++) {
      if (i!=0) {
        delete tree;
        tree = currentTree->clone();
      }
      makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i]);
      if (optimizeRates) 
        {
          computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, *tree, bestlogL);
        }
      else 
        {
          computeSpeciesTreeLikelihood(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, backupDupProba, backupLossProba, *tree);
        }
      if (logL+0.01<bestlogL) {
        betterTree = true;
        bestlogL =logL;
        deleteTreeProperties(*bestTree);
        delete bestTree;
        bestTree = tree->clone();  
        bestDupProba=duplicationProbabilities;
        bestLossProba=lossProbabilities;
        bestLossNumbers = lossNumbers;
        bestDuplicationNumbers = duplicationNumbers;
        bestBranchNumbers = branchNumbers;
        bestNum0Lineages = num0Lineages;
        bestNum1Lineages = num1Lineages;
        bestNum2Lineages = num2Lineages;
        bestIndex = index;
       std::cout << "SPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
       std::cout << TreeTools::treeToParenthesis(*tree, true)<< std::endl;
      }
    }
    if (betterTree) {
      logL = bestlogL; 
      duplicationProbabilities = bestDupProba;
      lossProbabilities = bestLossProba;
      backupDupProba=bestDupProba;
      backupLossProba=bestLossProba;
      numIterationsWithoutImprovement = 0;
      lossNumbers = bestLossNumbers;
      duplicationNumbers = bestDuplicationNumbers;
      branchNumbers = bestBranchNumbers;
      num0Lineages = bestNum0Lineages;
      num1Lineages = bestNum1Lineages;
      num2Lineages = bestNum2Lineages;
      deleteTreeProperties(*currentTree);
      delete currentTree;
      currentTree = bestTree->clone(); 
      breadthFirstreNumber (*currentTree);
     std::cout <<"SPRs: Improvement! : "<<numIterationsWithoutImprovement<< std::endl;
     std::cout << "\t\tServer: SPRs: new total Likelihood value "<<logL<< std::endl;
      //Send the new std::vectors to compute the new likelihood of the best tree
      broadcast(world, stop, server);
      broadcast(world, rearrange, server); 
      broadcast(world, lossProbabilities, server); 
      broadcast(world, duplicationProbabilities, server); 
      std::string currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, true);
      broadcast(world, currentSpeciesTree, server);
      //COMPUTATION IN CLIENTS
      index++;  
      bestIndex = index;
     std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
      logL = 0.0;
      std::vector<double> logLs;
      resetVector(duplicationNumbers);
      resetVector(lossNumbers);
      resetVector(branchNumbers);
      resetVector(num0Lineages);
      resetVector(num1Lineages);
      resetVector(num2Lineages);
      gather(world, logL, logLs, server);
      logL = VectorTools::sum(logLs);
      gather(world, duplicationNumbers, AllDuplications, server); 
      gather(world, lossNumbers, AllLosses, server);
      gather(world, branchNumbers, AllBranches, server);
      gather(world, num0Lineages, allNum0Lineages, server);
      gather(world, num1Lineages, allNum1Lineages, server);
      gather(world, num2Lineages, allNum2Lineages, server);
    }
    else {
      logL = bestlogL;  
      duplicationProbabilities = bestDupProba;
      lossProbabilities = bestLossProba;
      lossNumbers = bestLossNumbers;
      duplicationNumbers = bestDuplicationNumbers;
      branchNumbers = bestBranchNumbers;
      num0Lineages = bestNum0Lineages;
      num1Lineages = bestNum1Lineages;
      num2Lineages = bestNum2Lineages;
      delete currentTree;
      currentTree = bestTree->clone(); 
      numIterationsWithoutImprovement++;
     std::cout <<"SPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement<< std::endl;
    }
    deleteTreeProperties(*tree);
    delete tree;
  }
}
 


/************************************************************************
 * Tries all SPRs and all rerootings with average rates of duplication and loss, 
 * not branchwise rates. Only does SPRs at a given distance. 
 ************************************************************************/
void fastTryAllPossibleSPRsAndReRootings(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, int sprLimit, bool optimizeRates) {
  if (optimizeRates)
    {
     std::cout <<"Making SPRs and NNIs and optimizing duplication and loss rates."<< std::endl;
    }
  else 
    {
     std::cout <<"Making SPRs and NNIs but NOT optimizing duplication and loss rates."<< std::endl;
    }
  numIterationsWithoutImprovement=0;
  while (numIterationsWithoutImprovement<2*currentTree->getNumberOfNodes()) {
   std::cout <<"Dup Rates";
    VectorTools::print(duplicationProbabilities);
   std::cout <<"Loss Rates";
    VectorTools::print(lossProbabilities);
  
    fastTryAllPossibleSPRs(world, currentTree, bestTree, index, bestIndex,  stop, logL, bestlogL, lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, branchProbaOptimization, genomeMissing, sprLimit, optimizeRates);
   std::cout <<"Dup Rates";
    VectorTools::print(duplicationProbabilities);
   std::cout <<"Loss Rates";
    VectorTools::print(lossProbabilities);
       if (numIterationsWithoutImprovement>2*currentTree->getNumberOfNodes()) {	
      break;
    }
    
    
    
    fastTryAllPossibleReRootingsAndMakeBestOne(world, currentTree, bestTree, index, bestIndex,  stop, logL, bestlogL, lossNumbers, duplicationNumbers, branchNumbers, bestLossNumbers, bestDuplicationNumbers, bestBranchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages, num2Lineages, bestNum0Lineages, bestNum1Lineages, bestNum2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, averageDuplicationProbability, averageLossProbability, rearrange, numIterationsWithoutImprovement, server, branchProbaOptimization, genomeMissing, optimizeRates);
   std::cout << "after fastTryAllPossibleReRootingsAndMakeBestOne :currentTree: "<< std::endl;
	 std::cout << TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;
    
    
 }
 std::cout << "\n\n\t\t\tFirst fast step of SPRs and rerootings over.\n\n"<< std::endl;
  
}



/************************************************************************
 * Computes the likelihood of a species tree, only once. 
 ************************************************************************/


std::string computeSpeciesTreeLikelihood(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, std::vector <double> backupDupProba, std::vector <double> backupLossProba, TreeTemplate<Node> &tree) {
  broadcast(world, stop, server);
  broadcast(world, rearrange, server); 
  duplicationProbabilities = backupDupProba;
  lossProbabilities = backupLossProba;
  breadthFirstreNumber (tree, duplicationProbabilities, lossProbabilities);
  broadcast(world, lossProbabilities, server); 
  broadcast(world, duplicationProbabilities, server); 
  std::string currentSpeciesTree = TreeTools::treeToParenthesis(tree, true);
  broadcast(world, currentSpeciesTree, server);
  //COMPUTATION IN CLIENTS
  index++;  
 std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
  logL = 0.0;
  std::vector<double> logLs;
  resetVector(duplicationNumbers);
  resetVector(lossNumbers);
  resetVector(branchNumbers);
  resetVector(num0Lineages);
  resetVector(num1Lineages);
  resetVector(num2Lineages);
  gather(world, logL, logLs, server);
  logL = VectorTools::sum(logLs);
  gather(world, duplicationNumbers, AllDuplications, server); 
  gather(world, lossNumbers, AllLosses, server);
  gather(world, branchNumbers, AllBranches, server);
  gather(world, num0Lineages, allNum0Lineages, server);
  gather(world, num1Lineages, allNum1Lineages, server);
  gather(world, num2Lineages, allNum2Lineages, server);
  int temp = AllDuplications.size();
  for (int k =0; k<temp ; k++ ) {
    duplicationNumbers= duplicationNumbers +AllDuplications[k];	
    lossNumbers= lossNumbers+AllLosses[k];
    branchNumbers= branchNumbers+AllBranches[k];
    num0Lineages= num0Lineages+allNum0Lineages[k];
    num1Lineages= num1Lineages+allNum1Lineages[k];
    num2Lineages= num2Lineages+allNum2Lineages[k];
  }  
  return currentSpeciesTree;
}



/************************************************************************
 * Computes the likelihood of a species tree, optimizes the duplication and loss
 * rates, and updates the likelihood of the species tree. This way the likelihood
 * of this species tree is computed with adequate rates. 
 ************************************************************************/


void computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL) {
  //First we compute the likelihood with old duplication and loss rates
  computeDuplicationAndLossRatesForTheSpeciesTreeInitially (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, tree);
  std::string currentSpeciesTree = computeSpeciesTreeLikelihood(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages,num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, lossProbabilities, duplicationProbabilities, rearrange, server, branchProbaOptimization, genomeMissing, duplicationProbabilities, lossProbabilities, tree);
  double currentlogL = -UNLIKELY;
  int i=1;
  //Then we update duplication and loss rates based on the results of this first 
  //computation, until the likelihood stabilizes (roughly)
  while ((i<=1)&&(currentlogL-logL>logL-bestlogL)) { 
    currentlogL = logL;
    i++;
    computeDuplicationAndLossRatesForTheSpeciesTree (branchProbaOptimization, num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities, genomeMissing, tree);
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    std::vector<double> logLs;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
   std::cout<<std::setprecision(15);
   std::cout <<"LogL "<<logL<< std::endl;
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    //SHOULD WE USE ACCUMULATE ?
    int temp = AllDuplications.size();
    for (int k =0; k<temp ; k++ ) {
      duplicationNumbers= duplicationNumbers +AllDuplications[k];	
      lossNumbers= lossNumbers+AllLosses[k];
      branchNumbers= branchNumbers+AllBranches[k];
      num0Lineages= num0Lineages+allNum0Lineages[k];
      num1Lineages= num1Lineages+allNum1Lineages[k];
      num2Lineages= num2Lineages+allNum2Lineages[k];
    }  
    
  }
 std::cout << i<< " iterations of likelihood computation for optimizing duplication and loss rates have been done."<< std::endl;
}



/************************************************************************
 * Optimizes the duplication and loss rates with a numerical algorithm.
 * As we should start from overestimated values, we start by decreasing all 
 * values, and keep on trying to minimize these rates until the likelihood stops 
 * increasing. 
 ************************************************************************/
void numericalOptimizationOfDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL) {
  std::cout<<"Numerical optimization of Duplication And Loss expected numbers of events"<<std::endl;
  std::cout<<"Before optimization: Duplications"<<std::endl;
  VectorTools::print(duplicationProbabilities);  
  std::cout<<"Before optimization: Losses"<<std::endl;
  VectorTools::print(lossProbabilities);
  
  std::vector <double> formerDupRates = duplicationProbabilities;
  std::vector <double> formerLossRates = lossProbabilities;
  
  std::vector <double> newDupRates = duplicationProbabilities;
  std::vector <double> newLossRates = lossProbabilities;
  double currentlogL = -UNLIKELY;

  while (currentlogL > logL + 0.01) {
    currentlogL = logL;
    formerDupRates = newDupRates;
    formerLossRates = newLossRates;
    newDupRates = newDupRates * 0.5;
    newLossRates = newLossRates * 0.5;
    computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(world, index, stop, logL, lossNumbers, duplicationNumbers, branchNumbers, AllLosses, AllDuplications, AllBranches, num0Lineages, num1Lineages, num2Lineages, allNum0Lineages, allNum1Lineages, allNum2Lineages, newLossRates, newDupRates, rearrange, server, branchProbaOptimization, genomeMissing, tree, bestlogL);
  }
  
  //If we have improved the log-likelihood
  if (currentlogL < bestlogL) {
    duplicationProbabilities = formerDupRates;
    lossProbabilities = formerLossRates;
    bestlogL = currentlogL;
    std::cout<<"Better logLikelihood! logL: "<< currentlogL <<std::endl;
    std::cout<<"After optimization: Duplications"<<std::endl;
    VectorTools::print(duplicationProbabilities);  
    std::cout<<"After optimization: Losses"<<std::endl;
    VectorTools::print(lossProbabilities);
  }
  else {
    std::cout<<"No improvement in logLikelihood."<<std::endl;
  }
  
  return;
}














/**********************************CRAP******************************************/




/************************************************************************
 * Randomly chooses the modification to be done, either SPR or changeRoot. The probabilities of the two events are a function of expectedFrequencySPR and expectedFrequencyChangeRoot.
************************************************************************/
/*
void makeRandomModifications(TreeTemplate<Node> &tree, int expectedFrequencyNNI, int expectedFrequencySPR, int expectedFrequencyChangeRoot) {
  //breadthFirstreNumber (tree);
  int tot = expectedFrequencyNNI + expectedFrequencySPR + expectedFrequencyChangeRoot;
  int ran=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(tot);
  std::vector <int> allNodeIds;
  //TEMPORARY
  //ran = 0;
  if (ran < expectedFrequencyNNI) {//Make a NNI move
    //we benefit from the fact that the tree has been numbered with the lowest indices as closest to the root : we discard 0,1,2    
    for (int i = 3; i<tree.getNumberOfNodes(); i++) {
      allNodeIds.push_back(i);
    }
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size())+3;
    makeNNI(tree, ran);
  }
  else if (ran < expectedFrequencySPR ) { //Make a SPR move ! we do not want to get the root as one of the two nodes involved
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
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    int cutNodeId = allNodeIds[ran];
    //TEMPORARY
    //  cutNodeId = 2;
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(cutNodeId)));
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
    //TEMPORARY
    // ran2 = 0;
    int newBrotherId = allNodeIds[ran2];
    makeSPR(tree, cutNodeId, newBrotherId);
    }
  else { //Make a ChangeRoot move !
   std::cout << "CHANGE ROOT MOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< std::endl;
    allNodeIds = tree.getNodesId();
    ran = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allNodeIds.size());
    changeRoot(tree, allNodeIds[ran]);
  }
}




















  /*
    std::vector <int> firstHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(0)));
    std::vector <int> secondHalfNodeIds = TreeTemplateTools::getNodesId(*(N->getSon(1)));
    
    if (secondHalfNodeIds.size()<firstHalfNodeIds.size()) {
    allNodeIds = firstHalfNodeIds;
    allNodeIds.insert ( secondHalfNodeIds.end(), firstHalfNodeIds.begin(), firstHalfNodeIds.end() );
    //  allNodeIds = firstHalfNodeIds;
   // for (int i = 0 ; i< secondHalfNodeIds.size() ; i++ ) {
   // allNodeIds.push_back(secondHalfNodeIds[i]);
   // }
  
   }
   else { 
   allNodeIds = secondHalfNodeIds;
   allNodeIds.insert ( firstHalfNodeIds.end(), secondHalfNodeIds.end(),secondHalfNodeIds.end()); 
   //      allNodeIds = secondHalfNodeIds;
   //   for (int i = 0 ; i< firstHalfNodeIds.size() ; i++ ) {
   //   allNodeIds.push_back(firstHalfNodeIds[i]);
   //   }
    
   }
  */




 //TEST
    /*  std::cout << TreeTools::treeToParenthesis(*currentTree, true)<< std::endl;

    //Send the new std::vectors to compute the new likelihood of the best tree
    computeAverageDuplicationAndLossProbabilitiesForAllBranches (num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities);
    broadcast(world, stop, server);
    broadcast(world, rearrange, server); 
    broadcast(world, lossProbabilities, server); 
    broadcast(world, duplicationProbabilities, server); 
    currentSpeciesTree = TreeTools::treeToParenthesis(*currentTree, false);
    broadcast(world, currentSpeciesTree, server);
    //COMPUTATION IN CLIENTS
    index++;  
    bestIndex = index;
   std::cout <<"\t\tNumber of species trees tried : "<<index<< std::endl;
    logL = 0.0;
    resetVector(duplicationNumbers);
    resetVector(lossNumbers);
    resetVector(branchNumbers);
    resetVector(num0Lineages);
    resetVector(num1Lineages);
    resetVector(num2Lineages);
    gather(world, logL, logLs, server);
    logL = VectorTools::sum(logLs);
    gather(world, duplicationNumbers, AllDuplications, server); 
    gather(world, lossNumbers, AllLosses, server);
    gather(world, branchNumbers, AllBranches, server);
    gather(world, num0Lineages, allNum0Lineages, server);
    gather(world, num1Lineages, allNum1Lineages, server);
    gather(world, num2Lineages, allNum2Lineages, server);
    */
    //END TEST


