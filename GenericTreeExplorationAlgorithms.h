/*
 *  GenericTreeExplorationAlgorithms.h
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 17/12/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */
#ifndef _GENERICTREEEXPLORATIONALGORITHM_H_
#define _GENERICTREEEXPLORATIONALGORITHM_H_

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Node.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>

#include "ReconciliationTools.h"


using namespace bpp;
using namespace std;


void changeRoot(TreeTemplate<Node> &tree, int newOutGroup);
void makeSPR(TreeTemplate<Node> &tree, int cutNodeId, int newBrotherId, bool verbose = true);
void makeNNI(TreeTemplate<Node> &tree, int nodeId);
void buildVectorOfRegraftingNodes(TreeTemplate<Node> &tree, int nodeForSPR, std::vector <int> & nodeIdsToRegraft);
std::vector<int> getRemainingNeighbors(const Node * node1, const Node * node2);
void getRemainingNeighborsUntilDistance(const Node * node1, const Node * node2, int distance, int d, std::vector <int> & neighbors);
void getNeighboringNodesIdLimitedDistance (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds);
void buildVectorOfRegraftingNodesLimitedDistance(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft);
void makeDeterministicModifications(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForSPR, int & nodeForRooting);
void makeDeterministicNNIsAndRootChangesOnly(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForRooting);
bool checkChangeHasNotBeenDone(TreeTemplate<Node> &tree, TreeTemplate<Node> *bestTree, int & nodeForNNI, 
                               int & nodeForRooting, std::vector < double >  &NNILks, 
                               std::vector < double >  &rootLks);
double refineGeneTreeDLOnly (TreeTemplate<Node> * spTree, 
                             TreeTemplate<Node> * geneTree, 
                             std::map<std::string, std::string > seqSp,
                             std::map<std::string, int > spID,
                             std::vector< double> &lossExpectedNumbers, 
                             std::vector < double> &duplicationExpectedNumbers, 
                             int & MLindex, 
                             std::vector <int> &num0lineages, 
                             std::vector <int> &num1lineages, 
                             std::vector <int> &num2lineages, 
                             std::set <int> &nodesToTryInNNISearch);

#endif 