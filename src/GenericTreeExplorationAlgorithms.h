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
//#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>


//From Seqlib:
#include <Bpp/Seq/Alphabet/DNA.h>


#include "ReconciliationTools.h"


using namespace bpp;
using namespace std;


void changeRoot(TreeTemplate<Node> &tree, int newOutGroup);
std::vector<Node*>  makeSPR(TreeTemplate<Node> &tree, 
                            int cutNodeId, int newBrotherId, 
                            bool verbose = true, 
                            bool returnNodesToUpdate = false);
void makeMuffatoSPR(TreeTemplate<Node> &tree, 
             Node* cutNode, 
             Node* newFather, 
             bool verbose = true) ;

void makeNNI(TreeTemplate<Node> &tree, int nodeId);
void buildVectorOfRegraftingNodes(TreeTemplate<Node> &tree, int nodeForSPR, std::vector <int> & nodeIdsToRegraft);
std::vector<int> getRemainingNeighbors(const Node * node1, const Node * node2);
void getRemainingNeighborsUntilDistance(const Node * node1, const Node * node2, int distance, int d, std::vector <int> & neighbors);
void getNeighboringNodesIdLimitedDistance (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds);
//void getRemainingNeighborsUntilDistanceLowerNodes(TreeTemplate<Node> &tree, const Node * node1, const Node * node2, int distance, int d, std::vector <int> & neighbors);
void getNeighboringNodesIdLimitedDistanceLowerNodes (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds);
void buildVectorOfRegraftingNodesLimitedDistance(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft);
void buildVectorOfRegraftingNodesLimitedDistanceLowerNodes(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft);
void makeDeterministicModifications(TreeTemplate<Node> &tree, size_t & nodeForNNI, size_t & nodeForSPR, size_t & nodeForRooting);
void makeDeterministicNNIsAndRootChangesOnly(TreeTemplate<Node> &tree, size_t & nodeForNNI, size_t & nodeForRooting, const bool fixedOutgroupSpecies_) ;
bool checkChangeHasNotBeenDone(TreeTemplate<Node> &tree, TreeTemplate<Node> *bestTree, size_t & nodeForNNI, 
                               size_t & nodeForRooting, std::vector < double >  &NNILks, 
                               std::vector < double >  &rootLks);
void dropLeaves(TreeTemplate<Node> & tree, const std::vector<string> &spToDrop);
Tree* MRP(const vector<Tree*>& vecTr);
void rootTreeWithOutgroup (TreeTemplate<Node> &tree, 
						   const std::vector<std::string> outgroupTaxa) throw ( TreeException );
bool isTreeRootedWithOutgroup (const TreeTemplate<Node> &tree, const std::vector<std::string> outgroupTaxa) ;
#endif 