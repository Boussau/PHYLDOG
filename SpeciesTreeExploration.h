//
// File: SpeciesTreeExploration.h
// Created by: Bastien Boussau
// Created on: Tu March 26 2008
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _SPECIESTREEEXPLORATION_H_
#define _SPECIESTREEEXPLORATION_H_
// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/Node.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/RandomTools.h>
#include <NumCalc/NumConstants.h>


//#include <algorithm>

//From the BOOST library !!
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;

#include "ReconciliationTools.h"

using namespace bpp;

void changeRoot(TreeTemplate<Node> &tree, int newOutGroup);
void makeSPR(TreeTemplate<Node> &tree, int cutNodeId, int newBrotherId);
void makeNNI(TreeTemplate<Node> &tree, int nodeId);
void buildVectorOfRegraftingNodes(TreeTemplate<Node> &tree, int nodeForSPR, std::vector <int> & nodeIdsToRegraft);
std::vector<int> getRemainingNeighbors(const Node * node1, const Node * node2);
void getRemainingNeighborsUntilDistance(const Node * node1, const Node * node2, int distance, int d, std::vector <int> & neighbors);
void getNeighboringNodesIdLimitedDistance (TreeTemplate<Node> &tree, int nodeId, int distance, std::vector <int> & neighboringNodeIds);
void buildVectorOfRegraftingNodesLimitedDistance(TreeTemplate<Node> &tree, int nodeForSPR, int distance, std::vector <int> & nodeIdsToRegraft);

void tryAllPossibleSPRsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree,int nodeForSPR, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchProbaOptimization, std::map <std::string, int> genomeMissing);

void tryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchProbaOptimization, std::map <std::string, int> genomeMissing);
void tryAndMakeNNIs(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchProbaOptimization, std::map <std::string, int> genomeMissing);
void makeRandomModifications(TreeTemplate<Node> &tree, int expectedFrequencyNNI, int expectedFrequencySPR, int expectedFrequencyChangeRoot);
void makeDeterministicModifications(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForSPR, int & nodeForRooting);
void makeDeterministicNNIsAndRootChangesOnly(TreeTemplate<Node> &tree, int & nodeForNNI, int & nodeForRooting);

void localOptimizationWithNNIsAndReRootings(const mpi::communicator& world, TreeTemplate<Node> *tree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool &stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, int & nodeForNNI, int & nodeForRooting, std::string & branchProbaOptimization, std::map <std::string, int> genomeMissing, std::vector < double > & NNILks);
void optimizeOnlyDuplicationAndLossRates(const mpi::communicator& world, TreeTemplate<Node> *tree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool &stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, int & nodeForNNI, int & nodeForRooting, std::string & branchProbaOptimization, std::map <std::string, int> genomeMissing);
void fastTryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map <std::string, int> genomeMissing);
void fastTryAllPossibleSPRs(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map <std::string, int> genomeMissing, int sprLimit);
void fastTryAllPossibleSPRsAndReRootings(const mpi::communicator& world, TreeTemplate<Node> *currentTree, TreeTemplate<Node> *bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, double averageDuplicationProbability, double averageLossProbability, bool rearrange, int &numIterationsWithoutImprovement, int server, std::string &branchProbaOptimization, std::map <std::string, int> genomeMissing, int sprLimit, bool optimizeRates);
std::string computeSpeciesTreeLikelihood(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map <std::string, int> genomeMissing, std::vector <double> backupDupProba, std::vector <double> backupLossProba, TreeTemplate<Node> &tree);
void computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map <std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL);
void numericalOptimizationOfDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossProbabilities, std::vector<double> &duplicationProbabilities, bool rearrange, int server, std::string &branchProbaOptimization, std::map < std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL);

#endif 
