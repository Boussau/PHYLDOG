/* This file contains various functions useful for reconciliations, such as reconciliation computation, printing of trees with integer indexes, search of a root with reconciliation...*/

#ifndef _RECONCILIATIONTOOLS_H_
#define _RECONCILIATIONTOOLS_H_

#include <queue>
#include <set>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/Newick.h>


// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/NumConstants.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/Clonable.h>
#include <Utils/Number.h>
#include <Utils/BppString.h>
#include <Utils/KeyvalTools.h>

#include "mpi.h" 
using namespace bpp;

const std::string SPECIESID="SPECIESID";
const std::string EVENT="EVENT";
const std::string LOSSES="LOSSES";
const std::string DUPLICATIONS="DUPLICATIONS";
const std::string EVENTSPROBA="EVENTSPROBA";
const std::string LOWLIK="LOWLIK";
const std::string NUMGENES="NUMGENES";
const std::string NUMLINEAGES="NUMLINEAGES";

const double UNLIKELY=-100000000000000000000.0;
const double SMALLPROBA=0.0000000001;
const double BIGPROBA=0.9999999999;
const int MAXFILENAMESIZE = 500;
const int MAXSPECIESTREESIZE = 10000; //size of the species tree, in number of CHARs, as it is written in Newick format

void assignArbitraryBranchLengths(TreeTemplate<Node> & tree);
void reNumber (TreeTemplate<Node> & tree, Node * noeud, int & index);
void reNumber (TreeTemplate<Node> & tree);
std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree);
std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities);
void printVector(std::vector<int> & v);
void resetLossesAndDuplications(TreeTemplate<Node> & tree, /*std::vector <int> &lossNumbers, */std::vector <double> &lossProbabilities, /*std::vector <int> &duplicationNumbers, */std::vector <double> &duplicationProbabilities);
void resetLossesDuplicationsSpeciations(TreeTemplate<Node> & tree, std::vector <int> &lossNumbers, std::vector <double> &lossProbabilities, std::vector <int> &duplicationNumbers, std::vector <double> &duplicationProbabilities, std::vector <int> &branchNumbers);
void resetLossesDuplicationsSpeciationsForGivenNodes(TreeTemplate<Node> & tree, std::vector <int> & lossNumbers, std::vector <double> & lossProbabilities, std::vector <int> & duplicationNumbers, std::vector <double> & duplicationProbabilities, std::vector <int> & branchNumbers, std::vector <int> nodesToUpdate, std::map <int, std::vector<int> > & geneNodeIdToLosses, std::map <int, int > & geneNodeIdToDuplications, std::map <int, std::vector<int> > & geneNodeIdToSpeciations);
void resetVector(std::vector<int> & v);
void resetVector(std::vector<double> & v);
void resetSpeciesIds (TreeTemplate<Node> & tree);
void resetSpeciesIdsForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds); 
void resetSpeciesIdsForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate);
void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::vector<int> &branchNumbers, std::map <int,int> & geneNodeIdToDuplications, std::map <int, std::vector <int> > & geneNodeIdToLosses, std::map <int, std::vector <int> > & geneNodeIdToSpeciations) ;
void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::map <int,int> &geneNodeIdToDuplications, std::map <int, std::vector <int> > &geneNodeIdToLosses, std::map <int, std::vector <int> > &geneNodeIdToSpeciations) ;
void computeDuplicationAndLossProbabilities (int i, int j, int k, double & lossProbability, double & duplicationProbability);
void computeDuplicationAndLossProbabilitiesForAllBranches (std::vector <int> numOGenes, std::vector <int> num1Genes, std::vector <int> num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities);
void computeAverageDuplicationAndLossProbabilitiesForAllBranches (std::vector <int> numOGenes, std::vector <int> num1Genes, std::vector <int> num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities);
double computeBranchProbability (double duplicationProbability, double lossProbability, int numberOfLineages);
double computeLogBranchProbability (double duplicationProbability, double lossProbability, int numberOfLineages);
double computeBranchProbabilityAtRoot (double duplicationProbability, double lossProbability, int numberOfLineages);
double computeLogBranchProbabilityAtRoot (double duplicationProbability, double lossProbability, int numberOfLineages);
void computeScenarioScore (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::vector<int> &branchNumbers, std::map <int, std::vector <int> > &geneNodeIdToSpeciations, std::vector<double>duplicationProbabilities, std::vector <double> lossProbabilities, std::vector <int> &num0lineages, std::vector <int> &num1lineages, std::vector <int> &num2lineages);
std::string nodeToParenthesisWithIntNodeValues(const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName) throw (NodeNotFoundException);
std::string treeToParenthesisWithIntNodeValues(const Tree & tree, bool bootstrap, const std::string & propertyName);
std::string nodeToParenthesisWithDoubleNodeValues(const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName) throw (NodeNotFoundException);
std::string treeToParenthesisWithDoubleNodeValues(const Tree & tree, bool bootstrap, const std::string & propertyName);
void setLossesAndDuplications(TreeTemplate<Node> & tree, 
                              std::vector <int> &lossNumbers, 
                              std::vector <int> &duplicationNumbers);
void assignNumLineagesOnSpeciesTree(TreeTemplate<Node> & tree, 
                                    std::vector <int> &num0Lineages, 
                                    std::vector <int> &num1Lineages, 
                                    std::vector <int> &num2Lineages);
double makeReconciliationAtGivenRoot (TreeTemplate<Node> * tree, 
                                      TreeTemplate<Node> * geneTree, 
                                      std::map<std::string, std::string > seqSp, 
                                      std::vector< double> lossProbabilities, 
                                      std::vector < double> duplicationProbabilities, 
                                      int MLindex);
double findMLReconciliation (TreeTemplate<Node> * spTree, 
                             TreeTemplate<Node> * geneTreeSafe, 
                             std::map<std::string, std::string > seqSp, 
                             std::vector<int> & lossNumbers, 
                             std::vector< double> lossProbabilities, 
                             std::vector< int> & duplicationNumbers, 
                             std::vector < double> duplicationProbabilities, 
                             int & MLindex, 
                             std::vector<int> &branchNumbers, 
                             int speciesIdLimitForRootPosition, 
                             int heuristicsLevel, 
                             std::vector <int> &num0lineages, 
                             std::vector <int> &num1lineages, 
                             std::vector <int> &num2lineages, 
                             std::set <int> &nodesToTryInNNISearch);
int assignSpeciesIdToLeaf(Node * node,  
                          const std::map<std::string, 
                          std::string > & seqSp, 
                          const std::map<std::string, int > & spID);
void recoverLosses(Node & node, 
                   int & a, const int & b, int & olda, 
                   const TreeTemplate<Node> & tree, 
                   double & likelihoodCell, 
                   const std::vector< double> & lossRates, 
                   const std::vector< double> & duplicationRates);
void recoverLossesWithDuplication(const Node & nodeA, 
                                  const int &a, 
                                  const int &olda, 
                                  const TreeTemplate<Node> & tree,
                                  double & likelihoodCell, 
                                  const std::vector< double> & lossRates, 
                                  const std::vector< double> & duplicationRates);
double computeConditionalLikelihoodAndAssignSpId(TreeTemplate<Node> & tree,
                                                 std::vector <Node *> sons,  
                                                 double & rootLikelihood, 
                                                 double & son0Likelihood,
                                                 double & son1Likelihood,
                                                 const std::vector< double> & lossRates, 
                                                 const std::vector< double> & duplicationRates, 
                                                 int & rootSpId,
                                                 const int & son0SpId,
                                                 const int & son1SpId,
                                                 int & rootDupData,
                                                 int & son0DupData,
                                                 int & son1DupData,
                                                 bool atRoot);
double computeSubtreeLikelihoodPostorder(TreeTemplate<Node> & spTree, 
                                         TreeTemplate<Node> & geneTree, 
                                         Node * node, 
                                         const std::map<std::string, std::string > & seqSp, 
                                         const std::map<std::string, int > & spID, 
                                         std::vector <std::vector<double> > & likelihoodData, 
                                         const std::vector< double> & lossRates, 
                                         const std::vector < double> & duplicationRates, 
                                         std::vector <std::vector<int> > & speciesIDs, 
                                         std::vector <std::vector<int> > & dupData);
void computeRootingLikelihood(TreeTemplate<Node> & spTree, 
                              Node * node, 
                              std::vector <std::vector<double> > & likelihoodData, 
                              const std::vector< double> & lossRates, 
                              const std::vector < double> & duplicationRates, 
                              std::vector <std::vector<int> > & speciesIDs, 
                              std::vector <std::vector<int> > & dupData, 
                              int sonNumber, 
                              std::map <double, Node*> & LksToNodes);
void computeSubtreeLikelihoodPreorder(TreeTemplate<Node> & spTree, 
                                      TreeTemplate<Node> & geneTree, 
                                      Node * node, 
                                      const std::map<std::string, std::string > & seqSp, 
                                      const std::map<std::string, int > & spID, 
                                      std::vector <std::vector<double> > & likelihoodData, 
                                      const std::vector< double> & lossRates, 
                                      const std::vector < double> & duplicationRates, 
                                      std::vector <std::vector<int> > & speciesIDs, 
                                      std::vector <std::vector<int> > & dupData,
                                      int sonNumber, 
                                      std::map <double, Node*> & LksToNodes);
void recoverLossesAndLineages(Node & node, int & a, const int & b, int & olda, 
                              int & a0, 
                              const TreeTemplate<Node> & tree, 
                              int & dupData, 
                              std::vector<int> &num0lineages, 
                              std::vector<int> &num1lineages);
void recoverLossesAndLineagesWithDuplication(const Node & nodeA, 
                                             const int &a, 
                                             const int &olda, 
                                             const TreeTemplate<Node> & tree, 
                                             std::vector <int> &num0lineages);
void computeNumbersOfLineagesInASubtree(TreeTemplate<Node> & tree,
                                        std::vector <Node *> sons,  
                                        int & rootSpId,
                                        const int & son0SpId,
                                        const int & son1SpId,
                                        int & rootDupData,
                                        int & son0DupData,
                                        int & son1DupData,
                                        bool atRoot, 
                                        std::vector <int> &num0lineages, 
                                        std::vector <int> &num1lineages, 
                                        std::vector <int> &num2lineages, 
                                        std::set <int> &branchesWithDuplications);
void computeNumbersOfLineagesFromRoot(TreeTemplate<Node> * spTree, 
                                      TreeTemplate<Node> * geneTree, 
                                      Node * node, 
                                      std::map<std::string, std::string > seqSp,
                                      std::map<std::string, int > spID,
                                      std::vector <int> &num0lineages, 
                                      std::vector <int> &num1lineages, 
                                      std::vector <int> &num2lineages, 
                                      std::vector <std::vector<int> > & speciesIDs,
                                      std::vector <std::vector<int> > & dupData, 
                                      std::set <int> & branchesWithDuplications);
double findMLReconciliationDR (TreeTemplate<Node> * spTree, 
                               TreeTemplate<Node> * geneTree, 
                               std::map<std::string, std::string > seqSp,
                               std::map<std::string, int > spID,
/*vector<int> & lossNumbers,*/ 
                               std::vector< double> lossRates, 
/*vector< int> & duplicationNumbers,*/ 
                               std::vector < double> duplicationRates, 
                               int & MLindex, 
/*vector<int> &branchNumbers, int speciesIdLimitForRootPosition, int heuristicsLevel,*/ 
                               std::vector <int> &num0lineages, 
                               std::vector <int> &num1lineages, 
                               std::vector <int> &num2lineages, 
                               std::set <int> &nodesToTryInNNISearch);
double computeAverageLossProportionOnCompletelySequencedLineages(std::vector <int> & num0lineages, 
                                    std::vector <int> & num1lineages, 
                                    std::vector <int> & num2lineages, 
                                    std::map <std::string, int> & genomeMissing, 
                                    TreeTemplate<Node> & tree);
/*void alterLineageCountsWithCoverages(std::vector <int> & num0lineages, 
                                     std::vector <int> & num1lineages, 
                                     std::vector <int> & num2lineages, 
                                     std::map <std::string, int> & genomeMissing, 
                                     TreeTemplate<Node> & tree);*/
void alterLineageCountsWithCoverages(std::vector <int> & num0lineages, 
                                     std::vector <int> & num1lineages, 
                                     std::vector <int> & num2lineages, 
                                     std::map <std::string, int> & genomeMissing, 
                                     TreeTemplate<Node> & tree, 
                                     bool average);
/*void alterLineageCountsWithCoveragesAverage(std::vector <int> & num0lineages, 
                                            std::vector <int> & num1lineages, 
                                            std::vector <int> & num2lineages, 
                                            std::map <std::string, int> & genomeMissing, 
                                            TreeTemplate<Node> & tree);*/
void alterLineageCountsWithCoveragesInitially(std::vector <int> & num0lineages, 
                                              std::vector <int> & num1lineages, 
                                              std::vector <int> & num2lineages, 
                                              std::map <std::string, int> & genomeMissing, 
                                              TreeTemplate<Node> & tree);
void resetLineageCounts(std::vector <int> & num0lineages, 
                        std::vector <int> & num1lineages, 
                        std::vector <int> & num2lineages);
void deleteSubtreeProperties(Node &node);
void deleteTreeProperties(TreeTemplate<Node> & tree);
std::map <std::string, int> computeSpeciesNamesToIdsMap (TreeTemplate<Node> & tree);
void computeDuplicationAndLossRatesForTheSpeciesTree (std::string &branchProbaOptimization, 
                                                      std::vector <int> & num0Lineages, 
                                                      std::vector <int> & num1Lineages, 
                                                      std::vector <int> & num2Lineages, 
                                                      std::vector<double> & lossProbabilities, 
                                                      std::vector<double> & duplicationProbabilities, 
                                                      std::map <std::string, int> & genomeMissing, 
                                                      TreeTemplate<Node> & tree);
void computeDuplicationAndLossRatesForTheSpeciesTreeInitially (std::string &branchProbaOptimization, 
                                                               std::vector <int> & num0Lineages, 
                                                               std::vector <int> & num1Lineages, 
                                                               std::vector <int> & num2Lineages, 
                                                               std::vector<double> & lossProbabilities, 
                                                               std::vector<double> & duplicationProbabilities, 
                                                               std::map <std::string, int> & genomeMissing, 
                                                               TreeTemplate<Node> & tree);
void removeLeaf(TreeTemplate<Node> & tree, std::string toRemove);
void cleanVectorOfOptions (std::vector<std::string> & listOptions, bool sizeConstraint);
bool sortMaxFunction (std::pair <std::string, double> i, std::pair <std::string, double> j);
bool sortMinFunction (std::pair <std::vector<std::string>, double> i, std::pair <std::vector<std::string>, double> j);
void generateListOfOptionsPerClient(std::vector <std::string> listOptions, 
                                    int size, 
                                    std::vector <std::vector<std::string> > &listOfOptionsPerClient, 
                                    std::vector <int> &numbersOfGenesPerClient);
std::string removeComments(
                           const std::string & s,
                           const std::string & begin,
                           const std::string & end);

#endif  //_RECONCILIATIONTOOLS_H_

