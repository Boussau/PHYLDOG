/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a bioinformatics computer program whose purpose is to
simultaneously build gene and species trees when gene families have
undergone duplications and losses. It can analyze thousands of gene
families in dozens of genomes simultaneously, and was presented in
an article in Genome Research. Trees and parameters are estimated
in the maximum likelihood framework, by maximizing theprobability
of alignments given the species tree, the gene trees and the parameters
of duplication and loss.

This software is governed by the CeCILL license under French law and
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
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Node.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>


//#include <algorithm>

//From the BOOST library !!
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;

#include "ReconciliationTools.h"
#include "COALTools.h"

#include "GenericTreeExplorationAlgorithms.h"


using namespace bpp;

void localOptimizationWithNNIsAndReRootings(const mpi::communicator& world, 
                                            TreeTemplate<Node> *& tree, 
                                            TreeTemplate<Node> *& bestTree, 
                                            unsigned int &index, 
                                            unsigned int &bestIndex,  
                                            bool &stop, 
                                            int timeLimit, 
                                            double &logL, 
                                            double &bestlogL, 
                                            std::vector<int> &num0Lineages, 
                                            std::vector<int> &num1Lineages, 
                                            std::vector<int> &num2Lineages, 
                                            std::vector<int> &bestNum0Lineages, 
                                            std::vector<int> &bestNum1Lineages, 
                                            std::vector<int> &bestNum2Lineages, 
                                            std::vector<std::vector<int> > &allNum0Lineages, 
                                            std::vector<std::vector<int> > &allNum1Lineages, 
                                            std::vector<std::vector<int> > &allNum2Lineages, 
                                            std::vector<double> &lossExpectedNumbers, 
                                            std::vector<double> &duplicationExpectedNumbers, 
                                            std::vector<unsigned int> &num12Lineages, 
                                            std::vector<unsigned int> &num22Lineages,
                                            std::vector<unsigned int> &bestNum12Lineages, 
                                            std::vector<unsigned int> &bestNum22Lineages, 
                                            std::vector<double> &coalBls,
                                            std::string &reconciliationModel,
                                            bool rearrange, 
                                            unsigned int &numIterationsWithoutImprovement, 
                                            unsigned int server, 
                                            size_t & nodeForNNI, 
                                            size_t & nodeForRooting, 
                                            map <string, double > treesToLogLk,
                                            std::string & branchExpectedNumbersOptimization, 
                                            std::map <std::string, int> genomeMissing, 
                                            std::vector < double > & NNILks, 
                                            std::vector<double> &rootLks, unsigned int currentStep, 
											const bool fixedOutgroupSpecies_, 
											const std::vector < std::string > outgroupSpecies_);
void optimizeOnlyDuplicationAndLossRates(const mpi::communicator& world, 
                                         TreeTemplate<Node> *& tree, 
                                         TreeTemplate<Node> *& bestTree, 
                                         unsigned int &index, 
                                         unsigned int &bestIndex,  
                                         bool &stop,
                                         int timeLimit,
                                         double &logL, 
                                         double &bestlogL, 
                                         std::vector<int> &num0Lineages, 
                                         std::vector<int> &num1Lineages, 
                                         std::vector<int> &num2Lineages, 
                                         std::vector<int> &bestNum0Lineages, 
                                         std::vector<int> &bestNum1Lineages, 
                                         std::vector<int> &bestNum2Lineages, 
                                         std::vector<std::vector<int> > &allNum0Lineages, 
                                         std::vector<std::vector<int> > &allNum1Lineages, 
                                         std::vector<std::vector<int> > &allNum2Lineages, 
                                         std::vector<double> &lossExpectedNumbers, 
                                         std::vector<double> &duplicationExpectedNumbers, 
                                         std::vector<unsigned int> &num12Lineages, 
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange, 
                                         unsigned int &numIterationsWithoutImprovement, 
                                         unsigned int server, 
                                         std::string & branchExpectedNumbersOptimization, 
                                         std::map <std::string, int> genomeMissing,
                                         unsigned int currentStep);
void fastTryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, 
                                                TreeTemplate<Node> *& currentTree, 
                                                TreeTemplate<Node> *& bestTree, 
                                                unsigned int &index, 
                                                unsigned int &bestIndex,  
                                                bool stop, 
                                                int timeLimit,
                                                double &logL, 
                                                double &bestlogL, 
                                                std::vector<int> &num0Lineages, 
                                                std::vector<int> &num1Lineages, 
                                                std::vector<int> &num2Lineages, 
                                                std::vector<int> &bestNum0Lineages, 
                                                std::vector<int> &bestNum1Lineages, 
                                                std::vector<int> &bestNum2Lineages, 
                                                std::vector<std::vector<int> > &allNum0Lineages, 
                                                std::vector<std::vector<int> > &allNum1Lineages, 
                                                std::vector<std::vector<int> > &allNum2Lineages, 
                                                std::vector<double> &lossExpectedNumbers, 
                                                std::vector<double> &duplicationExpectedNumbers, 
                                                /*double averageDuplicationProbability, 
                                                double averageLossProbability,*/ 
                                                bool rearrange, 
                                                unsigned int &numIterationsWithoutImprovement, 
                                                unsigned int server, 
                                                std::string &branchExpectedNumbersOptimization, 
                                                std::map <std::string, int> genomeMissing, 
                                                bool optimizeRates, unsigned int currentStep);
void fastTryAllPossibleSPRs(const mpi::communicator& world, 
                            TreeTemplate<Node> *& currentTree, 
                            TreeTemplate<Node> *& bestTree, 
                            unsigned int &index, unsigned int &bestIndex,  
                            bool stop, 
                            int timeLimit, 
                            double &logL, 
                            double &bestlogL, 
                            std::vector<int> &num0Lineages, 
                            std::vector<int> &num1Lineages, 
                            std::vector<int> &num2Lineages, 
                            std::vector<int> &bestNum0Lineages, 
                            std::vector<int> &bestNum1Lineages, 
                            std::vector<int> &bestNum2Lineages, 
                            std::vector<std::vector<int> > &allNum0Lineages, 
                            std::vector<std::vector<int> > &allNum1Lineages, 
                            std::vector<std::vector<int> > &allNum2Lineages, 
                            std::vector<double> &lossExpectedNumbers, 
                            std::vector<double> &duplicationExpectedNumbers, 
                            std::vector<unsigned int> &num12Lineages, 
                            std::vector<unsigned int> &num22Lineages,
                            std::vector<unsigned int> &bestNum12Lineages, 
                            std::vector<unsigned int> &bestNum22Lineages, 
                            std::vector<double> &coalBls,
                            std::string &reconciliationModel,
                            bool rearrange, 
                            unsigned int &numIterationsWithoutImprovement, 
                            unsigned int server, 
                            std::string &branchExpectedNumbersOptimization, 
                            std::map <std::string, int> genomeMissing, 
                            int sprLimit, 
                            bool optimizeRates, unsigned int currentStep, 
							const bool fixedOutgroupSpecies_, 
							const std::vector < std::string > outgroupSpecies_);
void fastTryAllPossibleSPRsAndReRootings(const mpi::communicator& world, 
                                         TreeTemplate<Node> *& currentTree, 
                                         TreeTemplate<Node> *& bestTree, 
                                         unsigned int &index, 
                                         unsigned int &bestIndex,  
                                         bool stop, 
                                         int timeLimit,
                                         double &logL, 
                                         double &bestlogL, 
                                         std::vector<int> &num0Lineages, 
                                         std::vector<int> &num1Lineages, 
                                         std::vector<int> &num2Lineages, 
                                         std::vector<int> &bestNum0Lineages, 
                                         std::vector<int> &bestNum1Lineages, 
                                         std::vector<int> &bestNum2Lineages, 
                                         std::vector<std::vector<int> > &allNum0Lineages, 
                                         std::vector<std::vector<int> > &allNum1Lineages, 
                                         std::vector<std::vector<int> > &allNum2Lineages, 
                                         std::vector<double> &lossExpectedNumbers, 
                                         std::vector<double> &duplicationExpectedNumbers, 
                                         std::vector<unsigned int> &num12Lineages, 
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<unsigned int> &bestNum12Lineages, 
                                         std::vector<unsigned int> &bestNum22Lineages, 
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange, 
                                         unsigned int &numIterationsWithoutImprovement, 
                                         unsigned int server, 
                                         std::string &branchExpectedNumbersOptimization, 
                                         std::map <std::string, int> genomeMissing, 
                                         int sprLimit, 
                                         bool optimizeRates, unsigned int currentStep, 
										 const bool fixedOutgroupSpecies_, 
										 const std::vector < std::string > outgroupSpecies_);
void broadcastsAllInformation(const mpi::communicator& world, 
                              unsigned int server, bool stop, bool &rearrange, 
                              std::vector<double> &lossExpectedNumbers, 
                              std::vector<double> &duplicationExpectedNumbers, 
                              std::vector<double> &coalBls,
                              std::string & currentSpeciesTree, 
                              unsigned int &currentStep, 
                              std::string &reconciliationModel);
void broadcastsAllInformationButStop(const mpi::communicator& world, unsigned int server, 
                                     bool &rearrange, 
                                     std::vector<double> &lossExpectedNumbers, 
                                     std::vector<double> &duplicationExpectedNumbers, 
                                     std::vector<double> &coalBls, 
                                     std::string &currentSpeciesTree, 
                                     unsigned int &currentStep, 
                                     std::string &reconciliationModel);
std::string computeSpeciesTreeLikelihood(const mpi::communicator& world, 
                                         unsigned int &index, 
                                         bool stop, 
                                         double &logL, 
                                         std::vector<int> &num0Lineages, 
                                         std::vector<int> &num1Lineages, 
                                         std::vector<int> &num2Lineages, 
                                         std::vector<std::vector<int> > &allNum0Lineages, 
                                         std::vector<std::vector<int> > &allNum1Lineages, 
                                         std::vector<std::vector<int> > &allNum2Lineages, 
                                         std::vector<double> &lossExpectedNumbers, 
                                         std::vector<double> &duplicationExpectedNumbers, 
                                         std::vector<unsigned int> &num12Lineages, 
                                         std::vector<unsigned int> &num22Lineages,
                                         std::vector<double> &coalBls,
                                         std::string &reconciliationModel,
                                         bool rearrange, 
                                         unsigned int server, 
                                         std::string &branchExpectedNumbersOptimization, 
                                         std::map <std::string, int> genomeMissing, 
                                         TreeTemplate<Node> &tree, 
                                         unsigned int currentStep);
std::string computeSpeciesTreeLikelihoodWithGivenStringSpeciesTree(const mpi::communicator& world, 
                                                                   unsigned int &index, 
                                                                   bool stop, 
                                                                   double &logL, 
                                                                   std::vector<int> &num0Lineages, 
                                                                   std::vector<int> &num1Lineages, 
                                                                   std::vector<int> &num2Lineages, 
                                                                   std::vector< std::vector<int> > &allNum0Lineages, 
                                                                   std::vector< std::vector<int> > &allNum1Lineages, 
                                                                   std::vector< std::vector<int> > &allNum2Lineages, 
                                                                   std::vector<double> &lossExpectedNumbers, 
                                                                   std::vector<double> &duplicationExpectedNumbers, 
                                                                   std::vector<unsigned int> &num12Lineages, 
                                                                   std::vector<unsigned int> &num22Lineages,
                                                                   std::vector<double> &coalBls,
                                                                   std::string &reconciliationModel,
                                                                   bool rearrange, 
                                                                   unsigned int server, 
                                                                   std::string &branchExpectedNumbersOptimization, 
                                                                   std::map < std::string, int> genomeMissing, 
                                                                   TreeTemplate<Node> &tree, 
                                                                   std::string currentSpeciesTree, 
                                                                   bool firstTime, 
                                                                   unsigned int currentStep);
void computeSpeciesTreeLikelihoodWhileOptimizingDuplicationAndLossRates(const mpi::communicator& world, 
                                                                        unsigned int &index, bool stop, double &logL, 
/*std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, */ 
                                                                        std::vector<int> &num0Lineages, 
                                                                        std::vector<int> &num1Lineages, 
                                                                        std::vector<int> &num2Lineages, 
                                                                        std::vector<std::vector<int> > &allNum0Lineages, 
                                                                        std::vector<std::vector<int> > &allNum1Lineages, 
                                                                        std::vector<std::vector<int> > &allNum2Lineages, 
                                                                        std::vector<double> &lossExpectedNumbers, 
                                                                        std::vector<double> &duplicationExpectedNumbers, 
                                                                        std::vector<unsigned int> &num12Lineages, 
                                                                        std::vector<unsigned int> &num22Lineages,
                                                                        std::vector<double> &coalBls,
                                                                        std::string &reconciliationModel,
                                                                        bool rearrange, unsigned int server, 
                                                                        std::string &branchExpectedNumbersOptimization, 
                                                                        std::map <std::string, int> genomeMissing, 
                                                                        TreeTemplate<Node> &tree, double & bestlogL, 
                                                                        unsigned int currentStep);
void firstCommunicationsServerClient (const mpi::communicator & world, 
                                      unsigned int & server, std::vector <unsigned int>  & numbersOfGenesPerClient, 
                                      unsigned int & assignedNumberOfGenes, 
                                      std::vector<std::string> & assignedFilenames, 
                                      std::vector <std::vector<std::string> > & listOfOptionsPerClient, 
                                      bool & optimizeSpeciesTreeTopology, int & SpeciesNodeNumber, 
                                      std::vector <double> & lossExpectedNumbers, 
                                      std::vector <double> & duplicationExpectedNumbers, 
                                      std::vector <int> & num0Lineages, 
                                      std::vector <int> & num1Lineages, 
                                      std::vector <int> & num2Lineages, 
                                      std::vector <unsigned int> & num12Lineages, 
                                      std::vector <unsigned int> & num22Lineages,
                                      std::vector < double >& coalBls,
                                      std::string & currentSpeciesTree);
/******************************************************************************/
// This function runs the second communication between the server and the clients.
// The clients send back how many families they deal with, 
// and the server can then decide whether to stop or not.
/******************************************************************************/
void secondCommunicationsServerClient (const mpi::communicator & world , 
                                       unsigned int & server, 
                                       unsigned int & whoami, 
                                       unsigned int & finalNumberOfGeneFamilies, 
                                       vector<unsigned int> & finalNumbersOfGeneFamilies, 
                                       std::vector <double> & lossExpectedNumbers, 
                                       std::vector <double> & duplicationExpectedNumbers,
                                       std::vector<double> &coalBls,
                                       std::string & currentSpeciesTree);
/******************************************************************************/
// This function runs the communication between the server and the clients for building a MRP species tree.
// The clients send back gene trees with one gene per species, 
// and the server can then build the MRP species tree.
/******************************************************************************/
void mrpCommunicationsServerClient (const mpi::communicator & world, 
                                    unsigned int & server, 
                                    unsigned int & whoami, 
                                    string & trees1PerSpecies, 
                                    std::vector<string> & allTrees1PerSpecies);
/******************************************************************************/
// This function runs the communication between the server and the clients for counting gene families after filtering.
/******************************************************************************/

void numberOfFilteredFamiliesCommunicationsServerClient (const mpi::communicator & world, unsigned int & server, 
                                                         unsigned int & whoami, unsigned int & numberOfGeneFamilies);
/******************************************************************************/
// This function runs the communication between the server and the clients to stop the program.
// The server send stop and the index of the best tree found.
/******************************************************************************/
void lastCommunicationsServerClient (const mpi::communicator & world, 
                                     unsigned int & server, 
                                     bool & stop, 
                                     unsigned int & bestIndex);
/************************************************************************
 * Gathers information from clients. 
 ************************************************************************/
void gathersInformationFromClients (const mpi::communicator & world, 
                                    unsigned int & server,
                                    unsigned int &whoami, 
                                    double &logL, 
                                    std::vector<int> &num0Lineages, 
                                    std::vector<int> &num1Lineages, 
                                    std::vector<int> &num2Lineages, 
                                    std::vector< std::vector<int> > &allNum0Lineages, 
                                    std::vector< std::vector<int> > &allNum1Lineages, 
                                    std::vector< std::vector<int> > &allNum2Lineages, 
                                    std::vector<unsigned int> &num12Lineages,
                                    std::vector<unsigned int> &num22Lineages, 
                                    std::string &reconciliationModel );
void inputNNIAndRootLks(std::vector <double> & NNILks, 
                        std::vector <double> & rootLks, 
                        std::map<std::string, std::string> & params, 
                        std::string & suffix);
void outputNNIAndRootLks(std::vector <double> & NNILks, 
                         std::vector <double> & rootLks, 
                         std::map<std::string, std::string> & params, 
                         std::string & suffix);



/*void makeRandomModifications(TreeTemplate<Node> &tree, int expectedFrequencyNNI, int expectedFrequencySPR, int expectedFrequencyChangeRoot);

*/


/*void numericalOptimizationOfDuplicationAndLossRates(const mpi::communicator& world, int &index, bool stop, double &logL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector< std::vector<int> > AllLosses, std::vector< std::vector<int> > AllDuplications, std::vector< std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector< std::vector<int> > &allNum0Lineages, std::vector< std::vector<int> > &allNum1Lineages, std::vector< std::vector<int> > &allNum2Lineages, std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers, bool rearrange, int server, std::string &branchExpectedNumbersOptimization, std::map < std::string, int> genomeMissing, TreeTemplate<Node> &tree, double & bestlogL);*/
#endif 



/*void tryAllPossibleSPRsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *& currentTree, TreeTemplate<Node> *& bestTree,int nodeForSPR, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers,  bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchExpectedNumbersOptimization, std::map <std::string, int> genomeMissing);
void tryAllPossibleReRootingsAndMakeBestOne(const mpi::communicator& world, TreeTemplate<Node> *& currentTree, TreeTemplate<Node> *& bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers,  bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchExpectedNumbersOptimization, std::map <std::string, int> genomeMissing);
void tryAndMakeNNIs(const mpi::communicator& world, TreeTemplate<Node> *& currentTree, TreeTemplate<Node> *& bestTree, int &index, int &bestIndex,  bool stop, double &logL, double &bestlogL, std::vector<int> &lossNumbers, std::vector<int> &duplicationNumbers, std::vector<int> &branchNumbers, std::vector<int> &bestLossNumbers, std::vector<int> &bestDuplicationNumbers, std::vector<int> &bestBranchNumbers, std::vector<std::vector<int> > AllLosses, std::vector<std::vector<int> > AllDuplications, std::vector<std::vector<int> > AllBranches, std::vector<int> &num0Lineages, std::vector<int> &num1Lineages, std::vector<int> &num2Lineages, std::vector<int> &bestNum0Lineages, std::vector<int> &bestNum1Lineages, std::vector<int> &bestNum2Lineages, std::vector<std::vector<int> > &allNum0Lineages, std::vector<std::vector<int> > &allNum1Lineages, std::vector<std::vector<int> > &allNum2Lineages, std::vector<double> &lossExpectedNumbers, std::vector<double> &duplicationExpectedNumbers,  bool rearrange, int &numIterationsWithoutImprovement, int server, std::string & branchExpectedNumbersOptimization, std::map <std::string, int> genomeMissing);*/


