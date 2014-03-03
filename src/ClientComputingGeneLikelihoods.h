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

#ifndef _ClientComputingGeneLikelihoods_H_
#define _ClientComputingGeneLikelihoods_H_

// From NumCalc:
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/NumConstants.h>

// From core
#include <Bpp/Text/StringTokenizer.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Io/Nhx.h>


//From the BOOST library 
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>

#include "ReconciliationTools.h"
#include "COALTools.h"
#include "SpeciesTreeExploration.h"

//#include "GeneTreeAlgorithms.h"


#include "GeneTreeLikelihood.h"
#include "DLGeneTreeLikelihood.h"
#include "COALGeneTreeLikelihood.h"




namespace mpi = boost::mpi;


namespace bpp
{
    
  class ClientComputingGeneLikelihoods 
  {
  private:
    mpi::communicator  world_; 
    unsigned int  server_; 
    unsigned int  rank_;
    std::map<std::string, std::string>  params_; 
    std::vector<std::string>  assignedFilenames_; 
    unsigned int  numDeletedFamilies_; 
    TreeTemplate<Node> * geneTree_; 
    TreeTemplate<Node>* tree_; 
    std::vector <int> num0Lineages_;
    std::vector <int>  num1Lineages_; 
    std::vector <int>  num2Lineages_; 
    std::vector <std::vector<int> >  allNum0Lineages_; 
    std::vector <std::vector<int> >  allNum1Lineages_; 
    std::vector <std::vector<int> >  allNum2Lineages_; 
    std::vector <unsigned int> num12Lineages_; 
    std::vector <unsigned int> num22Lineages_; 
    std::vector <double>  lossExpectedNumbers_; 
    std::vector <double>  duplicationExpectedNumbers_;
    std::vector < std::vector < std::vector < std::vector <unsigned int> > > > coalCounts_;
    std::vector <double>  coalBls_;
    std::vector <std::vector<unsigned int> >  allNum12Lineages_; 
    std::vector <std::vector<unsigned int> >  allNum22Lineages_; 
    std::map <std::string, int>  spId_; 
    bool rearrange_;
    bool resetGeneTrees_;
    bool optimizeSpeciesTreeTopology_;
    bool recordGeneTrees_;
    bool stop_;
    unsigned int currentStep_;
    int speciesIdLimitForRootPosition_; 
    int heuristicsLevel_; 
    int MLindex_; 
    double logL_;
    int startRecordingTreesFrom_;
    std::vector <double>  allLogLs_; 
    std::vector <GeneTreeLikelihood *>  treeLikelihoods_;
    std::vector <std::map<std::string, std::string> >  allParams_; 
     std::vector <std::map<std::string, std::string> >  allParamsBackup_;
    std::vector <Alphabet *>  allAlphabets_; 
    std::vector <VectorSiteContainer *>  allDatasets_; 
    std::vector <SubstitutionModel *>  allModels_; 
    std::vector <DiscreteDistribution *>  allDistributions_; 
    std::vector <TreeTemplate<Node> *>  allGeneTrees_; 
    std::vector <TreeTemplate<Node> *>  allUnrootedGeneTrees_; 
    std::vector < std::map <std::string, std::string> > allSeqSps_; 
    std::vector<unsigned int> allSprLimitGeneTree_;
    std::vector <std::vector <std::string> > reconciledTrees_;
    std::vector <std::vector <std::string> > duplicationTrees_;
    std::vector <std::vector <std::string> > lossTrees_;
    string reconciliationModel_;
    string currentSpeciesTree_;
    
  public:   
//Simple constructor
ClientComputingGeneLikelihoods(const mpi::communicator& world, 
    unsigned int & server,
    unsigned int & rank,
    std::map<std::string, std::string> & params) :
    world_(world), server_(server), 
    rank_(rank), params_(params),
    assignedFilenames_(), 
    numDeletedFamilies_ (0), 
    geneTree_(0),
    tree_ (0),
    num0Lineages_(), 
    num1Lineages_(), 
    num2Lineages_(), 
    allNum0Lineages_(), 
    allNum1Lineages_(), 
    allNum2Lineages_(),
    num12Lineages_(), 
    num22Lineages_(),
    lossExpectedNumbers_(),
    duplicationExpectedNumbers_(),
    coalCounts_(),
    coalBls_(), 
    allNum12Lineages_(),
    allNum22Lineages_(),
    spId_(),
    rearrange_(false),
    resetGeneTrees_(false),
    optimizeSpeciesTreeTopology_(false),
    recordGeneTrees_(false),
    stop_(false),
    currentStep_(0),
    speciesIdLimitForRootPosition_(0),
    heuristicsLevel_(0),
    MLindex_(0),
    logL_(0),
    startRecordingTreesFrom_(0),
    allLogLs_(),
    treeLikelihoods_(),
    allParams_(),
    allParamsBackup_(),
    allAlphabets_(),
    allDatasets_(),
    allModels_(),
    allDistributions_(),
    allGeneTrees_(),
    allUnrootedGeneTrees_(),
    allSeqSps_(),
    allSprLimitGeneTree_(),
    reconciledTrees_(),
    duplicationTrees_(),
    lossTrees_(),
    reconciliationModel_("DL"), 
    currentSpeciesTree_("")
    {
      parseOptions();
      
    }
    
    //Copy constructor
    ClientComputingGeneLikelihoods(const ClientComputingGeneLikelihoods& c) :
    world_(c.world_), server_(c.server_), 
    rank_(c.rank_), params_(c.params_),
    assignedFilenames_(c.assignedFilenames_), 
    numDeletedFamilies_ (c.numDeletedFamilies_), 
    geneTree_(c.geneTree_),
    tree_(c.tree_),
    num0Lineages_(c.num0Lineages_), 
    num1Lineages_(c.num1Lineages_), 
    num2Lineages_(c.num2Lineages_), 
    allNum0Lineages_(c.allNum0Lineages_), 
    allNum1Lineages_(c.allNum1Lineages_), 
    allNum2Lineages_(c.allNum2Lineages_),
    num12Lineages_(c.num12Lineages_), num22Lineages_(c.num22Lineages_),
    lossExpectedNumbers_(c.lossExpectedNumbers_),
    duplicationExpectedNumbers_(c.duplicationExpectedNumbers_),
    coalCounts_(c.coalCounts_),
    coalBls_(c.coalBls_), 
    allNum12Lineages_(c.allNum12Lineages_),
    allNum22Lineages_(c.allNum22Lineages_),
    spId_(c.spId_),
    rearrange_(c.rearrange_),
    resetGeneTrees_(c.resetGeneTrees_),
    optimizeSpeciesTreeTopology_(c.optimizeSpeciesTreeTopology_),
    recordGeneTrees_(c.recordGeneTrees_),
    stop_(c.stop_),
    currentStep_(c.currentStep_),
    speciesIdLimitForRootPosition_(c.speciesIdLimitForRootPosition_),
    heuristicsLevel_(c.heuristicsLevel_),
    MLindex_(c.MLindex_),
    logL_(c.logL_),
    startRecordingTreesFrom_(c.startRecordingTreesFrom_),
    allLogLs_(c.allLogLs_),
    treeLikelihoods_(c.treeLikelihoods_),
    allParams_(c.allParams_),
    allParamsBackup_(c.allParamsBackup_),
    allAlphabets_(c.allAlphabets_),
    allDatasets_(c.allDatasets_),
    allModels_(c.allModels_),
    allDistributions_(c.allDistributions_),
    allGeneTrees_(c.allGeneTrees_),
    allUnrootedGeneTrees_(c.allUnrootedGeneTrees_),
    allSeqSps_(c.allSeqSps_),
    allSprLimitGeneTree_(c.allSprLimitGeneTree_),
    reconciledTrees_(c.reconciledTrees_),
    duplicationTrees_(c.duplicationTrees_),
    lossTrees_(c.lossTrees_),
    reconciliationModel_(c.reconciliationModel_),
    currentSpeciesTree_(c.currentSpeciesTree_)
    {}
    
    //= operator
    ClientComputingGeneLikelihoods& operator=(const ClientComputingGeneLikelihoods& c)
    {
      world_ = c.world_; 
      server_ = c.server_; 
      rank_ = c.rank_; 
      params_ = c.params_;
      assignedFilenames_ = c.assignedFilenames_; 
      numDeletedFamilies_  = c.numDeletedFamilies_; 
      geneTree_ = c.geneTree_;
      tree_ = c.tree_;
      num0Lineages_ = c.num0Lineages_; 
      num1Lineages_ = c.num1Lineages_; 
      num2Lineages_ = c.num2Lineages_; 
      allNum0Lineages_ = c.allNum0Lineages_; 
      allNum1Lineages_ = c.allNum1Lineages_; 
      allNum2Lineages_ = c.allNum2Lineages_;
      num12Lineages_ = c.num12Lineages_; num22Lineages_ = c.num22Lineages_;
      lossExpectedNumbers_ = c.lossExpectedNumbers_;
      duplicationExpectedNumbers_ = c.duplicationExpectedNumbers_;
      coalCounts_ = c.coalCounts_;
      coalBls_ = c.coalBls_; 
      allNum12Lineages_ = c.allNum12Lineages_;
      allNum22Lineages_ = c.allNum22Lineages_;
      spId_ = c.spId_;
      rearrange_ = c.rearrange_;
      resetGeneTrees_ = c.resetGeneTrees_;
      optimizeSpeciesTreeTopology_ = c.optimizeSpeciesTreeTopology_;
      recordGeneTrees_ = c.recordGeneTrees_;
      stop_ = c.stop_;
      currentStep_ = c.currentStep_;
      speciesIdLimitForRootPosition_ = c.speciesIdLimitForRootPosition_;
      heuristicsLevel_ = c.heuristicsLevel_;
      MLindex_ = c.MLindex_;
      logL_ = c.logL_;
      startRecordingTreesFrom_ = c.startRecordingTreesFrom_;
      allLogLs_ = c.allLogLs_;
      treeLikelihoods_ = c.treeLikelihoods_;
      allParams_ = c.allParams_;
      allParamsBackup_ = c.allParamsBackup_;
      allAlphabets_ = c.allAlphabets_;
      allDatasets_ = c.allDatasets_;
      allModels_ = c.allModels_;
      allDistributions_ = c.allDistributions_;
      allGeneTrees_ = c.allGeneTrees_;
      allUnrootedGeneTrees_ = c.allUnrootedGeneTrees_;
      allSeqSps_ = c.allSeqSps_;
      allSprLimitGeneTree_ = c.allSprLimitGeneTree_;
      reconciledTrees_ = c.reconciledTrees_;
      duplicationTrees_ = c.duplicationTrees_;
      lossTrees_ = c.lossTrees_;
      reconciliationModel_ = c.reconciliationModel_;
      currentSpeciesTree_ = c.currentSpeciesTree_;
      return *this;
    }
    
    //Destructor
    virtual ~ClientComputingGeneLikelihoods() 
    {			
      if (geneTree_) delete geneTree_;
      if (tree_) delete tree_;
      for (unsigned int i = 0 ; i< assignedFilenames_.size()-numDeletedFamilies_ ; i++) 
      {       
	if (allAlphabets_[i])
	  delete allAlphabets_[i];
	if (allDatasets_[i])
	  delete allDatasets_[i];
	if (allModels_[i])
	  delete allModels_[i];
	if (allDistributions_[i])
	  delete allDistributions_[i];
	if (allGeneTrees_[i])
	  delete allGeneTrees_[i];
	if (treeLikelihoods_[i])
	  delete treeLikelihoods_[i];
	if (allUnrootedGeneTrees_[i])
	  delete allUnrootedGeneTrees_[i];  
      }
    }
    
    //Clone function
    ClientComputingGeneLikelihoods* clone() const { return new ClientComputingGeneLikelihoods(*this); }
    
  public:
    void parseOptions() ;
    
    void parseAssignedGeneFamilies() ;
    
    void MLSearch();
    
    void outputGeneTrees ( unsigned int & bestIndex );
    
    //Get the logL of the species tree according to the gene families handled by the client
    double getValue() const throw (Exception) { return logL_; }
    
    
    
    
  };
  
  
  
  
  
  
  
  
  
}


#endif //_ClientComputingGeneLikelihoods_H_
