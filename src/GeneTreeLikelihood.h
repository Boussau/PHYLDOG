/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a bioinformatics computer program whose purpose is to
simultaneously build gene and species trees when gene families have
undergone duplications and losses. It can analyze thousands of gene
families in dozens of genomes simultaneously, and was presented in
an article in Genome Research. Trees and parameters are estimated
in the maximum likelihood framework, by maximizing the probability
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


#ifndef GeneTreeLikelihood_h
#define GeneTreeLikelihood_h

#include <exception>

#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Nhx.h>


// From NumCalc:
//#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/AutoParameter.h>

#include "FastRHomogeneousTreeLikelihood.h"
#include "ReconciliationTools.h"
#include "GeneTreeAlgorithms.h"
//#include "mpi.h" 



// From core
#include <Bpp/Text/StringTokenizer.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>


//From the BOOST library 
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/communicator.hpp>


namespace mpi = boost::mpi;


using namespace bpp;

/*namespace bpp 
{
  */  

    /**
     * @brief This class adds support for reconciliation to a species tree to the NNIHomogeneousTreeLikelihood class.
     */
    class GeneTreeLikelihood
    {
    protected:
        NNIHomogeneousTreeLikelihood * nniLk_;
        //  TreeTemplate<Node> * _tree;
        TreeTemplate<Node> * spTree_;
        TreeTemplate<Node> * rootedTree_;
        TreeTemplate<Node> * geneTreeWithSpNames_;
        std::map <std::string, std::string> seqSp_; //link between sequence and species
        std::map <std::string, int> spId_;
        std::set <int> nodesToTryInNNISearch_;
        double scenarioLikelihood_;
        //  mutable double _sequenceLikelihood;
        int MLindex_;
        bool rootOptimization_;
        mutable std::set <int> tentativeNodesToTryInNNISearch_;
        mutable int tentativeMLindex_;
        mutable double tentativeScenarioLikelihood_;
        mutable int totalIterations_;
        mutable int counter_;
        mutable std::vector <int> listOfPreviousRoots_;
        int _speciesIdLimitForRootPosition_;
        int heuristicsLevel_;
        mutable bool optimizeSequenceLikelihood_;
        mutable bool optimizeReconciliationLikelihood_;
        mutable bool considerSequenceLikelihood_;
        unsigned int sprLimit_;
	std::map <std::string, std::string > params_;
        unsigned int sprLimitGeneTree_;
	
    public:
        
        GeneTreeLikelihood(); 
        
	/**
	 * @brief Build a new DLGeneTreeLikelihood object.
	 *
	 * @param params The parameters to parse.
	 * @param spTree The species tree
	 * @throw Exception if an error occured.
	 */
	GeneTreeLikelihood(std::string file , map<string, string> params, TreeTemplate<Node> & spTree ) throw (exception);
	
	
        /**
         * @brief Build a new ReconciliationTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param spTree The species tree
         * @param rootedTree rooted version of the gene tree
         * @param seqSp link between sequence and species names
         * @param spId link between species name and species ID
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception if an error occured.
         */
        GeneTreeLikelihood(
                             const Tree & tree,
                             SubstitutionModel * model,
                             DiscreteDistribution * rDist,
                             TreeTemplate<Node> & spTree,  
                             TreeTemplate<Node> & rootedTree, 
                             TreeTemplate<Node> & geneTreeWithSpNames,
                             const std::map <std::string, std::string> seqSp,
                             std::map <std::string,int> spId,
                             int speciesIdLimitForRootPosition,
                             int heuristicsLevel,
                             int & MLindex, 
                             bool checkRooted = true,
                             bool verbose = false,
                             bool rootOptimization = false, 
                             bool considerSequenceLikelihood = true, 
                             unsigned int sprLimit = 2)
        throw (Exception);
        
        /**
         * @brief Build a new ReconciliationTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param data Sequences to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param spTree The species tree
         * @param rootedTree rooted version of the gene tree
         * @param seqSp link between sequence and species names
         * @param spId link between species name and species ID
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position     
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception if an error occured.
         */
        GeneTreeLikelihood(
                             const Tree & tree,
                             const SiteContainer & data,
                             SubstitutionModel * model,
                             DiscreteDistribution * rDist,
                             TreeTemplate<Node> & spTree,  
                             TreeTemplate<Node> & rootedTree,  
                             TreeTemplate<Node> & geneTreeWithSpNames,
                             const std::map <std::string, std::string> seqSp,
                             std::map <std::string,int> spId,
                             int speciesIdLimitForRootPosition,  
                             int heuristicsLevel,
                             int & MLindex, 
                             bool checkRooted = true,
                             bool verbose = false, 
                             bool rootOptimization = false, 
                             bool considerSequenceLikelihood = true, 
                             unsigned int sprLimit = 2)
        throw (Exception);
        
        
        /**
         * @brief Copy constructor.
         */ 
        GeneTreeLikelihood(const GeneTreeLikelihood & lik);
        
        GeneTreeLikelihood & operator=(const GeneTreeLikelihood & lik);
        
        virtual ~GeneTreeLikelihood() {};
        
        
        
#ifndef NO_VIRTUAL_COV
        GeneTreeLikelihood*
#else
        Clonable*
#endif
        clone() const { return new GeneTreeLikelihood(*this); }
                
        double getScenarioLikelihood() const throw (Exception) { return scenarioLikelihood_; }
        
        void setSpTree(TreeTemplate<Node> & spTree) { if (spTree_) delete spTree_; spTree_ = spTree.clone(); }
        
        void setSpId(std::map <std::string, int> & spId) {spId_ = spId;}
 
        ParameterList getParameters() {return nniLk_->getParameters();}
        
        TreeTemplate<Node> & getSpTree() const {return *spTree_;}
        
        TreeTemplate<Node> & getRootedTree() const {return *rootedTree_;}
        
        TreeTemplate<Node> & getGeneTreeWithSpNames() const {return *geneTreeWithSpNames_;}
        
        std::map <std::string, std::string> getSeqSp() {return seqSp_;}

        void OptimizeSequenceLikelihood(bool yesOrNo) const  {
            optimizeSequenceLikelihood_ = yesOrNo;
        }
        
        void OptimizeReconciliationLikelihood(bool yesOrNo) const {
            optimizeReconciliationLikelihood_ = yesOrNo;
        }
        
        NNIHomogeneousTreeLikelihood* getSequenceLikelihoodObject() const {
	  return nniLk_;
	}
        
        unsigned int getSprLimitGeneTree() const {
	 return sprLimitGeneTree_; 
	}
	  
        
        void optimizeNumericalParameters(map<string, string> params) {
            
            int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
            bool backupOpt = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
            params[ std::string("optimization.max_number_f_eval")] = 100;
            params[ std::string("optimization.topology")] = "false";
            PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(nniLk_), nniLk_->getParameters(), params, "", true, false);
            params[ std::string("optimization.max_number_f_eval")] = backup;
            params[ std::string("optimization.topology")] = backupOpt;
            
        };

		bool isInitialized() {
			return nniLk_->isInitialized();
		}
        
        
    };
    
    
//} //end of namespace bpp.



#endif
