/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a computer program whose purpose is to simultaneously build 
gene and species trees when gene families have undergone duplications and 
losses. It can analyze thousands of gene families in dozens of genomes 
simultaneously, and was presented in an article in Genome Research. Trees and 
parameters are estimated in the maximum likelihood framework, by maximizing 
theprobability of alignments given the species tree, the gene trees and the 
parameters of duplication and loss.

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


#ifndef COALGeneTreeLikelihood_h
#define COALGeneTreeLikelihood_h


#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>

#include "ReconciliationTools.h"
#include "COALTools.h"

#include "GeneTreeAlgorithms.h"
#include "GeneTreeLikelihood.h"
#include "DLGeneTreeLikelihood.h"

#include "mpi.h" 

using namespace bpp;


/*namespace bpp 
{
  */  
    /**
     * @brief This class adds support for coalescence-based reconciliation to a species tree to the NNIHomogeneousTreeLikelihood class.
     */
    class COALGeneTreeLikelihood:
    public GeneTreeLikelihood
    {
        //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
        std::vector < std::vector < std::vector < std::vector< unsigned int > > > > coalCounts_;
        mutable std::vector < std::vector < std::vector < std::vector < unsigned int > > > > tentativeCoalCounts_;

        //coalBl: length of a branch of the species tree, in coalescent units (1 coalescent unit = N generations)
        std::vector < double > coalBl_;
        
        //num12Lineages_ and num22Lineages_: counts of these particular patterns for each branch of the species tree.
        std::vector< unsigned int > num12Lineages_;
        std::vector< unsigned int > num22Lineages_;

        
    public:
      
        
	/**
	* @brief Build a new COALGeneTreeLikelihood object.
	*
	* @param params The parameters to parse.	 
	* @param spTree The species tree
	* @throw Exception in an error occured.
	*/
	    COALGeneTreeLikelihood(std::string file, 
				    map<string, string> params, 
				    TreeTemplate<Node> & spTree) throw (exception) ;
	    
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
         * @param coalCounts vector to store coalescent numbers per branch
         * @param coalBl vector to give number of coalescent units per branch of the species tree
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        COALGeneTreeLikelihood(
                               const Tree & tree,
                               SubstitutionModel * model,
                               DiscreteDistribution * rDist,
                               TreeTemplate<Node> & spTree,  
                               TreeTemplate<Node> & rootedTree, 
                               TreeTemplate<Node> & geneTreeWithSpNames,
                               const std::map <std::string, std::string> seqSp,
                               std::map <std::string,int> spId,
                               std::vector < std::vector < std::vector < std::vector<unsigned int> > > > coalCounts,
                               std::vector < double > coalBl,
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
         * @param coalCounts vector to store coalescent numbers per branch
         * @param coalBl vector to give number of coalescent units per branch of the species tree
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position     
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        COALGeneTreeLikelihood(
                               const Tree & tree,
                               const SiteContainer & data,
                               SubstitutionModel * model,
                               DiscreteDistribution * rDist,
                               TreeTemplate<Node> & spTree,  
                               TreeTemplate<Node> & rootedTree,  
                               TreeTemplate<Node> & geneTreeWithSpNames,
                               const std::map <std::string, std::string> seqSp,
                               std::map <std::string,int> spId,
                               std::vector < std::vector < std::vector < std::vector <unsigned int> > > > coalCounts,
                               std::vector < double > coalBl,
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
        COALGeneTreeLikelihood(const COALGeneTreeLikelihood & lik);
        
        COALGeneTreeLikelihood & operator=(const COALGeneTreeLikelihood & lik);
        
        virtual ~COALGeneTreeLikelihood();
        
        
        
#ifndef NO_VIRTUAL_COV
        COALGeneTreeLikelihood*
#else
        Clonable*
#endif
        clone() const { return new COALGeneTreeLikelihood(*this); }
        
        void initParameters();
        void resetMLindex() ;
        /**
         * @name The NNISearchable interface.
         *
         * Current implementation:
         * When testing a particular NNI, only the branch length of the parent node is optimized (and roughly).
         * All other parameters (substitution model, rate distribution and other branch length are kept at there current value.
         * When performing a NNI, only the topology change is performed.
         * This is up to the user to re-initialize the underlying likelihood data to match the new topology.
         * Usually, this is achieved by calling the topologyChangePerformed() method, which call the reInit() method of the LikelihoodData object.
         * @{
         */
        
        //double getLikelihood() const;
        
        double getLogLikelihood() const;
        
//         void computeSequenceLikelihood();
        
        void computeReconciliationLikelihood();
        
//         void computeTreeLikelihood();
        
        double getValue() const throw (Exception);
        
//         void fireParameterChanged(const ParameterList & params);
        
        double getTopologyValue() const throw (Exception) { return getValue(); } 
        
        double getScenarioLikelihood() const throw (Exception) { return scenarioLikelihood_; }
        
        void setSpTree(TreeTemplate<Node> & spTree) { if (spTree_) delete spTree_; spTree_ = spTree.clone(); }
        
        void setSpId(std::map <std::string, int> & spId) {spId_ = spId;}
        
        double testNNI(int nodeId) const throw (NodeException);
        
        void doNNI(int nodeId) throw (NodeException);
        
        std::vector < std::vector < std::vector<std::vector< unsigned int > > > > getCoalCounts() const;
        
        void computeNumLineagesFromCoalCounts () ;
        
        std::vector< unsigned int > getNum12Lineages() const;
        
        std::vector< unsigned int > getNum22Lineages() const;

        std::vector <double> getCoalBranchLengths() const;
        
//         ParameterList getParameters() {return nniLk_->getParameters();}
        
        TreeTemplate<Node> & getSpTree() const {return *spTree_;}
        
        TreeTemplate<Node> & getRootedTree() const {return *rootedTree_;}
        
        TreeTemplate<Node> & getGeneTreeWithSpNames() const {return *geneTreeWithSpNames_;}
        
        std::map <std::string, std::string> getSeqSp() {return seqSp_;}
        
        void setCoalBranchLengths (std::vector < double > coalBl);
        
        int getRootNodeindex();
        
        //void resetSequenceLikelihood();
        
        double getSequenceLikelihood();
        
        void OptimizeSequenceLikelihood(bool yesOrNo) const  {
            optimizeSequenceLikelihood_ = yesOrNo;
        }
        
        void OptimizeReconciliationLikelihood(bool yesOrNo) const {
            optimizeReconciliationLikelihood_ = yesOrNo;
        }
        
        void initialize();
        
        void print() const;
        
        /************************************************************************
         * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
         * and executes the ones with the highest likelihood. 
         ************************************************************************/

        void refineGeneTreeSPRsFast(map<string, string> params);
        
        
        
        /************************************************************************
         * Tries all NNIs, and accepts NNIs that improve the likelihood as soon as
         * they have been tried.
         ************************************************************************/
        void refineGeneTreeNNIs(map<string, string> params, unsigned int verbose = 0);
        
        /************************************************************************
         * Tells if the gene family is single copy (1 gene per sp)
         ************************************************************************/
        bool isSingleCopy();
        
        
    };
    
    
//} //end of namespace bpp.



#endif
