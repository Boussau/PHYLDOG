//
// File: DLGeneTreeLikelihood.h
// Created by: Bastien Boussau 
// Created on: Tue October 04 14:16 2011
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


#ifndef DLGeneTreeLikelihood_h
#define DLGeneTreeLikelihood_h


#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Likelihood.all>



// From NumCalc:
//#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/AutoParameter.h>

#include "FastRHomogeneousTreeLikelihood.h"
#include "ReconciliationTools.h"
#include "GeneTreeAlgorithms.h"
#include "GeneTreeLikelihood.h"

//#include "mpi.h" 


using namespace bpp;

/*namespace bpp 
{*/
    
    /**
     * @brief This class adds support for reconciliation to a species tree to the NNIHomogeneousTreeLikelihood class.
     */
    class DLGeneTreeLikelihood:
    public GeneTreeLikelihood
    {
       // NNIHomogeneousTreeLikelihood * nniLk_;
        std::vector <int> _duplicationNumbers;
        std::vector <int> _lossNumbers;
        std::vector <int>  _branchNumbers;        
        mutable std::vector <double> _duplicationProbabilities;
        mutable std::vector <double> _lossProbabilities; 
        std::vector <int> _num0Lineages;
        std::vector <int> _num1Lineages;
        std::vector <int> _num2Lineages;
        mutable std::vector <int> _tentativeDuplicationNumbers;
        mutable std::vector <int> _tentativeLossNumbers; 
        mutable std::vector <int> _tentativeBranchNumbers; 
        mutable std::vector <int> _tentativeNum0Lineages;
        mutable std::vector <int> _tentativeNum1Lineages; 
        mutable std::vector <int> _tentativeNum2Lineages;
        mutable bool _DLStartingGeneTree;
        
    public:
        /**
         * @brief Build a new DLGeneTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param spTree The species tree
         * @param rootedTree rooted version of the gene tree
         * @param seqSp link between sequence and species names
         * @param spId link between species name and species ID
         * @param lossNumbers vector to store loss numbers per branch
         * @param lossProbabilities vector to store expected numbers of losses per branch
         * @param duplicationNumbers vector to store duplication numbers per branch
         * @param duplicationProbabilities vector to store expected numbers of duplications per branch
         * @param branchNumbers vector to store branch numbers in the tree
         * @param num0Lineages vectors to store numbers of branches ending with a loss
         * @param num1Lineages vectors to store numbers of branches ending with 1 gene
         * @param num2Lineages vectors to store numbers of branches ending with 2 genes
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        DLGeneTreeLikelihood(
                                     const Tree & tree,
                                     SubstitutionModel * model,
                                     DiscreteDistribution * rDist,
                                     TreeTemplate<Node> & spTree,  
                                     TreeTemplate<Node> & rootedTree, 
                                     TreeTemplate<Node> & geneTreeWithSpNames,
                                     const std::map <std::string, std::string> seqSp,
                                     std::map <std::string,int> spId,
                                     std::vector <double> & lossProbabilities, 
                                     std::vector <double> & duplicationProbabilities, 
                                     std::vector <int> & num0Lineages,
                                     std::vector <int> & num1Lineages,
                                     std::vector <int> & num2Lineages, 
                                     int speciesIdLimitForRootPosition,
                                     int heuristicsLevel,
                                     int & MLindex, 
                                     bool checkRooted = true,
                                     bool verbose = false,
                                     bool rootOptimization = false, 
                                     bool considerSequenceLikelihood = true, 
                                     bool DLStartingGeneTree = false, 
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
         * @param lossNumbers vector to store loss numbers per branch
         * @param lossProbabilities vector to store expected numbers of losses per branch
         * @param duplicationNumbers vector to store duplication numbers per branch
         * @param duplicationProbabilities vector to store expected numbers of duplications per branch
         * @param branchNumbers vector to store branch numbers in the tree
         * @param num0Lineages vectors to store numbers of branches ending with a loss
         * @param num1Lineages vectors to store numbers of branches ending with 1 gene
         * @param num2Lineages vectors to store numbers of branches ending with 2 genes
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position     
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        DLGeneTreeLikelihood(
                                     const Tree & tree,
                                     const SiteContainer & data,
                                     SubstitutionModel * model,
                                     DiscreteDistribution * rDist,
                                     TreeTemplate<Node> & spTree,  
                                     TreeTemplate<Node> & rootedTree,  
                                     TreeTemplate<Node> & geneTreeWithSpNames,
                                     const std::map <std::string, std::string> seqSp,
                                     std::map <std::string,int> spId,
                                     //std::vector <int> & lossNumbers, 
                                     std::vector <double> & lossProbabilities, 
                                     //std::vector <int> & duplicationNumbers, 
                                     std::vector <double> & duplicationProbabilities, 
                                     //std::vector <int> & branchNumbers, 
                                     std::vector <int> & num0Lineages,
                                     std::vector <int> & num1Lineages,
                                     std::vector <int> & num2Lineages,  
                                     int speciesIdLimitForRootPosition,  
                                     int heuristicsLevel,
                                     int & MLindex, 
                                     bool checkRooted = true,
                                     bool verbose = false, 
                                     bool rootOptimization = false, 
                                     bool considerSequenceLikelihood = true, 
                                     bool DLStartingGeneTree = false, 
                                     unsigned int sprLimit = 2)
        throw (Exception);
        
        /**
         * @brief Copy constructor.
         */ 
        DLGeneTreeLikelihood(const DLGeneTreeLikelihood & lik);
        
        DLGeneTreeLikelihood & operator=(const DLGeneTreeLikelihood & lik);
        
        virtual ~DLGeneTreeLikelihood();
        
        
        
#ifndef NO_VIRTUAL_COV
        DLGeneTreeLikelihood*
#else
        Clonable*
#endif
        clone() const { return new DLGeneTreeLikelihood(*this); }
        
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
        
        void computeSequenceLikelihood();
        
        void computeReconciliationLikelihood();
        
        void computeTreeLikelihood();
        
        double getValue() const throw (Exception);
        
        void fireParameterChanged(const ParameterList & params);
        
        double getTopologyValue() const throw (Exception) { return getValue(); } 
        
        double getScenarioLikelihood() const throw (Exception) { return _scenarioLikelihood; }
        
        void setSpTree(TreeTemplate<Node> & spTree) { if (_spTree) delete _spTree; _spTree = spTree.clone(); }
        
        void setSpId(std::map <std::string, int> & spId) {_spId = spId;}
        
        double testNNI(int nodeId) const throw (NodeException);
        
        void doNNI(int nodeId) throw (NodeException);
        
        std::vector <int> getDuplicationNumbers();
        std::vector <int> getLossNumbers();
        std::vector <int> getBranchNumbers();
        std::vector <int> get0LineagesNumbers() const;
        std::vector <int> get1LineagesNumbers() const;
        std::vector <int> get2LineagesNumbers() const;
        
        ParameterList getParameters() {return nniLk_->getParameters();}
        
        TreeTemplate<Node> & getSpTree() const {return *_spTree;}
        
        TreeTemplate<Node> & getRootedTree() const {return *_rootedTree;}
        
        TreeTemplate<Node> & getGeneTreeWithSpNames() const {return *_geneTreeWithSpNames;}
        
        std::map <std::string, std::string> getSeqSp() {return _seqSp;}
        
        void setExpectedNumbers(std::vector <double> duplicationProbabilities, std::vector <double> lossProbabilities);
        
        int getRootNodeindex();
        
        //void resetSequenceLikelihood();
        
        double getSequenceLikelihood();
        
        void OptimizeSequenceLikelihood(bool yesOrNo) const  {
            _optimizeSequenceLikelihood = yesOrNo;
        }
        
        void OptimizeReconciliationLikelihood(bool yesOrNo) const {
            _optimizeReconciliationLikelihood = yesOrNo;
        }
        
        void optimizeNumericalParameters(map<string, string> params) {

            int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
            bool backupOpt = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
            params[ std::string("optimization.max_number_f_eval")] = 100;
            params[ std::string("optimization.topology")] = "false";
            PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(nniLk_), nniLk_->getParameters(), params, "", true, false);
            params[ std::string("optimization.max_number_f_eval")] = backup;
            params[ std::string("optimization.topology")] = backupOpt;

            

            /* auto_ptr<BackupListener> backupListener;
            unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, false);
            OptimizationTools::optimizeNumericalParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(nniLk_), nniLk_->getParameters(), backupListener.get(), nstep, );*/
        };
        
        void initialize();
        
        void print() const;
        
        /************************************************************************
         * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
         * and executes the ones with the highest likelihood. 
         ************************************************************************/
        void refineGeneTreeSPRs(map<string, string> params);
        
        /************************************************************************
         * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
         * and executes the ones with the highest likelihood. 
         * To do all this as fast as possible, we optimize only a few branch lengths on the SPR tree, 
         * and we use a simple recursion for that.
         ************************************************************************/
        void refineGeneTreeSPRsFast(map<string, string> params);
		void refineGeneTreeSPRsFast2(map<string, string> params);
		void refineGeneTreeSPRsFast3 (map<string, string> params);

        //Not up to date anymore, do not use without checking the code.
        void refineGeneTreeSPRs2(map<string, string> params);


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
