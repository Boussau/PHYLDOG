//
// File: GeneTreeLikelihood.h
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


#ifndef GeneTreeLikelihood_h
#define GeneTreeLikelihood_h


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
        TreeTemplate<Node> * _spTree;
        TreeTemplate<Node> * _rootedTree;
        TreeTemplate<Node> * _geneTreeWithSpNames;
        const std::map <std::string, std::string> _seqSp;
        std::map <std::string, int> _spId;
        std::set <int> _nodesToTryInNNISearch;
        double _scenarioLikelihood;
        //  mutable double _sequenceLikelihood;
        int _MLindex;
        bool _rootOptimization;
        mutable std::set <int> _tentativeNodesToTryInNNISearch;
        mutable int _tentativeMLindex;
        mutable double _tentativeScenarioLikelihood;
        mutable int _totalIterations;
        mutable int _counter;
        mutable std::vector <int> _listOfPreviousRoots;
        int _speciesIdLimitForRootPosition;
        int _heuristicsLevel;
        mutable bool _optimizeSequenceLikelihood;
        mutable bool _optimizeReconciliationLikelihood;
        mutable bool _considerSequenceLikelihood;
        unsigned int sprLimit_;
        
    public:
        
        GeneTreeLikelihood(); 
        
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
                
        double getScenarioLikelihood() const throw (Exception) { return _scenarioLikelihood; }
        
        void setSpTree(TreeTemplate<Node> & spTree) { if (_spTree) delete _spTree; _spTree = spTree.clone(); }
        
        void setSpId(std::map <std::string, int> & spId) {_spId = spId;}
 
        ParameterList getParameters() {return nniLk_->getParameters();}
        
        TreeTemplate<Node> & getSpTree() const {return *_spTree;}
        
        TreeTemplate<Node> & getRootedTree() const {return *_rootedTree;}
        
        TreeTemplate<Node> & getGeneTreeWithSpNames() const {return *_geneTreeWithSpNames;}
        
        std::map <std::string, std::string> getSeqSp() {return _seqSp;}

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
            
        };

        
        
    };
    
    
//} //end of namespace bpp.



#endif
