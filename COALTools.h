//
//  CoalTools.h
//  phyldog
//
//  Created by Bastien Boussau on 07/03/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#ifndef CoalTools_h
#define CoalTools_h
#include "ReconciliationTools.h"
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>


/*****************************************************************************
 * This function performs a postorder tree traversal in order to find 
 * fill up vectors of counts of coalescence events for rootings. 
 * When followed by the preorder tree traversal function, 
 * vectors of counts for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalCountsPostorder(TreeTemplate<Node> & spTree, 
                                 TreeTemplate<Node> & geneTree, 
                                 Node * node, 
                                  std::map<std::string, std::string > & seqSp, 
                                  std::map<std::string, int > & spID, 
                                 std::vector < std::vector< std::vector< std::vector<unsigned int > > > > & coalCounts,
                                 std::vector <std::vector<unsigned int> > & speciesIDs);


/*****************************************************************************
 * Utilitary functions to initialize vectors of counts at leaves.
 ****************************************************************************/
void initializeCountVectors(std::vector< std::vector< std::vector<unsigned int> > > & vec) ;
void initializeCountVector(std::vector<std::vector<unsigned int> >  & vec) ;

/*****************************************************************************
 * Utilitary functions to increment vectors of counts.
 ****************************************************************************/
void incrementOutCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos);
void incrementInCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos);


/*****************************************************************************
 * Computes a vector of counts from two son vectors, and assigns species ID to the
 * father node.
 ****************************************************************************/

void computeCoalCountsFromSons (TreeTemplate<Node> & tree, std::vector <Node *> sons, 
                                unsigned int & rootSpId, 
                                const unsigned int & son0SpId,
                                const unsigned int & son1SpId,
                                std::vector< std::vector<unsigned int> > & coalCountsFather,
                                std::vector< std::vector<unsigned int> > & coalCountsSon0,
                                std::vector< std::vector<unsigned int> > & coalCountsSon1);


/*****************************************************************************
 * This function recovers ILS by comparing a subtree in a gene tree to
 * a species tree.
 ****************************************************************************/
void recoverILS(Node & node, int & a, int & olda, 
                std::vector <std::vector <unsigned int> > &vec);


/*****************************************************************************
 * Computes the likelihood using our coalescence model, 
 * given a vector of vector giving, for each branch of the species tree,
 * the number of incoming lineages, and the number of outgoing lineages.
 * Formula from Degnan and Salter (2005), Evolution 59(1), pp. 24-37.
 * 3 versions of the function:
 * - working for lots of gene families at once
 * - working for one gene family
 * - working for one branch of one gene family
 ****************************************************************************/

double computeCoalLikelihood (std::vector < std::vector<std::vector<unsigned int> > > vec, std::vector < double > CoalBl ) ;
double computeCoalLikelihood (std::vector < std::vector<unsigned int> > vec, std::vector < double > CoalBl ) ;
double computeCoalLikelihood ( std::vector<unsigned int>  vec, double CoalBl ) ;


/*****************************************************************************
 * This function performs a preorder tree traversal in order to fill vectors of counts. 
 * When used after the postorder tree traversal function, counts for all rootings are computed.
 * coalCounts contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalCountsPreorder(TreeTemplate<Node> & spTree, 
                                      TreeTemplate<Node> & geneTree, 
                                      Node * node, 
                                      const std::map<std::string, std::string > & seqSp, 
                                      const std::map<std::string, int > & spID, 
                                      std::vector < std::vector< std::vector< std::vector< unsigned int > > > > & coalCounts, 
                                      std::vector<double> & bls, 
                                      std::vector <std::vector<unsigned int> > & speciesIDs, 
                                      int sonNumber, 
                                      std::map <double, Node*> & LksToNodes);


/*****************************************************************************
 * This function computes the Coalescent counts of a rooting. 
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 * 
 ****************************************************************************/
void computeRootingCoalCounts(TreeTemplate<Node> & spTree, 
                              Node * node, 
                              std::vector < std::vector< std::vector< std::vector< unsigned int > > > > & coalCounts, 
                              const std::vector< double> & bls, 
                              std::vector <std::vector<unsigned int> > & speciesIDs, 
                              int sonNumber, 
                              std::map <double, Node*> & LksToNodes) ;


/*****************************************************************************
 * This class maximizes the Coalescent likelihood for a branch, by optimizing
 * the branch length in coalescent units.
 ****************************************************************************/

class CoalBranchLikelihood :
public Function,
public AbstractParametrizable
{
protected:
    std::vector <std::vector <unsigned int> > vec_;
    std::vector <std::vector <unsigned int> > compressedVec_;
    std::vector <double> lks_;
    std::map <string, unsigned int> patternToWeights_; 
    double lnL_;
    
public:
    CoalBranchLikelihood(const std::vector <std::vector <unsigned int> >& vec) :
    AbstractParametrizable(""),
    vec_(vec), compressedVec_(0), lnL_(log(0.)), patternToWeights_(), lks_(0)
    {
        Parameter p("BrLen", 1, 0);
        addParameter_(p);
    }
    
    CoalBranchLikelihood(const CoalBranchLikelihood& bl) :
    AbstractParametrizable(bl),
    vec_(bl.vec_), compressedVec_(bl.compressedVec_), lnL_(bl.lnL_),
    patternToWeights_(bl.patternToWeights_), lks_(0)
    {}
    
    CoalBranchLikelihood& operator=(const CoalBranchLikelihood& bl)
    {
        AbstractParametrizable::operator=(bl);
        vec_ = bl.vec_;
        compressedVec_ = bl.compressedVec_;
        lnL_ = bl.lnL_;
        patternToWeights_ = bl.patternToWeights_;
        lks_ = bl.lks_;
        return *this;
    }
    
    virtual ~CoalBranchLikelihood() {}
    
    CoalBranchLikelihood* clone() const { return new CoalBranchLikelihood(*this); }
    
public:
    double initModel();
    
    /**
     * @warning No checking on alphabet size or number of rate classes is performed,
     * use with care!
     */
/*
 void initLikelihoods(const VVVdouble *array1, const VVVdouble *array2)
    {
        _array1 = array1;
        _array2 = array2;
    }
    
    void resetLikelihoods()
    {
        _array1 = 0;
        _array2 = 0;
    }
    */
    void setParameters(const ParameterList &parameters)
    throw (ParameterNotFoundException, ConstraintException)
    {
        setParametersValues(parameters);
    }
    
    double getValue() const throw (Exception) { return lnL_; }
    
    void fireParameterChanged(const ParameterList & parameters)
    {
       // computeAllTransitionProbabilities();
        computeLogLikelihood();
    }
    
protected:
   // void computeAllTransitionProbabilities();
    void computeLogLikelihood();
};



#endif
