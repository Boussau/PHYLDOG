//
//  AugmentedTreeLikelihood.h
//  ReconcileDuplications.proj
//
//  Created by Bastien Boussau on 15/11/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#ifndef AUGMENTEDTREELIKELIHOOD_H
#define AUGMENTEDTREELIKELIHOOD_H

/*
 Formulas according to Lartillot 2006:
 Likelihood of a mapping for a site s along a branch of length l and a rate r and equilibrium frequencies pi (Equation 32):
 p(mapping(s)|l, r, pi) = exp(-r*l) * (rl)^n/(n!) \product_{k=1}^n pi(sigma_k)
 where n is the total number of substitutions along the branch for the site s, and sigma_k is the identity of the character at step k of the substitution mapping.
 For R rates r_1, r_2... r_R:
 p(mapping(s)|l, R, pi) = \sum_{r=r_1}^{r=r_R} p(mapping(s)|l, r, pi)
 For S sites, R rates, one branch: 
 p(mapping|l, R, pi) = \product_S p(mapping(s)|l, R, pi)
 
 Total tree likelihood for a site s, a tree t, a rate r, a vector of branch lengths l and a vector of equilibrium frequencies pi on an alphabet of size O (Equation 34):
 p(mapping|t, r, l, pi) = (\product_{a=1}^{a=O} pi(a)^(w_a)) (\product_s r^n) (\product_j l_j^v_j) (\product_{ij} exp(-rl_j)/(n_{ij}!))
 
 where w_a is the number of substitutions to a + the number of times state a is seen at the root of the tree, 
 n is the total number of substitutions on the tree for site s
 l_j is the branch length of branch j, and v_j is the total number of substitutions for all sites on branch j
 n_ij is the number of substitutions for site i and branch j.
 
 If we average over rates of evolution (full augmented likelihood of the tree):
 p(mapping|t, R, l, pi) = (\product_{a=1}^{a=O} pi(a)^(w_a)) (\product_j l_j^v_j) [ \sum_{r=r_1}^{r=r_R} (\product_s r^n) (\product_{ij} exp(-rl_j)/(n_{ij}!)) ]

 
 */

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Mapping.all>

namespace bpp 
{
    class AugmentedTreeLikelihood :
        public Function,
        public AbstractParametrizable
    {
    protected:
        ProbabilisticSubstitutionMapping mapping_;
        DiscreteDistribution* rDist_;
        vector <int> numberOfSubstitutionsPerBranch_;
        vector <int> numberOfSubstitutionsPerSite_;
        vector<double> eqFreqs_;
        vector<double> bls_;
        unsigned int nbClasses_;
        unsigned int nbStates_;
        double logProductBlNumSubst_;
        double logOneOverFactorialNIJ_;
        double lnL_;

    public:
        AugmentedTreeLikelihood(ProbabilisticSubstitutionMapping& mapping, vector<double> eqFreqs, DiscreteDistribution* rDist, ParameterList &parameters) :
        AbstractParametrizable(""),
        eqFreqs_(eqFreqs), bls_(0), rDist_(rDist), nbClasses_(0),
        lnL_(log(0.)), mapping_(mapping)
        {
            addParameters_(parameters);
        }
        
        AugmentedTreeLikelihood(const AugmentedTreeLikelihood& al) :
        AbstractParametrizable(al),
        eqFreqs_(al.eqFreqs_), bls_(al.bls_), rDist_(al.rDist_),
        lnL_(al.lnL_), mapping_(al.mapping_)
        {}
        
        AugmentedTreeLikelihood& operator=(const AugmentedTreeLikelihood& al)
        {
            AbstractParametrizable::operator=(al);
            eqFreqs_ = al.eqFreqs_;
            rDist_ = al.rDist_,
            bls_ = al.bls_;
            nbClasses_ = al.nbClasses_;
            lnL_ = al.lnL_;
            mapping_ = al.mapping_;
            return *this;
        }
        
        virtual ~AugmentedTreeLikelihood() {}
        
        AugmentedTreeLikelihood* clone() const { return new AugmentedTreeLikelihood(*this); }
        
    public:
        void initModel();
                
        void setParameters(const ParameterList &parameters)
        throw (ParameterNotFoundException, ConstraintException)
        {
            setParametersValues(parameters);
        }
        
        double getValue() const throw (Exception) { return lnL_; }
        
        void fireParameterChanged(const ParameterList & parameters)
        {
            std::cout <<"HEHEH "<<std::endl;
            updateBls();
            std::cout <<"HEHEH 2"<<std::endl;

            computeLogLikelihood();
            std::cout <<"HEHEH 3"<<std::endl;

        }
        
    protected:
        void updateBls();
        void computeLogLikelihood();
        
    };
    
}

#endif //AUGMENTEDTREELIKELIHOOD_H
