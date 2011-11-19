//
//  AugmentedTreeLikelihood.cpp
//  ReconcileDuplications.proj
//
//  Created by Bastien Boussau on 15/11/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#include "AugmentedTreeLikelihood.h"

using namespace bpp;
using namespace std;


// From the STL:
#include <iostream>



/*******************************************************************************/

void AugmentedTreeLikelihood::initModel()
{
    std::cout <<"initModel 0"<<std::endl;
    nbStates_ = eqFreqs_.size();
    nbClasses_  = rDist_->getNumberOfCategories();

    double temp = 0;
    logOneOverFactorialNIJ_ = 0;

    numberOfSubstitutionsPerSite_.resize(mapping_.getNumberOfSites(), 0);
    numberOfSubstitutionsPerBranch_.resize(mapping_.getNumberOfBranches(), 0);

    bls_.resize(mapping_.getNumberOfBranches(), 0);
    
    vector <int > nodeIds = mapping_.getTree().getNodesId();
    nodeIds.pop_back(); //remove root id.

   // std::cout <<"initModel 4: "<< mapping_.getNumberOfSites() <<" ; "<<  mapping_.getNumberOfBranches()<<std::endl;

    for (unsigned int i = 0 ; i < mapping_.getNumberOfSites() ; i++) {
        for (unsigned int j = 0 ; j < nodeIds.size() ; j++) {
          /*  std::cout <<"initModel 5: "<< mapping_(2,2, 1)<<std::endl;
            std::cout <<" again initModel 5: "<< mapping_(2,2, 2)<<std::endl;

            std::cout <<"initModel 5: "<< mapping_.getNodeIndex(j)<<std::endl;
            
            mapping_.getNumberOfSubstitutions(mapping_.getNodeIndex(j), i);
            std::cout <<"initModel 501"<<std::endl;
            VectorTools::sum(mapping_.getNumberOfSubstitutions(mapping_.getNodeIndex(j), i));
            std::cout <<"initModel 502"<<std::endl;*/

            temp = VectorTools::sum(mapping_.getNumberOfSubstitutions(nodeIds[j], i));
            numberOfSubstitutionsPerSite_[i] += temp;
            numberOfSubstitutionsPerBranch_[nodeIds[j]] += temp;
            logOneOverFactorialNIJ_ -= log(NumTools::fact(floor( temp + 0.5 )));
        }
    }
    fireParameterChanged(getParameters());
    return;
}

/*******************************************************************************/

void AugmentedTreeLikelihood::computeLogLikelihood()
{
    //We compute a "fake" likelihood, with only the parts for which branch lengths matter
    //So we drop: (\product_{a=1}^{a=O} pi(a)^(w_a)) 
    //Similarly, we consider the rates-across-sites distribution to be fixed.
    //For the moment, we do not use the rates-across-sites distribution in the formula; 
    //we might not need it as the mapping has been done averaging over the distribution.
    lnL_ = logProductBlNumSubst_ + logOneOverFactorialNIJ_;
   // for (unsigned int i = 0 ; i < mapping_.getNumberOfSites() ; i++) {
      //  for (unsigned int r = 0 ; r < nbClasses_ ; r++) {
        //    lnL += log(rDist_->getCategory(r)) * numberOfSubstitutionsPerSite_[i] + ;            
     //   }
  //  }
    double sumBranches = 0;
    for (unsigned int j = 0 ; j < mapping_.getNumberOfBranches() ; j++) {
        sumBranches -= bls_[j];
    }
    lnL_ +=  sumBranches * mapping_.getNumberOfSites();

   // p(mapping|t, R, l, pi) = (\product_j l_j^v_j) [ \sum_{r=r_1}^{r=r_R} (\product_s r^n) (\product_{ij} exp(-rl_j)/(n_{ij}!)) ]
    lnL_ = -lnL_;
    
}

/*******************************************************************************/

void AugmentedTreeLikelihood::updateBls() {
    logProductBlNumSubst_ = 0;
    for (unsigned int j = 0 ; j < mapping_.getNumberOfBranches() ; j++) {
        bls_[j] = getParameterValue("BrLen" + TextTools::toString(j));
        logProductBlNumSubst_ += numberOfSubstitutionsPerBranch_[j] * log(bls_[j]);
    }
}
