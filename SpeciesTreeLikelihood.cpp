/*
 *  SpeciesTreeLikelihood.cpp
 *  
 *
 *  Created by boussau on 18/02/10.
 *  Copyright 2010 UC Berkeley. All rights reserved.
 *
 */

#include "SpeciesTreeLikelihood.h"
#include "SpeciesTreeExploration.h"

using namespace bpp;

using namespace std;

/*******************************************************************************/

void SpeciesTreeLikelihood::updateDuplicationAndLossRates() {
  double d = getParameterValue("coefDup");
  duplicationProbabilities_ = backupDuplicationProbabilities_ * d;
  double l = getParameterValue("coefLoss");
  lossProbabilities_ = backupLossProbabilities_ * d;
}


/*******************************************************************************/
void SpeciesTreeLikelihood::computeLogLikelihood() {
  computeSpeciesTreeLikelihood(world_, index_,  stop_, logL_, num0Lineages_, 
                               num1Lineages_, num2Lineages_, 
                               allNum0Lineages_, allNum1Lineages_, 
                               allNum2Lineages_, lossProbabilities_, 
                               duplicationProbabilities_, rearrange_, server_, 
                               branchProbaOptimization_, genomeMissing_, *tree_);
}
/*******************************************************************************/

