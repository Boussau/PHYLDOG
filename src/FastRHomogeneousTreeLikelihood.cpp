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

#include "FastRHomogeneousTreeLikelihood.h"
#include <Bpp/Phyl/PatternTools.h>

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

FastRHomogeneousTreeLikelihood::FastRHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_(usePatterns);
}

/******************************************************************************/

FastRHomogeneousTreeLikelihood::FastRHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_ = new DRASRTreeLikelihoodData(
    tree_,
    rateDistribution_->getNumberOfCategories(),
    usePatterns);
}

/******************************************************************************/

FastRHomogeneousTreeLikelihood::FastRHomogeneousTreeLikelihood(
  const FastRHomogeneousTreeLikelihood& lik) :
  AbstractHomogeneousTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
}

/******************************************************************************/

FastRHomogeneousTreeLikelihood& FastRHomogeneousTreeLikelihood::operator=(
  const FastRHomogeneousTreeLikelihood& lik)
{
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if (likelihoodData_) delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

FastRHomogeneousTreeLikelihood::~FastRHomogeneousTreeLikelihood()
{
  delete likelihoodData_;
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  if (data_) delete data_;
  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());
  if (verbose_) ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *model_);
  if (verbose_) ApplicationTools::displayTaskDone();

  nbSites_ = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_ = likelihoodData_->getNumberOfStates();

  if (verbose_) ApplicationTools::displayResult("Number of distinct sites",
                                                TextTools::toString(nbDistinctSites_));
  initialized_ = false;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  for (size_t i = 0; i < nbSites_; i++)
  {
    l *= getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; i++)
  {
    la[i] = getLogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return l;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  if (l < 0) l = 0; //May happen because of numerical errors.
  return log(l);
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    l += (*la)[i] * rootFreqs_[i];
  }
  return l;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    l += (*la)[i] * rootFreqs_[i];
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  return log(l);
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][state];
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][state]);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
      || model_->getParameters().getCommonParametersWith(params).size() > 0)
  {
    //Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if (params.size() > 0)
  {
    //We may save some computations:
    for (size_t i = 0; i < params.size(); i++)
    {
      string s = params[i].getName();
      if (s.substr(0, 5) == "BrLen")
      {
        //Branch length parameter:
        computeTransitionProbabilitiesForNode(nodes_[TextTools::to < size_t > (s.substr(5))]);
      }
    }
    rootFreqs_ = model_->getFrequencies();
  }

  computeTreeLikelihood();

  minusLogLik_ = -getLogLikelihood();

}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized()) throw Exception("FastRHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/

double FastRHomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
  size_t site,
  size_t rateClass) const
{
  double dl = 0;
  Vdouble* dla = &likelihoodData_->getDLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    dl += (*dla)[i] * rootFreqs_[i];
  }
  return dl;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getDLikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    dl += getDLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return dl;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getDLogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++)
  {
    dl += getDLogLikelihoodForASite(i);
  }
  return dl;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("FastRHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<FastRHomogeneousTreeLikelihood*>(this)->computeTreeDLikelihood(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* dpxy__son = &dpxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* dpxy__son_c = &(*dpxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* dpxy__son_c_x = &(*dpxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*dpxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood(father);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == NULL) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _dLikelihoods_son = &likelihoodData_->getDLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _dLikelihoods_son_i = &(*_dLikelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _dLikelihoods_son_i_c = &(*_dLikelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_dLikelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeDLikelihood(father);
}

/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/

double FastRHomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
  size_t site,
  size_t rateClass) const
{
  double d2l = 0;
  Vdouble* d2la = &likelihoodData_->getD2LikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    d2l += (*d2la)[i] * rootFreqs_[i];
  }
  return d2l;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getD2LikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    d2l += getD2LikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return d2l;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++)
  {
    dl += getD2LogLikelihoodForASite(i);
  }
  return dl;
}

/******************************************************************************/

double FastRHomogeneousTreeLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("FastRHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<FastRHomogeneousTreeLikelihood*>(this)->computeTreeD2Likelihood(variable);
  return -getD2LogLikelihood();
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* d2pxy__son = &d2pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* d2pxy__son_c = &(*d2pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* d2pxy__son_c_x = &(*d2pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*d2pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == NULL) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _d2Likelihoods_son = &likelihoodData_->getD2LikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _d2Likelihoods_son_i = &(*_d2Likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _d2Likelihoods_son_i_c = &(*_d2Likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_c_x)[y] * (*_d2Likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihood(tree_->getRootNode());
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
    if (node->isLeaf()) return;
   
    if ( (! (node->hasNodeProperty("toComp") && ( (dynamic_cast<const BppString *>(node->getNodeProperty("toComp") ))->toSTL() == "N" ) ) ) )
    {

        size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();
        size_t nbNodes = node->getNumberOfSons();
        
        // Must reset the likelihood array first (i.e. set all of them to 1):
        VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
# ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < nbSites; i++)
        {
            //For each site in the sequence,
            VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
            for (size_t c = 0; c < nbClasses_; c++)
            {
                //For each rate classe,
                Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                for (size_t x = 0; x < nbStates_; x++)
                {
                    //For each initial state,
                    (*_likelihoods_node_i_c)[x] = 1.;
                }
            }
        }
        
        for (size_t l = 0; l < nbNodes; l++)
        {
            //For each son node,
            
            const Node* son = node->getSon(l);
            
            computeSubtreeLikelihood(son); //Recursive method:
            VVVdouble* pxy__son = &pxy_[son->getId()];
            vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
            VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            
# ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < nbSites; i++)
            {
                //For each site in the sequence,
                VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
                VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
                for (size_t c = 0; c < nbClasses_; c++)
                {
                    //For each rate classe,
                    Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                    VVdouble* pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++)
                    {
                        //For each initial state,
                        Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
                        double likelihood = 0;
                        for (size_t y = 0; y < nbStates_; y++)
                        {
                           /* if ( ( (node->hasNodeProperty("toComp") && ( (dynamic_cast<const BppString *>(node->getNodeProperty("toComp") ))->toSTL() == "N" ) ) ) )
                            {
                            std::cout << "(*_likelihoods_son_i_c)[y]: "<< (*_likelihoods_son_i_c)[y]<<std::endl;
                            }*/
                            likelihood += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_likelihoods_node_i_c)[x] *= likelihood;
                    }
                }
            }
        }
    }
}

/******************************************************************************/

void FastRHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

/*******************************************************************************/


void FastRHomogeneousTreeLikelihood::initializeLikelihoodData(const Node* node)
{
    if (node->isLeaf()) return;
        
        size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();
        size_t nbNodes = node->getNumberOfSons();
        
        // Must reset the likelihood array first (i.e. set all of them to 1):
        VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
        for (size_t i = 0; i < nbSites; i++)
        {
            //For each site in the sequence,
            VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
            for (size_t c = 0; c < nbClasses_; c++)
            {
                //For each rate classe,
                Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                for (size_t x = 0; x < nbStates_; x++)
                {
                    //For each initial state,
                    (*_likelihoods_node_i_c)[x] = 1.;
                }
            }
        }
        
        for (size_t l = 0; l < nbNodes; l++)
        {
            //For each son node,
            
            const Node* son = node->getSon(l);
            
            initializeLikelihoodData(son); //Recursive method:
            
            VVVdouble* pxy__son = &pxy_[son->getId()];
            vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
            VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            
            for (size_t i = 0; i < nbSites; i++)
            {
                //For each site in the sequence,
                VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
                VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
                for (size_t c = 0; c < nbClasses_; c++)
                {
                    //For each rate classe,
                    Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                    VVdouble* pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++)
                    {
                        //For each initial state,
                        Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
                        double likelihood = 0;
                        for (size_t y = 0; y < nbStates_; y++)
                        {
                            likelihood += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_likelihoods_node_i_c)[x] *= likelihood;
                    }
                }
            }
        }
}

/*******************************************************************************/


void FastRHomogeneousTreeLikelihood::initializeLikelihoodData() {
    initializeLikelihoodData(tree_->getRootNode());
    minusLogLik_ = -getLogLikelihood();
}

