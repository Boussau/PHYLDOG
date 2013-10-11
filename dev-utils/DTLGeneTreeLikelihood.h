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
#ifndef _DTLGENETREELIKELIHOOD_H_
#define _DTLGENETREELIKELIHOOD_H_

// From the STL:
#include <iostream>
#include <iomanip>
#include <algorithm>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/HomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
//#include <Phyl/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/ClockTreeLikelihood.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>


// From NumCalc:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From Utils:
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/Text/KeyvalTools.h>


//#include <Phyl/SAHomogeneousTreeLikelihood.h>
#include "ReconciliationTools.h"
#include "ReconciliationTreeLikelihood.h"
#include "SpeciesTreeExploration.h"
#include "SpeciesTreeLikelihood.h"
#include "GeneTreeAlgorithms.h"



//From the BOOST library 
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

#include "DTL.h"

namespace bpp 
{


  class DTLGeneTreeLikelihood 
  {
  private:
  //Parameter list
  std::map<std::string, std::string> params_;
  //Species tree (used to compute DTL score)
  Species_tree * DTLScoringTree_;
  TreeTemplate<Node> * spTree_;
  //Method to compute the DT score:
  std::string scoringMethod_;
  //Gene tree
  TreeTemplate<Node> * geneTree_;
  TreeTemplate<Node> * unrootedGeneTree_;
  std::string strGeneTree_;
  vector <string> strGeneTrees_;
  //Objects useful to compute the sequence likelihood
  Alphabet * alphabet_ ;
  VectorSiteContainer * sites_;
  SubstitutionModel*    model_ ;
  DiscreteDistribution* rDist_ ;
  //std::map to store the link between sequence and species.
  std::map<std::string, std::string> seqSp_;
  //DTL rates (in the order D, T, L):
  double delta_;
  double tau_;
  double lambda_;
  //Whether we should stop tree search
  bool stop_;
  //Index of the tree under study, and index of the most likely tree
  int index_;
  int bestIndex_;  
  //logLk and best logLk
  double DTLlogL_;
  double bestDTLlogL_;
  double sequencelogL_;
  double bestSequencelogL_;
  //Whether we should rearrange the gene trees
  bool rearrange_;
  //Number of iterations of the search algorithm without improvement
  int numIterationsWithoutImprovement_;
  //How far can we regraft subtrees when doing a spr
  int sprLimit_;
  
  
  public:
  //Simple constructor
  DTLGeneTreeLikelihood(std::map<std::string, std::string> & params) :
  params_(params)
  {
  parseOptions();
  }
  
  //Copy constructor should be done. Ignore for now.
  /*  DTLGeneTreeLikelihood(const DTLGeneTreeLikelihood& gtl):
  geneTree_(gtl.geneTree_), index_(gtl.index_),
  bestIndex_(gtl.bestIndex_), logL_(gtl.logL_), 
  bestlogL_(gtl.bestlogL_), {
  }*/
  
  //= operator should be done. Ignore for now.
  /*  DTLGeneTreeLikelihood& operator=(const DTLGeneTreeLikelihood& gtl):
   geneTree_(gtl.geneTree_), index_(gtl.index_),
   bestIndex_(gtl.bestIndex_), logL_(gtl.logL_), 
   bestlogL_(gtl.bestlogL_), {
   }*/
  
  //Destructor
  virtual ~DTLGeneTreeLikelihood() 
  {	
    if (DTLScoringTree_)
      delete DTLScoringTree_;
    if (geneTree_)
      delete geneTree_;
    if (unrootedGeneTree_)
      delete unrootedGeneTree_;
  }
  
  //Clone function should be done. Ignore for now.
  //  DTLGeneTreeLikelihood* clone() const { return new DTLGeneTreeLikelihood(*this); }

  //Get the logL of the species tree
  double getValue() const throw (Exception) { return DTLlogL_ + sequencelogL_; }
  
  //Computes the DTL likelihood of the current unrootedGeneTree_.
  pair<double, string> computeDTLLikelihood ();
  
  //Computes the sequence likelihood of the current unrootedGeneTree_.
  double optimizeSequenceLikelihood ();
  
  //Initializes various fields in the species tree
  void initialize();
  
  //Does a ML search for the best gene tree using SPRs and NNIs
  void MLSearch();
  
  //Outputs the gene tree
  void printGeneTree();
  
  protected:
  //Computes the loglk of the gene tree
  void computeLogLikelihood();
  //Parses the options and builds the SpeciesTreeLikelihoodObject
  void parseOptions();
  
  
  };

}

#endif //_DTLGENETREELIKELIHOOD_H_

