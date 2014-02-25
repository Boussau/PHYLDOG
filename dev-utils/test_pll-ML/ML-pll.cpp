//
// File: bppML.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
//

/*
 *   Copyright or © or Copr. Bio++ Development Team
 * 
 *   This software is a computer program whose purpose is to estimate
 *   phylogenies and evolutionary parameters from a dataset according to
 *   the maximum likelihood principle.
 * 
 *   This software is governed by the CeCILL  license under French law and
 *   abiding by the rules of distribution of free software.  You can  use,
 *   modify and/ or redistribute the software under the terms of the CeCILL
 *   license as circulated by CEA, CNRS and INRIA at the following URL
 *   "http://www.cecill.info".
 * 
 *   As a counterpart to the access to the source code and  rights to copy,
 *   modify and redistribute granted by the license, users are provided only
 *   with a limited warranty  and the software's author,  the holder of the
 *   economic rights,  and the successive licensors  have only  limited
 *   liability.
 * 
 *   In this respect, the user's attention is drawn to the risks associated
 *   with loading,  using,  modifying and/or developing or reproducing the
 *   software by the user in light of its specific status of free software,
 *   that may mean  that it is complicated to manipulate,  and  that  also
 *   therefore means  that it is reserved for developers  and  experienced
 *   professionals having in-depth computer knowledge. Users are therefore
 *   encouraged to load and test the software's suitability as regards their
 *   requirements in conditions enabling the security of their systems and/or
 *   data to be ensured and,  more generally, to use and operate it in the
 *   same conditions as regards security.
 * 
 *   The fact that you are presently reading this means that you have had
 *   knowledge of the CeCILL license and that you accept its terms.
 */

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

extern "C" {
  #include <pll/pll.h>
}

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppml parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  // cout << "******************************************************************" << endl;
  // cout << "*       Bio++-PLL Maximum Likelihood Computation                 *" << endl;
  // cout << "* TESTING PURPOSE ONLY                                           *" << endl;
  // cout << "******************************************************************" << endl;
  // cout << endl;
  
  if (args == 1)
  {
    help();
    return 0;
  }
  
  try
  {
    BppApplication bppml(args, argv, "BppML");
    bppml.startTimer();
    
    
    // PLL
    /* Set the PLL instance attributes */
    pllInstanceAttr attr;
    attr.rateHetModel     = PLL_GAMMA;
    attr.fastScaling      = PLL_TRUE;
    attr.saveMemory       = PLL_FALSE;
    attr.useRecom         = PLL_FALSE;
    attr.randomNumberSeed = 0xDEADBEEF;
    attr.numberOfThreads  = 8;            /* This only affects the pthreads version */
    
    
    
    pllInstance * tr = pllCreateInstance (&attr);
    
    
    /* Parse a PHYLIP/FASTA file */
    pllAlignmentData * alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, bppml.getParam("input.sequence.file").c_str());
    if (!alignmentData)
    {
      std::cerr << "PLL: Error while parsing " << bppml.getParam("input.sequence.file") << std::endl;
      return (EXIT_FAILURE);
    }
    
    pllNewickTree * newick = pllNewickParseFile(bppml.getParam("input.tree.file").c_str());
    if (!newick)
    {
      std::cerr << "PLL: Error while parsing newick file" << std::endl;
      return (EXIT_FAILURE);
    }
    if (!pllValidateNewick (newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
    {
      cerr << "Invalid phylogenetic tree" << endl;
      cerr << errno << endl;
      return (EXIT_FAILURE);
    }
    
    /* Parse the partitions file into a partition queue structure */
    pllQueue * partitionInfo = pllPartitionParse (bppml.getParam("pll.sites.partition").c_str());
    
    /* Validate the partitions */
    if (!pllPartitionsValidate (partitionInfo, alignmentData))
    {
      cerr << "Error: Partitions do not cover all sites\n";
      return (EXIT_FAILURE);
    }
    // cout << "PLL: Commit the partitions and build a partitions structure" << std::endl ;
    
    /* Commit the partitions and build a partitions structure */
    partitionList * partitions = pllPartitionsCommit (partitionInfo, alignmentData);
    
    // cout << "PLL: We don't need the the intermedia partition queue structure anymore" << std::endl ;
    /* We don't need the the intermedia partition queue structure anymore */
    pllQueuePartitionsDestroy (&partitionInfo);
    
    // cout << "PLL: eliminate duplicate sites from the alignment and update weights vector" << std::endl ;
    /* eliminate duplicate sites from the alignment and update weights vector */
    pllAlignmentRemoveDups (alignmentData, partitions);
    
    // cout << "PLL: Set the topology of the PLL tree from a parsed newick tree" << std::endl;
    /* Set the topology of the PLL tree from a parsed newick tree */
    pllTreeInitTopologyNewick (tr, newick, PLL_FALSE);
    
    // cout << "PLL: Connect the alignment and partition structure with the tree structure" << std::endl ;
    /* Connect the alignment and partition structure with the tree structure */
    if (!pllLoadAlignment (tr, alignmentData, partitions, PLL_DEEP_COPY))
    {
      std::cerr << "Incompatible tree/alignment combination" << std::endl;
      return (EXIT_FAILURE);
    }
    
    // PLL END
    
    for(unsigned int iteration = 0; iteration <1000; iteration++){
    
    /*
     *  //Check initial likelihood:
     *  double logL = tl->getValue();*/
    
    // cout << "PLL: initializing model" << std::endl ;
    pllInitModel(tr, partitions, alignmentData);
    
    
    //using median instead of mean
    tr->useMedian = PLL_FALSE;
    double alpha = 1.0;
    
    pllSetFixedAlpha(alpha, 0, partitions, tr);
    
    
    // cout << "For this experiment, we use " << (tr->useMedian ? "median" : "mean") << "." << endl;
    
    
    // cout << "PLL: Log-likelihood of topology, without any optimization: " << tr->likelihood << std::endl ;
    
    // cout << "FIXING FREQUENCIES" << endl;
    double frequencies[] = {.25, .25, .25, .25};
    pllSetFixedBaseFrequencies(frequencies, 4, 0, partitions, tr);
    
    pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    
    // cout << "PLL: Log-likelihood of topology, without any optimization: " << tr->likelihood << std::endl ;
    
    
    
    // cout << "\n- Gamma parameters" << endl;
    alpha = pllGetAlpha (partitions, 0);
    // cout << "alpha = " << alpha << endl;
    
    double gammaRates[4];
    pllGetGammaRates(partitions, 0, gammaRates);
    // cout << "rates = " << gammaRates[0] << " ; " << gammaRates[1] << " ; " << gammaRates[2] << " ; " << gammaRates[3]   << endl;
    
    
    double baseFreq[4];
    pllGetBaseFrequencies(tr, partitions, 0, baseFreq);
    // cout << "\n Base frequencies\n = " << baseFreq[0] << " ; " << baseFreq[1] << " ; " << baseFreq[2] << " ; " << baseFreq[3]   << endl;
    
    
    // cout << "\n- Subst Matrix (exchangeability):" << endl;
    double subsMatrix[4];
    pllGetSubstitutionMatrix(tr, partitions, 0, subsMatrix);
    // cout << "SubsMatrix = "<< subsMatrix[0] << " ; " << subsMatrix[1] << " ; " << subsMatrix[2] << " ; " << subsMatrix[3]   << endl;
    
    
    //Optimizing br. Length
      // cout << "PLL: optimizing branch length..." << std::endl ;
      pllTreeEvaluate (tr, partitions, 64);
      // cout << "PLL: Log-likelihood of topology " << tr->likelihood << std::endl ;
    
    
//     // Displaying resulting tree
//     Tree2String(tr->tree_string, tr, partitions, tr->start->back, true, true, 0, 0, 0, true, 0,0);
//     
//     // cout << "Tree: \n" << tr->tree_string << endl;
    
    }
    bppml.done();
  }
  catch (exception& e)
  {
    // cout << e.what() << endl;
    return 1;
  }
  
  return 0;
}

