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

#include "../src/LikelihoodEvaluator.h"


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
    
    LikelihoodEvaluator le(bppml.getParams());
    le.initialize();
    
    cout << "Computed Likelihood = " << le.getLogLikelihood() << endl;
    
//     // Displaying resulting tree
//     Tree2String(tr->tree_string, tr, partitions, tr->start->back, true, true, 0, 0, 0, true, 0,0);
//     
//     // cout << "Tree: \n" << tr->tree_string << endl;
    
    
    bppml.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  
  return 0;
}

