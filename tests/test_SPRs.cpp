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


using namespace std;

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

#include <Bpp/Numeric/Random/RandomTools.h>

#include "../src/GenericTreeExplorationAlgorithms.h"


using namespace bpp;

/******************************************************************************/


// performs tons of SPRs on a file and sees if the number of node changes

int main(int argc, char** argv)
{
  
  if (argc !=2)
  {
    cout << "Usage " << argv[0] << " treeFile.fasta" << endl; 
    return 0;
  }
  
  try
  {
    // makespr
    Newick nw;
    TreeTemplate<Node>* tree = nw.read(string(argv[1]));
    TreeTemplate<Node>* treeOrig = tree->clone();
    set<int> origNodesID;
    vector<int> origNodesIDVect = tree->getNodesId();
    origNodesID.insert(origNodesIDVect.begin(),origNodesIDVect.end());
    
    for(unsigned int count = 0; count != 100000; count++){
      unsigned int a,b;
      do
      {
        a = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<int>(tree->getNumberOfNodes()-1);
        b = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<int>(tree->getNumberOfNodes()-1);
      } while (!(tree->getNode(a)->hasFather() && tree->getNode(b)->hasFather()));
      makeSPR(*tree,a,b,false);
      set<int> currNodesID;
      vector<int> currNidV = tree->getNodesId();
      currNodesID.insert(currNidV.begin(),currNidV.end());
      if (origNodesID != currNodesID)
      {
        cout << "Error: lost node when doing " << a << " -> " << b <<endl;
        cout << "Had " << origNodesID.size() << " nodes, now has " << currNodesID.size() << endl;
        delete tree;
        tree = treeOrig->clone();
      }
    }
    
    
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  
  return 0;
}

