//
// File: bppML.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team

   This software is a computer program whose purpose is to estimate
   phylogenies and evolutionary parameters from a dataset according to
   the maximum likelihood principle.

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
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Maximum Likelihood Computation, version 1.6.0      *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. 29/01/13 *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Guéguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    BppApplication bppml(args, argv, "BppML");
    bppml.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    auto_ptr<GeneticCode> gCode;

    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams());

    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false);
    delete allSites;

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

    // Get the initial tree
    Tree* tree = 0;
      tree = PhylogeneticsApplicationTools::getTree(bppml.getParams()); //looks for input.tree.file
      ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allow to check file existence before running optimization!
    PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());

    DiscreteRatesAcrossSitesTreeLikelihood* tl;

    bool checkTree    = ApplicationTools::getBooleanParameter("input.tree.check_root", bppml.getParams(), true, "", true, false);
    bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, false);
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, false);

    SubstitutionModel*    model    = 0;
    SubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist    = 0;

      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, bppml.getParams());
      SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());

        string compression = ApplicationTools::getStringParameter("likelihood.recursion_simple.compression", bppml.getParams(), "recursive", "", true, false);
        ApplicationTools::displayResult("Likelihood data compression", compression);
        if (compression == "simple")
          if (dynamic_cast<MixedSubstitutionModel*>(model))
            tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, false);
          else
            tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, false);

        else if (compression == "recursive")
          if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
            tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, true);
          else
            tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, true);

        else throw Exception("Unknown likelihood data compression method: " + compression);

      tl->initialize();

    delete tree;

    //Check initial likelihood:
    double logL = tl->getValue();

    // Write parameters to screen:
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
      tree = new TreeTemplate<Node>(tl->getTree());
    
      //optimizing branch lenth
      auto_ptr<BackupListener> backupListener;
      // TODO: phyldog original: 0.1, try 1 or 10 to speed up, and get closer to PLL results?
      double tolerance = .1;
      unsigned int tlEvalMax = 1000000;
      OutputStream* messageHandler = 0 ; 

      unsigned int numEval = OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (tl), tl->getBranchLengthsParameters(),backupListener.get(), tolerance, tlEvalMax, messageHandler, messageHandler, 0);
      
      
      
      
      // Write resulting tree:
      PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
    
ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));

cout << "\n\nMORE PARAMETERS" << endl;

cout << "getRateDistributionParameters" << endl;
tl->getRateDistributionParameters().printParameters(cout);

cout << "getSubstitutionModelParameters" << endl;
tl->getSubstitutionModelParameters().printParameters(cout);


cout << "\n\nGAMMA LAW POSTERIOR PARAMS" << endl;
  // Getting posterior rate class distribution:
DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl);
ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
ApplicationTools::displayMessage("\n");




    delete alphabet;
    delete sites;
    if (model) delete model;
    if (modelSet) delete modelSet;
    delete rDist;
    delete tl;
    delete tree;
    bppml.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

