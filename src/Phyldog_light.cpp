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

#include <iostream>
#include <iomanip>
#include <algorithm>
/*
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;
*/

#include "GeneTreeAlgorithms.h"



/*
This program takes a species tree file, sequence alignments and files describing the relations between sequences and species.
To compile it, use the Makefile !
*/



/******************************************************************************/

TreeTemplate<Node> *  getTreeOrTreesFromOptions(map <string,string> params, vector<Tree*> trees, bool cont) {

  TreeTemplate<Node> *  rootedTree = 00;
  // Get the initial gene tree
  std::string geneTree_File =ApplicationTools::getStringParameter ( "gene.tree.file",params,"none" );
  if ( geneTree_File=="none" ) {
    std::cout << "\n\nNo Gene tree was provided. The option gene.tree.file must be filled with the path of a valid tree file. \n\n" << std::endl;
  }
  if ( !FileTools::fileExists ( geneTree_File ) ) {
    std::cerr << "Error: geneTree_File "<< geneTree_File <<" not found. The option gene.tree.file must be filled with the path of a valid tree file." << std::endl;
    cont = false;
    return rootedTree;
  }
  if ( rootedTree ) {
    delete rootedTree;
    rootedTree =0;
  }
  IMultiTree* treeReader;
  treeReader = new Newick(true);
  treeReader->read(geneTree_File, trees);
  delete treeReader;

  rootedTree = new TreeTemplate<Node> (*(trees[0]->clone() ) );

  if ( !rootedTree->isRooted() ) {
    rootedTree->newOutGroup ( 0 );
  }

  ApplicationTools::displayResult ( "Gene Tree file", geneTree_File );
  ApplicationTools::displayResult ( "Number of leaves", TextTools::toString ( rootedTree->getNumberOfLeaves() ) );
  cont = true;
  return rootedTree;
}


/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "phyldog parameter1_name=parameter1_value parameter2_name=parameter2_value"   ).endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message << "Example of some options: ").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "species.tree.file                    | path to a rooted species tree, Newick format" ).endLine();
  (*ApplicationTools::message << "branch.expected.numbers.optimization | average, branchwise, average_then_branchwise or no: how we optimize duplication and loss parameters").endLine();
  //(*ApplicationTools::message << "init.gene.tree=user # Starting gene tree(s). Could be "user" or "bionj". "user" requires that at least one user-input tree is given with the "gene.tree.file" option, whereas the option "bionj" makes phyldog use this algorithm to create a starting gene tree.").endLine();
  (*ApplicationTools::message << "gene.tree.file                       | path to a gene tree, Newick format").endLine();
  (*ApplicationTools::message << "output.tree.file                     | path where to write the end gene tree").endLine();
  (*ApplicationTools::message << "taxaseq.file=file.link               | file giving the link between species and sequence names (more on this below)").endLine();
  (*ApplicationTools::message << "input.sequence.file=file.fasta       | file giving the input sequence alignment for the gene family").endLine();
  (*ApplicationTools::message << "alphabet=DNA                         | could also be 'RNA', 'protein', or 'Codon'.").endLine();


  (*ApplicationTools::message << "input.sequence.format=Fasta          | format of the sequence alignment. Could be Fasta, Phylip, Clustal, Mase, Nexus... Please see the bppsuite help for more details.").endLine();

  (*ApplicationTools::message << "output.file=basename                 | base name used to output the reconstructed gene trees. \n\tThe program will output basename_reconciled.tree, in NHX format, where duplication and speciation nodes are annotated, with the tag 'Ev=D' or 'Ev=S' respectively. \n\tIt will also output basename_events.txt, where events of duplication and loss are written, along with the species ID information. The format is one event per line, with an event described as: event(SpeciesID, 'FamilyName', duplication|loss). \n\tIt will also output basename_orthologs.txt, where orthologs and paralogs are written. The format is one orthology/paralogy relationship per line, with first the family name, the type of relationship (Orthology or paralogy) then a list of genes, '<===>' and the other series of genes that are in relationship to the first ones.").endLine();

  (*ApplicationTools::message << "output.numbered.tree.file=file       |  file where the species tree topology is saved, annotated with node indices. No file is produced if this option is not provided.").endLine();

  (*ApplicationTools::message << "input.sequence.sites_to_use=all      |  tells whether we should use all sites in the alignment or not. Could be 'all', 'nogap', or 'complete'. Please see the bppsuite help for more details.").endLine();

  (*ApplicationTools::message << "input.sequence.max_gap_allowed=100%  |  maximum number of gaps tolerated for including a site in the analysis.").endLine();
  (*ApplicationTools::message << "genome.coverage.file                 | file giving the percent coverage of the genomes used").endLine();
  (*ApplicationTools::message << "spr.limit                            | integer giving the breadth of SPR movements, in number of nodes. 0.1* number of nodes in the gene tree might be OK.").endLine();
  (*ApplicationTools::message << "model=GTR                            | options of the model used. Should match the alphabet. Starting values may be provided, for instance: model=GTR(a=1.17322, b=0.27717, c=0.279888, d=0.41831, e=0.344783, initFreqs=observed, initFreqs.observedPseudoCount=1). Please see the bppsuite help for more details. ").endLine();
  (*ApplicationTools::message << "rate_distribution=Gamma(n=4)         | Rate heterogeneity among sites option. Here we assume a gamma law with 4 categories. If we add a category of invariants: ate_distribution=Invariant(dist=Gamma(n=4,alpha=1.0), p=0.1)").endLine();

  //######## Finally, optimization options ########

  (*ApplicationTools::message << "optimization.topology=yes | if we choose to optimize the gene tree topology, 'no' otherwise.").endLine();
  //  (*ApplicationTools::message << "optimization.ignore_parameter=parameters | Parameters that we choose not to optimize these 10 parameters in order to save computing time, as we have provided reasonable input values. However, in cases where good input values are not available, it may be wise to leave this field empty and optimize these parameters.").endLine();



  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}


/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
//////////////////////////////////////////////////////MAIN/////////////////////////////////////////////////
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


int main(int args, char ** argv)
{
  WHEREAMI( __FILE__ , __LINE__ );
  if(args == 1)
  {
    help();
    exit(0);
  }


  try {
    ApplicationTools::startTimer();
    //All processors parse the main options
    std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);
    double scenarioLikelihood_ = 0.0;

    // getting the species tree
    std::string spTreeFile =ApplicationTools::getStringParameter("species.tree.file",params,"none");
    ApplicationTools::displayResult("Species Tree file", spTreeFile);
    Newick newick(true);
    TreeTemplate<Node> * rootedSpeciesTree_ = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(rootedSpeciesTree_->getNumberOfLeaves()));
    breadthFirstreNumber ( *rootedSpeciesTree_ );

    // Getting the output file name
    std::string outputFile =ApplicationTools::getStringParameter("output.file",params,"none", "", false, 1);



    // Now we start the gene family side
    string alnFile = ApplicationTools::getStringParameter ( "input.sequence.file",params,"none" );
    LikelihoodEvaluator * levaluator_;
    levaluator_ = new LikelihoodEvaluator(params);


    //TODO: dirty cont to eliminate
    bool cont = true;

    std::map<std::string, std::deque<std::string> > spSeq;
    std::map <std::string, std::string> seqSp;
    getCorrespondanceSequenceSpeciesFromOptions(params, cont, seqSp, spSeq );

    if (!cont)
    throw(Exception("Unable to load this family"));
    removeUselessSequencesFromAlignment( rootedSpeciesTree_, levaluator_->getSites(), cont , spSeq, alnFile) ;

    TreeTemplate<Node> * rootedTree_;
    if (cont) {
      /****************************************************************************
      * Then we need to get the file containing the gene tree,
      * or build the gene tree.
      *****************************************************************************/
      rootedTree_ = getTreeFromOptions(params, levaluator_->getAlphabet(), levaluator_->getSites(), levaluator_->getSubstitutionModel(), levaluator_->getRateDistribution(), cont);
    }

    if (cont) {
      // set the levaluator tree to the modified one
      levaluator_->setTree(rootedTree_);
    }

    levaluator_->initialize();

    if (!rootedTree_->isRooted()) {
      std::cout <<"gene tree is not rooted!!!\n "<< TreeTemplateTools::treeToParenthesis(*rootedTree_) <<std::endl;
    }

    size_t sprLimitGeneTree = ApplicationTools::getIntParameter("SPR.limit", params, 3, "", false, false);



    // getting the gene tree(s), in case there are several.
    //Here we check whether there are several trees inside the gene tree file (if init.gene.tree was user).
    //If there are several trees, then we test each of them, and we will choose the tree
    //with the highest total likelihood as a starting point.
    vector<Tree*> trees;
    std::string geneTree_File =ApplicationTools::getStringParameter ( "gene.tree.file",params,"none" );
    IMultiTree* treeReader;
    treeReader = new Newick(true);
    treeReader->read(geneTree_File, trees);
    delete treeReader;

    std::map <std::string, int>  spId = computeSpeciesNamesToIdsMap ( *rootedSpeciesTree_ );

    // initialize vectors
    size_t numNodes = rootedSpeciesTree_->getNumberOfNodes();
    std::vector<double> lossExpectedNumbers (numNodes, 0.3141514799);
    std::vector<double> duplicationExpectedNumbers (numNodes, 0.0016105081);
    std::vector<int> num0Lineages (numNodes,0);
    std::vector<int> num1Lineages (numNodes,0);
    std::vector<int> num2Lineages (numNodes,0);
    int MLindex  ;
    double bestSequenceLogL = UNLIKELY;
    double scenarioLikelihood = UNLIKELY;
    std::set<int> nodesToTryInNNISearch;
    if (trees.size() > 1) {
      //We have several trees, we need to choose which one is the best!
      double bestSequenceLogL;
      findBestGeneTreeAmongSeveralCandidates(trees, rootedTree_, bestSequenceLogL, scenarioLikelihood, rootedSpeciesTree_, seqSp, spId, lossExpectedNumbers,
      duplicationExpectedNumbers, MLindex,
      num0Lineages, num1Lineages,
      num2Lineages, nodesToTryInNNISearch, true, levaluator_);
    }
    else {
      scenarioLikelihood = findMLReconciliationDR (rootedSpeciesTree_, rootedTree_,
        seqSp, spId, lossExpectedNumbers,
        duplicationExpectedNumbers, MLindex,
        num0Lineages, num1Lineages,
        num2Lineages, nodesToTryInNNISearch, false);
      }
      std::cout << "\n\t\tTotal initial logLikelihood value: "<< levaluator_->getLogLikelihood() +  scenarioLikelihood << "\n\t\tSequence loglk: "<< levaluator_->getLogLikelihood()<<" and scenario loglk: "<< scenarioLikelihood<< std::endl;

      //Rooting the tree properly:
      vector<Node*> nodes = rootedTree_->getNodes();
      for (unsigned int j = 0 ; j < nodes.size() ; j++) {

        if (nodes[j]->hasNodeProperty("outgroupNode")) {

          if (rootedTree_->getRootNode() == nodes[j]) {
            if (j < nodes.size()-1)
            {
              rootedTree_->rootAt(nodes[nodes.size()-1]);
            }
            else {
              rootedTree_->rootAt(nodes[nodes.size()-2]);
            }
          };
          rootedTree_->newOutGroup( nodes[j] );
          break;
        }
      }
      MLindex = -1;

      // Now we have the best rooted gene tree
      // We can start gene tree exploration to improve the gene tree.
      string temp ;
      temp = outputFile+"_reconciled.tree";

      refineGeneTreeWithSPRsFast2 (params, rootedSpeciesTree_, rootedTree_,
        seqSp, spId,
        lossExpectedNumbers,
        duplicationExpectedNumbers,
        MLindex,
        num0Lineages, num1Lineages, num2Lineages, nodesToTryInNNISearch,
        sprLimitGeneTree, levaluator_, scenarioLikelihood, temp);

      std::cout << "\n\t\tTotal final logLikelihood value: "<< levaluator_->getLogLikelihood() +  scenarioLikelihood << "\n\t\tSequence loglk: "<< levaluator_->getLogLikelihood()<<" and scenario loglk: "<< scenarioLikelihood<< std::endl;


      // Now, outputting the gene tree
      temp = outputFile+"_reconciled.tree";
      writeReconciledGeneTreeToFile ( params, rootedTree_->clone(), rootedSpeciesTree_, seqSp, temp ) ;
      temp = outputFile+"_events.txt";
      outputNumbersOfEventsToFile( params, rootedTree_->clone(), rootedSpeciesTree_, seqSp, outputFile, temp );
      temp =outputFile+"_orthologs.txt";
      outputOrthologousAndParalogousGenesToFile( params, rootedTree_->clone(), rootedSpeciesTree_, seqSp, outputFile, temp ) ;

      std::cout << "PHYLDOG_LIGHT's done. Bye." << std::endl;
      ApplicationTools::displayTime("Total execution time:");

    }
    catch(std::exception & e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return (0);
  }
