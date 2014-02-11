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


#ifndef LikelihoodWrapper_hpp
#define LikelihoodWrapper_hpp

#include<string>

extern "C" {
#include <pll/pll.h>
}



Class LikelihoodWrapper {
private:
  std::string treeFile;
  std::string alignmentFile;
  
  
  /** @name Data and commands for PLL
  *  PLL specific objects and routine calls.
  */
  ///@{
  
  /**
  * PLL attributes. Example:
  * attr.rateHetModel     = PLL_GAMMA;
    attr.fastScaling      = PLL_FALSE;
    attr.saveMemory       = PLL_FALSE;
    attr.useRecom         = PLL_FALSE;
    attr.randomNumberSeed = 0xDEADBEEF;
  */
  pllInstanceAttr attr_PLL;
  
  /**
  The PLL tree.
  */
  pllInstance * tr_PLL;
  
  /**
  The PLL alignment.
  */
  pllAlignmentData * alignmentData_PLL;
  
  /**
  The PLL newick data.
  */
  pllNewickTree * newick_PLL;
  
  /**
  The PLL partition information.
  */  
  pllQueue * partitionInfo_PLL;
    
  /**
  Loads the PLL tree.
  */
  loadPLLtree();
  ///@}
  
  
  /** @name Translation data
  *  Data used for conversion purposes BPP <-> PLL
  */
  ///@{
  
  /**
  Alignment strictly formatted for PLL
  */
  string fastaForPLL;
  
  /**
  Newick strictly formatted for PLL
  */
  string newickForPLL;
  
  /**
  Defines the real sequences names to simplified ones for PLL
  */
  std::map<string,string> realToStrict;
  
  /**
  Defines simplified names for PLL to real sequence ones
  */
  std::map<string,string> strictToReal;

  
  /**
  Loads the PLL tree.
  */
  
  ///@}
  
  
  
public:
  
  LikelihoodWrapper(std::string treeFile, std::string alignmentFile);
  
  double getLikelihood();
  
  
};

#else

Class LikelihoodWrapper;

#endif