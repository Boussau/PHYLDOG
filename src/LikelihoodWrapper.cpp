


extern "C" {
#include <pll/pll.h>
}

#include<iostream>


#include "LikelihoodWrapper.h"


using namespace std;


LikelihoodWrapper::initializePLL(){
  
  // PLL
  /* Set the PLL instance attributes */
  pllInstanceAttr attr;
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_TRUE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = 8;            /* This only affects the pthreads version */
  
  tr_PLL = pllCreateInstance (&attr);
}


LikelihoodWrapper::LikelihoodWrapper(string treeFile, string alignmentFile):
  treeFile(treeFile),
  alignmentFile(alignmentFile)
{
  // For PLL
  initializePLL();
}


LikelihoodWrapper::loadPLLalignment(char* path)
{
  /* Parse a PHYLIP/FASTA file */
  pllAlignmentData * alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, path);
  if (!alignmentData)
  {
    cerr << "PLL: Error while parsing " << path << std::endl;
    return (EXIT_FAILURE);
  }
}


LikelihoodWrapper::loadPLLtree(char* path)
{
  newick_PLL = pllNewickParseFile(path);
  if (!newick_PLL)
  {
    cerr << "PLL: Error while parsing newick file" << std::endl;
    return (EXIT_FAILURE);
  }
  if (!pllValidateNewick (newick_PLL))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    cerr << "Invalid phylogenetic tree" << endl;
    cerr << errno << endl;
    return (EXIT_FAILURE);
  }
}

LikelihoodWrapper::loadPLLpartitions(char* path)
{
  /* Parse the partitions file into a partition queue structure */
  partitionInfo_PLL = pllPartitionParse (path);
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (partitionInfo_PLL, alignmentData_PLL))
  {
    cerr << "Error: Partitions do not cover all sites\n";
    return (EXIT_FAILURE);
  }
  
  /* Commit the partitions and build a partitions structure */
  partitions_PLL = pllPartitionsCommit (partitionInfo_PLL, alignmentData_PLL);
  
  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo_PLL);
  
  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (alignmentData_PLL, partitions_PLL);
}

