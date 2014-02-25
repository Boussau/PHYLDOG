#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll/pll.h>

int main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  pllInstance * tr;
  pllNewickTree * newick;
  partitionList * partitions;
  pllQueue * partitionInfo;
  pllInstanceAttr attr;
  double alpha = 1.0;

#ifdef _FINE_GRAIN_MPI
  pllInitMPI (&argc, &argv);
#endif

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [fasta-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Set the PLL instance attributes */
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_TRUE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = 8;            /* This only affects the pthreads version */

  /* Create a PLL tree */
  tr = pllCreateInstance (&attr);

  /* Parse a PHYLIP file */
  alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, argv[1]);


  if (!alignmentData)
   {
     fprintf (stderr, "Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }
 
  /* Parse a NEWICK file */
  newick = pllNewickParseFile (argv[2]);

  if (!newick)
   {
     fprintf (stderr, "Error while parsing newick file %s\n", argv[2]);
     return (EXIT_FAILURE);
   }
  if (!pllValidateNewick (newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
   {
     fprintf (stderr, "Invalid phylogenetic tree\n");
     printf ("%d\n", errno);
     //return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  partitionInfo = pllPartitionParse (argv[3]);
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (partitionInfo, alignmentData))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* Commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);

  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (alignmentData, partitions);

  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (tr, newick, PLL_FALSE);

  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels); */

  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (tr, alignmentData, partitions, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  
  /* Initialize the model. Note that this function will also perform a full
     tree traversal and evaluate the likelihood of the tree. Therefore, you
     have the guarantee that tr->likelihood the valid likelihood */

  pllInitModel(tr, partitions, alignmentData);
  printf ("Log-likelihood of topology: %f\n", tr->likelihood);

  //pllOptimizeModelParameters(tr, partitions, 0.1);

  printf ("Setting Alpha to %f..\n", alpha);
  pllSetFixedAlpha(alpha, 0, partitions, tr);
  double f[4] = { 0.25, 0.25, 0.25, 0.25 };
  printf ("Setting frequencies to 0.25 each..\n");
  pllSetFixedBaseFrequencies(f, 4, 0, partitions, tr);

  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

  printf ("Log-likelihood of topology: %f\n", tr->likelihood);
  
  return (1);
}
