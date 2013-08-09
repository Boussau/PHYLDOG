#include <iostream>
#include "DTL.h"
using namespace std;
using namespace bpp;

int main(int argc, char ** argv)
{
  
  //string tree_file_name=argv[1];
  //ifstream tree_file (tree_file_name.c_str());
  //string tree;
  //getline (tree_file,tree);

  // I made a function that computes the DTL prob. for a set of trees given as a vector of newvick strings..
  vector <string> trees;


  //TO see how this works let's look at a simple species tree: 
  string tree="(((a:2,b:2):1,(c:1,d:1):1):1,e:1);";

  // .. but you can look at a bigger one as well, the species tree for 36 cyanobacteria ..
  //tree="(((SYNR3:0.08883,(((PROMM:0.00369,PROM3:0.00596)1.000:0.04156,((PRMAR1:0.05936,PROM4:0.05350)1.000:0.02083,((PROMT:0.00141,PROM1:0.00237)1.000:0.07407,((PROM9:0.01473,(PROMS:0.00750,(PROM0:0.00638,PROM2:0.01282)0.235:0.00224)0.999:0.00686)1.000:0.01910,(PROMP:0.01726,PROM5:0.01678)1.000:0.02574)1.000:0.12115)1.000:0.02174)1.000:0.04000)0.991:0.01529,((SYNPX:0.02959,(SYNS9:0.03639,SYNSC:0.02024)0.747:0.00840)1.000:0.02860,(SYNPW:0.02383,SYNS3:0.04709)0.915:0.01443)0.998:0.01757)1.000:0.05878)1.000:0.21997,(SYNE7:0.00014,SYNP6:0.00093)1.000:0.10878):0.06110,(((ACAM1:0.14206,(THEEB:0.13968,CYAP4:0.09867)0.734:0.02489)1.000:0.04263,((SYNP2:0.14042,((SYNY3:0.11974,(CYAA5:0.07248,CYAP8:0.07026)1.000:0.03407)0.942:0.01982,(MICAN:0.11279,CYAP7:0.08960)0.995:0.02367)1.000:0.03038)1.000:0.06996,(TRIEI:0.16578,(NOSP7:0.04794,(ANASP:0.00449,ANAVT:0.00429)1.000:0.03152)1.000:0.09843)0.954:0.02102)0.147:1.01643)0.988:0.02429,(GLVIO1:0.27430,(SYNJA:0.02794,SYNJB:0.02719)1.000:0.16100):0.08719)0.994:0.00322)11.000;";

  TreeTemplate<Node> * S = TreeTemplateTools::parenthesisToTree(tree);

  // .. we are just going to use the species tree as a gene tree ..
  string gtree="(((((a%1:2,b%2:2):1,(c%3:1,d%4:1):1):1,e%4:1):1,(a%5:2,b%6:2):1):1,c%7:1);";
  //string gtree="(((((a:2,b:2):1,(c:1,d:1):1):1,e:1):1,(a:2,b:2):1):1,c:1);";
  // .. but you could look at also an actual gene tree ..
  //gtree="(((ACAM1:0.12373,ACAM1:0.30013)0.969:0.13815,((SYNP2:0.66664,CYAP4:0.20048)0.505:0.03227,(((SYNY3:0.33669,(ANAVT:0.27812,(NOSP7:0.18865,(CYAA5:0.14079,TRIEI:0.07957)0.947:0.07931)0.605:0.05217)0.216:0.02901)0.058:0.01155,((MICAN:0.60525,(CYAP8:0.40781,(ANASP:0.38404,NOSP7:0.42938)0.94:0.17483)0.523:0.16004)1:0.93419,((ACAM1:0.37346,((ANASP:0.06184,ANAVT:0.03934)0.987:0.11348,(CYAP7:0.26180,TRIEI:0.22849)0.305:0.04874)0.799:0.17081)0.994:0.47097,(ACAM1:0.25626,(CYAP8:0.20323,(ANASP:0.01480,ANAVT:0.02393)1:0.23794)0.78:0.10720)1:1.12944)0.828:0.25448)0.998:0.53287)0.429:0.01465,(MICAN:0.27194,CYAP7:0.22574)0.457:0.03627)0.868:0.06418)0.911:0.07587):0.00815,(((SYNS3:0.17097,SYNPW:0.15676)0.894:0.05326,(SYNS9:0.14098,(SYNPX:0.05225,SYNSC:0.09008)0.869:0.06600)0.559:0.05007)1:0.39930,(SYNP6:0.00000,SYNE7:0.00000)1:0.28627):0.15494)0.909;";
  
  
  //.. some DTL rates ..
  scalar_type delta=0.2;
  scalar_type tau=0.1;
  scalar_type lambda=0.7;


  // To get some trees I take all NNIs on the rooted gene tree (including the ones around the root) .. 
  trees=all_NNIs(gtree);
  // ..unroot these trees, as the DTL calculation works on unrooted trees, and gives a rooted ML reconciliation
  for (vector <string> ::iterator it=trees.begin();it!=trees.end();it++)
    {
      tree_type * G=TreeTemplateTools::parenthesisToTree((*it));
      G->unroot();
      (*it)=TreeTemplateTools::treeToParenthesis(*G);
    }
  
  // We construct a species tree object that will calculate the probability of G-s in the vector trees summed over all roots and reconciliations:
  Species_tree * sum_tree = new Species_tree(S);
  //init_treewise sets up different pieces of the calcultion .. 
  sum_tree->init_treewise(S,delta,tau,lambda,"sum");

  //NB: if the input species tree is not ultra metric init implicitly converts the S tree to an ultrametric tree with an _arbitrary_ time order.
  // This should not be relied upon, but should work..
  cout << " The S tree with branch lengths in units of time such that the tree has height one.." << endl << endl;
  sum_tree->strip_virtual_leaves();
  sum_tree->write(&cout);

  
  cout << " Time orders, which also serve as branch IDs";
  cout << can_order(sum_tree->tree) << endl;
  cout << " All IDs:" << endl;
  for (int i = 0;i<(sum_tree->N_slices+sum_tree->tree->getNumberOfLeaves());i++)
    cout << i <<" "<< sum_tree->branchname[i]<<endl;

  //LL_treewise does the actual calculation and returns a vector of LLs and newvick strings
  cout << endl <<"The probability of G-s in the vector trees summed over all roots and reconciliations:" << endl << endl;
  vector < pair<scalar_type,string> > LLs=sum_tree->LL_treewise(trees);

  for ( int i=0;i<trees.size();i++)
    cout << LLs[i].first <<" "<< LLs[i].second << endl; 

  sum_tree->init_treewise(S,delta,tau,lambda,"root");

  vector < pair<scalar_type,string> > rLLs=sum_tree->LL_treewise(trees);
  cout << endl <<"The probability of G-s in the vector trees summed over all reconciliations for the ML root :" << endl << endl;
  for ( int i=0;i<trees.size();i++)
    cout << rLLs[i].first <<" "<< rLLs[i].second << endl; 

  // Alternatively we can calculate the likelihood of the ML reconciliation .
  cout << endl <<" The probability of G-s for the ML root and reconciliation:" << endl << endl;
  Species_tree * max_tree = new Species_tree(S);
  max_tree->init_treewise(S,delta,tau,lambda,"tree");

  vector < pair<scalar_type,string> > mLLs=max_tree->LL_treewise(trees);


  //Reconciliations
  for ( int i=0;i<trees.size();i++)
    cout << mLLs[i].first <<" "<< mLLs[i].second << endl; 

  // Changing the origination probability can be done by passing an extra arg. to LL_treewise, this sets P_orig._root = C x P_orig._outside_root
  cout << endl <<"The probability of G-s for the ML root and reconciliation, with the probability of originating at the root being 1e3 more probable:" << endl << endl;
  mLLs.clear();
  max_tree->init_treewise(S,delta,tau,lambda,"tree");
  mLLs=max_tree->LL_treewise(trees,1e3);
  
  for ( int i=0;i<trees.size();i++)
    cout << mLLs[i].first <<" "<< mLLs[i].second << endl; 

  // Changing the DTL rates ..
  cout << endl <<"The probability of G-s for the ML root and reconciliation, with zero transfer rate:" <<endl << endl;

  //Species_tree * noT_max_tree = new Species_tree(S);
  max_tree->init_treewise(S,delta,0,lambda,"tree");

  mLLs.clear();
  mLLs=max_tree->LL_treewise(trees,1);
  for ( int i=0;i<trees.size();i++)
    cout << mLLs[i].first <<" "<< mLLs[i].second << endl; 

  // To illustrate the format of the reconciliations:
  // -19.2316 (e:1,((b&SL@3|3|7:2,(c:1,d:1)S@4|4|13:1)S@2|2|4:1,a&SL@2|2|3&SL@3|3|7:2)D@2|2|3:1)S@1|1|0;
  // Remember to ignore things after the pipes "|" except for the first two numbers if the event is a T!
  // Remember that the gene tree is rooted!
  // Reading from left to right:
  // b&SL@3|3|7 should be read b&SL@3 and it means that the gene linage leading to b underwent a speciation S at species tree node 3, but the only the lineage leading to b survived to present
  // (c:1,d:1)S@4 the gene tree linage leading to the node (c,d) split in a speciation at node 4
  // S@2 same..
  // a&SL@2|2|3&SL@3|3|7 read a & SL@2 & SL@3 the gene linage leading to a underwent a series of speciation at 2 and 3, but only the linage leading to a  survived to present
  // D@2 the lineage leading to the gene tree node ((b,(c,d)),a) split as a result of duplication somewhere above the species tree node with ID 2
  // S@1 the gene tree originated somewhere above the root of the species tree (ID 1).
  // The same tree if transfer is allowed:
  // -17.6017 ((c:1,d:1)S@4|4|13:0.5,(b:2,(a:2,e:2)T@5|5|9|18:1)S@3|3|8:0.5)S@2|2|4&SL@1|1|1;
  // reading chronologically:
  // - G root S@2 & SL@1; the G tree linage leading to the root of G originated above the root (this is because we gave a prior of 1e3 for originating at the root)
  // there was a loss after S@1 before the linage split at S@2 as a results of a speciation
  // - both descendent linages under went uneventful vertical inheritane spliting at S@3 and S@4 as a result  of a speciations.
  // - the linages resulting from S@4 we see today as genes in the genomes of c and d.
  // - the linages descending from the split at S@3 we obsorve in three genomes, a and b as we expect from vertical inheritance   
  // and e as a result of a transfer T@5|5|9|18 read as T@5|5|9, meaning a transfer in time slice 5 from above node 5 to above node 9 
  // (to be clear T@5|8|9 is trf in time slice 5 from 8 to 9)
  return 1;
}
