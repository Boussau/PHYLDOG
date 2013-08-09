#include <iostream>
#include "DTL.h"
using namespace std;
using namespace bpp;


int main(int argc, char ** argv)
{ 

  string tree_file_name=argv[1];
  ifstream tree_file (tree_file_name.c_str());
  string tree;
  getline (tree_file,tree);

  // --Initialization--  
  // A species tree instance is created:
  TreeTemplate<Node> * S = TreeTemplateTools::parenthesisToTree(tree);
  Species_tree * sim_tree = new Species_tree(S);
  // To output the random branch-wise rates we are generating we define the stream "sim_tree->event_stream"  
  ofstream rates_file(((string)argv[1]+".rates").c_str());    
  sim_tree->event_stream=&rates_file;
  // The simulation has to be initialized with the length of the stem above the root (where D&L may occur)
  // - the S tree has unit height _below_ the root, i.e. if, as below we specify a stem length of 0.1 we have a tree with total height of 1.1 
  // - branch-wise D,T,L and origination rates are drawn from a gamma distribution given mean and variance 
  // - the T_ef factors are branchpair-wise multiplicative factors that have mean 1 and are drawn from a gamma distribution with the specified s.d.
  // the transfer rate from branche e to branch f is tau_f X T_ef , where tau_f is the rate for the receiving branch and T_ef is the modifier factor for transfer from e to f.
  // setting any of the s.d.-s to zero gives homogenous branch-wise rates (D,T,L and origination) and branchpair-wise homogenous rates.  
  sim_tree->init_sim(0.01, 0.2, 0.01, 0.7, 0.0001,1e-7,1e-7,0.,0.,0.01);
  // stem length, D rate mean, T rate mean, L rate mean,  D rate standard deviation, T rate s.d., L rate s.d., Origination rate , Origination rate s.d., s.d of T_ef factor) 

  rates_file.close();

  // --Simulation--  
  // There are two ways to output trees:

  // i) one can optionally define the stream "sim_tree->event_stream", each tree is given a sequential id #SBG 1 - #SBG number of trees
  // the format of the file is: tree id \newline Newick tree \newline .. 
  ofstream trees_file(((string)argv[1]+".trees").c_str());    
  sim_tree->event_stream=&trees_file;

  // ii) one must always pass a string vector that will contain the Newick trees coming from the simulation 
  // emit_G_trees( number of trees , vector of trees to be filled)
  vector <string> trees;
  sim_tree->emit_G_trees(1000,&trees);  

  trees_file.close();
  for (vector<string>::iterator it=trees.begin();it!=trees.end();it++)
    cout << (*it);

    
  trees.clear();

  // --misc.--  
  // the mapping between node id-s, branch id-s and leaf and node names can be out put using the legend file:
  ofstream legend_file(((string)argv[1]+".leg").c_str());    
  legend_file << sim_tree->legend();
  legend_file.close();

}
//taxa={"PROM2":1,"PROMP":1,"PROM1":1,"PRMAR1":1,"PROMM":1,"SYNS3":1,"SYNPW":1,"SYNSC":1,"SYNPX":1,"SYNS9":1,"SYNR3":1,"SYNE7":1,"TRIEI":1,"NOSP7":1,"ANASP":1,"SYNP2":1,"SYNY3":1,"CYAA5":1,"CYAP8":1,"CYAP7":1,"MICAN":1,"ACAM1":1,"THEEB":1,"CYAP4":1,"GLVIO1":1,"SYNJB":1,"SYNJA":1,"GUITH":2,"ARATH":2,"OSTTA":2}
