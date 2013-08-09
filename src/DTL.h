
#include <iostream>
#include <set>

//bio++
//#include <Bpp/Numeric/random>
#include <Bpp/Numeric/Random.all>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
//#include <Bpp/Phyl/iotree>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTools.h>

// BOOST ublas
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/mpi.hpp>
//#include <boost/fusion/tuple.hpp>
//#include <boost/fusion/support.hpp>
//#include <boost/fusion/sequence.hpp>
//#include <boost/fusion/mpl.hpp>
//#include <boost/fusion/iterator.hpp>
//#include <boost/fusion/functional.hpp>
//#include <boost/fusion/container.hpp>
//#include <boost/fusion/algorithm.hpp>
//#include <boost/fusion/support/detail/access.hpp>
//#include <boost/mpl/int.hpp>
//#include <boost/fusion/adapted/adt/detail/adapt_base.hpp>
//#include <boost/fusion/support/detail/access.hpp> 

// BOOST ublas blas (& lapack, atlas etc.) bindings, nonstandard
// cf. http://svn.boost.org/svn/boost/sandbox/numeric_bindings/boost/numeric/bindings/
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/blas/blas.hpp>



typedef long double scalar_type;

typedef boost::numeric::ublas::vector<short> int_vector_type;

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::vector<long double> long_vector_type;
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> array_type;
typedef bpp::TreeTemplate<bpp::Node> tree_type;
typedef bpp::Node node_type;

// the Species tree class can emit a gene tree and compute it's likelihood given rates of D,T and L events
// construction of the both the simulation and the likelihood are based on Tofig h et al.

class Species_tree
{ 
 public:
  bool MPI;
  std::map <std::string,scalar_type> code2tree_lls;
  std::vector <std::string> gcodes;
  std::vector <scalar_type> glls;

  
  std::vector<scalar_type> tree_lls;
  std::vector<std::string> tree_codes;
  int n_client_trees;
  /* Makes uderlying tree clocklike, returns height */
  scalar_type node_makeclocklike_height(node_type* nd);
  std::string panj_string;
  //*****************************************************
  //************** PUBLIC ** PUBLIC ** PUBLIC ***********
  //*****************************************************
  std::vector<std::string> in_trees;
  std::map<std::string,std::string> codes;
  std::vector<std::string> recs;
  std::vector<scalar_type> peas;
  scalar_type tmp_pea;
  scalar_type tmp_pea_sum;

  std::map<int,int> branch_2_x,below_branch0,below_branch1;
  std::map<int,scalar_type> branch_length;

  
  std::map <std::string,scalar_type> LL_cache;

  scalar_type T_avg,D_avg,L_avg,lambda_avg,delta_avg,tau_avg,omega_avg;
  scalar_type T_c,D_c,L_c,beta,w,intp;
  std::string mode;
  std::string gene_token;
  std::string REF_Sstring;
  std:: map< int ,std::map < int ,array_type> > slice_2_slice;
  std::map< int, std::map< int , std::map<int, std::map< int, array_type  > > > >  edge_2_edge;
  std::map< int , scalar_type > edges_in_slice;
  std::map < int , std::map <int , std::map <int , int > > > flat_map;
  
  //vector_type Qe_s_t[N_slices][ddc];
  std:: map< int ,std::map < int ,vector_type> >  Qe;
  int N_slices;
  bool branchwise;
  std::map <int,int> lin_slice;
  std::map <int,int> lin_si;
  std::map <int,int> lin_edge;
  
  scalar_type omega_root;
  scalar_type zerotree;
  scalar_type onetree;
  scalar_type llonetree;
  scalar_type local_max;
  scalar_type tree_samples;
  tree_type * tree, * delta_tree, * tau_tree, * lambda_tree, * real_tree ,* S_tree, *G, *deco_tree;
  node_type * root, * delta_root, * tau_root, * lambda_root, * G_root, * real_root;
  std::map<int,scalar_type> branchwise_delta;
  std::map<int,scalar_type> branchwise_tau;
  std::map<int,scalar_type> branchwise_lambda;
  std::map<int,scalar_type> branchwise_omega;
  std::map<int,scalar_type> mem_branchwise_delta;
  std::map<int,scalar_type> mem_branchwise_tau;
  std::map<int,scalar_type> mem_branchwise_lambda;
  std::map<int,scalar_type> mem_branchwise_omega;
  std::map<int,scalar_type> mem_branch_count;
  std::map<int,scalar_type> mmem_branch_count;

  std::map<int,scalar_type> mem_branch_sum;
  std::map<int,scalar_type> mmem_branch_sum;

  std::map<int,scalar_type> mem_branch_fam_count;
  
  std::map<int,int > sim_Ds_count;
  std::map<int,int > sim_Ts_count;
  std::map<int,int > sim_Ls_count;
  int sim_x;
  std::map<int,scalar_type > mem_sim_O_profile;
  std::map<int,scalar_type > mem_rec_O_profile;
  //std::map<int,scalar_type > sim_O_profile;
  std::vector<scalar_type> sim_O_profile;

  std::map<int,scalar_type > rec_O_profile;

  std::map<int,int > rec_Ds_count;
  std::map<int,int > rec_Ts_count;
  std::map<int,int > rec_Ls_count;


  std::map<int,scalar_type> sim_D_count;
  std::map<int,scalar_type> sim_T_count;
  std::map<int,scalar_type> sim_L_count;
  std::map<int,scalar_type> sim_O_count;
  //std::vector<scalar_type> sim_D_count;
  //std::vector<scalar_type> sim_T_count;
  //std::vector<scalar_type> sim_L_count;
  //std::vector<scalar_type> sim_O_count;

  std::map<int,scalar_type> sim_branch_count;
  std::map<int,scalar_type> sim_branch_sum;
  std::map<int,scalar_type> sim_branch_genome_size;


  std::map<int,scalar_type> mem_D_count;
  std::map<int,scalar_type> mem_T_count;
  std::map<int,scalar_type> mem_L_count;
  std::map<int,scalar_type> mem_O_count;

  std::map<int,scalar_type> mem_branch_genome_size;
  std::map<int,scalar_type> mmem_branch_genome_size;

  std::map<int,scalar_type> mmem_D_count;
  std::map<int,scalar_type> mmem_T_count;
  std::map<int,scalar_type> mmem_L_count;
  std::map<int,scalar_type> mmem_O_count;

  std::map<std::string,scalar_type> namewise_delta;
  std::map<std::string,scalar_type> namewise_tau;
  std::map<std::string,scalar_type> namewise_lambda;
  std::map<std::string,scalar_type> namewise_omega;
  std::map<int,std::string> branchname;
  std::map<std::string,int> branchi;
  std::map<std::string,int> tried_spr;
  
  std::map<std::string,std::string> rate_trees;

  long_vector_type p_omega;
  long_vector_type tmp_omega,tmpp_omega;
  std::map <int,scalar_type> tmp_branch_delta,tmp_branch_tau,tmp_branch_lambda,tmp_branch_omega,count_branch_omega,count_branch_D,count_branch_T,count_branch_L;
  std::map <scalar_type,std::vector<std::pair<int,std::pair<node_type*,node_type*> > > > summary;
  scalar_type tree_norm;
  scalar_type tmp_ll;
  std::string tmp_Sstring;
  
  boost::mpi::communicator  world;
  vector_type branch_delta;
  vector_type branch_tau;
  vector_type branch_lambda;
  vector_type branch_omega;
  std::map <int, std::map <int, scalar_type> > Tef_var;

  bool infer_mode;

  int event_count,D_count,T_count,L_count,TL_count;
  scalar_type max_P,max_L;

  std::ostream * event_stream;
  std::ofstream nullstream ;

  // resolutions for init L
  int dc;
  int ddc;
  int dc_mod;
  // ..
  scalar_type unrooted_register(std::vector <std::string> tree_strings,scalar_type treeN=1e20);
  scalar_type unrooted_L(std::string tree_string,std::string mode="N2");

  scalar_type unrooted_run(scalar_type unseen=0,std::string mode="N2");
  void DPstep_L_N3(std::pair <node_type*,node_type*> dnode ,std::pair <node_type*,node_type*> left_dnode,std::pair <node_type*,node_type*> right_dnode);
  void DPstep_L_N2(std::pair <node_type*,node_type*> dnode ,std::pair <node_type*,node_type*> left_dnode,std::pair <node_type*,node_type*> right_dnode);
  void DPstep_m_N2(std::pair <node_type*,node_type*> dnode ,std::pair <node_type*,node_type*> left_dnode,std::pair <node_type*,node_type*> right_dnode);

  void backstep_m_N2(std::pair <node_type*,node_type*> dnode ,std::pair <node_type*,node_type*> left_dnode,std::pair <node_type*,node_type*> right_dnode);
  void backtrace(std::pair <node_type*,node_type*> root_dnode,tree_type* tree,int in_o_x);
  void backtrace(std::pair <node_type*,node_type*> dnode,int x);
  void name_leaf(node_type * leaf, int x, int leaf_x=-1); 

  void tracecount(std::pair <node_type*,node_type*> root_dnode,tree_type* tree,int in_o_x);
  void tracecount(std::pair <node_type*,node_type*> dnode,int x);
  void tracecount(node_type * leaf, int x, int leaf_x=-1); 
  void record_event(int e_type, int at_x,int xc1,int xc2=-1,int o=1);

  void tracestream(std::pair <node_type*,node_type*> root_dnode,tree_type* tree, int x=-1);
  void tracestream(std::pair <node_type*,node_type*> dnode,int x);
  void tracestream(node_type * leaf, int x, int leaf_x=-1); 
  void remember_rates();
  void recall_rates();
  void print_rates();


  std::vector <scalar_type> O_counts;
  std::vector <scalar_type> D_counts;
  std::vector <scalar_type> T_counts;
  std::vector< std::vector< scalar_type> > Ttf;
  std::vector <scalar_type> branch_sum;
  std::vector <scalar_type> branch_count;
  std::vector <scalar_type> branch_genome_size;
  std::vector <scalar_type> branch_zero;
  std::vector <scalar_type> branch_one;
  std::vector <scalar_type> branch_two;
  std::vector<scalar_type> branch_from_count;
  std::vector<scalar_type> from_count;

  scalar_type DTLclock_height(int i);

  int x_down_f(int x);
  boost::numeric::ublas::vector<int> x_down;
  int node_down_f(int x);
  boost::numeric::ublas::vector<int> node_down;

  std::map <int,int> node_down_cache;

  std::string legend();
  void x_DFS();
  void x_DFS(node_type*);
  void Qef_m_N2(node_type* leaf);
  void Qef_bck_N2(node_type* leaf);

  void Qef_N2(node_type* leaf);


  void init_sim(scalar_type stem_len, scalar_type delta_mean, scalar_type tau_mean, scalar_type lambda_mean, scalar_type delta_sd, scalar_type tau_sd, scalar_type lambda_sd, scalar_type omega_mean, scalar_type omega_sd, scalar_type pair_sd, scalar_type out_tau=0.);
  void init_L(scalar_type delta_in,scalar_type tau_in,scalar_type lambda_in);
  void init_L_improved(scalar_type delta_in,scalar_type tau_in,scalar_type lambda_in,scalar_type stem_omega=0.,scalar_type out_tau=0.,std::string mode="bck");
  void init_L_improved(scalar_type stem_omega=0.1);
  void reset(scalar_type delta=0.01,scalar_type tau=0.01,scalar_type lambda=0.01,scalar_type stem_omega=0.1);  
  void init_x();
  void init_parsimony(scalar_type delta_in,scalar_type tau_in,scalar_type lambda_in);
  // simulate a gene tree stored in global var tree with root pointed to by global var root
  void emit_G_trees(int N_trees ,std::vector <std::string> * trees = NULL,std::vector <std::string> * codes = NULL,std::vector <std::string> * stream = NULL);
  bool homogenous_model;
  scalar_type nclasses;

  std::map< int,std::map<int,int> > slice_2_branches;  
  std::map< int,std::map<int,int> > branch_2_slices;  
  //count most on for each ... 
  scalar_type count_events(std::vector<std::string> recs,bool homogenous =false,bool subT=true);
  void count_labels(std::vector<std::string> trees);

  void apparent_rates(scalar_type pea_count,bool homogenous,bool subT );
  void count_events(node_type * head);

  void emit_G_tree(scalar_type delta_in, scalar_type tau_in, scalar_type lambda_in,scalar_type omega_in=0 ,scalar_type eta_in=0 ,bool remember_slices=false);
  void origination(scalar_type omega);
  // name nodes using above
  void name_internal(node_type * node, std::string pre="");
  // strip label leaves
  void strip_label_leaves();
  // strip virtual leaves 
  void strip_virtual_leaves();
  void G_strip_virtual_leaves(tree_type * G_tree,bool add_xL=false);
  void add_xL(node_type* node,node_type*child);
  std::string stream_xL(tree_type * G_tree);
  std::string stream_xL(node_type * node);
  
  void init_treewise(tree_type * S, scalar_type delta, scalar_type tau, scalar_type lambda,std::string calc_mode="Nest");
  std::vector < std::pair<scalar_type,std::string> > LL_treewise(std::vector<std::string> 
								 trees,scalar_type root_p=1);

  // strip S prime type nodes from tree
  void strip(node_type * node=NULL);
  void strip(node_type * node,bool add_xL);

  void G_strip(node_type * node,tree_type * G_tree);

  void defaults();

  // construct slices without emitting G tree
  void construct_slices();

  void construct_random_tree(int stop); 
  //constructors
  //OMFG
  ~Species_tree();
  // construct from random tree 
  Species_tree(int leaves);

  // construct from input tree
  Species_tree(tree_type * S);
  void reconstruct(tree_type * S);
  void add_stem();

  //newick format output
  void write(std::ostream * out_stream);
  void write(tree_type * out_tree,std::ostream * out_stream);
  void canonic_order();
  void canonic_order(node_type * node);

  void make_rate_table(bool reinit);
  void make_rate_table(node_type * tree_node,node_type * delta_node,node_type * tau_node,node_type * lambda_node);
  std::vector<scalar_type> DTL_lk(tree_type* S,std::vector<std::string> trees,std::map<node_type*,scalar_type> expD,std::map<node_type*,scalar_type> expT,std::map<node_type*,scalar_type> expL,std::string output="no_sum_orig");


  std::map <node_type *,int > flat_node_map;
  std::map <int,node_type *> inverse_flat_node_map;
  int lin_size;
  int ur_S_N_leaves;
  std::string outstream_file_name;

//*****************************************************
//************ PRIVATE ** PRIVATE ** PRIVATE **********
//*****************************************************
  std::vector<tree_type *> trees;

 private:  

  //unrooted maps are global
  std::map <std::string,int> ur_sigma;
  std::vector< std::vector< std::pair<node_type*,node_type*> > > run_vec;
  std::vector< std::vector< std::pair<node_type*,node_type*> > > trees_roots;
  std::vector<tree_type *> gin_trees;
  std::map < std::pair <node_type*,node_type*>, std ::pair <node_type*,node_type*> > root_dnodes;


  std::map < std::pair<node_type *,node_type *>,  long_vector_type *> ur_a;
  std::map < std::pair<node_type *,node_type *>,  int_vector_type *> max_event_name;
  std::map < std::pair<node_type *,node_type *>,  std::map <int , std::vector<int> > > max_hidden_events;
  std::map < std::pair<node_type *,node_type *>,  int_vector_type *> max_event_x;
  std::map < std::pair<node_type *,node_type *>,  int_vector_type *> max_event_xc1;
  std::map < std::pair<node_type *,node_type *>,  int_vector_type *> max_event_xc2;

  std::map < std::pair<int,int>,  int_vector_type *> cherry_max_event_name;
  std::map < std::pair<int,int>, std::map <int , std::vector<int> > > cherry_max_hidden_events;
  std::map < int,std::map< int , int > >  leaf_hidden_event_name;
  std::map < int,std::map < int , int > > leaf_hidden_event_x;
  std::map < std::pair<int,int>,  int_vector_type *> cherry_max_event_x;
  std::map < std::pair<int,int>,  int_vector_type *> cherry_max_event_xc1;
  std::map < std::pair<int,int>,  int_vector_type *> cherry_max_event_xc2;

  std::map < std::pair<int,int>,  int_vector_type *> root_cherry_max_event_name;
  std::map < std::pair<int,int>, std::map <int , std::vector<int> > > root_cherry_max_hidden_events;
  std::map < std::pair<int,int>,  int_vector_type *> root_cherry_max_event_x;
  std::map < std::pair<int,int>,  int_vector_type *> root_cherry_max_event_xc1;
  std::map < std::pair<int,int>,  int_vector_type *> root_cherry_max_event_xc2;


  std::map < std::pair<node_type *,node_type *>,std::pair<node_type *,node_type *> > dnode_c1,dnode_c2;

  std::map <std::pair<node_type *,node_type *>, scalar_type> root_max;
  std::map <std::pair<node_type *,node_type *>, scalar_type> root_max_s;
  std::map <std::pair<node_type *,node_type *>, int> root_max_x;
  std::map < std::pair<int,int>, scalar_type> cherry_root_max;
  std::map < std::pair<int,int>, scalar_type> cherry_root_max_s;
  std::map < std::pair<int,int>, int> cherry_root_max_x;


  std::map < std::pair<int,int>, long_vector_type *> cherry_a;
  std::map < std::pair<int,int>, int > cherry_reg;

  std::map < std::pair<node_type *,node_type *>,std::pair<int,int> > cherries;

  std::map<std::pair<node_type *,node_type *>, std::vector<int> > ur_sum_is;

  // speciation 1 or internal edge 0
  //int type_of_x[lin_size];
  std::map<int,int> type_of_x;
  // for type 1 xp and xpp
  // for type 0 ..
  std::map<int,int> below_son0;
  std::map<int,int> below_son1;
  std::map<int,int> advance_from;
  std::map<int,int> advance_until;

  std::map<int,int> below_from;
  std::map<int,int> below_until;
  //int below_from[lin_size];
  //int below_until[lin_size];

  // flat maps for DP
  std::map <int,vector_type> flat_Qef;
  std::map <int,scalar_type> Qegp;
  std::map <int,vector_type> p_flat_Qef;
  std::map <int,scalar_type> delta_dt;
  std::map <int,scalar_type> tau_dt;

  std::map <int, std::map <int, node_type * > > node_map;

  

  // global rates of D,T and L
  scalar_type delta,tau,lambda,eta;
  // rates of D,T and L per branch 
  std::map< std::string , scalar_type > delta_table, tau_table, lambda_table;
  std::map<int , scalar_type > id_delta_table, id_tau_table, id_lambda_table;
  std::map<int , int> id_slice_table;
  // node with id
  std::map<int , node_type*> id_table;
  // length of edge with head with id
  std::map<int , scalar_type> l_table;
  // number of dt long subslices on edge with head with id
  std::map<int , scalar_type> N_table;
  bool branchwise_rates;

  // current time for simulation with t=0 at root
  scalar_type t;
  // time ordered candidates
  std::map<scalar_type, node_type * > candidates;
  std::map<scalar_type, std::string > candidate_names;

  // edges in current time slice given by head nodes
  std::vector<node_type * > current_slice;
  std::vector< std::vector< node_type * > > time_slices;
  std::vector< std::vector <int> > id_time_slices;
  // some nonessential behaviour modifiers
  bool marker_nodes,strip_Sp_nodes;
  bpp::Newick newick;

  //bool isorisnot(int slice,int si,int edge, int N_edges)
  //{return (si==0 && (edge==N_edges-1 || slice == 0)) || (si>0);}

  //#########################################################
  //################### Borrowed from HGT_simul #############
  //#########################################################

  /* Multiply branch lengths by f at nd and underlying nodes */
  void scale_subtree(node_type* nd, scalar_type f);

  /* Returns the height of node nd (assumingly from a clock-like tree) */
  scalar_type node_height(node_type* nd);

  /* Makes uderlying tree clocklike, returns height */
  //scalar_type node_makeclocklike_height(node_type* nd);

  //#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  //################### Borrowed from HGT_simul #############
  //#########################################################

  //#########################################################
  //################### Some general S tree operatinos ######
  //#########################################################

  // "safe" distance to father
  scalar_type distance_to_father(node_type * node);

  // add marker nodes to display internal node names in newick format output 
  void add_marker_nodes(node_type * node);
  void add_marker_node(node_type * node,std::string name);
  scalar_type distance_to_root(node_type * node, node_type * root);

  // rename leaves 0,1,2,..#leaves-1 to keep internal node names short
  void rename_leaves(node_type * node);
  // get names of all leaves in subtree below node
  std::string get_leaf_names(node_type * node);

  // construct name of nodes based on leaves in subtree below them
  std::string get_name_in_S(node_type * node, std::string pre="");

  scalar_type get_x_in_S(node_type * node);


  //#########################################################
  //################### Some general S tree operations ######
  //#########################################################

  //#########################################################
  //################### S & S' tree construction ############
  //#########################################################

  // look for node that happened first after current time t  
  void lookforCandidates(node_type* node);

  // find node defining next time slice given t time has passed
  void findCandidates();

  void add_leaf(node_type * to, scalar_type d, std::string name);

  void add_speciation(node_type * leaf);

  // insert a node at Delta_t below tail of edge leading to the father of head
  void insert_node(node_type * head, scalar_type Delta_t, std::string name);

  void remove_leaf(node_type * node);

  // remove a S prime type node from tree
  void remove_node(node_type * node);
  void remove_node(node_type * node,bool add_xL);

  // test if a node is an S prime type node
  bool isSpNode(node_type * node);

  //#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  //################### S & S' tree construction ############
  //#########################################################

  //#########################################################
  //################### DTL events and necc. operations #####
  //#########################################################

  // prune subtree below degree one node node
  void prune_at(node_type* node);

  void prune(node_type* node);
  
  node_type * get_S_node(node_type * node);

  // graft subtree below from to below event node to
  void graft(node_type * from, node_type * to,  node_type * root);

  //D event Delta_t below start of time slice on branch with head node
  void D(node_type * node, scalar_type Delta_t);

  // T event Delta_t below start of time slice, 
  // transfer from branch with head from to branch with head to 
  // corresponds to grafting on the branch with head from 
  // the subtree that branch with head to defines 
  // (transfer to "to" from "from" means subtree from "to" must be grafted to "from")
  void T(node_type * from, node_type * to,scalar_type Delta_t);
  void TL(node_type * from, node_type * to,scalar_type Delta_t);

  // L event Delta_t below start of time slice on branch with head node
  void L(node_type * node, scalar_type Delta_t);

  // construct time slice defined by node
  // i.e. construct a std::vector of nodes that are heads of branches that are contemporary with node
  void construct_slice(node_type * node);

  // filter time slice to get valid T donors
  std::vector< node_type * > filter_slice(node_type * to);
  
  // simulate next event in current_slice
  scalar_type simulate_slice( scalar_type l_j, node_type * node_j, bool  no_events=false);

  //#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  //################### DTL events and necc. operations #####
  //#########################################################
  void init();

  // 
  scalar_type delta_cost,tau_cost,lambda_cost;

  
  //END Species_tree    
};

std::vector<std::string > all_NNIs(std::string Sstring);

scalar_type LL(Species_tree * infer_tree, tree_type * S, scalar_type delta, scalar_type tau, scalar_type lambda, scalar_type beta=1.,std::string mode="max",scalar_type stem_omega=0.,scalar_type out_tau=0.);


scalar_type oLL(Species_tree * infer_tree, tree_type * S, scalar_type delta, scalar_type tau, scalar_type lambda, scalar_type beta=1.,std::string mode="max",scalar_type stem_omega=0.,scalar_type out_tau=0.);

double GS_alpha(Species_tree * infer_tree,tree_type * S,double rho,double height, double beta=1,std::string mode="max", double stem_l=0,double stem_o=0);
double GS_rho(Species_tree * infer_tree,tree_type * S,double alpha,double height, double beta=1,std::string mode="max", double stem_l=0,double stem_o=0);
double GS_height(Species_tree * infer_tree,tree_type * S,double alpha,double rho,double height, double beta=1,std::string mode="max", double stem_l=0,double stem_o=0);
double GS_omega(Species_tree * infer_tree,tree_type * S,double alpha,double rho,double height, double beta=1,std::string mode="max",double stem_o=0);

std::pair <scalar_type,scalar_type> observe_BD(Species_tree * infer_tree, tree_type * S, std::vector<std::string> trees,scalar_type delta,scalar_type tau,scalar_type lambda);
std::pair < scalar_type,std::string> ML_step(Species_tree * infer_tree, std::string Sstring,std::string mode="max",std::string move_mode = "joint",bool greedy=false);
std::pair < scalar_type,std::string> last_step(Species_tree * infer_tree, std::string Sstring,std::string mode="max",int supp=0,bool greedy=false);

std::pair < scalar_type,std::string> ML_nni_step(Species_tree * infer_tree, std::string Sstring,std::string mode="max",bool greedy="false");

std::string clock_root(Species_tree * infer_tree, std::string String,std::string mode);


tree_type * random_step(tree_type * S,std::string mode);
tree_type * random_nni_step(tree_type * S);

std::string randomize_time_order(std::string Sstring,int steps=1000,std::string="order");
std::string randomize_nni(std::string Sstring,int steps=1000);
std::string get_name(node_type * node);


std::string print_time_order(std::string Sstring);
std::string can_order(tree_type * S);
std::string DTL_time_order(Species_tree * infer_tree, std::string Sstring);
scalar_type   rate_estimate(tree_type * S,std::vector <std::string> trees,Species_tree * infer_tree, bool homogenous);
std::string DTL_clock(std::string S,Species_tree * infer_tree, std::string mode="Hest");
std::string canonic_time_order(tree_type * T);
void Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = " ");

std::pair < tree_type *,std::vector<std::string> > restrict(tree_type * S,std::vector<std::string> sample_species,std::vector<std::string> sample_trees);

std::vector<scalar_type> gamma_class(std::vector<scalar_type>  vals, int nclasses );

std::string PANJ(std::string Sstring, std::vector<std::string> treestrings);

unsigned int good_seed();
