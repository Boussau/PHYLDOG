#include "DTL.h"

using namespace std;
using namespace bpp;

//*****************************************************
//************** PUBLIC ** PUBLIC ** PUBLIC ***********
//*****************************************************


 
Species_tree::~Species_tree()
{
  for (map< int ,std::map < int ,array_type> >::iterator it=slice_2_slice.begin();it!=slice_2_slice.end();it++)
    {
      for (std::map < int ,array_type>::iterator iit=(*it).second.begin();iit!=(*it).second.end();iit++)
	(*iit).second.resize(0,0);
      (*it).second.clear();
    }
  slice_2_slice.clear();
  //std::map< int, std::map< int , std::map<int, std::map< int, array_type  > > > >  edge_2_edge;
  for (std::map< int, std::map< int , std::map<int, std::map< int, array_type  > > > >::iterator it=edge_2_edge.begin();it!=edge_2_edge.end();it++)
    {
      for (std::map< int , std::map<int, std::map< int, array_type  > > >::iterator iit=(*it).second.begin();iit!=(*it).second.end();iit++)
	{
	  for (std::map<int, std::map< int, array_type  > > ::iterator iiit=(*iit).second.begin();iiit!=(*iit).second.end();iiit++)
	    {
	      for (std::map< int, array_type  > ::iterator iiiit=(*iiit).second.begin();iiiit!=(*iiit).second.end();iiiit++)
		{
		  (*iiiit).second.resize(0,0);		  
		}
	      (*iiit).second.clear();
	    }
	  (*iit).second.clear();
	}
      (*it).second.clear();
    }
  edge_2_edge.clear();

  edges_in_slice.clear();
  //std::map < int , std::map <int , std::map <int , int > > > flat_map;
  for (std::map < int , std::map <int , std::map <int , int > > > ::iterator it=flat_map.begin();it!=flat_map.end();it++)
    {
      for ( std::map <int , std::map <int , int > >::iterator iit=(*it).second.begin();iit!=(*it).second.end();iit++)
	{
	  (*iit).second.clear();
	}
      (*it).second.clear();        
    }
  flat_map.clear();

  //std:: map< int ,std::map < int ,vector_type> >  Qe;
  for (map< int ,std::map < int ,vector_type> >::iterator it=Qe.begin();it!=Qe.end();it++)
    {
      for (std::map < int ,vector_type>::iterator iit=(*it).second.begin();iit!=(*it).second.end();iit++)
	(*iit).second.resize(0);
      (*it).second.clear();
    }
  Qe.clear();
  //tree_type * tree, * delta_tree, * tau_tree, * lambda_tree, * real_tree ,* S_tree, *G, *deco_tree;
  //node_type * root, * delta_root, * tau_root, * lambda_root, * G_root, * real_root;
  ur_sigma.clear();
  for (std::map < std::pair<node_type *,node_type *>, long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)  
    delete (*it).second;
  ur_a.clear();
  for (  std::map < std::pair<int,int>, long_vector_type *>::iterator it= cherry_a.begin();it!= cherry_a.end();it++)  
    delete (*it).second;
  cherry_a.clear();
  
  cherries.clear();
  type_of_x.clear();
  below_son0.clear();
  below_son1.clear();
  advance_from.clear();
  advance_until.clear();
  below_from.clear();
  below_until.clear();
  flat_Qef.clear();
  p_flat_Qef.clear();
  delta_dt.clear();
  tau_dt.clear();
  flat_node_map.clear();
  inverse_flat_node_map.clear();
  lin_slice.clear();
  lin_si.clear();
  lin_edge.clear();
  delta_table.clear(); tau_table.clear(); lambda_table.clear();
  id_delta_table.clear(); id_tau_table.clear(); id_lambda_table.clear();
  id_slice_table.clear();
  id_table.clear();
  l_table.clear();
  N_table.clear();
  candidates.clear();
  candidate_names.clear();
  current_slice.clear();
  //std::vector< std::vector< node_type * > > time_slices;
  //std::vector< std::vector <int> > id_time_slices;

  time_slices.clear();
  id_time_slices.clear();
  
}

  /*

Species_tree::~Species_tree()
{
  std:: map< int ,std::map < int ,array_type> > slice_2_slice;
  std::map< int, std::map< int , std::map<int, std::map< int, array_type  > > > >  edge_2_edge;
  std::map< int , scalar_type > edges_in_slice;
  std::map < int , std::map <int , std::map <int , int > > > flat_map;
  std:: map< int ,std::map < int ,vector_type> >  Qe;
  tree_type * tree, * delta_tree, * tau_tree, * lambda_tree, * real_tree ,* S_tree, *G, *deco_tree;
  node_type * root, * delta_root, * tau_root, * lambda_root, * G_root, * real_root;
  std::map <node_type *,int> ur_sigma;
  std::map < std::pair<node_type *,node_type *>, vector_type> ur_a,ur_am,ur_pam,ur_m_head,ur_m_tail;
  std::map < std::pair<int,int>, vector_type> cherry_a;
  std::map < std::pair<node_type *,node_type *>,std::pair<int,int> > cherries;
  std::map<std::pair<node_type *,node_type *>, std::vector<int> > ur_sum_is;
  std::map<int,int> type_of_x;
  std::map<int,int> below_son0;
  std::map<int,int> below_son1;
  std::map<int,int> advance_from;
  std::map<int,int> advance_until;
  std::map<int,int> below_from;
  std::map<int,int> below_until;
  std::map <int,vector_type> flat_Qef;
  std::map <int,vector_type> p_flat_Qef;
  std::map <int,scalar_type> delta_dt;
  std::map <int,scalar_type> tau_dt;
  std::map <node_type *,int > flat_node_map;
  std::map <int,node_type *> inverse_flat_node_map;
  std::map <int,int> lin_slice;
  std::map <int,int> lin_si;
  std::map <int,int> lin_edge;
  scalar_type delta,tau,lambda,eta;
  std::map< std::string , scalar_type > delta_table, tau_table, lambda_table;
  std::map<int , scalar_type > id_delta_table, id_tau_table, id_lambda_table;
  std::map<int , int> id_slice_table;
  std::map<int , node_type*> id_table;
  std::map<int , scalar_type> l_table;
  std::map<int , scalar_type> N_table;
  std::map<scalar_type, node_type * > candidates;
  std::map<scalar_type, std::string > candidate_names;
  std::vector<node_type * > current_slice;
  std::vector< std::vector< node_type * > > time_slices;
  std::vector< std::vector <int> > id_time_slices;
  bpp::Newick newick;
  
}
  */


// name nodes using above
void Species_tree::name_internal(Node * node, string pre)
{ 
  node->setName(get_name_in_S(node,pre));
  vector<Node * > sons = node->getSons();
  for (int i = 0;i < node->getNumberOfSons();i++)
    name_internal(sons[i]);
} 

// strip label leaves 
void Species_tree::strip_label_leaves()
{
  vector <Node * > leaves = TreeTemplateTools::getLeaves(*root);
  vector<Node *>::iterator it;	
  for ( it = leaves.begin(); it!=leaves.end(); it++ )
    {
      if ((*it)!=root && (*it)->getName()[0]=='*' && (*it)->getDistanceToFather() == 0)
	remove_leaf((*it));
    }
}

// strip virtual leaves 
void Species_tree::strip_virtual_leaves()
{
  bool striped=false;
  while (!striped)
    {
      vector <Node * > leaves = TreeTemplateTools::getLeaves(*root);
      vector<Node *>::iterator it;	
      if (tree->getNumberOfNodes()==1)
	break;
      for ( it = leaves.begin(); it!=leaves.end(); it++ )
	{
	  if ((*it)!=root && (*it)->getName()[0]=='*' )
	    remove_leaf((*it));
	}
      striped=true;
      for ( it = leaves.begin(); it!=leaves.end(); it++ )
      	if ((*it)!=root && (*it)->getName()[0]=='*' )
	  striped=false;
    }
  while ( root->getNumberOfSons()==1 )
    {
      Node * old_root = root;
      root=root->getSon(0);
      tree->rootAt(root);
      remove_leaf(old_root);
    }
  strip();
}


// strip S prime type nodes from tree
void Species_tree::strip(Node * node)
{
  if (node==NULL)
    node=root;
  vector<Node * > sons = node->getSons();
  for( int i = 0; i < node->getNumberOfSons(); i++ )
    strip(sons[i]);
  // removes Sp nodes
  if (root != node && !(node->isLeaf()) && isSpNode(node))    
    remove_node(node);
  
  return;
}


// construct slices without emitting G tree
void Species_tree::construct_slices()
{
  marker_nodes = false;
  if (strip_Sp_nodes)
    {
      strip_Sp_nodes=false;
      emit_G_tree(0.,0.,0.,0.,0.,true);	
      strip_Sp_nodes=true;
    }
  else
    emit_G_tree(0.,0.,0.,0.,0.,true);	
  return;
}
void Species_tree::origination(scalar_type omega)
{
  vector <Node *> nodes = tree->getNodes();
  scalar_type time_sum = 0.;
  for (vector <Node *>::iterator it=nodes.begin(); it!=nodes.end(); it++)
    if ((*it)!=root)
      time_sum+=(*it)->getDistanceToFather();
  scalar_type expected_number_of_originations = time_sum*omega;
  //Poission(lambda,0) = exp(-lambda)
  //if one or more origination events occured 
  if( exp(-expected_number_of_originations) < RandomTools::giveRandomNumberBetweenZeroAndEntry(1) )
    {
      //we chose one at random
      scalar_type  tmp=RandomTools::giveRandomNumberBetweenZeroAndEntry(time_sum);
      time_sum=0.;
      for (vector <Node *>::iterator it=nodes.begin(); it!=nodes.end(); it++)
	if ((*it)!=root)
	  {
	    time_sum+=(*it)->getDistanceToFather();
	    if (time_sum>tmp)
	      {
		(*event_stream) << "#event# " << "origination" << " at " << (*it)->getName() << endl;

		root = TreeTemplateTools::cloneSubtree<Node>(*(*it)); 
		TreeTemplate<Node>  new_tree(root); 
		tree = new_tree.clone();
		root = tree->getRootNode();
		break;
	      }
	  }
    }
  return;
}

//simulate a gene tree stored in global var tree with root pointed to by global var root
void Species_tree::emit_G_tree(scalar_type delta_in, scalar_type tau_in, scalar_type lambda_in,scalar_type omega_in, scalar_type eta_in,bool remember_slices)
{
  event_count=0;
  D_count=0;
  T_count=0;
  L_count=0;
  TL_count=0;

  //cout << endl<<endl;
  scalar_type until,step_t;
  Node * node;
  if (!remember_slices) init();
  //global rates of D,T and L
  delta=delta_in;
  tau=tau_in;
  lambda=lambda_in;
  eta=eta_in;
  t=0;
  // simulate origination 
  if (omega_in>0)
    origination(omega_in);
  if (tree->getNumberOfNodes()<2)
    return;
  // find node defining next time slice 
  candidate_names.clear();
  findCandidates();
  node = candidates.begin()->second;
  until = candidates.begin()->first;
  construct_slice(node);
  if (remember_slices)
    time_slices.push_back(current_slice);	  
  candidate_names.erase(candidate_names.begin());
  
  // keep going until there are candidates
  while (1)
    {
      if (current_slice.size()==0)
	break;
      // simulate next event in slice
      step_t = simulate_slice(until-t, node, remember_slices);	      
      //cout << "#at" << t <<  " " << until << " " << step_t << " " << node->getName()<< endl;
      if (step_t == until-t)
	{
	  //cout << "pass slice" << endl;
	  if (candidate_names.size())
	    {
	      // find node defining next time slice 
	      findCandidates();
	      if (candidates.size())
		{
		  node = candidates.begin()->second;
		  until = candidates.begin()->first;
		  // construct the time slice defined by node
		  construct_slice(node);	  
		  if (remember_slices)
		    time_slices.push_back(current_slice);	  		  
		}
	      if (candidate_names.size())		
		candidate_names.erase(candidate_names.begin());
	      else
		break;
	    }
	  else
	    break;
	}
      t += step_t;	
    }
  // strip away time S prime nodes corresponding to time slices 
  if (strip_Sp_nodes)
    strip(root);
  // add marker nodes 
  if (marker_nodes)
    add_marker_nodes(root);
  return;
}

void Species_tree::defaults()
{
  // time resolution for gene tree likelihood
  // marker nodes to include internal node names in newick format
  branchwise=0;
  gene_token="%";
  marker_nodes=true;
  strip_Sp_nodes=true;
  branchwise_rates=false;
  beta=1.;					
  mode="max";
  nclasses=-1;
  nullstream.open("/dev/null"); 
  if (!event_stream)
    event_stream = &nullstream;
  return;
}

void Species_tree::construct_random_tree(int stop)
{
  vector <Node * > leaves = TreeTemplateTools::getLeaves(*root);
  if (leaves.size()<stop)
    {
      scalar_type Delta_t=RandomTools::randExponential(1./leaves.size());
      for ( vector<Node*>::iterator leaf = leaves.begin();  leaf != leaves.end(); leaf++)
	(*leaf)->setDistanceToFather((*leaf)->getDistanceToFather()+Delta_t);
      vector <Node *> tos;
      tos.push_back(root);
      RandomTools::getSample(leaves,tos);
      Node * to = tos[0];
      add_speciation(to);

      construct_random_tree(stop);
    }
  else
    {
      for ( vector<Node*>::iterator leaf = leaves.begin();  leaf != leaves.end(); leaf++)
	(*leaf)->setDistanceToFather((*leaf)->getDistanceToFather()+1./leaves.size());
    }
  return;
}
  
//constructors
// construct from random tree 
Species_tree::Species_tree(int leaves)
{
  defaults();
  // this should be replaced by a proper model
  root = new Node();  
  tree = new TreeTemplate<Node>(root);
  root->setName("root");
  add_speciation(root);
    
  construct_random_tree(leaves);
  rename_leaves(root);
  //write(&cout);
  S_tree = tree->clone(); 
  name_internal(root);

  //  add_stem();

}    
// construct from input tree
Species_tree::Species_tree(TreeTemplate<Node> * S)
{
  defaults();
  node_makeclocklike_height((S->getRootNode()));    
  tree = S->clone();
  root = tree->getRootNode();
  S_tree = tree->clone(); 
  name_internal(root);

  //  add_stem();
  return;
}

void Species_tree::reconstruct(TreeTemplate<Node> * S)
{

  vector <Node*> nodes=S->getNodes();
  for (vector<Node * > :: iterator it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      if ((*it)->getDistanceToFather()==0)
	(*it)->setDistanceToFather(0.001);
  
  bool tmp=branchwise;
  defaults();
  branchwise=tmp;
  node_makeclocklike_height((S->getRootNode()));    
  delete tree;
  tree = S->clone();
  root = tree->getRootNode();
  S_tree = tree->clone(); 
  name_internal(root);
  //canonic_order();
  //  add_stem();
  nodes.clear();
  return;
}


void Species_tree::add_stem()
{ 
  rename_leaves(root);

   
  real_tree = tree->clone();     
  real_root = root;
  scalar_type h = TreeTemplateTools::getHeight(*root);
  root = new Node();
  root -> setName("R");
  root->addSon(real_root);
  //real_root->addSon(root);
  //root->setFather(real_root);
  // stem length
  real_root->setDistanceToFather(0.1);
  tree->resetNodesId();

  if (false)
    {
      Node * out_group = new Node();
      out_group -> setName("O");
      root -> addSon(out_group);
      //out_group->setFather(root);
      out_group->setDistanceToFather(h+0.2);
    }

  //tree->rootAt(root->getId());
  
  S_tree=tree->clone();
  //rename_leaves();
  //
  name_internal(root);
    
}

//newick format output
void Species_tree::write(TreeTemplate<Node> * out_tree,ostream * out_stream)
{
  if(out_tree->getNumberOfNodes()>1)
    (*out_stream) << TreeTemplateTools::treeToParenthesis((*out_tree),false,"ID");
  else
    (*out_stream) << "()" << endl;
}
void Species_tree::write(ostream * out_stream)
{
  if(tree->getNumberOfNodes()>1)
    (*out_stream) << TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
  else
    if (tree->getNumberOfNodes()==1 && tree->getNodes()[0]->getName()[0]!='*')
      (*out_stream) << "("<< tree->getNodes()[0]->getName() <<")" << endl;
    else
      (*out_stream) << "("<< ")" << endl;
  
}
void Species_tree::canonic_order()
{
  canonic_order(root);
}
void Species_tree::canonic_order(Node * node)
{
  map <string,Node*> order_map;
  vector<Node*>sons=node->getSons();
  if (sons.size()==2)
    {
      string name0=get_name_in_S(sons[0]);
      string name1=get_name_in_S(sons[1]);
      node->removeSons();
      
      if ( ( name0 )>( name1 ) )
	{
	  node->addSon(sons[0]);
	  node->addSon(sons[1]);
	}
      else
	{
	  node->addSon(sons[1]);
	  node->addSon(sons[0]);	  
	}
      vector<Node*>sons=node->getSons();
      canonic_order(sons[0]);
      canonic_order(sons[1]);	
    }
}


//*****************************************************
//************ PRIVATE ** PRIVATE ** PRIVATE **********
//*****************************************************


//#########################################################
//################### Borrowed from HGT_simul #############
//#########################################################

/* Multiply branch lengths by f at nd and underlying nodes */

void Species_tree::scale_subtree(Node* nd, scalar_type f){

  scalar_type l=nd->getDistanceToFather();
  nd->setDistanceToFather(l*f);
  for(int i=0;i<nd->getNumberOfSons(); i++)
    scale_subtree(nd->getSon(i), f);
}

/* Returns the height of node nd (assumingly from a clock-like tree) */

scalar_type Species_tree::node_height(Node* nd){

  int nsons=nd->getNumberOfSons();
  if(nsons==0) return 0.;

  return node_height(nd->getSon(0)) + nd->getSon(0)->getDistanceToFather();
}

/* Makes uderlying tree clocklike, returns height */

scalar_type Species_tree::node_makeclocklike_height(Node* nd){

  scalar_type height=0.;
  int i, nsons;

  nsons=nd->getNumberOfSons();
  if(nsons==0) return 0.;

  vector<scalar_type> und_l(nsons);
  for(i=0;i<nsons;i++){
    Node* son=nd->getSon(i);
    und_l[i]=node_makeclocklike_height(son);
    und_l[i]+=son->getDistanceToFather();
    height+=und_l[i];
  }
  height/=nsons;

  for(i=0;i<nsons;i++)
    scale_subtree(nd->getSon(i), height/und_l[i]);

  return height;
}

//#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
//################### Borrowed from HGT_simul #############
//#########################################################

//#########################################################
//################### Some general S tree operatinos ######
//#########################################################

// "safe" distance to father
scalar_type Species_tree::distance_to_father(Node * node)
{
  if (node == root)
    return 0;
  else
    return node->getDistanceToFather();
}

// add marker nodes to display internal node names in newick format output 
void Species_tree::add_marker_nodes(Node * node)
{
  vector<Node * > sons = node->getSons();
  for (int i = 0;i < node->getNumberOfSons();i++)
    add_marker_nodes(sons[i]);
  if (!node->isLeaf() && node != root )
    {
      Node * tmp = new Node();
      Node * child = node->getSon(0);
      node->removeSon(child);
      tmp->setName(node->getName());
      // for display 
      node->addSon(tmp);
      node->addSon(child);
      tmp->setDistanceToFather(0.);

    }
  return;
}
void Species_tree::add_marker_node(Node * node,string name)
{
  if (!node->isLeaf() && node != root )
    {
      Node * tmp = new Node();
      Node * child = node->getSon(0);
      node->removeSon(child);
      tmp->setName(name);
      // for display 
      node->addSon(tmp);
      node->addSon(child);
      tmp->setDistanceToFather(0.);

    }
  return;


}


scalar_type Species_tree::distance_to_root(Node * node,Node * root)
{
  return TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*root,*node);
}

// rename leaves 0,1,2,..#leaves-1 to keep internal node names short
void Species_tree::rename_leaves(Node* node)
{
  int i=0;
  vector<Node * > leaves = TreeTemplateTools::getLeaves(* node);
  vector<Node *>::iterator it;
	
  for ( it = leaves.begin(); it!=leaves.end(); it++ )
    {
      stringstream out;
      out << i; 
      (*it)->setName(out.str());
      i+=1;
    }
}




// get names of all leaves in subtree below node
string Species_tree::get_leaf_names(Node * node)
{
  if (node->isLeaf())
    return node->getName();
  string name;
  vector<Node * > leaves = TreeTemplateTools::getLeaves(*node);
  vector<string> leaf_names;
  vector<Node *>::iterator it;
  for ( it = leaves.begin(); it!=leaves.end(); it++ )
    leaf_names.push_back( (*it)->getName());
  sort(leaf_names.begin(),leaf_names.end());
  vector<string>::iterator st;
  for ( st = leaf_names.begin(); st!=leaf_names.end(); st++ )
    name += (*st) + "+";
  if (name.size())
    name.erase(name.end()-1);
  return name;
}

// construct name of nodes based on leaves in subtree below them
string Species_tree::get_name_in_S(Node * node, string pre)
{
  string name;
  if (node==root)
    name = pre+"*R.";
  else if (node->isLeaf())
    name ="";//name = node->getName();
  else      
    name = pre+"*S.";
  return name + get_leaf_names(node);
}

scalar_type Species_tree::get_x_in_S(Node * node)
{
  vector<Node * > leaves = TreeTemplateTools::getLeaves(*node);
  vector<Node *>::iterator it;
  scalar_type tmp=0;
  for ( it = leaves.begin(); it!=leaves.end(); it++ )
    tmp += atof(((*it)->getName()).c_str())+1;
  tmp /= float( leaves.size() ) ; 
  return tmp;
}
  
void Species_tree::make_rate_table(bool reinit)
{
  branchwise_rates=reinit;
  if (reinit)
    {  
      
      tree = S_tree->clone();
      root = tree->getRootNode();
      name_internal(root);
    }
  Node * tree_node=root, * delta_node=delta_root, * tau_node=tau_root, * lambda_node=lambda_root;
  delta_table[tree_node->getName()] = distance_to_father(delta_node);
  tau_table[tree_node->getName()] = distance_to_father(tau_node);
  lambda_table[tree_node->getName()] = distance_to_father(lambda_node);
  string name = get_name_in_S(tree_node);
  delta_table["*Sp."+name] = distance_to_father(delta_node);
  tau_table["*Sp."+name] = distance_to_father(tau_node);
  lambda_table["*Sp."+name] = distance_to_father(lambda_node);
  vector<Node * > tree_sons = tree_node->getSons(), delta_sons = delta_node->getSons(),  tau_sons = tau_node->getSons(),  lambda_sons = lambda_node->getSons();
  for (int i = 0;i < tree_node->getNumberOfSons();i++)
    make_rate_table(tree_sons[i],delta_sons[i],tau_sons[i],lambda_sons[i]);
  return;

}
void Species_tree::make_rate_table(Node * tree_node,Node * delta_node,Node * tau_node,Node * lambda_node)
{
  delta_table[tree_node->getName()] = distance_to_father(delta_node);
  tau_table[tree_node->getName()] = distance_to_father(tau_node);
  lambda_table[tree_node->getName()] = distance_to_father(lambda_node);
  string name = get_name_in_S(tree_node);
  delta_table["*Sp."+name] = distance_to_father(delta_node);
  tau_table["*Sp."+name] = distance_to_father(tau_node);
  lambda_table["*Sp."+name] = distance_to_father(lambda_node);

  vector<Node * > tree_sons = tree_node->getSons(), delta_sons = delta_node->getSons(),  tau_sons = tau_node->getSons(),  lambda_sons = lambda_node->getSons();
  for (int i = 0;i < tree_node->getNumberOfSons();i++)
    make_rate_table(tree_sons[i],delta_sons[i],tau_sons[i],lambda_sons[i]);
  return;
}

//#########################################################
//################### Some general S tree operations ######
//#########################################################

//#########################################################
//################### S & S' tree construction ############
//#########################################################

// look for node that happened first after current time t  
void Species_tree::lookforCandidates(Node* node)
{
  scalar_type d = distance_to_root(node,root);
      
  // if we are after t return
  if ( d - t > 1e-5  && !isSpNode(node))
    {
      candidates[d] = node;
      return;
    }
  else
    for( int i = 0; i < node->getNumberOfSons(); i++ )
      lookforCandidates(node->getSon(i));
  return;
}

// find node defining next time slice given t time has passed
void Species_tree::findCandidates()
{
  candidates.clear();
  if (candidate_names.size()==0)
    {
      vector<Node *> nodes = tree->getNodes();
      for (vector<Node *>::iterator it=nodes.begin(); it!=nodes.end(); it++)
	if (! ((*it)->isLeaf()) && !((*it)==root))
	  {
	    if ( candidate_names.count(distance_to_root((*it),root)) )		
	      {
		(*it)->setDistanceToFather((*it)->getDistanceToFather()+2e-5);
		candidate_names[distance_to_root((*it),root)]=(*it)->getName();		
	      }
	    else
	      candidate_names[distance_to_root((*it),root)]=(*it)->getName();
	  }
      vector<Node *> leaves = tree->getLeaves();
      for (vector<Node *>::iterator it=leaves.begin(); it!=leaves.end(); it++)
	if (!((*it)==root))
	  {
	    //candidate_names[distance_to_root((*it),root)]=(*it)->getName();
	    candidate_names[distance_to_root((*it),root)+RandomTools::giveRandomNumberBetweenZeroAndEntry(1e-10)]=(*it)->getName();

	    if (infer_mode)	      
	      break;
	  }

    }
  for (map<scalar_type, string>::iterator it=candidate_names.begin(); it!=candidate_names.end(); it++)
	{
	vector <Node *> tmp;
	//cout << (*it).second << endl;
	TreeTemplateTools::searchNodeWithName(*root,(*it).second,tmp);
	if (tmp.size()>0)
	  {	    
	    candidates[(*it).first] = tmp[0];
	    break;
	  }	    
	else
	  candidate_names.erase(it);       
	}
  return;
  // use map to implicitly time order candidates 
  candidates.clear();
  // construct by DFS down to first nodes later than t 
  //cout << root->getName() << endl;
  lookforCandidates(root);
}

void Species_tree::add_leaf(Node * to, scalar_type d,string name)
{
  Node * new_node = new Node();
  to->addSon(new_node);

  new_node->setDistanceToFather(d);
  new_node->setName(name);

  return;
}

void Species_tree::add_speciation(Node * leaf)
{
  Node * new_node_1 = new Node();
  leaf->addSon(new_node_1);
  new_node_1->setDistanceToFather(0.);
  stringstream out;
  out << new_node_1->getId(); 
  new_node_1->setName(out.str());

  Node * new_node_2 = new Node();
  leaf->addSon(new_node_2);
  new_node_2->setDistanceToFather(0.);
  out << new_node_2->getId(); 
  new_node_2->setName(out.str());
    

  return;
}


// insert a node at Delta_t below tail of edge leading to the father of head
void Species_tree::insert_node(Node * head, scalar_type Delta_t, string name)
{
  scalar_type edge_t = distance_to_father(head);  
  Node * new_node = new Node();
  Node * head_father = head->getFather();
  
  head_father->removeSon(head);
  // set new_node's father as the father of head
  head_father->addSon(new_node);  
  // set head as son of new_node
  new_node->addSon(head);
  // set distances
  new_node->setDistanceToFather(Delta_t);
  
  head->setDistanceToFather(edge_t-Delta_t);
  new_node->setName(name);
  return;
}
void Species_tree::remove_leaf(Node * node)
{
  Node * father =  node->getFather();
  father->removeSon(node);
  delete node;
}
// remove a S prime type node from tree
void Species_tree::remove_node(Node * node)
{
  vector<Node *> sons = node->getSons();
  for( int i = 0; i < node->getNumberOfSons(); i++ )
    if (sons[i]->getDistanceToFather()==0)
      {
	node->removeSon(sons[i]);
	delete sons[i];
      }
  Node * father =  node->getFather();  
  Node * child = node->getSon(0);
  scalar_type d = TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*father,*child);
  father->removeSon(node);
  delete node;

  //cout << " up " <<endl;
  father->addSon(child);
  //cout << " down " <<endl;
  child->setDistanceToFather(d);

  return;
}

// test if a node is an S prime type node
bool Species_tree::isSpNode(Node * node)
{
  int nsons = node->getNumberOfSons();
  if (nsons==1 && node!=root)
    return true;
  else if(marker_nodes && nsons==2  && node!=root)
    {
      vector<Node *> sons = node->getSons();
      for( int i = 0; i < nsons; i++ )
	if (sons[i]->getDistanceToFather()==0)
	  return true;
    }
  return false;
}

//#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
//################### S & S' tree construction ############
//#########################################################

//#########################################################
//################### DTL events and necc. operations #####
//#########################################################


// prune subtree below degree one node node
void Species_tree::prune_at(Node* node)
{ 
  prune(node);

  return;
}

void Species_tree::prune(Node* node)
{
  if(node->getNumberOfSons()==0) return;    
  int nsons=node->getNumberOfSons();    
  for(int i=0;i<nsons;i++)
    prune(node->getSon(i));
  while (node->getNumberOfSons()!=0)
    node->removeSon( (unsigned int) 0);
  return;
}
  
Node * Species_tree::get_S_node(Node * node)
{    
  Node * S_node;
  if ( isSpNode(node) )
    S_node=node->getSons()[0];
  else
    S_node=node;
  return S_node;
}

// graft subtree below from to below event node to
void Species_tree::graft(Node * from, Node * to, Node * root)
{
  // head of the donor branch S_node must be a speciation node 
  from = get_S_node(from);
  // to is always an event node
  // tail of donor brach with head S_node and event node both start at the start of the slice
  scalar_type add_t;
  add_t = distance_to_root(from,root) - distance_to_root(to,root);
  // we clone the subtree defined by S_node
  Node * subtree_clone = TreeTemplateTools::cloneSubtree<Node>(*from);	  
  // and graft it to the event node to
  //subtree_clone->setFather(to);
  to->addSon(subtree_clone);
  subtree_clone->setDistanceToFather(add_t);

  return;
}

//D event Delta_t below start of time slice on branch with head node
void Species_tree::D(Node * node, scalar_type Delta_t)
{
  string event_name = "*D."+node->getName();
  (*event_stream) << "#event# " << event_name << " at time " << t+Delta_t << " " << Delta_t <<endl;
  insert_node(node, Delta_t, event_name);
  Node * new_node = node->getFather();
  graft(node , new_node, root);
  for (int i=0;i<2;i++)
    if(new_node->getSon(i)!=node)
      current_slice.insert(current_slice.begin(),new_node->getSon(i));
  D_count++;
  event_count++;
  return;
  
}  

// T event Delta_t below start of time slice, 
// transfer from branch with head from to branch with head to 
// corresponds to grafting on the branch with head from 
// the subtree that branch with head to defines 
// (transfer to "to" from "from" means subtree from "to" must be grafted to "from")
void Species_tree::T(Node * from, Node * to, scalar_type Delta_t)
{
  string event_name = "*T."+to->getName()+"<-"+from->getName();
  (*event_stream) << "#event# " << event_name << " at time " << t+Delta_t << endl;
  insert_node(from, Delta_t, event_name);
  Node* new_node = from->getFather();
  graft(to,new_node, root);
  for (int i=0;i<2;i++)
    if(new_node->getSon(i)!=from)
      current_slice.insert(current_slice.begin(),new_node->getSon(i));
  T_count++;
  event_count++;
  return;
}  
// combined TL event
void Species_tree::TL(Node * from, Node * to, scalar_type Delta_t)
{
  string event_name = "*T."+to->getName()+"<-"+from->getName();
  (*event_stream) << "#event# " << event_name << " at time " << t+Delta_t << endl;
  insert_node(from, Delta_t, event_name);
  Node* new_node = from->getFather();
  graft(to,new_node, root);
  for (int i=0;i<2;i++)
    if(new_node->getSon(i)!=from)
      current_slice.insert(current_slice.begin(),new_node->getSon(i));
  L(to,Delta_t);  
  TL_count++;
  event_count++;
  return;
}  

// L event Delta_t below start of time slice on branch with head node
void Species_tree::L(Node * node, scalar_type Delta_t)
{
  string event_name = "*L."+node->getName();
  (*event_stream) << "#event# " << event_name << " at time " << t+Delta_t << " " << Delta_t <<endl;
  insert_node(node, Delta_t, event_name);
  Node* new_node = node->getFather();
  vector<Node * > :: iterator erase_pos;
  for (vector<Node * > :: iterator slice_iter = current_slice.begin(); slice_iter != current_slice.end(); slice_iter++)
    if ((*slice_iter)==node)
      erase_pos=slice_iter;
  current_slice.erase(erase_pos);
  prune_at(new_node);
  L_count++;
  event_count++;
  return;
}

// construct time slice defined by node
// i.e. construct a vector of nodes that are heads of branches that are contemporary with node
void Species_tree::construct_slice(Node * node)
{
  current_slice.clear();
  if (node->isLeaf())
    {
      vector<Node * > leaves = tree->getLeaves();    
      for (vector<Node * > :: iterator jt = leaves.begin(); jt != leaves.end(); jt++ )
	if ((*jt)!=root)
	  current_slice.push_back((*jt));
      leaves.clear();
      return;
    }

  //cout << node->getName() << endl;
  scalar_type d_node = distance_to_root(node,root);
  scalar_type d_Fnode = distance_to_root(node,root) - distance_to_father(node); 

  vector<Node * > nodes = tree->getNodes();    
  for (vector<Node * > :: iterator jt = nodes.begin(); jt != nodes.end(); jt++ )
    {
      Node * new_node = *jt;
      scalar_type d_new_node_to_root = distance_to_root(new_node,root);
      scalar_type d_father_of_new_node_to_root = 0; 
      if (new_node != root)
	d_father_of_new_node_to_root = distance_to_root(new_node->getFather(),root);

      // for all nodes in the future
      // that have fathers in the past	  
      if ( node != new_node &&  d_node <= d_new_node_to_root  && d_father_of_new_node_to_root < d_node && node != new_node->getFather() )
	{
	  // insert an Sp node above jt at t_i - t_f(j)
	  if ( abs(d_node - d_new_node_to_root) > 1e-5 && !node->isLeaf() && node->getName()[1]!='L')
	    {
	      //string name =  new_node->getName();
	      string name = get_name_in_S(new_node);
	      insert_node(new_node, d_node - d_father_of_new_node_to_root, "*Sp."+name);
	      new_node = new_node->getFather();
	    }	      
	  // insert this into the current slice
	  if (count(current_slice.begin(),current_slice.end(),new_node)==0)
	    current_slice.push_back(new_node);
	}
    } 
  if (count(current_slice.begin(),current_slice.end(),node)==0)
    current_slice.push_back(node);

  nodes.clear();
  return;
}

// filter time slice to get valid T donors
vector< Node * > Species_tree::filter_slice(Node * to)
{
  vector<Node *> donors;
  for (vector<Node * > :: iterator slice_iter = current_slice.begin(); slice_iter != current_slice.end(); slice_iter++)
    // T not possible from 1. self, 2. gene tree branch leading to the same species tree node 
    if (*slice_iter!=to && (*slice_iter)->getName() != to->getName())
      donors.push_back(*slice_iter);
  return donors;
}
  
// simulate next event in current_slice
scalar_type Species_tree::simulate_slice( scalar_type l_j, Node * node_j, bool  no_events)
{
  if (no_events)
    return l_j;
    
  scalar_type sigma_j=0.;
  // sum the rate of all possible events
  vector< Node * > donors;
	    
  for (vector<Node * > :: iterator slice_iter = current_slice.begin(); slice_iter != current_slice.end(); slice_iter++)
    {
      
      //D rate per gene tree branch
      sigma_j += delta * delta_table[(*slice_iter)->getName()];
      //T rate per gene tree branch
      //with unifrom probability of donors
      donors = filter_slice(*slice_iter);
      if (donors.size() && tau>0)
	sigma_j += tau * tau_table[(*slice_iter)->getName()];
      if (donors.size() && eta>0)
	sigma_j += eta;

      //L rate per gene tree branch
      sigma_j += lambda * lambda_table[(*slice_iter)->getName()];
      //cout << (*slice_iter)->getName() << " " << lambda_table[(*slice_iter)->getName()] <<" | ";
    }
  //  cout << endl ;

  // draw time to next event
  scalar_type Delta_t;
  if (sigma_j>0.)
    Delta_t=RandomTools::randExponential(1./sigma_j);
  else
    {
      return l_j;
    }
  // if an event happened 
  if (Delta_t < l_j )
    {
      // draw what happened (standard event based simulation scheme)
      scalar_type stop = RandomTools::giveRandomNumberBetweenZeroAndEntry(sigma_j);
      sigma_j=0;

      for (vector<Node * > :: iterator slice_iter = current_slice.begin(); slice_iter != current_slice.end(); slice_iter++)
	{
	  //D
	  sigma_j += delta * delta_table[(*slice_iter)->getName()];
	  if (sigma_j>stop)
	    {	
	      D(*slice_iter,Delta_t);
	      break;
	    }
	  donors = filter_slice(*slice_iter);
	  if (donors.size() && tau>0)
	    {
	      //T
	      sigma_j += tau * tau_table[(*slice_iter)->getName()];
	      if (sigma_j>stop)
		{
		  //draw where from
		  vector <Node *> froms;
		  froms.push_back(root);
		  RandomTools::getSample(donors,froms);
		  Node * from = froms[0];		  
		  T(from,*slice_iter,Delta_t);
		  break;
		}
	    }
	  if (donors.size() && eta>0)
	    {
	      //TL
	      sigma_j += eta;
	      if (sigma_j>stop)
		{
		  //draw where from
		  vector <Node *> froms;
		  froms.push_back(root);
		  RandomTools::getSample(donors,froms);
		  Node * from = froms[0];		  
		  TL(from,*slice_iter,Delta_t);
		  break;
		}
	    }

	  // L
	  sigma_j += lambda * lambda_table[(*slice_iter)->getName()];
	  if (sigma_j>stop)
	    {
	      L(*slice_iter,Delta_t);
	      break;
	    }

	}
      return Delta_t;
    }
  else
    {
      return l_j; 
    }
}

//#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
//################### DTL events and necc. operations #####
//#########################################################
void Species_tree::init()
{
  // the tree to make a gene tree or Sp tree
  delete tree;
  tree = S_tree->clone();
  root = tree->getRootNode();
  time_slices.clear();
  // name speciation nodes
  name_internal(root);

  // by default 
  /*
  if (!branchwise_rates)
    {
      // trees with D,T and L rates on branches  
      delta_tree = S_tree->clone();
      tau_tree = S_tree->clone();
      lambda_tree = S_tree->clone();
      
      delta_root = delta_tree->getRootNode();
      tau_root = tau_tree->getRootNode();
      lambda_root = lambda_tree->getRootNode();

      name_internal(delta_root);
      name_internal(tau_root);
      name_internal(lambda_root);
      
      TreeTemplateTools::setBranchLengths(*delta_root,1);	
      TreeTemplateTools::setBranchLengths(*tau_root,1);	
      TreeTemplateTools::setBranchLengths(*lambda_root,1);	
      make_rate_table(false);
    }
  */
  // rename leaves to keep internal node names short
  //if (marker_nodes)
  //  rename_leaves();
  // scale tree to unit (max) height
  TreeTemplateTools::scaleTree(*root,1./TreeTemplateTools::getHeight(*root));
  return;
}
