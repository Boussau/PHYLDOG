#include "DTL.h"
#include "DTL_mpi.h"
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Seq/DistanceMatrix.h>
namespace mpi = boost::mpi;
using namespace std;
using namespace bpp;


vector<string> all_NNIs(string Sstring)
{
  tree_type * S = TreeTemplateTools::parenthesisToTree(Sstring);
  vector<string> NNIs;
  vector < Node * > nodes = S->getNodes();
  NNIs.push_back(TreeTemplateTools::treeToParenthesis(*S));
  if (nodes.size()<4)
    return NNIs;
  Node * root=S->getRootNode();
  for (vector < Node *> :: iterator ni=nodes.begin();ni!=nodes.end(); ni++)
    if (!((*ni)==root) && !((*ni)->isLeaf()))
      {
	Node * node0=(*ni)->getSon(0);
	Node * node1=(*ni)->getSon(1);
	Node * father=(*ni)->getFather();
	Node * node2;
	if (father->getSon(0)==(*ni))
	  node2=father->getSon(1);
	else
	  node2=father->getSon(0);
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	//NNI 1.
	(*ni)->addSon(node0);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node1);
       
	//clone out
	TreeTemplate<Node> * nni_S=S->clone();
	NNIs.push_back(TreeTemplateTools::treeToParenthesis(*nni_S));
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	
	//NNI 2.
	(*ni)->addSon(node1);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node0);

	//clone out
	nni_S=S->clone();
	NNIs.push_back(TreeTemplateTools::treeToParenthesis(*nni_S));
	//disolve
	(*ni)->removeSons();
	father->removeSons();

	//restore
	(*ni)->addSon(node0);
	(*ni)->addSon(node1);
	father->addSon((*ni));
	father->addSon(node2);
      }
  
  return NNIs; 
}

pair < scalar_type,string> ML_nni_step(Species_tree * infer_tree,string Sstring,string mode,bool greedy)
{

  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int mpi_rank = world.rank();
  int size = world.size();
  broadcast(world, Sstring,server);
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);

  // get inital set up from branch lengts
  Node * root = S->getRootNode();  
  //vector<Node * > nodes_tmp = S->getNodes();    
  vector<Node * > nodes= S->getNodes();
  //RandomTools::getSample(nodes_tmp,nodes);
  map <scalar_type,Node * > times;
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if (!(*ni)->isLeaf())
      {
	scalar_type h=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*(*ni),*root);
	if (times.count(h))
	  h+=1e-20;      
	times[h]=(*ni);
      }
  // register ranks 
  map <Node *,int> node_2_rank;
  map <int,Node *> rank_2_node;
  map <int,int> id_2_rank;
  map <int,int> rank_2_id;
  vector <pair <int,int> > moves;
  int rank=0;
  for (map < scalar_type,Node *> :: iterator hi = times.begin(); hi != times.end(); hi++ )
    {
      rank++;
      rank_2_node[rank]=(*hi).second;
      node_2_rank[(*hi).second]=rank;
      rank_2_id[rank]=(*hi).second->getId();
      id_2_rank[(*hi).second->getId()]=rank;

      stringstream out;out<<rank;
      (*hi).second->setBranchProperty("ID",BppString(out.str()));
    }

  // set branch lengths according to ranks 
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - node_2_rank[(*ni)->getFather()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( node_2_rank[(*ni)] - node_2_rank[(*ni)->getFather()] );
      
  scalar_type max_dLL=-1e30;

  /*
  // find moves
  for (map <int,Node *> :: iterator nri = rank_2_node.begin(); nri != rank_2_node.end(); nri++ )
    if ((*nri).first!=1)
      {
	Node * node=(*nri).second;
	int r=(*nri).first;
	if (rank_2_node[r-1]->getSon(0)!=node && rank_2_node[r-1]->getSon(1)!=node)
	  {
	    pair <int,int> move;
	    move.first=r-1;
	    move.second=r;
	    moves.push_back(move);
	  }
      }
  //cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID");

  scalar_type oLL=LL_mpi(infer_tree,S,mode);
  //XX
  //cout << oLL << endl;
  TreeTemplate<Node> * max_dS;
  vector<pair<int,int> > good_moves;

  map <scalar_type, TreeTemplate<Node> * > move_trees;
  nodes = S->getNodes();    
  root=S->getRootNode();

  // preform moves
  for (vector <pair <int,int> > :: iterator mi = moves.begin(); mi != moves.end(); mi++ )
    {
      //TreeTemplate<Node> * dS=S->clone();
      //vector<Node * > dS_nodes = dS->getNodes();    
      //Node * dS_root=dS->getRootNode();

      id_2_rank[rank_2_id[(*mi).first]] = (*mi).second;
      id_2_rank[rank_2_id[(*mi).second]] = (*mi).first;

      for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
	if  ((*ni)->isLeaf() )
	  (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
	else if ((*ni)!=root)
	  (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
      
      
      scalar_type dLL = LL_mpi(infer_tree,S,mode)-oLL;
      //cout << (*mi).first << " <-> " << (*mi).second << endl;
      //cout << dLL <<endl;
      id_2_rank[rank_2_id[(*mi).first]] = (*mi).first;
      id_2_rank[rank_2_id[(*mi).second]] = (*mi).second;

      for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
	if  ((*ni)->isLeaf() )
	  (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
	else if ((*ni)!=root)
	  (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
      if (dLL>max_dLL)
	{
	  max_dLL=dLL;
	  //max_dS=dS;
	}

      if(dLL>0)
	good_moves.push_back((*mi));

            
    }
  
  //cout << max_dLL << endl;
  map<int,int> moved;
  for (vector <pair <int,int> > :: iterator mi = good_moves.begin(); mi != good_moves.end(); mi++ )
    {
      if (moved.count((*mi).first)==0 && moved.count((*mi).first)==0)
	{
	  //cout << (*mi).first << " <-G-> " << (*mi).second << endl;
      
	  id_2_rank[rank_2_id[(*mi).first]] = (*mi).second;
	  id_2_rank[rank_2_id[(*mi).second]] = (*mi).first;
	  moved[(*mi).first]=1;
	  moved[(*mi).second]=1;
	}
    }

  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
  */
  if (mpi_rank==server) cout << Sstring<<endl;
  if (mpi_rank==server) cout << TreeTemplateTools::treeToParenthesis(*S)<<endl;
  if (mpi_rank==server) cout << can_order(S)<<endl;

  scalar_type oLL=LL_mpi(infer_tree,can_order(S),mode);
  if (mpi_rank==server) cout << "oLL " << oLL << endl;

  scalar_type maxdLL=0;
  int max_id=-1;
  int max_father_id=-1;

  //TreeTemplate<Node> * maxdS;
  string maxdS;
  root=S->getRootNode();
  for (vector < Node *> :: iterator ni=nodes.begin();ni!=nodes.end(); ni++)
    if (!((*ni)==root) && !((*ni)->isLeaf()) )
      {
	Node * node0=(*ni)->getSon(0);
	Node * node1=(*ni)->getSon(1);
	Node * father=(*ni)->getFather();
	Node * node2;

	if (father->getSon(0)==(*ni))
	  node2=father->getSon(1);
	else
	  node2=father->getSon(0);
	
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	//NNI 1.
	(*ni)->addSon(node0);       
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node1);
	//clone out
	int old_rank=id_2_rank[(*ni)->getId()];
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second+=1;
	id_2_rank[(*ni)->getId()]=id_2_rank[(*ni)->getFather()->getId()]+1;

	vector<Node * > dS_nodes = S->getNodes();    
	Node * dS_root=S->getRootNode();
	
	for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	  if  ((*nii)->isLeaf() )
	    (*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	  else if ((*nii)!=dS_root)
	    (*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		
	//TreeTemplate<Node> * dS=S->clone();
	string nni_string=can_order(S);	
	scalar_type dLL = LL_mpi(infer_tree,nni_string,mode)-oLL;
	//cout << get_name(node0) << " ^ "  << get_name(node2) << " -- " << get_name(node1) <<endl;
  
	if (dLL>maxdLL)
	  {
	    maxdLL=dLL;
	    maxdS=nni_string;//S->clone();
	    max_id=(*ni)->getId();
	    max_father_id=father->getId();
	    if (mpi_rank==server) cout << dLL << " " << oLL << endl; 
	    if (mpi_rank==server) cout << nni_string << endl; 
	    if (greedy) break;
	  }

	//disolve
	(*ni)->removeSons();
	father->removeSons();
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<old_rank+1 && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second-=1;
	id_2_rank[(*ni)->getId()]=old_rank;

	//NNI 2.
	(*ni)->addSon(node1);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node0);	
	
	//clone out
	old_rank=id_2_rank[(*ni)->getId()];
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second+=1;
	id_2_rank[(*ni)->getId()]=id_2_rank[(*ni)->getFather()->getId()]+1;

	dS_nodes = S->getNodes();    
	dS_root=S->getRootNode();
	
	for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	  if  ((*nii)->isLeaf() )
	    (*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	  else if ((*nii)!=dS_root)
	    (*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		
	//dS
	nni_string=can_order(S);	
	dLL = LL_mpi(infer_tree,nni_string,mode)-oLL;
	//cout << dLL << "nni move"<<endl;
	//cout << get_name(node1) << " ^ "  << get_name(node2) << " -- " << get_name(node0) <<endl;
	if (dLL>maxdLL)
	  {
	    maxdLL=dLL;
	    maxdS=nni_string;//S->clone();;
	    max_id=(*ni)->getId();
	    max_father_id=father->getId();
	    if (mpi_rank==server) cout << dLL << " ."<< endl; 
	    if (mpi_rank==server) cout << nni_string << endl; 

	    if (greedy) break;
	  }
	
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<old_rank+1 && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second-=1;
	id_2_rank[(*ni)->getId()]=old_rank;
	
	//restore
	(*ni)->addSon(node0);
	(*ni)->addSon(node1);
	father->addSon((*ni));
	father->addSon(node2);
      }
  times.clear();
  nodes.clear();
  node_2_rank.clear();
  rank_2_node.clear();
  id_2_rank.clear();
  rank_2_id.clear();
  moves.clear();  

  pair < scalar_type,string> ML_step;

  delete S;
  ML_step.first=maxdLL;

  if (maxdLL>0)
    ML_step.second=maxdS;
  else
    ML_step.second=Sstring;
  return ML_step;
 /* if (maxdLL>0)
    {
      S=maxdS;
      max_dLL=maxdLL;
      int old_rank=id_2_rank[max_id];
      for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	if ((*ii).second<id_2_rank[max_id] && (*ii).second>id_2_rank[max_father_id])
	  (*ii).second+=1;
      id_2_rank[max_id]=id_2_rank[max_father_id]+1;      
    }
  nodes = S->getNodes();    
  root=S->getRootNode();

  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )    
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
  pair < scalar_type,string> ML_step;
  
  ML_step.first=max_dLL;
  ML_step.second=can_order(S);
  //ML_step.second=maxdS;

  return ML_step;
*/
}


TreeTemplate<Node> * random_nni_step(TreeTemplate<Node> * S)
{

  // get inital set up from branch lengts
  Node * root = S->getRootNode();  
  vector<Node * > nodes = S->getNodes();    
  map <scalar_type,Node * > times;
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if (!(*ni)->isLeaf())
      {
	scalar_type h=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*(*ni),*root);
	if (times.count(h))
	  h+=1e-20;      
	times[h]=(*ni);
      }
  // register ranks 
  map <Node *,int> node_2_rank;
  map <int,Node *> rank_2_node;
  map <int,int> id_2_rank;
  map <int,int> rank_2_id;
  vector <pair <int,int> > moves;
  int rank=0;
  for (map < scalar_type,Node *> :: iterator hi = times.begin(); hi != times.end(); hi++ )
    {
      rank++;
      rank_2_node[rank]=(*hi).second;
      node_2_rank[(*hi).second]=rank;
      rank_2_id[rank]=(*hi).second->getId();
      id_2_rank[(*hi).second->getId()]=rank;

      stringstream out;out<<rank;
      (*hi).second->setBranchProperty("ID",BppString(out.str()));
    }
  // set branch lengths according to ranks 
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - node_2_rank[(*ni)->getFather()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( node_2_rank[(*ni)] - node_2_rank[(*ni)->getFather()] );
      
  // find moves
  for (map <int,Node *> :: iterator nri = rank_2_node.begin(); nri != rank_2_node.end(); nri++ )
    if ((*nri).first!=1)
      {
	Node * node=(*nri).second;
	int r=(*nri).first;
	if (rank_2_node[r-1]->getSon(0)!=node && rank_2_node[r-1]->getSon(1)!=node)
	  {
	    pair <int,int> move;
	    move.first=r-1;
	    move.second=r;
	    //cout<<move.first << " " << move.second << endl;
	    moves.push_back(move);
	  }
      }

  vector <pair <int,int> > rmove;
  rmove.push_back(moves.front());
  RandomTools::getSample(moves,rmove);
  pair<int,int> move = rmove[0];		  
  
  id_2_rank[rank_2_id[move.first]] = move.second;
  id_2_rank[rank_2_id[move.second]] = move.first;
  
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
  
  root=S->getRootNode();

  for (vector < Node *> :: iterator ni=nodes.begin();ni!=nodes.end(); ni++)
    if (!((*ni)==root) && !((*ni)->isLeaf()) )
      {
	Node * node0=(*ni)->getSon(0);
	Node * node1=(*ni)->getSon(1);
	Node * father=(*ni)->getFather();
	Node * node2;

	if (father->getSon(0)==(*ni))
	  node2=father->getSon(1);
	else
	  node2=father->getSon(0);
	
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	//NNI 1.
	(*ni)->addSon(node0);       
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node1);
	//clone out
	int old_rank=id_2_rank[(*ni)->getId()];
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second+=1;
	id_2_rank[(*ni)->getId()]=id_2_rank[(*ni)->getFather()->getId()]+1;
	

	vector<Node * > dS_nodes = S->getNodes();    
	Node * dS_root=S->getRootNode();
	
	for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	  if  ((*nii)->isLeaf() )
	    (*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	  else if ((*nii)!=dS_root)
	    (*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		

	if (RandomTools::giveIntRandomNumberBetweenZeroAndEntry(400)==1)
	  break;
  
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<old_rank+1 && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second-=1;
	id_2_rank[(*ni)->getId()]=old_rank;

	//NNI 2.
	(*ni)->addSon(node1);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node0);
		
	//clone out
	old_rank=id_2_rank[(*ni)->getId()];
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second+=1;
	id_2_rank[(*ni)->getId()]=id_2_rank[(*ni)->getFather()->getId()]+1;
	
	dS_nodes = S->getNodes();    
	dS_root=S->getRootNode();
	
	for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	  if  ((*nii)->isLeaf() )
	    (*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	  else if ((*nii)!=dS_root)
	    (*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		

	if (RandomTools::giveIntRandomNumberBetweenZeroAndEntry(400)==1)
	  break;
	
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	(*ni)->removeSons();
	father->removeSons();
	for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	  if ((*ii).second<old_rank+1 && (*ii).second>id_2_rank[father->getId()])
	    (*ii).second-=1;
	id_2_rank[(*ni)->getId()]=old_rank;

	
	//restore
	(*ni)->addSon(node0);
	(*ni)->addSon(node1);
	father->addSon((*ni));
	father->addSon(node2);
      }

  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );


  times.clear();
  nodes.clear();
  node_2_rank.clear();
  rank_2_node.clear();
  id_2_rank.clear();
  rank_2_id.clear();
  //for (vector <pair <int,int> > :: iterator it=moves.begin();it!=moves.end();it++)
  //  (*it).clear();
  moves.clear();
  return S;  
}


string randomize_nni(string Sstring,int steps)
{
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);
  for (int i=0;i<steps;i++)
    {
      //cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID");      
      S=random_nni_step(S);  
    }
  string tmp=can_order(S);
  delete S;
  return tmp;
}


string get_name(Node * node)
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


string clock_root(Species_tree * infer_tree,string Sstring,string mode)
{

  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int mpi_rank = world.rank();
  int size = world.size();
  broadcast(world, Sstring,server);
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);

  map < Node*,int > do_01;
  map < Node*,int > do_12;
  map < Node*,int > do_02;
  map < pair <Node*,Node*>,int > dnodes;
  vector< pair <Node*,Node*> > seen;
  map < pair <node_type*,node_type*>, pair <node_type*,node_type*> > root_dnodes;
  map < pair <Node*,Node*>,int > roots;

  S->unroot();
  vector <Node *> S_leaves=S->getLeaves();
  vector<Node*>::iterator it,jt;
  vector <pair <Node*,Node*> > roots_vec;

  for (it = S_leaves.begin(); it!=S_leaves.end(); it++ )
    {
      Node * node = (*it);
      Node * head = node->getFather();
      if (S_leaves.size()==2)
	if (S_leaves[0]==(*it))
	  head=S_leaves[1];
	else
	  head=S_leaves[0];
      pair <Node*,Node*> tmp_pair(node,head);
      if (node->getId()>head->getId())
	{tmp_pair.first=head;tmp_pair.second=node;}
      roots[tmp_pair]++;		  
      if (roots[tmp_pair]==2)
	{
	  Node* a_root=new Node();
	  //a_root->setName(tmp_pair.first->getName() + " -0- " + tmp_pair.second->getName());
	  pair <Node*,Node*> root_dnode(a_root,a_root);
	  pair <Node*,Node*> left_pair(tmp_pair.first,tmp_pair.second);
	  pair <Node*,Node*> right_pair(tmp_pair.second,tmp_pair.first);
	  root_dnodes[root_dnode]=tmp_pair;
	  roots_vec.push_back(root_dnode);

	}
      vector <Node*> neigbours = head->getSons();
      if (head->hasFather())
	neigbours.push_back(head->getFather());
      int node_id;
      if (node==neigbours[0])
	{
	  node_id=0;
	  do_01[head]++;
	  do_02[head]++;	
	}
      else if (node==neigbours[1])
	{
	  node_id=1;
	  do_01[head]++;
	  do_12[head]++;	
	}
      else if (node==neigbours[2])
	{
	  node_id=2;       	
	  do_12[head]++;
	  do_02[head]++;	
	}
      pair <Node*,Node*> left_pair(node,head);
      pair <Node*,Node*> right_pair(node,head);     
      if ((node_id==0 || node_id==1) && do_01[head]==2)
	{
	  pair <Node*,Node*> head_pair(head,neigbours[2]);
	  if (node_id==0)
	    right_pair.first=neigbours[1];
	  else 
	    right_pair.first=neigbours[0];
	}
      if ((node_id==1 || node_id==2) && do_12[head]==2)
	{
	  pair <Node*,Node*> head_pair(head,neigbours[0]);
	  if (node_id==1)
	    right_pair.first=neigbours[2];
	  else 
	    right_pair.first=neigbours[1];
	}
      if ((node_id==0 || node_id==2) && do_02[head]==2)
	{
	  pair <Node*,Node*> head_pair(head,neigbours[1]);
	  if (node_id==0)
	    right_pair.first=neigbours[2];
	  else 
	    right_pair.first=neigbours[0];
	}
            
      for (jt=neigbours.begin();jt!=neigbours.end();jt++)
	if (node!=(*jt))
	  {
	    pair <Node*,Node*> dnode (head,(*jt));
	    dnodes[dnode]++;
	  }     
      neigbours.clear(); 
    }

  while(1)
    {
      for (map < pair <Node*,Node*>,int >::iterator dit=dnodes.begin(); dit!=dnodes.end(); dit++)
	if ((*dit).second==2)
	  {	    
	    Node * node = (*dit).first.first;
	    Node * head = (*dit).first.second;
	    pair <Node*,Node*> tmp_pair(node,head);
	    if (node->getId()>head->getId())
	      {tmp_pair.first=head;tmp_pair.second=node;}
	    roots[tmp_pair]++;	
	    if (roots[tmp_pair]==2)
	      {
		Node* a_root=new Node();
		//a_root->setName(tmp_pair.first->getName() + " -0- " + tmp_pair.second->getName());
		pair <Node*,Node*> root_dnode(a_root,a_root);
		pair <Node*,Node*> left_pair(tmp_pair.first,tmp_pair.second);
		pair <Node*,Node*> right_pair(tmp_pair.second,tmp_pair.first);
		root_dnodes[root_dnode]=tmp_pair;
		roots_vec.push_back(root_dnode);
	      }
	    vector <Node*> neigbours = head->getSons();
	    if (head->hasFather())
	      neigbours.push_back(head->getFather());
	    int node_id;
	    if (node==neigbours[0])
	      {
		node_id=0;
		do_01[head]++;
		do_02[head]++;	
	      }
	    else if (node==neigbours[1])
	      {
		node_id=1;
		do_01[head]++;
		do_12[head]++;	
	      }
	    else if (node==neigbours[2])
	      {
		node_id=2;       	
		do_12[head]++;
		do_02[head]++;	
	      }
	    pair <Node*,Node*> left_pair(node,head);
	    pair <Node*,Node*> right_pair(node,head);     
	    if ((node_id==0 || node_id==1) && do_01[head]==2)
	      {
		pair <Node*,Node*> head_pair(head,neigbours[2]);
		if (node_id==0)
		  right_pair.first=neigbours[1];
		else 
		  right_pair.first=neigbours[0];
	      }
	    if ((node_id==1 || node_id==2) && do_12[head]==2)
	      {
		pair <Node*,Node*> head_pair(head,neigbours[0]);
		if (node_id==1)
		  right_pair.first=neigbours[2];
		else 
		  right_pair.first=neigbours[1];
	      }
	    if ((node_id==0 || node_id==2) && do_02[head]==2)
	      {
		pair <Node*,Node*> head_pair(head,neigbours[1]);
		if (node_id==0)
		  right_pair.first=neigbours[2];
		else 
		  right_pair.first=neigbours[0];
	      }      
	    for (jt=neigbours.begin();jt!=neigbours.end();jt++)
	      if (node!=(*jt))
		{
		  pair <Node*,Node*> dnode (head,(*jt));
		  dnodes[dnode]++;
		}
	    seen.push_back((*dit).first);
	  }
      if (seen.size()==0)
	break;
      for (vector< pair <Node*,Node*> >::iterator sit=seen.begin();sit!=seen.end();sit++)
	dnodes.erase((*sit));
      seen.clear();
    }


  
  int root_i=0;

  string max_S;
  scalar_type max_ll=-1e30;
  for (map < pair <node_type*,node_type*>, pair <node_type*,node_type*> >::iterator rt=root_dnodes.begin();rt!=root_dnodes.end();rt++)
    {
      root_i+=1;
      Node* head=(*rt).second.first;
      Node* node=(*rt).second.second;
      vector <Node*> head_neigbours = head->getSons();
      if (head->hasFather())
	head_neigbours.push_back(head->getFather());
      vector <Node*> node_neigbours = node->getSons();
      if (head->hasFather())
	node_neigbours.push_back(node->getFather());
      int node_id=-1;
      if (node==head_neigbours[0])
	node_id=0;
      else if (node==head_neigbours[1])
	node_id=1;
      else if (node==head_neigbours[2])
	node_id=2;
      int did;
      Node * new_root = new Node();
      new_root->setName("R");
      if (( node_id==2 && head->hasFather() ) || head->isLeaf())
	{
	  did=1;
	  node->removeSon(head);
	  node->addSon(new_root);
	  new_root->addSon(head);
	  new_root->setDistanceToFather(1);
	  head->setDistanceToFather(1);

	}
      else
	{
	  did=0;
	  head->removeSon(node);
	  head->addSon(new_root);
	  new_root->addSon(node);
	  new_root->setDistanceToFather(1);
	  node->setDistanceToFather(1);

	}
      S->rootAt(new_root);
      //stringstream out;
      //out << root_i;
      //string fname=forest_file_name+".root"+out.str();
      //tree_stream.open(fname.c_str());
      vector <node_type*> tmp=S->getNodes();
      for (vector <node_type*>::iterator it=tmp.begin();it!=tmp.end();it++)
	if ((*it)->hasFather())
	{	 
	  scalar_type d=(*it)->getDistanceToFather();
	  (*it)->setDistanceToFather(1);
	}
      
      string clock_S =DTL_clock(TreeTemplateTools::treeToParenthesis(*S),infer_tree);
      scalar_type ll=LL_mpi(infer_tree,clock_S,mode);
      if (max_ll<ll)
	{
	  max_S=clock_S;
	  max_ll=ll;
	  }
      //tree_stream<<TreeTemplateTools::treeToParenthesis(*S);     
      //tree_stream.close();
    }
  return max_S;

}


string PANJ(string Sstring, vector<string> treestrings)
{  
  
  map < int , vector < scalar_type > > presence;
  map <string,int> map_pa;
  vector <string> names;
  vector <Node *> leaves = TreeTemplateTools::parenthesisToTree(Sstring)->getLeaves();
  int id=0;
  int N_leaves=0;
  for (vector<Node *>::iterator it=leaves.begin();it!=leaves.end();it++)
    {
      string name=(*it)->getName();
      if (name!="OUT_GROUP")
	{	  
	  names.push_back(name);
	  map_pa[name]=N_leaves;
	  for (vector <string>::iterator it=treestrings.begin();it!=treestrings.end();it++)
	    presence[N_leaves].push_back(0.);
	  N_leaves+=1;
	}
    }
  int tree_i=0;
  for (vector <string>::iterator it=treestrings.begin();it!=treestrings.end();it++)
    {
      vector <Node *> leaves =TreeTemplateTools::parenthesisToTree((*it))->getLeaves();
      for (vector<Node *>::iterator jt=leaves.begin();jt!=leaves.end();jt++)
	{
	  string name=(*jt)->getName();
	  if (name!="OUT_GROUP")
	    {	 
	      presence[map_pa[name]][tree_i]++;
	    }
	}
      tree_i++;
    }
  
  DistanceMatrix * d =new DistanceMatrix(names);

  for (int i=0;i<N_leaves;i++)
    for (int j=0;j<N_leaves;j++)
      (*d)(i,j)=0;
  
  for (int i=0;i<N_leaves;i++)
    for (int j=0;j<N_leaves;j++)
	for (int k=0;k<tree_i;k++)
	  //if(presence[i][k]>0 && presence[j][k]==0 || presence[j][k]>0 && presence[i][k]==0)
	  (*d)(i,j)+=((presence[j][k]-presence[i][k])*(presence[j][k]-presence[i][k]));//max(presence[j][k],presence[i][k])/(double)tree_i;	

  
  NeighborJoining  * bj =new NeighborJoining (*d,true,true);
  //BioNJ  * bj =new BioNJ (*d,true,true);
  bj->computeTree(true);
  tree_type * pabj_tree=bj->getTree();
  TreeTools::midpointRooting(*pabj_tree);
  
  vector <Node*> nodes= pabj_tree->getNodes();
  Node * root= pabj_tree->getRootNode();
  for (vector<Node * > :: iterator it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      if ((*it)->getDistanceToFather()==0)
	(*it)->setDistanceToFather(0.001);
  
  //node_makeclocklike_height(( pabj_tree->getRootNode()));    
  //pabj_tree->scaleTree(1./TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*root,*( pabj_tree->getLeaves()[0])));      
  

  string pabj_tree_string=TreeTemplateTools::treeToParenthesis(*(pabj_tree));
  return  pabj_tree_string;
  
  
}
