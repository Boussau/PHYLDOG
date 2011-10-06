#include "DTL.h"
#include "DTL_mpi.h"
using namespace std;
using namespace bpp;
namespace mpi = boost::mpi;

scalar_type LL(Species_tree * infer_tree,TreeTemplate<Node> * S, scalar_type delta, scalar_type tau, scalar_type lambda, scalar_type beta,string mode, scalar_type stem_omega,scalar_type out_tau)
{

  stringstream out;
  out << can_order(S);
  out << "|" << mode;out << "|" << delta;out << "|" << tau;out << "|" << lambda;out << "|" << stem_omega;
  scalar_type ll,ll2;
  if (infer_tree->LL_cache.count(out.str())==0 || 1) //fix!
    {
      infer_tree->intp=0.;
      infer_tree->w=0.;
      
      if (mode=="Hest")
	{
	  infer_tree->reset();
	  delta=0.01;
	  tau=0.01;
	  lambda=0.01;

	  infer_tree->branchwise=0;

	  ll= oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
	  infer_tree->branchwise=0;
	  infer_tree->count_events(infer_tree->recs,true,true);
	  infer_tree->recs.clear();infer_tree->peas.clear();
	  cout<<"#p "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	  cout<<"#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p "<<endl;
	  infer_tree->print_rates();		  		  	  
	  cout<<"#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p#p "<<endl;

	  int tmpc=1;
	  while (tmpc<40 )
	    {
	      tmpc+=1;
	      scalar_type tmp=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
	      cout<<"#bw "<< tmp  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	      ll2=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");
	      scalar_type tt;//=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate");

	      infer_tree->branchwise=1;	      
	      cout<<"#e "<< ll2  << " " << tt << " "<< infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	  
	      
	      //cout << ll2 << endl;
	      
	      if (ll2>ll  && abs(ll2-ll)>0)
		{
		  infer_tree->remember_rates();
		  infer_tree->count_events(infer_tree->recs,true,true);
		  infer_tree->recs.clear();infer_tree->peas.clear();

		  ll=ll2;
		}
	      else
		{
		  infer_tree->recall_rates();	
		  infer_tree->print_rates();		  		  	  
		  infer_tree->recs.clear();infer_tree->peas.clear();
		  break;
		}
	    }
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;	  
	  cout << ll << endl;
	}    
      else if (mode=="IHest")
	{
	  infer_tree->reset();
	  delta=0.01;
	  tau=0.01;
	  lambda=0.01;

	  infer_tree->branchwise=0;
	  infer_tree->reset();
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
	  //cout << ll << endl;
	  infer_tree->branchwise=1;
	  infer_tree->count_events(infer_tree->recs,false,true);
	  infer_tree->recs.clear(); infer_tree->peas.clear();
	  cout<<"#p "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
		  
	  int tmpc=1;
	  while (tmpc<40 )
	    {
	      tmpc++;
	      scalar_type tmp=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
	      cout<<"#bw "<< tmp  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	      ll2=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");
	      cout<<"#e "<< ll2  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;

	      infer_tree->branchwise=1;
	      if (ll2>ll  && abs(ll2-ll)>0.1)
		{
		  infer_tree->remember_rates();
		  infer_tree->count_events(infer_tree->recs,false,true);
		  infer_tree->recs.clear();infer_tree->peas.clear();
		  ll=ll2;
		}
	      else
		{
		  infer_tree->recall_rates();		  
		  infer_tree->print_rates();		 
		  infer_tree->recs.clear();infer_tree->peas.clear(); 		  
		  break;
		}
	      //cout<<"#l "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	      
	    }
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;	  
	  cout << ll << endl;

	}     
      else if (mode=="estimate")
	{
	  infer_tree->branchwise=0;
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"estimate") ;	  	  	 
	  //cout<<ll<<endl;
	  infer_tree->branchwise=1;
	  int tmpc=1;
	  while (tmpc<40)
	    {
	      tmpc++;
	      ll2=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");	  	  	 
	      //cout<<ll2<< " " << infer_tree->branchwise_omega[1] << endl;
	      if (ll2 > ll && abs(ll2-ll)>1)
		ll=ll2;
	      else
		break;
	    }
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;	  

	}
      else
	{
	  infer_tree->branchwise=1;
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;	  
	  //cout <<ll <<endl;
	}
      infer_tree->LL_cache[out.str()]=ll;
    }
  else
    {
      ll=infer_tree->LL_cache[out.str()];
    }
  cout << ll <<endl;
  return ll;//infer_tree->LL_cache[out.str()];

}
scalar_type oLL(Species_tree * infer_tree,TreeTemplate<Node> * S, scalar_type delta, scalar_type tau, scalar_type lambda, scalar_type beta,string mode, scalar_type stem_omega,scalar_type out_tau)
{

  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  scalar_type LL;
      
  infer_tree->mode=mode;
  infer_tree->reconstruct(S);
  infer_tree->init_x();
  infer_tree->beta=beta;  
  infer_tree->init_L_improved(delta,tau,lambda,0.1,0.01,mode);   
  LL=infer_tree->unrooted_run(0,mode);     
  return LL; 
}


pair <scalar_type,scalar_type> observe_BD(Species_tree * infer_tree, TreeTemplate<Node> * S, vector<string> trees,scalar_type delta,scalar_type tau,scalar_type lambda)
{
  //cout << delta << " " << tau <<" "<< lambda <<endl; 
  //cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID");

  infer_tree->reconstruct(S);
  infer_tree->init_L_improved(delta,tau,lambda); 

  scalar_type zerotree= infer_tree->zerotree;
  scalar_type onetree= infer_tree->onetree;

  //cout <<zerotree <<" " << onetree << endl;
  
  vector < Node * > S_leaves = S->getLeaves();
  map <string,scalar_type> leafnames;
  map <string,scalar_type> see0;
  map <string,scalar_type> seeavg;
  for (vector < Node *> :: iterator li=S_leaves.begin();li!=S_leaves.end(); li++)
    leafnames[(*li)->getName()]=1;
  scalar_type c=0;

  for (vector < string> :: iterator tss=trees.begin();tss!=trees.end(); tss++)
    {      
      tree_type * G_tree= TreeTemplateTools::parenthesisToTree((*tss));
      vector < Node * > G_leaves = G_tree->getLeaves();
      map <string,scalar_type> Gleafnames;
      for (vector < Node *> :: iterator li=G_leaves.begin();li!=G_leaves.end(); li++)
	Gleafnames[(*li)->getName()]+=1;
      for (map < string,scalar_type> :: iterator li=leafnames.begin();li!=leafnames.end(); li++)
	{
	  if (Gleafnames.count((*li).first)==0)
	    {
	      see0[(*li).first]+=1;	    
	    }
	  else
	    seeavg[(*li).first]+=Gleafnames[(*li).first];	    
	}
      c+=1;
    }
  c+=(onetree+zerotree)*trees.size();
  scalar_type cm=0;
  scalar_type ml=0;
  scalar_type md=0;

  for (vector < Node *> :: iterator li=S_leaves.begin();li!=S_leaves.end(); li++)
    {
      see0[(*li)->getName()]+=zerotree*trees.size()+onetree*trees.size()*((scalar_type)S_leaves.size()-1)/(scalar_type)S_leaves.size();
      seeavg[(*li)->getName()]+=onetree*trees.size()/(scalar_type)S_leaves.size();
      
      scalar_type f0=see0[(*li)->getName()]/c;
      scalar_type p=f0;
      scalar_type m=seeavg[(*li)->getName()]/c;
      scalar_type q=1.- (1.-f0)/m ;
      
      scalar_type obs_lambda=( p * log((q-1)/(p-1))/((p-q)*1));
      scalar_type obs_delta=q/p * obs_lambda;
      //cout <<  obs_lambda << " " <<  obs_delta << endl;
      cm+=1;
      ml+=obs_lambda;
      md+=obs_delta;
    }
  //cout << ml/cm << " " << md/cm << endl; 

  scalar_type obs_lambda=ml/cm;
  scalar_type obs_delta=md/cm;
  
  pair <scalar_type,scalar_type> obs;
  obs.first=obs_lambda;
  obs.second=obs_delta;
  return obs;
}

pair < scalar_type,string> ML_step(Species_tree * infer_tree,string Sstring,string mode, string move_mode, bool greedy)
{
  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int mpi_rank = world.rank();
  int size = world.size();
  broadcast(world, Sstring,server);
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);
 

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
  //cout << "WTF" <<endl;
  //for (map <int,int> :: iterator nri = id_2_rank.begin(); nri != id_2_rank.end(); nri++ )
  //  cout << (*nri).second << " ";
  //cout << endl;

  // set branch lengths according to ranks 
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - node_2_rank[(*ni)->getFather()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( node_2_rank[(*ni)] - node_2_rank[(*ni)->getFather()] );
      
  scalar_type max_dLL=-1e30;

  if (move_mode=="order" || move_mode=="joint")
    {
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
      scalar_type oLL=LL_mpi(infer_tree,can_order(S),mode);
      //cout << oLL << endl;

      TreeTemplate<Node> * max_dS;
      vector<pair<int,int> > good_moves;

      map <scalar_type, TreeTemplate<Node> * > move_trees;
      vector<Node * > nodes = S->getNodes();    
      Node * root=S->getRootNode();

      // preform moves
      for (vector <pair <int,int> > :: iterator mi = moves.begin(); mi != moves.end(); mi++ )
	{
	  //cout << (*mi).first << " <-> " << (*mi).second << endl;
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
      
      
	  scalar_type dLL = LL_mpi(infer_tree,can_order(S),mode)-oLL;
	  //cout << dLL << " "  << (*mi).first << " <-> " << (*mi).second << endl;

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

	  if(dLL>0)// && dLL>max_dLL)
	    {
	      good_moves.push_back((*mi));
	      if (infer_tree->world.rank()==0) 
		{
		  ofstream cout_stream(infer_tree->outstream_file_name.c_str(),ios::app);  
		  cout_stream << (*mi).first << " <-> " << (*mi).second << endl;
		  cout_stream.close();
		}

	      if (greedy) break;
	    }
	}
  
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
    }
  

  if (move_mode=="root" || move_mode=="joint")
    {
      scalar_type oLL=LL_mpi(infer_tree,can_order(S),mode);

      scalar_type maxdLL=0;
      int max_id=-1;
      TreeTemplate<Node> * maxdS;
      root=S->getRootNode();
      vector < Node *> sons=root->getSons();
      for (vector < Node *> :: iterator ni=sons.begin();ni!=sons.end(); ni++)
	if (!((*ni)->isLeaf()))
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
	      if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>1)
		(*ii).second+=1;
	    id_2_rank[(*ni)->getId()]=2;

	    vector<Node * > dS_nodes = S->getNodes();    
	    Node * dS_root=S->getRootNode();
	
	    for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	      if  ((*nii)->isLeaf() )
		(*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	      else if ((*nii)!=dS_root)
		(*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		
	    //TreeTemplate<Node> * dS=S->clone();
	    scalar_type dLL = LL_mpi(infer_tree,can_order(S),mode)-oLL;

	    //cout << dLL << "root move"<<endl;
  
	    if (dLL>maxdLL)
	      {
		maxdLL=dLL;
		maxdS=S->clone();
		max_id=(*ni)->getId();
	      }

	    //disolve
	    (*ni)->removeSons();
	    father->removeSons();
	    for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	      if ((*ii).second<old_rank+1 && (*ii).second>1)
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
	      if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>1)
		(*ii).second+=1;
	    id_2_rank[(*ni)->getId()]=2;

	    dS_nodes = S->getNodes();    
	    dS_root=S->getRootNode();
	
	    for (vector<Node * > :: iterator nii = dS_nodes.begin(); nii != dS_nodes.end(); nii++ )
	      if  ((*nii)->isLeaf() )
		(*nii)->setDistanceToFather( rank+1 - id_2_rank[(*nii)->getFather()->getId()] );
	      else if ((*nii)!=dS_root)
		(*nii)->setDistanceToFather( id_2_rank[(*nii)->getId()] - id_2_rank[(*nii)->getFather()->getId()] );		
	    //dS=S->clone();
	    dLL = LL_mpi(infer_tree,can_order(S),mode)-oLL;
	    //cout << dLL << "root move"<<endl;

	    if (dLL>maxdLL)
	      {
		maxdLL=dLL;
		maxdS=S->clone();
		max_id=(*ni)->getId();
	      }
	
	    //disolve
	    (*ni)->removeSons();
	    father->removeSons();
	    for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	      if ((*ii).second<old_rank+1 && (*ii).second>1)
		(*ii).second-=1;
	    id_2_rank[(*ni)->getId()]=old_rank;
	
	    //restore
	    (*ni)->addSon(node0);
	    (*ni)->addSon(node1);
	    father->addSon((*ni));
	    father->addSon(node2);
	  }
      if (maxdLL>0)
	{
	  S=maxdS;
	  max_dLL=maxdLL;
	  int old_rank=id_2_rank[max_id];
	  for (map<int,int> :: iterator ii = id_2_rank.begin(); ii != id_2_rank.end(); ii++ )
	    if ((*ii).second<id_2_rank[max_id] && (*ii).second>1)
	      (*ii).second+=1;
	  id_2_rank[max_id]=2;      

	}
      nodes = S->getNodes();    
      root=S->getRootNode();

      for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )    
	if  ((*ni)->isLeaf() )
	  (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
	else if ((*ni)!=root)
	  (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
    }
  pair < scalar_type,string> ML_step;
  
  ML_step.first=max_dLL;
  ML_step.second=can_order(S);
  //ML_step.second=maxdS;

  times.clear();
  nodes.clear();
  node_2_rank.clear();
  rank_2_node.clear();
  id_2_rank.clear();
  rank_2_id.clear();
  //for (vector <pair <int,int> > :: iterator it=moves.begin();it!=moves.end();it++)
  //  (*it).clear();
  moves.clear();  
  return ML_step; 
}

pair < scalar_type,string> last_step(Species_tree * infer_tree,string Sstring,string mode,int supp,bool greedy)
{

  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int mpi_rank = world.rank();
  int size = world.size();
  broadcast(world, Sstring,server);
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);

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
  //cout << "WTF" <<endl;
  //for (map <int,int> :: iterator nri = id_2_rank.begin(); nri != id_2_rank.end(); nri++ )
  //  cout << (*nri).second << " ";
  //cout << endl;

  // set branch lengths according to ranks 
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if  ((*ni)->isLeaf() )
      (*ni)->setDistanceToFather( rank+1 - node_2_rank[(*ni)->getFather()] );
    else if ((*ni)!=root)
      (*ni)->setDistanceToFather( node_2_rank[(*ni)] - node_2_rank[(*ni)->getFather()] );
      
  scalar_type max_dLL=-1e-30;

  scalar_type oLL=LL_mpi(infer_tree,can_order(S),mode);
  // XXOO
  if ( mpi_rank==server  && supp!=0)
    {
      string foutname=infer_tree->outstream_file_name+".fams";
      ofstream cout_stream(foutname.c_str(),ios::app);  
      cout << "PRINTING to" << foutname << endl;
      cout_stream << "#BEGIN " << oLL << " "<<can_order(S) << endl;
      cout_stream << "#WTF " << oLL << " "<<TreeTemplateTools::treeToParenthesis(*S) << endl;

      cout_stream << "#BTOMOVE" << endl;
      cout_stream << "#tree" << " " <<0<< " "<< can_order(S) << endl;      
      for ( int i=0;i<infer_tree->glls.size();i++)
	cout_stream << "#FAM " << infer_tree->glls[i] << " " << infer_tree->gcodes[i] <<endl;
      cout_stream << "#ETOMOVE" << endl;
      cout_stream.close();     
    }
  
  //cout << " oll " << oLL <<endl;
  //cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID");

  // find moves
  for (int ri=2;ri<S->getNumberOfLeaves();ri++)
    {
      Node * node_i=rank_2_node[ri];
      int rci,rfi;
      if (node_i->isLeaf() || (node_i->getSon(0)->isLeaf()&&node_i->getSon(1)->isLeaf()))
	rci=S->getNumberOfLeaves();
      else
	{
	  rci=min(node_2_rank[node_i->getSon(0)],node_2_rank[node_i->getSon(1)]);
	  if (node_2_rank[node_i->getSon(0)]==0)
	    rci=node_2_rank[node_i->getSon(1)];
	  if (node_2_rank[node_i->getSon(1)]==0)
	    rci=node_2_rank[node_i->getSon(0)];
	}
      rfi=node_2_rank[node_i->getFather()];
      //	cout << "> " << ri << " " << rci << " " << rfi << endl; 
      for (int rj=rfi+1;rj<rci;rj++)
	{
	  if (ri!=rj)
	    {	  
	      pair <int,int> move;
	      move.first=ri;
	      move.second=rj;
	      moves.push_back(move);	    
	    }
	  //cout << ri << " " << rj << endl;
	}
    }
  nodes = S->getNodes();    
  root= S->getRootNode();
  
  pair<int,int> max_move;
  for (vector <pair <int,int> > :: iterator mi = moves.begin(); mi != moves.end(); mi++ )
    {
      //cout << (*mi).first << " <-> " << (*mi).second << endl;
      id_2_rank[rank_2_id[(*mi).first]] = (*mi).second;
      //move up
      if ((*mi).second<(*mi).first)
	for (int ri=(*mi).second;ri<(*mi).first;ri++)
	  {
	    id_2_rank[rank_2_id[ri]] += 1;	  
	  }
      //move down            
      else
	for (int ri=(*mi).second;ri>(*mi).first;ri--)
	  {
	    id_2_rank[rank_2_id[ri]] -= 1;	  
	  }            


      for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
	if  ((*ni)->isLeaf() )
	  (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
	else if ((*ni)!=root)
	  (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
      scalar_type dLL =LL_mpi(infer_tree,can_order(S),mode)-oLL;
      if (mpi_rank==server && supp!=0)
	{
	  string foutname=infer_tree->outstream_file_name+".fams";
	  ofstream cout_stream(foutname.c_str(),ios::app);  
	  cout_stream << "#BTOMOVE" << endl;
	  cout_stream << "#tree" << " " <<dLL<< " "<< can_order(S) << endl;      
	  cout_stream << "#WTF " << dLL << " "<<TreeTemplateTools::treeToParenthesis(*S) << endl;
	  cout_stream << "#WTF " << (*mi).first << " <-> " << (*mi).second  << endl;


	  for ( int i=0;i<infer_tree->glls.size();i++)
	    cout_stream << "#FAM " << infer_tree->glls[i] << " " << infer_tree->gcodes[i] <<endl;
	  cout_stream << "#ETOMOVE" << endl;
	  cout_stream.close();     
	}
      
      //cout << dLL << " "  << (*mi).first << " <-> " << (*mi).second << endl;



      id_2_rank[rank_2_id[(*mi).first]] = (*mi).first;
      //move up
      if ((*mi).second<(*mi).first)
	for (int ri=(*mi).second;ri<(*mi).first;ri++)
	  {
	    id_2_rank[rank_2_id[ri]] -= 1;	  
	  }
      //move down            
      else
	for (int ri=(*mi).second;ri>(*mi).first;ri--)
	  {
	    id_2_rank[rank_2_id[ri]] += 1;	  
	  }            

      for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
	if  ((*ni)->isLeaf() )
	  (*ni)->setDistanceToFather( rank+1 - id_2_rank[(*ni)->getFather()->getId()] );
	else if ((*ni)!=root)
	  (*ni)->setDistanceToFather( id_2_rank[(*ni)->getId()] - id_2_rank[(*ni)->getFather()->getId()] );
      if (dLL>max_dLL)
	{
	  max_dLL=dLL;
	  //max_dS=dS;
	  max_move.first=(*mi).first;
	  max_move.second=(*mi).second;
	  if (dLL>0 && greedy) break;
	}
      
    }

  if (max_dLL>0)
    {
      id_2_rank[rank_2_id[max_move.first]] = max_move.second;
      //move up
      if (max_move.second<max_move.first)
	for (int ri=max_move.second;ri<max_move.first;ri++)
	  {
	    id_2_rank[rank_2_id[ri]] += 1;	  
	  }
      //move down            
      else
	for (int ri=max_move.second;ri>max_move.first;ri--)
	  {
	    id_2_rank[rank_2_id[ri]] -= 1;	  
	  }            
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

  times.clear();
  nodes.clear();
  node_2_rank.clear();
  rank_2_node.clear();
  id_2_rank.clear();
  rank_2_id.clear();
  //for (vector <pair <int,int> > :: iterator it=moves.begin();it!=moves.end();it++)
  //  (*it).clear();
  moves.clear();  
  if (mpi_rank==server)
    {
      string foutname=infer_tree->outstream_file_name+".fams";
      ofstream cout_stream(foutname.c_str(),ios::app);  
      cout << "#END " << can_order(S) << endl;
    }
  return ML_step;
  
  

}





TreeTemplate<Node> * random_step(TreeTemplate<Node> * S,string mode)
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
  

  if (mode=="root")
    {
      vector < Node *> sons=root->getSons();
      for (vector < Node *> :: iterator ni=sons.begin();ni!=sons.end(); ni++)
	if (!((*ni)->isLeaf()) )
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
	      if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>1)
		(*ii).second+=1;
	    id_2_rank[(*ni)->getId()]=2;
	

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
	      if ((*ii).second<old_rank+1 && (*ii).second>1)
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
	      if ((*ii).second<id_2_rank[(*ni)->getId()] && (*ii).second>1)
		(*ii).second+=1;
	    id_2_rank[(*ni)->getId()]=2;
	
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
	      if ((*ii).second<old_rank+1 && (*ii).second>1)
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
    }


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


string randomize_time_order(string Sstring,int steps,string mode)
{
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);
  for (int i=0;i<steps;i++)
    S=random_step(S,mode);  
  return can_order(S);
}
string can_order(TreeTemplate<Node> * S)
{
  stringstream out;

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


  return TreeTemplateTools::treeToParenthesis(*S,false,"ID");

}
string print_time_order(string Sstring)
{
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);
  stringstream out;

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
  //cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID");

  for (map <int,int> :: iterator nri = id_2_rank.begin(); nri != id_2_rank.end(); nri++ )
    {
      //cout << (*nri).second << " ";
      out << (*nri).second << " ";
    }
  //cout << endl;

  return out.str();

}

scalar_type GS_alpha(Species_tree * infer_tree,TreeTemplate<Node> * S,scalar_type rho,scalar_type height, scalar_type beta,string mode, scalar_type stem_l,scalar_type stem_o)
{

  //Golden Section search
  scalar_type GR=0.3819;
  scalar_type prec=0.1;
  scalar_type left=0.01;
  scalar_type middle=GR;
  scalar_type right=1.;
  scalar_type diff=2*prec;
  scalar_type alpha;
  alpha=middle;
  scalar_type LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
  scalar_type LL_x;
  while(prec<diff)
    {
      alpha= (right-middle)*GR+middle;       	  
      LL_x=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
      diff=abs(LL_x-LL_middle);
      if(LL_x>LL_middle)
	{
	  left=middle;
	  middle=alpha;
	  LL_middle=LL_x;
	}
      else
	{
	  right=alpha;
	  middle=(middle-left)*GR+left;
	  alpha=middle;
	  LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
	}
    }
  return alpha;
  
}
scalar_type GS_omega(Species_tree * infer_tree,TreeTemplate<Node> * S,scalar_type alpha,scalar_type rho, scalar_type height, scalar_type beta,string mode,scalar_type stem_o)
{

  //Golden Section search
  scalar_type GR=0.3819;
  scalar_type prec=0.1;
  scalar_type left=1e-10;
  scalar_type middle=GR;
  scalar_type right=1.;
  scalar_type diff=2*prec;
  
  scalar_type stem_l=middle;
  scalar_type LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
  scalar_type LL_x;
  while(prec<diff)
    {
      stem_l= (right-middle)*GR+middle;       	  
      LL_x=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
      diff=abs(LL_x-LL_middle);
      if(LL_x>LL_middle)
	{
	  left=middle;
	  middle=stem_l;
	  LL_middle=LL_x;
	}
      else
	{
	  right=stem_l;
	  middle=(middle-left)*GR+left;
	  stem_l=middle;
	  LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
	}
    }
  return stem_l;
  
}


scalar_type GS_rho(Species_tree * infer_tree,TreeTemplate<Node> * S,scalar_type alpha,scalar_type height, scalar_type beta,string mode, scalar_type stem_l,scalar_type stem_o)
{

  //Golden Section search
  scalar_type GR=0.3819;
  scalar_type prec=0.1;
  scalar_type left=0.01;
  scalar_type middle=GR;
  scalar_type right=1.;
  scalar_type diff=2*prec;
  scalar_type rho;
  rho=middle;
  scalar_type LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
  scalar_type LL_x;
  while(prec<diff)
    {
      rho= (right-middle)*GR+middle;       	  
      LL_x=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
      diff=abs(LL_x-LL_middle);
      if(LL_x>LL_middle)
	{
	  left=middle;
	  middle=rho;
	  LL_middle=LL_x;
	}
      else
	{
	  right=rho;
	  middle=(middle-left)*GR+left;
	  rho=middle;
	  LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
	}
      //cout << rho << endl;
    }
  return rho;
  
}

scalar_type GS_height(Species_tree * infer_tree,TreeTemplate<Node> * S,scalar_type alpha,scalar_type rho,scalar_type height, scalar_type beta,string mode, scalar_type stem_l,scalar_type stem_o)
{

  //Golden Section search
  scalar_type GR=0.3819;
  scalar_type prec=0.1;
  scalar_type left=0.1;
  scalar_type right=2.*height;
  scalar_type middle=GR*right;
  scalar_type diff=2*prec;
  height=middle;
  scalar_type LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
  scalar_type LL_x;
  while(prec<diff)
    {
      height= (right-middle)*GR+middle;       	  
      LL_x=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
      diff=abs(LL_x-LL_middle);
      if(LL_x>LL_middle)
	{
	  left=middle;
	  middle=height;
	  LL_middle=LL_x;
	}
      else
	{
	  right=height;
	  middle=(middle-left)*GR+left;
	  height=middle;
	  LL_middle=LL(infer_tree,S,alpha*rho*height,(1-alpha)*rho*height,(1-rho)*height,beta,mode,stem_l,stem_o);
	}
    }
  return height;
  
}



void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

string canonic_time_order(tree_type * T)
{
  // get inital set up from branch lengts
  Node * root = T->getRootNode();  
  vector<Node * > nodes = T->getNodes();    
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



  map <string,int> order_map;
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if (!(*ni)->isLeaf())
      {
	string name="";
  	vector <string> leaves= TreeTemplateTools::getLeavesNames(*(*ni));
	sort(leaves.begin(),leaves.end());
	for (vector <string>::iterator it=leaves.begin();it!=leaves.end();it++)
	  name+=(*it);
	order_map[name]= node_2_rank[(*ni)];
      }
  stringstream out;

  for (map < string,int> :: iterator hi = order_map.begin(); hi != order_map.end(); hi++ )
    if ((*hi).second>9)
      out<<" "<<(*hi).second;
    else
      out<<"  "<<(*hi).second;
  //out<<"\n";
  return out.str();
}


string DTL_time_order(Species_tree * infer_tree, string Sstring)
{
  const mpi::communicator  world = infer_tree->world;   
  int server = 0;  
  int mpi_rank = world.rank();
  int size = world.size();

  //  infer_tree->strip_virtual_leaves();
  //tree_type * S = infer_tree->tree;
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);

  Node * root=S->getRootNode();
      
  vector <Node*> nodes=S->getNodes();
  for (vector<Node * > :: iterator it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      if ((*it)->getDistanceToFather()==0)
	(*it)->setDistanceToFather(0.001);
  
  infer_tree->node_makeclocklike_height((S->getRootNode()));    
  S->scaleTree(1./TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*root,*(S->getLeaves()[0])));      
  double stem_len=0.1;
  if (true)
    {
      Node * new_root = new Node();  
      Node * out_group = new Node();  
      root->addSon(new_root);
      new_root->addSon(out_group);      
      out_group->setDistanceToFather(1.+stem_len);
      out_group -> setName("OUT_GROUP");
      S->rootAt(new_root);
      root=S->getRootNode();
      if (root->getSon(0)->isLeaf())
	root->getSon(1)->setDistanceToFather(stem_len);
      else
	root->getSon(0)->setDistanceToFather(stem_len);
    }
  //root->getSon(1)->setDistanceToFather(stem_len);
  
  infer_tree->name_internal(root);
  

  nodes = S->getNodes();    
  map <scalar_type,Node * > times;
  scalar_type root_h=infer_tree->DTLclock_height(1);
  for (vector<Node * > :: iterator ni = nodes.begin(); ni != nodes.end(); ni++ )
    if (!(*ni)->isLeaf() && (*ni)->getSons().size()>1)
      {
	scalar_type h=root_h-infer_tree->DTLclock_height(infer_tree->branchi[(*ni)->getName()]);
	//if (mpi_rank==server) cout <<" H "<< h <<endl;
	while (times.count(h))
	  h+=1e-4;      
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
    {
      if  ((*ni)->isLeaf() )
	(*ni)->setDistanceToFather( max(rank+1 - node_2_rank[(*ni)->getFather()],1) );
      else if ((*ni)!=root)
	(*ni)->setDistanceToFather( max(node_2_rank[(*ni)] - node_2_rank[(*ni)->getFather()],1) );
      //if (mpi_rank==server && (*ni)->hasDistanceToFather()) cout << " D " << (*ni)->getDistanceToFather() <<endl;	
    }

  node_2_rank.clear();
  rank_2_node.clear();
  id_2_rank.clear();
  rank_2_id.clear();
  nodes.clear();
  times.clear();
  moves.clear();
  string tmp=can_order(S);
  delete S;
  tmp=strip_out(tmp);
  if (mpi_rank==server) cout << " R " << tmp <<endl;
  return tmp; 
}


scalar_type   rate_estimate(tree_type * S,vector <string> trees,Species_tree * infer_tree, bool homogenous)
{
  scalar_type obs_lambda=0.1;
  scalar_type obs_delta=0.1;  
  pair<scalar_type,scalar_type> obs=observe_BD(infer_tree, S, trees, obs_delta/2.,obs_delta/2., obs_lambda);
  for (int i=0;i<3;i++)
    {
      obs_lambda=obs.first;
      obs_delta=obs.second;
      obs=observe_BD(infer_tree, S, trees, obs_delta*0.5,obs_delta*0.5, obs_lambda);
    }
  scalar_type rho=0.1;
  scalar_type delta=round(obs_delta*(rho)*10000.)/10000.;
  scalar_type tau=round(obs_delta*(1-rho)*10000.)/10000.;
  scalar_type lambda=round(obs_lambda*10000.)/10000.;
  
  obs_delta=max(rho,obs_delta);
  cout << obs_delta << " " << obs_lambda <<endl;
  string mode="intbck";
  infer_tree->w=0.;
  infer_tree->intp=0.;
  infer_tree->branchwise=0;  
  scalar_type ll= oLL(infer_tree,S,delta,tau,lambda,1,mode);

  infer_tree->count_events(infer_tree->recs,true,true);
  infer_tree->recs.clear();
  infer_tree->peas.clear();
  if (!homogenous)
    {
      ll=oLL(infer_tree,S,delta,tau,lambda,1,mode) ;
      infer_tree->count_events(infer_tree->recs,false,true);
      infer_tree->recs.clear(); 
      infer_tree->peas.clear();
      ll=oLL(infer_tree,S,infer_tree->delta_avg,infer_tree->tau_avg,infer_tree->lambda_avg,1,"noestimate") ;
    }
  else
    {
      infer_tree->branchwise=0;  
      ll=oLL(infer_tree,S,infer_tree->delta_avg,infer_tree->tau_avg,infer_tree->lambda_avg,1,"noestimate") ;
    }
  return ll;
}



string DTL_clock(string Sstring,Species_tree * infer_tree,string mode)
{
  const mpi::communicator  world = infer_tree->world;   
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  broadcast(world, Sstring,server);
  scalar_type ll,lll;
  if (rank==server) cout << "." << 0 <<endl;
  
  ll=LL_mpi(infer_tree,Sstring,mode);
  lll=ll;
  Sstring =DTL_time_order(infer_tree,Sstring);
  
  for (int i=0;i<10;i++)
    {
      if (rank==server) cout << "." << i+1 <<endl;
      ll=LL_mpi(infer_tree,Sstring,mode);
      if (abs(ll-lll)<1 || ll<lll)
	break;
      Sstring=DTL_time_order(infer_tree,Sstring);

      lll=ll;
    }
  if (rank==server) cout << "." <<endl;
  infer_tree->tmp_ll=lll;
  
  return Sstring;
}



pair < tree_type *,vector<string> > restrict(tree_type * S,vector<string> sample_species,vector<string> sample_trees)
{
  tree_type * tree = S->clone();
  vector <Node* > leaves=tree->getLeaves();
  for (vector <Node *>::iterator si=leaves.begin(); si!=leaves.end();si++)
    //not in sample
    if (find(sample_species.begin(), sample_species.end(), (*si)->getName())==sample_species.end())
      {
	Node * unsampled = (*si);	
	Node * father = (*si)->getFather();
	Node * root = tree->getRootNode();
	father->removeSon(unsampled);
	if (father!=root)
	  {
	    Node * fathersfather = father->getFather();
	    Node * son = father->getSon(0);
	    scalar_type s2f = son->getDistanceToFather();
	    scalar_type f2fsf = father->getDistanceToFather();
	    fathersfather->removeSon(father);
	    father->removeSon(son);
	    fathersfather->addSon(son);
	    delete father;
	    son->setDistanceToFather(s2f+f2fsf);
	  }
	else
	  {
	    father->removeSon(unsampled);	    
	    Node * old_root = root;
	    root=root->getSon(0);
	    tree->rootAt(root);
	    root->removeSon(old_root);
	  }
	tree->resetNodesId();
      }
  //cout << TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
  tree_type * S_unsample = tree->clone();

  vector <string> unsample_trees;       
  for (vector <string>::iterator ti=sample_trees.begin();ti!=sample_trees.end();ti++)
    {
      //restrict tree
      tree_type * tree = TreeTemplateTools::parenthesisToTree((*ti));
      vector <Node* > leaves=tree->getLeaves();
      Node * root = tree->getRootNode();
      int sampled_genes=0;
      for (vector <Node *>::iterator si=leaves.begin(); si!=leaves.end();si++)
	{
	  vector<string> tokens;
	  Tokenize((*si)->getName(),tokens,"&");
	  string name=tokens[0];	      
	  
	  //unsampled
	  if (find(sample_species.begin(), sample_species.end(),name)!=sample_species.end())
	    sampled_genes++;
	}
      if (sampled_genes>1)
	{
	  for (vector <Node *>::iterator si=leaves.begin(); si!=leaves.end();si++)
	    //unsampled
	    {								
	      vector<string> tokens;
	      Tokenize((*si)->getName(),tokens,"&");
	      string name=tokens[0];	      
	      
	    if (find(sample_species.begin(), sample_species.end(),name )==sample_species.end())
	      {
		Node * unsampled = (*si);	
		Node * father = unsampled->getFather();
		father->removeSon(unsampled);
		if (father!=root)
		  {
		    Node * fathersfather = father->getFather();
		    Node * son = father->getSon(0);
		    scalar_type s2f = son->getDistanceToFather();
		    scalar_type f2fsf = father->getDistanceToFather();
		    fathersfather->removeSon(father);
		    father->removeSon(son);
		    fathersfather->addSon(son);
		    delete father;
		    son->setDistanceToFather(s2f+f2fsf);
		  }
		else
		  {
		    Node * old_root = root;
		    root=root->getSon(0);
		    tree->rootAt(root);
		    root->removeSon(old_root);
		  }
		tree->resetNodesId();
	    

	      }
	}
	  unsample_trees.push_back(TreeTemplateTools::treeToParenthesis(*tree,false,"ID"));    
	  //cout << TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
	  //cout << tree->getNumberOfLeaves() << endl;
	}  
      delete tree;
	  
    }    
  pair  < tree_type *,vector<string> > return_val;
  return_val.first=S_unsample;
  return_val.second=unsample_trees;
  return return_val;  
}

unsigned int good_seed()
{
  long random_seed,random_seed_a, random_seed_b; 
  std::ifstream file ("/dev/urandom", std::ios::binary);
  if (file.is_open())
    {
      char * memblock;
      int size = sizeof(long);
      memblock = new char [size];
      file.read (memblock, size);
      file.close();
      random_seed_a = *reinterpret_cast<int*>(memblock);
      cout << random_seed_a <<endl;
      delete[] memblock;
    }// end if
  else
    {
      random_seed_a = 1;
    }
    random_seed_b = std::time(0);
    random_seed = random_seed_a xor random_seed_b;
    return random_seed;
} // end good_seed()
