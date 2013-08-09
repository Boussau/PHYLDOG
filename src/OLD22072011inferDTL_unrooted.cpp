#define  isorisnot( slice, si, edge, N_edges) (si==0 && (edge==N_edges-1 || slice == 0)) || (si>0)


#include "DTL.h"
#include <float.h>
using namespace std;
using namespace bpp;


scalar_type Species_tree::unrooted_run(scalar_type unseen, string mode)
{

  //########################## GLOBAL #########################
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  //long_vector_type tmp_omega,tmpp_omega;
  vector<long_vector_type> gather_tmpp_omega;
  vector<scalar_type> broadcast_branch_omega;
  vector<scalar_type> gather_lls;
  
  scalar_type broadcast_ll;
  scalar_type client_ll;
  tmpp_omega.clear();tmp_omega.clear();
  tmpp_omega.resize(lin_size);tmp_omega.resize(lin_size);
  scalar_type c=0;
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    broadcast_branch_omega.push_back(0); 

  if (mode=="intbck")
    {
      tmp_pea_sum=0.;
      O_counts.clear();  D_counts.clear();  T_counts.clear(); from_count.clear();
      for (vector< vector < scalar_type > >::iterator it=Ttf.begin();it!=Ttf.end();it++)
	(*it).clear();
      Ttf.clear();
      for (int i = 0; i<N_slices+tree->getNumberOfLeaves()-1; i++)
	{
	  vector <scalar_type> tmp;
	  for (int j = 0; j<N_slices+tree->getNumberOfLeaves()-1; j++)	  
	    tmp.push_back(0.);
	  Ttf.push_back(tmp);
	}
      
      branch_sum.clear();  branch_count.clear();  branch_zero.clear();branch_genome_size.clear();
      for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
	{
	  from_count.push_back(0); O_counts.push_back(0); D_counts.push_back(0); T_counts.push_back(0);
	  branch_sum.push_back(0); branch_count.push_back(0); branch_zero.push_back(0);branch_genome_size.push_back(0);
	}  
    }
  //########################## GLOBAL #########################

  scalar_type ll,ll2;
  ll=-1e30;
  ll2=-1e30;
  

  if (rank==server)
    {
	        
      //#!# done on SEVER
      for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
	{	 
	  broadcast_branch_omega[i]=branchwise_omega[i];
	}


      //## - ## ########## bcast p_omega - SERER ################## ## - ##
      broadcast(world,broadcast_branch_omega,server);
      //## - ## ########## bcast p_omega - SERVER ################## ## - ##      

      while ( (mode=="estimate" && (ll-ll2)>0.1) || c==0)
	{	  

	  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
	    {
	      branchwise_omega[i]=broadcast_branch_omega[i];	 	  
	      string namei=branchname[i];
	      namewise_omega[namei]=branchwise_omega[i];
	    }

	  scalar_type norm=0;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    if (type_of_x[x]!=1)
	      {
		p_omega(x)=broadcast_branch_omega[x_down(x)];//XX this is initially set to uniform	    
		norm+=p_omega(x);
	      }
	    else
	      p_omega(x)=0;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    p_omega(x)/=norm;	  
	  c+=1;
	  
	  // ## clients do calculations 
 	  if (mode=="estimate" || branchwise==false)
	    {
	      //## - ## ########## gather tmpp_omega - SERVER ################## ## - ##
	      gather(world,tmpp_omega,gather_tmpp_omega,server);
	      //## - ## ########## gather tmpp_omega - SERVER ################## ## - ##
	      
	      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		tmpp_omega(x)=0;
	      scalar_type renorm=0.;
	      for (int j=1;j<size;j++)	
		for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		  if (type_of_x[x]!=1)
		    {
		      renorm+=gather_tmpp_omega[j](x);
		      tmpp_omega(x)+=gather_tmpp_omega[j](x);
		    }
	      
	      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		tmpp_omega(x)/=renorm;	     
	      
	      //#!# only done on SEVER
	      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		if (type_of_x[x]!=1)
		  {
		    int slice=lin_slice[x];     
		      rec_O_profile[N_slices-slice]+=tmpp_omega(x)/edges_in_slice[slice];
		  }
	      //#!# only done on SEVER
	      
	      map <int,scalar_type> sumtmp;
	      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		if (type_of_x[x]!=1)		  
		  sumtmp[x_down(x)]+=tmpp_omega(x);
	      for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
		broadcast_branch_omega[i]=sumtmp[i];
	      //#!# done on SEVER
	      
	      //## - ## ########## bcast p_omega - SERER ################## ## - ##
	      broadcast(world,broadcast_branch_omega,server);
	      //## - ## ########## bcast p_omega - SERVER ################## ## - ##
	      sumtmp.clear();
	    }
	  //## - ## ########## gather ll for while loop - SERVER ################## ## - ##
	  client_ll=0.;
	  gather(world,client_ll,gather_lls,server);

	  broadcast_ll=0.;
	  for (int j=1;j<size;j++)	
	    {
	      broadcast_ll+=gather_lls[j];
	    }
	  broadcast(world,broadcast_ll,server);
	  //## - ## ########## gather ll for while loop - SERVER ################## ## - ##
	  ll2=ll;	  
	  ll=broadcast_ll;
	  //cout << ll << " " << ll2 << " r"<<rank<<endl;
	}
      gather_tmpp_omega.clear();
      broadcast_branch_omega.clear();
      gather_lls.clear();
      //CCC
      if (mode=="estimate")
	return ll;
      else
	return ll;
    }
  else
    {
      //counting reset
      cherries.clear();
      cherry_reg.clear();
      for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ri= run_vec.begin();ri!= run_vec.end();ri++)
	if (mode=="bck" || mode =="intbck" || mode=="tree" )
	  backstep_m_N2((*ri)[0],(*ri)[1],(*ri)[2]);
	else
	  DPstep_L_N2((*ri)[0],(*ri)[1],(*ri)[2]);
      scalar_type sumsumLL,root_x_sumsum;

      int i;
      //map <scalar_type,vector<pair<int,pair<node_type*,node_type*> > > > summary;
      //CCC
      for (map <scalar_type,vector<pair<int,pair<node_type*,node_type*> > > > :: reverse_iterator vsi= summary.rbegin();vsi!= summary.rend();vsi++)
	{
	  //for (vector<pair<int,pair<node_type*,node_type*> > >  :: reverse_iterator si= (*vsi).second.rbegin();si!= (*vsi).second.rend();si++)	     
	  //{
	  //   (*si).second.clear();
	  //   (*si).clear();
	  // }
	  (*vsi).second.clear();
	}
    
      summary.clear();

      //## - ## ########## bcast p_omega - CLIENT ################## ## - ##
      broadcast(world,broadcast_branch_omega,server);
      //## - ## ########## bcast p_omega - CLIENT ################## ## - ##
      
      //CCC
      while (  ((mode=="estimate" && (ll-ll2)>0.1) || c==0))
	{
	  i=0;
	  for (int j = 0;j<N_slices+tree->getNumberOfLeaves()-1;j++)
	    {
	      branchwise_omega[j]=broadcast_branch_omega[j];	 	  
	      string namei=branchname[j];
	      namewise_omega[namei]=branchwise_omega[j];
	    }
	  scalar_type norm=0;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    if (type_of_x[x]!=1)
	      {
		p_omega(x)=broadcast_branch_omega[x_down(x)];//XX this is initially set to uniform
		norm+=p_omega(x);
	      }
	    else
	      p_omega(x)=0;	  
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    p_omega(x)/=norm;	  

	  if (mode=="estimate")
	    rec_O_profile.clear();

	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    {tmp_omega(x)=0;tmpp_omega(x)=0;}
	  sumsumLL=0.,root_x_sumsum=0.;
	  //i=0;
	  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ts= trees_roots.begin(); ts!= trees_roots.end();ts++ )
	    {
	      summary.clear();      
	      c+=1;root_x_sumsum=0;
	      //scalar_type root_x_sumsum=0;
	      scalar_type norm=lin_size-ur_S_N_leaves,psum=0;
	      scalar_type root_norm=(*ts).size(),resum=0.;
	      for (std::vector< std::pair<node_type*,node_type*> > :: iterator ri= (*ts).begin();ri!= (*ts).end();ri++)	
		{
		  pair <Node*,Node*> root_dnode=(*ri);
		  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		    if (type_of_x[x]!=1)
		      {
			scalar_type tmp=0;
			if (cherries.count(root_dnode)==0)
			  tmp=(*ur_a[root_dnode])[x]/norm/root_norm*p_omega(x);
			else
			  tmp=(*cherry_a[cherries[root_dnode]])[x]/norm/root_norm*p_omega(x);
			psum+=tmp;
			tmp_omega(x)+=tmp;
			root_x_sumsum+=tmp;
		      }	  
		  if (mode=="intbck" || mode=="tree")
		    {	     
		      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
			if (type_of_x[x]!=1)
			  if (cherries.count(root_dnode)==0) {
			    pair<int,pair<node_type*,node_type*> > tmp;
			    tmp.first=x;
			    tmp.second=root_dnode;
			    scalar_type tmpp=(*ur_a[root_dnode])[x]/norm/root_norm*p_omega(x);
			    if (tmpp>0 )
			      {
				summary[tmpp].push_back(tmp);
				resum+=tmpp;
			      }}
			  else {
			    pair<int,pair<node_type*,node_type*> > tmp;
			    tmp.first=x;
			    tmp.second=root_dnode;
			    scalar_type tmpp=(*cherry_a[cherries[root_dnode]])[x]/norm/root_norm*p_omega(x);
			    if (tmpp>0 )
			      {
				summary[tmpp].push_back(tmp);
				resum+=tmpp;
			      }
			  }}}
	      scalar_type rootx_resum=resum;
	      if (mode=="intbck" || mode=="tree")
		{
		  int tmpc=0;
		  scalar_type resum_norm=0;
		  for (map <scalar_type,vector<pair<int,pair<node_type*,node_type*> > > > :: reverse_iterator vsi= summary.rbegin();vsi!= summary.rend();vsi++)
		    if (resum_norm<rootx_resum*(1.-intp) && tmpc<1e6 || tmpc==0)
		      {
			tmpc++;
			for (vector<pair<int,pair<node_type*,node_type*> > >  :: reverse_iterator si= (*vsi).second.rbegin();si!= (*vsi).second.rend();si++)	      
			  resum_norm+=(*vsi).first;

		      }
		    else
		      break;
		  tmpc=0;
		  resum=0;
		  for (map <scalar_type,vector<pair<int,pair<node_type*,node_type*> > > > :: reverse_iterator vsi= summary.rbegin();vsi!= summary.rend();vsi++)
		    if (resum<rootx_resum*(1.-intp) && tmpc<1e6 || tmpc==0)
		      {
			tmpc++;
			for (vector<pair<int,pair<node_type*,node_type*> > >  :: reverse_iterator si= (*vsi).second.rbegin();si!= (*vsi).second.rend();si++)	      
			  {
			    resum+=(*vsi).first;
			    tmp_pea=(*vsi).first/resum_norm;
			    tmp_pea_sum+=tmp_pea;
			    
			    tracecount((*si).second ,gin_trees[i], (*si).first);		 
			    
			    //stringstream out;
			      //event_stream=&out;
			      //tracestream((*si).second , gin_trees[i] , (*si).first);		  
			
			    //if (mode=="tree"  ||1)
			    // {
			    //	event_stream=&cout;	
			    //	backtrace((*si).second ,gin_trees[i], (*si).first);	   
			
			    // }
			      //recs.push_back(out.str().erase(out.str().size()-2));
			      //peas.push_back(1);//XX(*vsi).first/resum_norm);
			  }
		      }
		    else
		      break;
		}
	      root_x_sumsum=log(root_x_sumsum);
	      if (root_x_sumsum>0 || root_x_sumsum<-1e20 || root_x_sumsum!=root_x_sumsum)
		{
		  //cout << root_x_sumsum << endl;
		  //cout << " $$$$$$$$   "<<in_trees[i] <<endl;
		  root_x_sumsum=-1e5;
		}
	      sumsumLL+=root_x_sumsum;
	      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
		{
		  tmpp_omega(x)+=tmp_omega(x)/psum;//global normalization
		  tmp_omega(x)=0.;
		}
	      i++;
	    }

	  if (mode=="estimate" || branchwise==false)
	    {

	      //## - ## ########## gather tmpp_omega - CLIENT ################## ## - ##
	      gather(world,tmpp_omega,gather_tmpp_omega,server);
	      //## - ## ########## gather tmpp_omega - CLIENT ################## ## - ##
	      //server sums up omeg_p-s
	      //## - ## ########## gather tmpp_omega - CLIENT ################## ## - ##
	      broadcast(world,broadcast_branch_omega,server);
	      //## - ## ########## gather tmpp_omega - CLIENT ################## ## - ##

	    }
	  client_ll=sumsumLL;
	  //## - ## ########## gather ll for while loop - CLIENT ################## ## - ##
	  gather(world,client_ll,gather_lls,server);
	  // ## server sums lls

	  broadcast(world,broadcast_ll,server);
	  //## - ## ########## gather ll for while loop - CLIENT ################## ## - ##
	  ll2=ll;
	  ll=broadcast_ll;
	  //cout << ll << " " << ll2 << " cr "<<rank<<endl;

	}

      
      //
      gather_tmpp_omega.clear();
      broadcast_branch_omega.clear();
      gather_lls.clear();
      tmp_omega.clear();
      tmpp_omega.clear();
      ///CCC
      return ll; 
    }

}


scalar_type Species_tree::unrooted_register(vector < string> tree_strings, scalar_type treeN)
{
  
  cherries.clear();
  for (std::map < std::pair<node_type *,node_type *>,long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)  
    delete (*it).second;
  ur_a.clear();
  tree_samples=tree_strings.size();
  int count=0; 
  bool skip=0;
  for (vector < string> :: iterator tss=tree_strings.begin();tss!=tree_strings.end(); tss++)
    {						
      if (treeN<count)
	break;

      tree_type * G_tree= TreeTemplateTools::parenthesisToTree((*tss));
      vector <Node *> tmp = G_tree->getNodes();
	    
      int twos=0;
      for (vector<Node *>::iterator it = tmp.begin(); it!=tmp.end(); it++ )	
	if ((*it)->getNumberOfSons()>3)
	  skip=1;
	else if ((*it)->getNumberOfSons()>2)
	  twos+=1;
      if (twos>1)
	skip=1;
      if (skip)
	cout << "*"<<endl;

      if (!skip)
	{
	  count++;      
	  vector <Node *> leaves = G_tree->getLeaves();
	  for (vector<Node *>::iterator it = leaves.begin(); it!=leaves.end(); it++ )
	    {
	      vector <string> tokens;
	      string name = (*it)->getName();
	      Tokenize(name,tokens,"&");
	      (*it)->setName(tokens[0]);	      
	    }

	  vector <Node*>::iterator it,jt;	
	  name_internal(G_tree->getRootNode());
	  //string input_root_1 = G_tree->getRootNode()->getSon(0)->getName();
	  //string input_root_2 = G_tree->getRootNode()->getSon(1)->getName();
	  if (G_tree->isRooted())
	    G_tree->unroot(); 

	  vector < Node * > G_leaves = G_tree->getLeaves();
	  if ( G_leaves.size()>0)//2 || G_leaves.size()<40)
	    {

	      //vector < Node * > S_leaves = tree->getLeaves();
	      //ur_S_N_leaves = S_leaves.size();

   
	      map < pair <Node*,Node*>,int > dnodes;
	      map < pair <Node*,Node*>,int > roots;
	      vector <pair <Node*,Node*> > roots_vec;

	      map < Node*,int > do_01;
	      map < Node*,int > do_12;
	      map < Node*,int > do_02;
  
	      vector< pair <Node*,Node*> > seen;
	      for ( it = G_leaves.begin(); it!=G_leaves.end(); it++ )
		{
		  Node * node = (*it);
		  Node * head = node->getFather();
		  if (G_leaves.size()==2)
		    if (G_leaves[0]==(*it))
		      head=G_leaves[1];
		    else
		      head=G_leaves[0];
		  //cout << node->getName() <<"->"<< head->getName()<<endl;
	  
		  pair <Node*,Node*> tmp_pair(node,head);
		  if (node->getId()>head->getId())
		    {tmp_pair.first=head;tmp_pair.second=node;}
		  roots[tmp_pair]++;	

	  
		  if (roots[tmp_pair]==2)
		    {
		      Node* a_root=new Node();
		      a_root->setName(tmp_pair.first->getName() + " -0- " + tmp_pair.second->getName());
		      pair <Node*,Node*> root_dnode(a_root,a_root);
		      pair <Node*,Node*> left_pair(tmp_pair.first,tmp_pair.second);
		      pair <Node*,Node*> right_pair(tmp_pair.second,tmp_pair.first);
		      //cout << "oR:"<<tmp_pair.first->getName() << " -0- " << tmp_pair.second->getName() << endl; 
		      //if (mode=="N3")
		      //  DPstep_L_N3(root_dnode,left_pair,right_pair);
		      //else
		      //  DPstep_L_N2(root_dnode,left_pair,right_pair);
		      vector <pair <Node*,Node*> > pairs;
		      pairs.push_back(root_dnode);
		      pairs.push_back(left_pair);
		      pairs.push_back(right_pair);	  
		      run_vec.push_back(pairs);
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
		      //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
		      //if (mode=="N3")
		      //  DPstep_L_N3(head_pair,left_pair,right_pair);
		      //else 
		      //  DPstep_L_N2(head_pair,left_pair,right_pair);
		      vector <pair <Node*,Node*> > pairs;
		      pairs.push_back(head_pair);
		      pairs.push_back(left_pair);
		      pairs.push_back(right_pair);	  
		      run_vec.push_back(pairs);
		    }
		  if ((node_id==1 || node_id==2) && do_12[head]==2)
		    {
		      pair <Node*,Node*> head_pair(head,neigbours[0]);
		      if (node_id==1)
			right_pair.first=neigbours[2];
		      else 
			right_pair.first=neigbours[1];
		      //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
		      //if (mode=="N3")
		      //  DPstep_L_N3(head_pair,left_pair,right_pair);
		      //else 
		      //  DPstep_L_N2(head_pair,left_pair,right_pair);
		      vector <pair <Node*,Node*> > pairs;
		      pairs.push_back(head_pair);
		      pairs.push_back(left_pair);
		      pairs.push_back(right_pair);	  
		      run_vec.push_back(pairs);
		    }
		  if ((node_id==0 || node_id==2) && do_02[head]==2)
		    {
		      pair <Node*,Node*> head_pair(head,neigbours[1]);
		      if (node_id==0)
			right_pair.first=neigbours[2];
		      else 
			right_pair.first=neigbours[0];
		      //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
		      //if (mode=="N3")
		      //  DPstep_L_N3(head_pair,left_pair,right_pair);
		      //else 
		      //  DPstep_L_N2(head_pair,left_pair,right_pair);
		      vector <pair <Node*,Node*> > pairs;
		      pairs.push_back(head_pair);
		      pairs.push_back(left_pair);
		      pairs.push_back(right_pair);	  
		      run_vec.push_back(pairs);
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
			//cout << node->getName() <<"->"<< head->getName()<<endl;
	    
			pair <Node*,Node*> tmp_pair(node,head);
			if (node->getId()>head->getId())
			  {tmp_pair.first=head;tmp_pair.second=node;}
			roots[tmp_pair]++;	
			if (roots[tmp_pair]==2)
			  {
			    Node* a_root=new Node();
			    a_root->setName(tmp_pair.first->getName() + " -0- " + tmp_pair.second->getName());
			    pair <Node*,Node*> root_dnode(a_root,a_root);
			    pair <Node*,Node*> left_pair(tmp_pair.first,tmp_pair.second);
			    pair <Node*,Node*> right_pair(tmp_pair.second,tmp_pair.first);
			    //cout << "R:"<<tmp_pair.first->getName() << " -0- " << tmp_pair.second->getName() << endl; 
			    //if (mode=="N3")
			    //  DPstep_L_N3(root_dnode,left_pair,right_pair);
			    //else
			    //  DPstep_L_N2(root_dnode,left_pair,right_pair);
			    vector <pair <Node*,Node*> > pairs;
			    pairs.push_back(root_dnode);
			    pairs.push_back(left_pair);
			    pairs.push_back(right_pair);	  
			    run_vec.push_back(pairs);
			    root_dnodes[root_dnode]=tmp_pair;
			    roots_vec.push_back(root_dnode);
			  }
			//

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
			    //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
			    //if (mode=="N3")
			    //  DPstep_L_N3(head_pair,left_pair,right_pair);
			    //else 
			    //  DPstep_L_N2(head_pair,left_pair,right_pair);
			    vector <pair <Node*,Node*> > pairs;
			    pairs.push_back(head_pair);
			    pairs.push_back(left_pair);
			    pairs.push_back(right_pair);	  
			    run_vec.push_back(pairs);

			  }
			if ((node_id==1 || node_id==2) && do_12[head]==2)
			  {
			    pair <Node*,Node*> head_pair(head,neigbours[0]);
			    if (node_id==1)
			      right_pair.first=neigbours[2];
			    else 
			      right_pair.first=neigbours[1];
			    //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
			    //if (mode=="N3")
			    //  DPstep_L_N3(head_pair,left_pair,right_pair);
			    //else 
			    //  DPstep_L_N2(head_pair,left_pair,right_pair);
			    vector <pair <Node*,Node*> > pairs;
			    pairs.push_back(head_pair);
			    pairs.push_back(left_pair);
			    pairs.push_back(right_pair);	  
			    run_vec.push_back(pairs);

			  }
			if ((node_id==0 || node_id==2) && do_02[head]==2)
			  {
			    pair <Node*,Node*> head_pair(head,neigbours[1]);
			    if (node_id==0)
			      right_pair.first=neigbours[2];
			    else 
			      right_pair.first=neigbours[0];
			    //cout << head->getName() << " <= " << left_pair.first->getName() << "," << right_pair.first->getName() <<endl;
			    //if (mode=="N3")
			    //  DPstep_L_N3(head_pair,left_pair,right_pair);
			    //else 
			    //  DPstep_L_N2(head_pair,left_pair,right_pair);
			    vector <pair <Node*,Node*> > pairs;
			    pairs.push_back(head_pair);
			    pairs.push_back(left_pair);
			    pairs.push_back(right_pair);	  
			    run_vec.push_back(pairs);
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
	      trees_roots.push_back(roots_vec);

	      in_trees.push_back((*tss));
	      gin_trees.push_back(G_tree);

	    }
	}
    }
  return 1;

}





void Species_tree:: DPstep_L_N3(pair <Node*,Node*> dnode ,pair <Node*,Node*>  left_dnode,pair <Node*,Node*>  right_dnode)
{
  //cout << dnode.first->getName() <<"->"<< dnode.second->getName() <<endl;
  //cout << left_dnode.first->getName() <<"->"<< left_dnode.second->getName() <<endl;
  //cout << right_dnode.first->getName() <<"->"<< right_dnode.second->getName() <<endl;

  int leaf_1,leaf_2,xp,xpp,x,z;
  long_vector_type *ax;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x;
  //vector_type *ax,*amx,*vmx_head,*wmx_head,*vmx_tail,*wmx_tail,*fQef_xp,*fQef_xpp,*fQef_x;
  //vector <int> * sum_isx;

  //int tmp=0,max_x,min_x;
  //int below_xp,below_xpp,z;
  pair <Node*,Node*> child_1, child_2;
  long_vector_type *a_c1,*a_c2;
  bool is_leaf_1,is_leaf_2;
  scalar_type sev,sew,sfv,sfw,sum_tmp;
  //scalar_type max_sev,max_sew,max_sfv,max_sfw,tmp_1,tmp_2,max;
  //int max_sev_z,max_sew_z,max_sfv_z,max_sfw_z;

  // v and w
  child_1 = left_dnode;
  child_2 = right_dnode;
  pair <Node*,Node*> v_key(child_1.first,dnode.second);
  pair <Node*,Node*> w_key(child_2.first,dnode.second);

  if (ur_a.count(dnode)==0)
    ur_a[dnode] = new long_vector_type(lin_size);
  else
    ur_a[dnode]->resize(lin_size);
  ax = ur_a[dnode];
  
  //ur_a[dnode].resize(lin_size);
  //ax = & ur_a[dnode];

  is_leaf_1=child_1.first->isLeaf();
  if (is_leaf_1)
    {
      leaf_1 = ur_sigma[child_1.first->getName()];
    }
  is_leaf_2=child_2.first->isLeaf();
  if (is_leaf_2)
    {
      leaf_2 = ur_sigma[child_2.first->getName()];
    }
  // a(z,v) and a(z,w)
  //a_c1 = & ur_a[child_1];
  //a_c2 = & ur_a[child_2];
  a_c1 = ur_a[child_1];
  a_c2 = ur_a[child_2];


  for (x=lin_size-ur_S_N_leaves;x>=0;x--)
    {
      //cout << x << " ";
      if (type_of_x[x] == 1)
	{
	  // e and f

	  xp = below_son0[x];
	  xpp = below_son1[x];		   
	  fQef_x = & flat_Qef[x];
	  if (is_leaf_1 || is_leaf_2)
	    {
	      fQef_xp = & flat_Qef[xp];
	      fQef_xpp = & flat_Qef[xpp];
	    }

	  // a speciation
	  // ############## S #############
		  
	  // place v in e and f
	  if (is_leaf_1)
	    {
	      sev = (*fQef_xp)( leaf_1 - xp);
	      sfv = (*fQef_xpp)( leaf_1 - xpp);
	    }
	  else		    
	    {
	      sev=(*a_c1)(xp);
	      sfv=(*a_c1)(xpp);
	    }
	  // place w in e and f
	  if (is_leaf_2)
	    {
	      sew = (*fQef_xp)( leaf_2 - xp);
	      sfw = (*fQef_xpp)( leaf_2 - xpp);		      
	    }
	  else	
	    {
	      sew=(*a_c2)(xp);
	      sfw=(*a_c2)(xpp);
	    }
	  sum_tmp=sev*sfw + sfv*sew;

	  //cout << " S at " <<  dnode.first->getName() << " -> " << dnode.second->getName() <<endl; 
	  //cout << "     " <<inverse_flat_node_map[x]->getName()   << "+= " << inverse_flat_node_map[x]->getName() << " " << (*ax)[x ]<< endl;
	  //OK now we add the improvment

	  //for (z = advance_from[x]; z<=advance_until[x]; z++)
	  //{
	  //sum_tmp+=(*fQef_x)(z-x)*(*ax)(z);	      
	  //BORG
	  sum_tmp+=(*fQef_x)(xp-x)*(*ax)(xp);	      
	  //cout << "     " <<inverse_flat_node_map[x]->getName()   << "+= " << inverse_flat_node_map[xp]->getName() << " " << (*fQef_x)(xp-x) << " " << (*ax)[xp] << " s:" << sum_tmp << endl;
	  //cout << "     " <<x << "+= " <<z << " " <<  endl;

	  sum_tmp+=(*fQef_x)(xpp-x)*(*ax)(xpp);	      
	  //cout << "     " <<inverse_flat_node_map[x]->getName()   << "+= " << inverse_flat_node_map[xpp]->getName() << " " << (*fQef_x)(xpp-x) << " " << (*ax)[xpp] << " s:" << sum_tmp << endl;
	  //cout << "     " <<x << "+= " <<z << " " <<  endl;
	  //  }
	  (*ax)(x)=sum_tmp;
	  //cout << "."<<endl;
	  //###########################

	  // ############## S #############
	}
      else 
	{
	  sum_tmp = 0;
	  //xp = below_son0[x];//x;// + ( below_until[x] - below_from[x]);
	  xp=x;
	  //if (is_leaf_1 || is_leaf_2)
	  if (is_leaf_1 || is_leaf_2)
	    fQef_xp = & flat_Qef[xp];

	  for (xpp = below_from[x]; xpp<below_until[x]; xpp++)
	    //for (xpp = advance_from[x]; xpp<advance_until[x]; xpp++)
	    {
	      
	      if (xpp == xp)
		// a duplication
		// ############## D #############			
		{
		  
		  //below_xp = below_until[xp]+1;
		  // place v in e and f
		  if (is_leaf_1)
		    {
		      sev = (*fQef_xp)( leaf_1 - xp);
		    }
		  else
		    sev=(*a_c1)(xp);
		  if (is_leaf_2)
		    {
		      sew = (*fQef_xp)( leaf_2 - xp);
		    }		  
		  else		    
		    sew=(*a_c2)(xp);
		  sum_tmp +=  delta_dt[x] * sev*sew;;

		  // ############## D #############
		}
	      else 
		//a transfer
		// ############## T #############	

		{
		  fQef_xpp = & flat_Qef[xpp];
		  //below_xp = below_until[xp]+1;
		  // place v in e and f
		  if (is_leaf_1)
		    {
		      sev = (*fQef_xp)( leaf_1 - xp);
		      sfv = (*fQef_xpp)( leaf_1 - xpp);			  
		    }
		  else		    
		    {
		      sev=(*a_c1)(xp);
		      sfv=(*a_c1)(xpp);
		    }
		  // place w in e and f
		  if (is_leaf_2)
		    {
		      sew = (*fQef_xp)( leaf_2 - xp);
		      //cout << leaf_2 <<  " " << xpp << " " << x << " " << ur_S_N_leaves << " " << lin_size << endl;
		      sfw = (*fQef_xpp)( leaf_2 - xpp);		      			  
		    }		
		  else		    		      
		    {
		      sew=(*a_c2)(xp);
		      sfw=(*a_c2)(xpp);
		    }
		  sum_tmp += tau_dt[xpp]*(sev*sfw+sfv*sew);		  
		}
	      
	      // ############## T #############
	    }
	
	  //OK now we add the improvment
	  fQef_x = & flat_Qef[x];
	  //cout << " DT at " <<  dnode.first->getName() << " -> " << dnode.second->getName() <<endl; 
	  //cout << "     " <<inverse_flat_node_map[x]->getName()   << "+= " << inverse_flat_node_map[x]->getName() << " " << sum_tmp << endl;
	  //cout << "     " <<x << "= " <<endl;
	  //cout << advance_from[x] << " " <<  advance_until[x] << endl;
	  //OK now we add the improvment

	  //????NO ROOT SUM???
	  for (z = advance_from[x]; z<=advance_until[x]; z++)
	    {
	      if(type_of_x[advance_from[x]]==1 && z!=below_son0[advance_from[x]] &&z!=below_son1[advance_from[x]])
		sum_tmp+=(*fQef_x)(z-x)*(*ax)(z);	      
	      //cout << "     " <<inverse_flat_node_map[x]->getName()   << "+= " << inverse_flat_node_map[z]->getName() << " " << (*fQef_x)(z-x) << " " << (*ax)[z] << " s:" << sum_tmp << endl;
	      //cout << "     " <<x << "+= " <<z << " " <<  endl;

	    }
	  (*ax)(x) = sum_tmp;
	  
	  //cout << "."<<endl;
	  //###########################
	}
    }
  
}



void Species_tree:: DPstep_L_N2(pair <Node*,Node*> dnode ,pair <Node*,Node*>  left_dnode,pair <Node*,Node*>  right_dnode)
{
  pair <Node*,Node*> child_1, child_2;
  child_1 = left_dnode;
  child_2 = right_dnode;
  bool is_leaf_1,is_leaf_2;
  is_leaf_1=child_1.first->isLeaf();
  is_leaf_2=child_2.first->isLeaf();

  scalar_type max_D_c=0;
  scalar_type max_T_c=0;

  int leaf_1,leaf_2,xp,xpp,x,z,y;
  long_vector_type *ax;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x;
  long_vector_type *a_c1,*a_c2;
  scalar_type sev,sew,sfv,sfw,sum_tmp,nantmp,tv_sum_tmp,tw_sum_tmp,T_sum,D_sum,scalar_tmp,tmp_Qxy,tmp_ay;

  // v and w
  //pair <Node*,Node*> v_key(child_1.first,dnode.second);
  //pair <Node*,Node*> w_key(child_2.first,dnode.second);


  if (is_leaf_1)
    {
      leaf_1 = ur_sigma[child_1.first->getName()];
    }
  if (is_leaf_2)
    {
      leaf_2 = ur_sigma[child_2.first->getName()];
    }
  
  // cherry optimization .. not really worth the trouble ~10%
  if (is_leaf_1 && is_leaf_2)
    {       
      pair<int,int> cherry;
      if (leaf_1>leaf_2)
	{
	  cherry.first=leaf_1;
	  cherry.second=leaf_2;
	}
      else
	{
	  cherry.first=leaf_2;
	  cherry.second=leaf_1;
	}      

      // we already saw a cherry 
      if (cherry_reg.count(cherry))
	{
	  cherries[dnode]=cherry;
	  //cout << cherries.size() <<endl;
	  return;
	}
      // first time we see a cherry, we will store it's a vector in cherry_a 
      else
	{
	  if (cherry_a.count(cherry)==0)
	    cherry_a[cherry] = new long_vector_type(lin_size);
	  else
	    cherry_a[cherry]->resize(lin_size);	  	  
	  ax = cherry_a[cherry];

	  cherries[dnode]=cherry;
	  cherry_reg[cherry]=1;
	}
    }
  else
    {
      
      if (ur_a.count(dnode)==0)
	ur_a[dnode] = new long_vector_type(lin_size);
      else
	ur_a[dnode]->resize(lin_size);      
      ax = ur_a[dnode];

      if (cherries.count(child_1))
	{
	  pair<int,int> cherry=cherries[child_1];
	  a_c1 = cherry_a[cherry];
	}
      else
	a_c1 = ur_a[child_1];

      if (cherries.count(child_2))
	{
	  pair<int,int> cherry=cherries[child_2];
	  a_c2 = cherry_a[cherry];
	}
      else
	a_c2 = ur_a[child_2];

    }
   
  //pre sum for transfers
  tv_sum_tmp=0;
  tw_sum_tmp=0;
  //D_sum=0;
  //T_sum=0;
  for (xpp = lin_size-ur_S_N_leaves; xpp>lin_size-2*ur_S_N_leaves+1; xpp--)
    {	      
      //half a transfer
      // ############## T #############		     
      fQef_xpp = & flat_Qef[xpp];
      //y=below_son0[xpp];
      scalar_tmp=tau_dt[xpp];
      //tmp_ay=(*ax)(y);
      
      if (is_leaf_1)
	sfv = (*fQef_xpp)( leaf_1 - xpp);			  
      else		    
	sfv=(*a_c1)(xpp);
      if (is_leaf_2)
	sfw = (*fQef_xpp)( leaf_2 - xpp);		      			  
      else		    		      
	sfw=(*a_c2)(xpp);
      nantmp=scalar_tmp*sfv;
      if (nantmp==nantmp)
	tv_sum_tmp +=nantmp;
      nantmp=scalar_tmp*sfw;
      if (nantmp==nantmp)       
	tw_sum_tmp +=nantmp;
      //T_sum+=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
      //D_sum+= tau_dt[y] * tmp_ay;		  
      // ############## T #############
    }
  for (x=lin_size-ur_S_N_leaves;x>=0;x--)
    {
      (*ax)(x)=0;
      if (type_of_x[x] == 1)
	{
	  // e and f
	  xp = below_son0[x];
	  xpp = below_son1[x];		   
	  if (is_leaf_1 || is_leaf_2)
	    {
	      fQef_xp = & flat_Qef[xp];
	      fQef_xpp = & flat_Qef[xpp];
	    }
	  // a speciation
	  fQef_x = & flat_Qef[x];
	  // ############## S #############		  
	  // place v in e and f
	  if (is_leaf_1)
	    {
	      sev = (*fQef_xp)( leaf_1 - xp);
	      sfv = (*fQef_xpp)( leaf_1 - xpp);
	    }
	  else		    
	    {
	      sev=(*a_c1)(xp);
	      sfv=(*a_c1)(xpp);
	    }
	  // place w in e and f
	  if (is_leaf_2)
	    {
	      sew = (*fQef_xp)( leaf_2 - xp);
	      sfw = (*fQef_xpp)( leaf_2 - xpp);		      
	    }
	  else	
	    {
	      sew=(*a_c2)(xp);
	      sfw=(*a_c2)(xpp);

	    }
	  nantmp=Qegp[xp]*Qegp[xpp]*(sev*sfw + sfv*sew);
	  //nantmp=(sev*sfw + sfv*sew);
	  if (nantmp==nantmp)
	    sum_tmp=nantmp;
	  else
	    sum_tmp=0;
	  // ############## S #############
	  // ############## SL #############
	  nantmp=(*fQef_x)(xp-x)*(*ax)(xp);	      
	  if (nantmp==nantmp)
	    sum_tmp+=nantmp;
	  nantmp=(*fQef_x)(xpp-x)*(*ax)(xpp);	   
	  if (nantmp==nantmp)
	    sum_tmp+=nantmp;

	  // ############## SL #############
	  
	  (*ax)(x)=sum_tmp;

	  //pre sum for transfers
	  D_sum=0;
	  T_sum=0;
	  tv_sum_tmp=0;
	  tw_sum_tmp=0;
	  for (xpp = advance_from[x]; xpp<=advance_until[x]; xpp++)
	    {	      	      
	      //half a transfer
	      // ############## T #############		     
	      fQef_xpp = & flat_Qef[xpp];
	      y=below_son0[xpp];
	      scalar_tmp=tau_dt[xpp];
	      tmp_ay=(*ax)(y);

	      if (is_leaf_1)
		sfv = (*fQef_xpp)( leaf_1 - xpp);			  
	      else		    
		sfv=(*a_c1)(xpp);
	      if (is_leaf_2)
		sfw = (*fQef_xpp)( leaf_2 - xpp);		      			  
	      else		    		      
		sfw=(*a_c2)(xpp);

	      nantmp=scalar_tmp*sfv;		  
	      if (nantmp==nantmp)
		tv_sum_tmp +=nantmp;		  
		
	      nantmp=scalar_tmp*sfw;
	      if (nantmp==nantmp)
		tw_sum_tmp +=nantmp;		  
	      
	      nantmp=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	      if (nantmp==nantmp)
		T_sum+=nantmp;
	      nantmp= tau_dt[y] * tmp_ay;		  
	      if (nantmp==nantmp)
		D_sum+=nantmp;

	      // ############## T #############
	    }
	}
      else 
	{
	  sum_tmp = 0;
	  if (is_leaf_1 || is_leaf_2)
	    fQef_x = & flat_Qef[x];
	  // a duplication
	  // ############## D #############			
	  // place v in e and f
	  if (is_leaf_1)
	    sev = (*fQef_x)( leaf_1 - x);
	  else
	    sev=(*a_c1)(x);
	  if (is_leaf_2)
	    sew = (*fQef_x)( leaf_2 - x);
	  else		    
	    sew=(*a_c2)(x);

	  nantmp=  delta_dt[x] * sev*sew;	   
	  if (nantmp==nantmp)
	    sum_tmp+=nantmp;

	  // ############## D #############

	  // ############## T #############	  
	  nantmp= tv_sum_tmp*sew+tw_sum_tmp*sev-2*(tau_dt[x]*sew*sev);	 
	  if (nantmp==nantmp)
	    sum_tmp+=nantmp;

	  // ############## T #############	  

	  // ############## 0&TL #############
	  y=below_son0[x];
	  if (y>0)
	    {
	      fQef_x = & flat_Qef[x];

	      tmp_Qxy=(*fQef_x)(y-x);
	      tmp_ay=(*ax)(y);
	      scalar_tmp=tmp_Qxy*tmp_ay;

	      nantmp=scalar_tmp;	      
	      if (nantmp==nantmp)
		sum_tmp+=nantmp;

	      nantmp=Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum -  tau_dt[x] * scalar_tmp);	      
	      if (nantmp==nantmp)
		sum_tmp+=nantmp;

	      nantmp=Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum - scalar_tmp * tau_dt[y]) ; 		
	      if (nantmp==nantmp)
		sum_tmp+=nantmp;

	    }
	  // ############## 0&TL #############	  
	  (*ax)(x) = sum_tmp;

	}
    }
}


void Species_tree:: DPstep_m_N2(pair <Node*,Node*> dnode ,pair <Node*,Node*>  left_dnode,pair <Node*,Node*>  right_dnode)
{
  pair <Node*,Node*> child_1, child_2;
  child_1 = left_dnode;
  child_2 = right_dnode;

  bool is_leaf_1,is_leaf_2;
  is_leaf_1=child_1.first->isLeaf();
  is_leaf_2=child_2.first->isLeaf();

  scalar_type max_root=0;
  int max_x;

  scalar_type max_D_c=0;
  scalar_type max_T_c=0;

  int leaf_1,leaf_2,xp,xpp,x,z,y;
  long_vector_type *ax;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x;

  long_vector_type *a_c1,*a_c2;
  scalar_type sev,sew,sfv,sfw,sum_tmp,tv_sum_tmp,tw_sum_tmp,T_sum,D_sum,scalar_tmp,tmp_Qxy,tmp_ay;
  scalar_type max_sev,max_sew,max_sfv,max_sfw,tmp_1,tmp_2,max_T_sum,max_D_sum,max_tv_sum_tmp,max_tw_sum_tmp,max_sum_tmp;
  scalar_type nd2_tv_sum_tmp,nd2_tw_sum_tmp,nd2_T_sum,nd2_D_sum;
  int T_y,D_y,tv_x,tw_x;

  // v and w
  //pair <Node*,Node*> v_key(child_1.first,dnode.second);
  //pair <Node*,Node*> w_key(child_2.first,dnode.second);


  if (is_leaf_1)
    {
      leaf_1 = ur_sigma[child_1.first->getName()];
    }
  if (is_leaf_2)
    {
      leaf_2 = ur_sigma[child_2.first->getName()];
    }
  
  // cherry optimization .. not really worth the trouble ~10%
  if (is_leaf_1 && is_leaf_2)
    {       
      pair<int,int> cherry;
      if (leaf_1>leaf_2)
	{
	  cherry.first=leaf_1;
	  cherry.second=leaf_2;
	}
      else
	{
	  cherry.first=leaf_2;
	  cherry.second=leaf_1;
	}      
      // we already saw a cherry 
      if (cherry_reg.count(cherry))
	{
	  cherries[dnode]=cherry;	  
	  return;
	}
      // first time we see a cherry, we will store it's a vector in cherry_a 
      else
	{
	  if (cherry_a.count(cherry)==0)
	    {
	      //cout << "CoH"; 
	      cherry_a[cherry] = new long_vector_type(lin_size);
	    }
	  else
	    {
	      cherry_a[cherry]->resize(lin_size);
	    }
	  ax = cherry_a[cherry];
	  cherries[dnode]=cherry;
	  cherry_reg[cherry]=1;
	}
    }
  else
    {
      if (ur_a.count(dnode)==0)
	{
	  //cout << "dod"; 
	  ur_a[dnode] = new long_vector_type(lin_size);
	}
      else
	ur_a[dnode]->resize(lin_size);
      ax = ur_a[dnode];
      if (cherries.count(child_1))
	{
	  pair<int,int> cherry=cherries[child_1];
	  a_c1 = cherry_a[cherry];
	}
      else
	a_c1 = ur_a[child_1];
      if (cherries.count(child_2))
	{
	  pair<int,int> cherry=cherries[child_2];
	  a_c2 = cherry_a[cherry];
	}
      else
	a_c2 = ur_a[child_2];
    }
  
  tv_sum_tmp=0;
  tw_sum_tmp=0;
  max_tv_sum_tmp=0;
  max_tw_sum_tmp=0;

  for (xpp = lin_size-ur_S_N_leaves; xpp>lin_size-2*ur_S_N_leaves+1; xpp--){	      
    //half a transfer
    // ############## T #############		     
    fQef_xpp = & flat_Qef[xpp];
    scalar_tmp=tau_dt[xpp];
    
    if (is_leaf_1)
      sfv = (*fQef_xpp)( leaf_1 - xpp);			  
    else		    
      sfv=(*a_c1)(xpp);
    if (is_leaf_2)
      sfw = (*fQef_xpp)( leaf_2 - xpp);		      			  
    else		    		      
      sfw = (*a_c2)(xpp);
    tv_sum_tmp = scalar_tmp*sfv;
    if (max_tv_sum_tmp<tv_sum_tmp) 
      {
	nd2_tv_sum_tmp=max_tv_sum_tmp;
	tv_x=xpp;
	max_tv_sum_tmp=tv_sum_tmp;
      }
    //tv_sum_tmp +=scalar_tmp*sfv;
    tw_sum_tmp = scalar_tmp*sfw;
    if (max_tw_sum_tmp<tw_sum_tmp) 
      {
	nd2_tw_sum_tmp=max_tw_sum_tmp;
	tw_x=xpp;
	max_tw_sum_tmp=tw_sum_tmp;		  
      }
    //tw_sum_tmp +=scalar_tmp*sfw;
    // ############## T #############
  }
  scalar_type norm=lin_size-ur_S_N_leaves;
  for (x=lin_size-ur_S_N_leaves;x>=0;x--){
    (*ax)(x)=0;//xxx
    if (type_of_x[x] == 1){
      // e and f
      xp = below_son0[x];
      xpp = below_son1[x];		   
      if (is_leaf_1 || is_leaf_2){
	fQef_xp = & flat_Qef[xp];
	fQef_xpp = & flat_Qef[xpp];}
      // a speciation
      // ############## S #############		  
      // place v in e and f
      if (is_leaf_1){
	sev = (*fQef_xp)( leaf_1 - xp);
	sfv = (*fQef_xpp)( leaf_1 - xpp);}
      else{
	sev=(*a_c1)(xp);
	sfv=(*a_c1)(xpp);}
      // place w in e and f
      if (is_leaf_2){
	sew = (*fQef_xp)( leaf_2 - xp);
	sfw = (*fQef_xpp)( leaf_2 - xpp);}
      else{
	sew=(*a_c2)(xp);
	sfw=(*a_c2)(xpp);}

      max_sum_tmp=max(sev*sfw , sfv*sew);
      //sum_tmp=sev*sfw + sfv*sew;

      // ############## S #############
      fQef_x = & flat_Qef[x];
      // ############## SL #############
      sum_tmp=(*fQef_x)(xp-x)*(*ax)(xp);	      
      if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;
      //sum_tmp+=(*fQef_x)(xp-x)*(*ax)(xp);	   
      sum_tmp=(*fQef_x)(xpp-x)*(*ax)(xpp);	      
      if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;   
      //sum_tmp+=(*fQef_x)(xpp-x)*(*ax)(xpp);	   
      // ############## SL #############
      (*ax)(x)=max_sum_tmp;
      //(*ax)(x)=sum_tmp;
      //pre sum for transfers
      D_sum=0;
      T_sum=0;
      tv_sum_tmp=0;
      tw_sum_tmp=0;
      max_D_sum=0;
      max_T_sum=0;
      max_tv_sum_tmp=0;
      max_tw_sum_tmp=0;

      for (xpp = advance_from[x]; xpp<=advance_until[x]; xpp++){	      	      
	//half a transfer
	// ############## T #############		     
	fQef_xpp = & flat_Qef[xpp];
	y=below_son0[xpp];
	scalar_tmp=tau_dt[xpp];
	tmp_ay=(*ax)(y);
	  
	if (is_leaf_1)
	  sfv = (*fQef_xpp)( leaf_1 - xpp);			  
	else		    
	  sfv=(*a_c1)(xpp);
	if (is_leaf_2)
	  sfw = (*fQef_xpp)( leaf_2 - xpp);		      			  
	else		    		      
	  sfw=(*a_c2)(xpp);

	tv_sum_tmp =scalar_tmp*sfv;
	if (max_tv_sum_tmp<tv_sum_tmp) 
	  {
	    nd2_tv_sum_tmp=max_tv_sum_tmp;
	    tv_x=xpp;
	    max_tv_sum_tmp=tv_sum_tmp;
	  }
	//tv_sum_tmp +=scalar_tmp*sfv;
	tw_sum_tmp =scalar_tmp*sfw;
	if (max_tw_sum_tmp<tw_sum_tmp)
	  {
	    nd2_tw_sum_tmp=max_tw_sum_tmp;
	    tw_x=xpp;
	    max_tw_sum_tmp=tw_sum_tmp;		  
	  }
	//tw_sum_tmp +=scalar_tmp*sfw;
	T_sum=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	if (max_T_sum<T_sum) max_T_sum=T_sum;       
	//T_sum+=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	D_sum= tau_dt[y] * tmp_ay;	
	if (max_D_sum<D_sum) max_D_sum=D_sum;       
	//D_sum+= tau_dt[y] * tmp_ay;		  
	// ############## T #############
      }}
    else {
      sum_tmp = 0;
      max_sum_tmp=0;
      if (is_leaf_1 || is_leaf_2)
	fQef_x = & flat_Qef[x];
      // a duplication
      // ############## D #############			
      // place v in e and f
      if (is_leaf_1)
	sev = (*fQef_x)( leaf_1 - x);
      else
	sev=(*a_c1)(x);
      if (is_leaf_2)
	sew = (*fQef_x)( leaf_2 - x);
      else		    
	sew=(*a_c2)(x);
      sum_tmp =  delta_dt[x] * sev*sew;	   
      if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;
      //sum_tmp +=  delta_dt[x] * sev*sew;	   
      // ############## D #############
	
      // ############## T #############	  
      //sum_tmp = max_tv_sum_tmp*sew+max_tw_sum_tmp*sev-2*(tau_dt[x]*sew*sev);	 
      if (x==tv_x && x==tw_x)
	sum_tmp = max(nd2_tv_sum_tmp*sew , nd2_tw_sum_tmp*sev);
      else if (x==tv_x)
	sum_tmp = max(nd2_tv_sum_tmp*sew , max_tw_sum_tmp*sev);
      else if (x==tw_x )
	sum_tmp = max(max_tv_sum_tmp*sew , nd2_tw_sum_tmp*sev);
      else
	sum_tmp = max(max_tv_sum_tmp*sew , max_tw_sum_tmp*sev);
      //sum_tmp = max(max_tv_sum_tmp*sew , max_tw_sum_tmp*sev) - 1*(tau_dt[x]*sew*sev);	 
      if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;
      //sum_tmp += tv_sum_tmp*sew+tw_sum_tmp*sev-2*(tau_dt[x]*sew*sev);	 
      // ############## T #############	  
	
      // ############## 0&TL #############
      y=below_son0[x];
      if (y>0) {
	fQef_x = & flat_Qef[x];
	  
	tmp_Qxy=(*fQef_x)(y-x);
	tmp_ay=(*ax)(y);
	scalar_tmp=tmp_Qxy*tmp_ay;
	sum_tmp=scalar_tmp;	      
	if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;
	//sum_tmp+=scalar_tmp;	      
	sum_tmp= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (max_T_sum);	      
	if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;      
	//sum_tmp+= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum -  tau_dt[x] * scalar_tmp);	      
	sum_tmp= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * max_D_sum );     
	if (max_sum_tmp<sum_tmp) max_sum_tmp=sum_tmp;
	//sum_tmp+= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum - scalar_tmp * tau_dt[y]) ;
      }
      // ############## 0&TL #############	  
      (*ax)(x) = max_sum_tmp;
    }
  }
 
}

void Species_tree:: Qef_m_N2(Node* leaf)
{
  scalar_type max_sum_tmp,sum_tmp,scalar_tmp,tmp_ay,tmp_Qxy,T_sum,D_sum,max_T_sum,max_D_sum;
  int xp,xpp,x,y;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x,*fQef_y;
  int T_sum_x,D_sum_x;

  int leaf_x = ur_sigma[leaf->getName()];

  for (x=lin_size-ur_S_N_leaves;x>=0;x--){
    //cout << "xxx" << x << endl;
    fQef_x = & flat_Qef[x];
    max_sum_tmp=0;
    if (type_of_x[x] == 1){
      // e and f
      xp = below_son0[x];
      xpp = below_son1[x];		   
      fQef_xp = & flat_Qef[xp];
      fQef_xpp = & flat_Qef[xpp];
      // ############## SL #############
      sum_tmp=(*fQef_x)(xp-x)*(*fQef_xp)( leaf_x - xp);	      
      if (max_sum_tmp<sum_tmp) 
	{
	  //..
	  max_sum_tmp=sum_tmp;
	}
      //sum_tmp+=(*fQef_x)(xp-x)*(*ax)(xp);	   
      sum_tmp=(*fQef_x)(xpp-x)*(*fQef_xpp)( leaf_x - xpp);	      

      if (max_sum_tmp<sum_tmp) 
	{
	  //..
	  max_sum_tmp=sum_tmp;   
	}
      //sum_tmp+=(*fQef_x)(xpp-x)*(*ax)(xpp);	   
      // ############## SL #############
      //cout << leaf_x << " " << x<< " " << xp << " " << (*fQef_xp)( leaf_x - xp) << " " << (*fQef_x)(xp-x) <<endl;
      //cout << leaf_x << " " << x<< " " << xpp << " " << (*fQef_xpp)( leaf_x - xpp) << " " << (*fQef_x)(xpp-x) <<endl;

      //cout << "s " << leaf_x << " " << x << endl;
      (*fQef_x)( leaf_x - x )=max_sum_tmp;
      
      //pre sum for transfers
      D_sum=0;
      T_sum=0;
      max_D_sum=0;
      max_T_sum=0;
      for (xpp = advance_from[x]; xpp<=advance_until[x]; xpp++){	      	      
	//half a transfer
	// ############## T #############		     
	fQef_xpp = & flat_Qef[xpp];
	y=below_son0[xpp];

	scalar_tmp=tau_dt[xpp];

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);	  

	T_sum=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	if (max_T_sum<T_sum) 
	  {
	    //..
	    max_T_sum=T_sum;       
	    T_sum_x=y;
	  }
	//T_sum+=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	D_sum= tau_dt[y] * tmp_ay;	
	if (max_D_sum<D_sum) 
	  {
	    //..
	    max_D_sum=D_sum;       
	    D_sum_x=y;
	  }
	//D_sum+= tau_dt[y] * tmp_ay;		  
	// ############## T #############
      }}
    else {
      // ############## 0&TL #############	  	
      y=below_son0[x];
      if (y>0) {
	tmp_Qxy=(*fQef_x)(y-x);

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);

	//cout <<" y " << leaf_x << " " << x << " " << y << " " <<tmp_Qxy << " " <<tmp_ay <<endl;

	scalar_tmp=tmp_Qxy*tmp_ay;
	sum_tmp=scalar_tmp;	 
	if (max_sum_tmp<sum_tmp) 
	  {
	    max_sum_tmp=sum_tmp;
	  }
	//sum_tmp+=scalar_tmp;	      
	sum_tmp= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (max_T_sum);	      
	if (max_sum_tmp<sum_tmp) 
	  {
	    //..
	    max_sum_tmp=sum_tmp;      
	  }
	//sum_tmp+= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum -  tau_dt[x] * scalar_tmp);	      
	sum_tmp= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * max_D_sum );     
	if (max_sum_tmp<sum_tmp) 
	  {
	    //..
	    max_sum_tmp=sum_tmp;
	  }
	//sum_tmp+= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum - scalar_tmp * tau_dt[y]) ;
      }
      else
	{
	  max_sum_tmp=flat_Qef[x](leaf_x-x);
	}
      // ############## 0&TL #############
      //cout << "s " << leaf_x << " " << x << endl;
      
      (*fQef_x)( leaf_x - x )=max_sum_tmp;
      //cout << leaf_x << " " << x << endl;
    }
  
  }
}

void Species_tree:: Qef_bck_N2(Node* leaf)
{
  scalar_type max_sum_tmp,sum_tmp,scalar_tmp,tmp_ay,tmp_Qxy,T_sum,D_sum,max_T_sum,max_D_sum;
  int xp,xpp,x,y;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x,*fQef_y;
  int T_sum_x,D_sum_x;

  int leaf_x = ur_sigma[leaf->getName()];
  leaf_hidden_event_name[leaf_x][leaf_x]=1;
  leaf_hidden_event_x[leaf_x][leaf_x]=leaf_x;

  for (x=lin_size-ur_S_N_leaves;x>=0;x--){
    //cout << "xxx" << x << endl;
    fQef_x = & flat_Qef[x];
    max_sum_tmp=0;
    if (type_of_x[x] == 1){
      // e and f
      xp = below_son0[x];
      xpp = below_son1[x];		   
      fQef_xp = & flat_Qef[xp];
      fQef_xpp = & flat_Qef[xpp];
      // ############## SL #############
      sum_tmp=(*fQef_x)(xp-x)*(*fQef_xp)( leaf_x - xp);	      
      if (max_sum_tmp<sum_tmp) 
	{
	  //..
	  max_sum_tmp=sum_tmp;
	  leaf_hidden_event_name[leaf_x][x]=2;
	  leaf_hidden_event_x[leaf_x][x]=xp;
	}
      //sum_tmp+=(*fQef_x)(xp-x)*(*ax)(xp);	   
      sum_tmp=(*fQef_x)(xpp-x)*(*fQef_xpp)( leaf_x - xpp);	      

      if (max_sum_tmp<sum_tmp) 
	{
	  //..
	  max_sum_tmp=sum_tmp;   
	  leaf_hidden_event_name[leaf_x][x]=2;
	  leaf_hidden_event_x[leaf_x][x]=xpp;
	}
      //sum_tmp+=(*fQef_x)(xpp-x)*(*ax)(xpp);	   
      // ############## SL #############
      //cout << leaf_x << " " << x<< " " << xp << " " << (*fQef_xp)( leaf_x - xp) << " " << (*fQef_x)(xp-x) <<endl;
      //cout << leaf_x << " " << x<< " " << xpp << " " << (*fQef_xpp)( leaf_x - xpp) << " " << (*fQef_x)(xpp-x) <<endl;

      //cout << "s " << leaf_x << " " << x << endl;
      (*fQef_x)( leaf_x - x )=max_sum_tmp;
      
      //pre sum for transfers
      D_sum=0;
      T_sum=0;
      max_D_sum=0;
      max_T_sum=0;
      for (xpp = advance_from[x]; xpp<=advance_until[x]; xpp++){	      	      
	//half a transfer
	// ############## T #############		     
	fQef_xpp = & flat_Qef[xpp];
	y=below_son0[xpp];

	scalar_tmp=tau_dt[xpp];

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);	  

	T_sum=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	if (max_T_sum<T_sum) 
	  {
	    //..
	    max_T_sum=T_sum;       
	    T_sum_x=y;
	  }
	//T_sum+=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	D_sum= tau_dt[y] * tmp_ay;	
	if (max_D_sum<D_sum) 
	  {
	    //..
	    max_D_sum=D_sum;       
	    D_sum_x=y;
	  }
	//D_sum+= tau_dt[y] * tmp_ay;		  
	// ############## T #############
      }}
    else {
      // ############## 0&TL #############	  	
      y=below_son0[x];
      if (y>0) {
	tmp_Qxy=(*fQef_x)(y-x);

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);

	//cout <<" y " << leaf_x << " " << x << " " << y << " " <<tmp_Qxy << " " <<tmp_ay <<endl;

	scalar_tmp=tmp_Qxy*tmp_ay;
	sum_tmp=scalar_tmp;	 
	if (max_sum_tmp<sum_tmp) 
	  {
	    max_sum_tmp=sum_tmp;
	    leaf_hidden_event_name[leaf_x][x]=leaf_hidden_event_name[leaf_x][y];
	    leaf_hidden_event_x[leaf_x][x]=leaf_hidden_event_x[leaf_x][y];
	  }
	//sum_tmp+=scalar_tmp;	      
	sum_tmp= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (max_T_sum);	      
	if (max_sum_tmp<sum_tmp) 
	  {
	    //..
	    max_sum_tmp=sum_tmp;      
	    leaf_hidden_event_name[leaf_x][x]=5;
	    leaf_hidden_event_x[leaf_x][x]=T_sum_x;
	  }
	//sum_tmp+= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum -  tau_dt[x] * scalar_tmp);	      
	sum_tmp= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * max_D_sum );     
	if (max_sum_tmp<sum_tmp) 
	  {
	    //..
	    max_sum_tmp=sum_tmp;
	    leaf_hidden_event_name[leaf_x][x]=5;
	    leaf_hidden_event_x[leaf_x][x]=D_sum_x;
	  }
	//sum_tmp+= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum - scalar_tmp * tau_dt[y]) ;
      }
      else
	{
	  max_sum_tmp=flat_Qef[x](leaf_x-x);
	  leaf_hidden_event_name[leaf_x][x]=leaf_hidden_event_name[leaf_x][leaf_x];
	  leaf_hidden_event_x[leaf_x][x]=leaf_hidden_event_x[leaf_x][leaf_x];	  
	}
      // ############## 0&TL #############
      //cout << "s " << leaf_x << " " << x << endl;
      
      (*fQef_x)( leaf_x - x )=max_sum_tmp;
      //cout << leaf_x << " " << x << endl;
    }
  
  }
}


void Species_tree:: Qef_N2(Node* leaf)
{
  scalar_type max_sum_tmp,sum_tmp,scalar_tmp,tmp_ay,tmp_Qxy,T_sum,D_sum,max_T_sum,max_D_sum;
  int xp,xpp,x,y;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x,*fQef_y;

  int leaf_x = ur_sigma[leaf->getName()];

  for (x=lin_size-ur_S_N_leaves;x>=0;x--){
    //cout << "xxx" << x << endl;
    fQef_x = & flat_Qef[x];
    sum_tmp=0;
    if (type_of_x[x] == 1){
      // e and f
      xp = below_son0[x];
      xpp = below_son1[x];		   
      fQef_xp = & flat_Qef[xp];
      fQef_xpp = & flat_Qef[xpp];
      // ############## SL #############
      sum_tmp+=(*fQef_x)(xp-x)*(*fQef_xp)( leaf_x - xp);	      
      sum_tmp+=(*fQef_x)(xpp-x)*(*fQef_xpp)( leaf_x - xpp);	      
      // ############## SL #############
      //cout << leaf_x << " " << x<< " " << xp << " " << (*fQef_xp)( leaf_x - xp) << " " << (*fQef_x)(xp-x) <<endl;
      //cout << leaf_x << " " << x<< " " << xpp << " " << (*fQef_xpp)( leaf_x - xpp) << " " << (*fQef_x)(xpp-x) <<endl;

      //cout << "s " << leaf_x << " " << x << endl;
      (*fQef_x)( leaf_x - x )=sum_tmp;
      
      //pre sum for transfers
      D_sum=0;
      T_sum=0;

      for (xpp = advance_from[x]; xpp<=advance_until[x]; xpp++){	      	      
	//half a transfer
	// ############## T #############		     
	fQef_xpp = & flat_Qef[xpp];
	y=below_son0[xpp];

	scalar_tmp=tau_dt[xpp];

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);	  

	T_sum+=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	D_sum+= tau_dt[y] * tmp_ay;		  
	// ############## T #############
      }}
    else {
      // ############## 0&TL #############	  	
      y=below_son0[x];
      if (y>0) {
	tmp_Qxy=(*fQef_x)(y-x);

	fQef_y = & flat_Qef[y];
	tmp_ay=(*fQef_y)(leaf_x-y);

	//cout <<" y " << leaf_x << " " << x << " " << y << " " <<tmp_Qxy << " " <<tmp_ay <<endl;

	scalar_tmp=tmp_Qxy*tmp_ay;
	sum_tmp+=scalar_tmp;	      
	sum_tmp+= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum);	      
	sum_tmp+= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum );     
      }
      else
	sum_tmp+=flat_Qef[x](leaf_x-x);
      
      // ############## 0&TL #############
      //cout << "s " << leaf_x << " " << x << endl;
      
      (*fQef_x)( leaf_x - x )=sum_tmp;
      //cout << leaf_x << " " << x << endl;
    }

    
  }
}

vector<scalar_type> Species_tree::DTL_lk(tree_type* S,vector<string> trees,map<node_type*,scalar_type> expD,map<node_type*,scalar_type> expT,map<node_type*,scalar_type> expL,string output)
{
  vector<scalar_type> return_vec;
  vector<node_type*> nodes=S->getNodes();  
  node_type * root=S->getRootNode();
  for (vector<node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
    {
      string name=get_leaf_names((*nt));      
      if (!(*nt)->isLeaf())
	{
	  if ((*nt)!=root)       
	    name="*S."+name;
	  else
	    name="*R."+name;
	}
      scalar_type zero=0.;
      namewise_delta[name]=max(expD[(*nt)],zero);
      namewise_tau[name]=max(expT[(*nt)],zero);
      namewise_lambda[name]=max(expL[(*nt)],zero);
      namewise_omega[name]=1;      	
    }
  init_x();
  p_omega.resize(lin_size);
  branchwise=1;
  reconstruct(S);
  beta=1;  
  init_L_improved(0.01,0.01,0.01,0.1,0.01,"noestimate");   
  //tree_type * tmp=tree->clone();
  //cout << print_time_order(tmp);
  unrooted_register(trees);
  
  cherries.clear();
  cherry_reg.clear();
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ri= run_vec.begin();ri!= run_vec.end();ri++)
    DPstep_L_N2((*ri)[0],(*ri)[1],(*ri)[2]);
  scalar_type norm=lin_size-ur_S_N_leaves;
  scalar_type xLL,rootsum_LL=0;
  scalar_type xsum_rootLL,max_xsum_rootLL=0;
  scalar_type xmax_rootLL,max_xmax_rootLL=0;
  int max_x;
  int tree_i=0;
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ts= trees_roots.begin(); ts!= trees_roots.end();ts++ )
    {
      
      scalar_type root_norm=(*ts).size();
      max_xsum_rootLL=0;
      max_xmax_rootLL=0;
      rootsum_LL=0;
      for (std::vector< std::pair<node_type*,node_type*> > :: iterator ri= (*ts).begin();ri!= (*ts).end();ri++)	
	{
	  pair <Node*,Node*> root_dnode=(*ri);
	  xmax_rootLL=0;
	  xsum_rootLL=0;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    {	      
	      if (cherries.count(root_dnode)==0)
		xLL=(*ur_a[root_dnode])[x];
	      else
		xLL=(*cherry_a[cherries[root_dnode]])[x];	  
	      xsum_rootLL+=xLL/norm;
	      if (xLL>xmax_rootLL){ xmax_rootLL=xLL;max_x=x;}	      
	    }
	  if (xmax_rootLL>max_xmax_rootLL) max_xmax_rootLL=xmax_rootLL;	      
	  if (xsum_rootLL>max_xsum_rootLL) max_xsum_rootLL=xsum_rootLL;	      
	  rootsum_LL+=xsum_rootLL/root_norm;	  
	}
      if (output=="sum_orig")
	return_vec.push_back(log(max_xsum_rootLL));
      else 
	return_vec.push_back(log(max_xmax_rootLL/norm));	
      if (max_xsum_rootLL==0)
	cout << in_trees[tree_i] << " " << max_xsum_rootLL <<endl;
      //cout << tree_i << " " << max_xmax_rootLL/norm/root_norm << " " << max_xsum_rootLL/root_norm << " " << rootsum_LL <<" "<<root_norm <<" "<<norm <<" "<< 	branchname[x_down(max_x)] <<endl;
      tree_i+=1;
    }
  cherries.clear();
  for (std::map < std::pair<node_type *,node_type *>,long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)  
    delete (*it).second;
  ur_a.clear();
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ts= trees_roots.begin(); ts!= trees_roots.end();ts++ )
    (*ts).clear();
  trees_roots.clear();
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ri= run_vec.begin();ri!= run_vec.end();ri++)
    (*ri).clear();
  run_vec.clear();

  return return_vec;

}








void Species_tree::init_treewise(tree_type * S, scalar_type delta_in, scalar_type tau_in, scalar_type lambda_in,string calc_mode)
{
  if (delta_in<1e-10)
    delta_in=1e-10;
  if (tau_in<1e-10)
    tau_in=1e-10;
  if (lambda_in<1e-10)
    lambda_in=1e-10;

  mode=calc_mode;
  beta=1;  
  reconstruct(S);
  init_x();
  reset();
  branchwise=0;
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      branchwise_delta[i]=delta_in;
      branchwise_tau[i]=tau_in;
      branchwise_lambda[i]=lambda_in;
      string name=branchname[i];
      scalar_type t=branch_length[i];
      namewise_delta[name]=branchwise_delta[i]*t;
      namewise_tau[name]=branchwise_tau[i]*t;
      namewise_lambda[name]=branchwise_lambda[i]*t;
      namewise_omega[name]=branchwise_omega[i];
    }

  mode=calc_mode;
  init_L_improved(delta_in,tau_in,lambda_in,0.1,0.01,mode);   
}

vector < pair<scalar_type,string> > Species_tree::LL_treewise(vector<string>  trees,scalar_type root_p)
{
  vector < pair<scalar_type,string> > return_v;
  for (int i=0;i<trees.size();i++)
    {
      stringstream out;out<<i;
      codes[trees[i]]=out.str();  
    }  
  unrooted_register(trees);

  //computation

  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ri= run_vec.begin();ri!= run_vec.end();ri++)
    if (mode=="tree")
      backstep_m_N2((*ri)[0],(*ri)[1],(*ri)[2]);
    else
      DPstep_L_N2((*ri)[0],(*ri)[1],(*ri)[2]);
  //sum LL

  scalar_type norm=0;
  map <int,scalar_type> p_omega;
  
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
    if (type_of_x[x]!=1)
      if (x_down(x)==1)
	norm+=root_p;
      else
	norm+=1;
  
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
    if (type_of_x[x]!=1)
      if (x_down(x)==1)
	p_omega[x]=root_p/norm;
      else
	p_omega[x]=1/norm;

  int tree_i=0; //index of the gene tree we're considering.  
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ts= trees_roots.begin(); ts!= trees_roots.end();ts++ )
    {
    
    scalar_type psum=0;
    scalar_type root_norm=(*ts).size();
    pair <Node*,Node*> max_root_dnode;
    pair<Node*,Node*> max_rootsum_dnode;
    
    int max_root_x=-1;
    scalar_type max_root_L=-1;
    scalar_type max_root_sum_L=-1;
    
    for (std::vector< std::pair<node_type*,node_type*> > :: iterator ri= (*ts).begin();ri!= (*ts).end();ri++)	
      {
      pair <Node*,Node*> root_dnode=(*ri);
      scalar_type root_sum_L=0;
      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
        if (type_of_x[x]!=1)
          {
          scalar_type tmp=0;
          if (cherries.count(root_dnode)==0)
            tmp=(*ur_a[root_dnode])[x]/norm/root_norm*p_omega[x];
          else
            tmp=(*cherry_a[cherries[root_dnode]])[x]/norm/root_norm*p_omega[x];
          if (max_root_L<tmp)
            {
            max_root_L=tmp;
            max_root_x=x;
            max_root_dnode.first=root_dnode.first;
            max_root_dnode.second=root_dnode.second;
            }
          psum+=tmp;
          root_sum_L+=tmp;
          }	  	  
      if (max_root_sum_L<root_sum_L)
        {
	      max_root_sum_L=root_sum_L;
	      max_rootsum_dnode=root_dnodes[root_dnode];
        }
      }
    pair<scalar_type,string> return_pair;
    if (mode=="tree") // DTL likelihood using ML reconciliations and ML root
      {	  
        backtrace(max_root_dnode,gin_trees[tree_i], max_root_x);	  
        return_pair.first=log(max_root_L);
        return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    else if (mode=="root") //DTL likelihood summed over reconciliations, using ML root
      {
      return_pair.first=log(max_root_sum_L);
      
      
      Node* head=max_rootsum_dnode.first;
      Node* node=max_rootsum_dnode.second;
      vector <Node*> head_neigbours = head->getSons();
      if (head->hasFather())
        head_neigbours.push_back(head->getFather());
      vector <Node*> node_neigbours = node->getSons();
      if (node->hasFather())
        node_neigbours.push_back(node->getFather());
      int node_id=-1;
      
      if (node==head_neigbours[0])
        node_id=0;
      else if (node==head_neigbours[1])
        node_id=1;
      else if (node==head_neigbours[2])
        node_id=2;
      Node * new_root = new Node();
      new_root->setName("R");
      if (( node_id==2 && head->hasFather() ) || head->isLeaf())
        {
	      if (!head->hasDistanceToFather())
          head->setDistanceToFather(0.1);
	      double d=head->getDistanceToFather();
	      node->removeSon(head);
	      node->addSon(new_root);
	      new_root->addSon(head);
	      new_root->setDistanceToFather(d/2.);
	      head->setDistanceToFather(d/2.);
        }
      else
        {
	      if (!node->hasDistanceToFather())
          node->setDistanceToFather(0.1);
	      
	      double d=node->getDistanceToFather();
	      head->removeSon(node);
	      head->addSon(new_root);
	      new_root->addSon(node);
	      new_root->setDistanceToFather(d/2.);
	      node->setDistanceToFather(d/2.);
        }     
      gin_trees[tree_i]->rootAt(new_root);
      return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    else //DTL likelihood summed over roots and reconciliations.
      {
      return_pair.first=log(psum);
      return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    
    return_v.push_back(return_pair);
    tree_i++;  
    }

  
  // delocate
  cherries.clear();
  cherry_reg.clear();
  leaf_hidden_event_name.clear();
  leaf_hidden_event_x.clear();
  max_hidden_events.clear();
  max_event_x.clear();
  max_event_xc1.clear();
  max_event_xc2.clear();
  max_event_name.clear();

  for (std::vector< std::vector< std::pair<node_type*,node_type*> > >::iterator it=run_vec.begin();it!=run_vec.end();it++)
    (*it).clear();
  run_vec.clear();  


  //for (std::vector< std::pair<node_type*,node_type*> > ::iterator jt=it.begin();jt!=it.end();jt++)
  for (std::map < std::pair<node_type *,node_type *>,  long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)
    delete ur_a[(*it).first];  
  ur_a.clear();

  for (std::vector< std::vector< std::pair<node_type*,node_type*> > > ::iterator it=trees_roots.begin();it!=trees_roots.end();it++)
    (*it).clear();
  trees_roots.clear();

  
  for ( std::vector<tree_type *>::iterator it= gin_trees.begin();it!=gin_trees.end();it++)
    delete (*it);
  gin_trees.clear();
  in_trees.clear();
  return return_v;
}





/**************************************************************************
 * This function roots a tree between two nodes.
 **************************************************************************/

void rootBetweenTwoNodes(pair<node_type*,node_type*> max_rootsum_dnode, tree_type * tree) {
  Node* head=max_rootsum_dnode.first;
  Node* node=max_rootsum_dnode.second;
  vector <Node*> head_neigbours = head->getSons();
  if (head->hasFather())
    head_neigbours.push_back(head->getFather());
  vector <Node*> node_neigbours = node->getSons();
  if (node->hasFather())
    node_neigbours.push_back(node->getFather());
  int node_id=-1;
  
  if (node==head_neigbours[0])
    node_id=0;
  else if (node==head_neigbours[1])
    node_id=1;
  else if (node==head_neigbours[2])
    node_id=2;
  Node * new_root = new Node();
  new_root->setName("R");
  if (( node_id==2 && head->hasFather() ) || head->isLeaf())
    {
    if (!head->hasDistanceToFather())
      head->setDistanceToFather(0.1);
    double d=head->getDistanceToFather();
    node->removeSon(head);
    node->addSon(new_root);
    new_root->addSon(head);
    new_root->setDistanceToFather(d/2.);
    head->setDistanceToFather(d/2.);
    }
  else
    {
    if (!node->hasDistanceToFather())
      node->setDistanceToFather(0.1);
    
    double d=node->getDistanceToFather();
    head->removeSon(node);
    head->addSon(new_root);
    new_root->addSon(node);
    new_root->setDistanceToFather(d/2.);
    node->setDistanceToFather(d/2.);
    }     
  tree->rootAt(new_root);
  

}

/**************************************************************************
 * This function roots a gene tree between two nodes, using their node Ids.
 **************************************************************************/
void rootBetweenTwoNodesWithIds(pair<node_type*,node_type*> max_rootsum_dnode, tree_type * tree) {
  pair <Node*,Node*> pairOfNodes;
  pairOfNodes.first = tree->getNode(max_rootsum_dnode.first->getId());
  pairOfNodes.second = tree->getNode(max_rootsum_dnode.second->getId());
  rootBetweenTwoNodes(pairOfNodes,tree);
}



/**************************************************************************
 * This function does the same as LL_treewise above, but takes as input gene trees, 
 * with gene names at the leaves, not species names. Gene trees should have 
 * Node Ids numbered as usual with Bio++. Running resetNodesId() on these
 * gene trees is probably safe.
 **************************************************************************/



vector < pair<scalar_type,string> > Species_tree::LL_treewiseWithSeqSpTranslation(vector<tree_type * >  trees, scalar_type root_p, vector <std::map<std::string, std::string > > seqSp)
{
  vector < pair<scalar_type,string> > return_v;
  
  vector<string>  translatedTrees ;
  
  for (int i=0;i<trees.size();i++)
    {
    translatedTrees.push_back( geneTreeToParenthesisWithSpeciesNames(trees[i], seqSp[i]) );
    stringstream out;out<<i;
    codes[translatedTrees[i]]=out.str();  
    }  
  unrooted_register(translatedTrees);
  
  //computation
  
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ri= run_vec.begin();ri!= run_vec.end();ri++)
    if (mode=="tree")
      backstep_m_N2((*ri)[0],(*ri)[1],(*ri)[2]);
    else
      DPstep_L_N2((*ri)[0],(*ri)[1],(*ri)[2]);
  //sum LL
  
  scalar_type norm=0;
  map <int,scalar_type> p_omega;
  
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
    if (type_of_x[x]!=1)
      if (x_down(x)==1)
        norm+=root_p;
      else
        norm+=1;
  
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
    if (type_of_x[x]!=1)
      if (x_down(x)==1)
        p_omega[x]=root_p/norm;
      else
        p_omega[x]=1/norm;
  
  int tree_i=0; //index of the gene tree we're considering.  
  for ( std::vector< std::vector< std::pair<node_type*,node_type*> > > :: iterator ts= trees_roots.begin(); ts!= trees_roots.end();ts++ )
    {
    
    scalar_type psum=0;
    scalar_type root_norm=(*ts).size();
    pair <Node*,Node*> max_root_dnode;
    pair<Node*,Node*> max_rootsum_dnode;
    
    int max_root_x=-1;
    scalar_type max_root_L=-1;
    scalar_type max_root_sum_L=-1;
    
    for (std::vector< std::pair<node_type*,node_type*> > :: iterator ri= (*ts).begin();ri!= (*ts).end();ri++)	
      {
      pair <Node*,Node*> root_dnode=(*ri);
      scalar_type root_sum_L=0;
      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
        if (type_of_x[x]!=1)
          {
          scalar_type tmp=0;
          if (cherries.count(root_dnode)==0)
            tmp=(*ur_a[root_dnode])[x]/norm/root_norm*p_omega[x];
          else
            tmp=(*cherry_a[cherries[root_dnode]])[x]/norm/root_norm*p_omega[x];
          if (max_root_L<tmp)
            {
            max_root_L=tmp;
            max_root_x=x;
            max_root_dnode.first=root_dnode.first;
            max_root_dnode.second=root_dnode.second;
            }
          psum+=tmp;
          root_sum_L+=tmp;
          }	  	  
      if (max_root_sum_L<root_sum_L)
        {
	      max_root_sum_L=root_sum_L;
	      max_rootsum_dnode=root_dnodes[root_dnode];
        }
      }
    pair<scalar_type,string> return_pair;
    if (mode=="tree") // DTL likelihood using ML reconciliations and ML root
      {	  
        backtrace(max_root_dnode,gin_trees[tree_i], max_root_x);	  
        return_pair.first=log(max_root_L);
        return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    else if (mode=="root") //DTL likelihood summed over reconciliations, using ML root
      {
      return_pair.first=log(max_root_sum_L);
      rootBetweenTwoNodesWithIds(max_rootsum_dnode, trees[tree_i]); //We appropriately root the gene tree with gene names 
     /* 
      Node* head=max_rootsum_dnode.first;
      Node* node=max_rootsum_dnode.second;
      vector <Node*> head_neigbours = head->getSons();
      if (head->hasFather())
        head_neigbours.push_back(head->getFather());
      vector <Node*> node_neigbours = node->getSons();
      if (node->hasFather())
        node_neigbours.push_back(node->getFather());
      int node_id=-1;
      
      if (node==head_neigbours[0])
        node_id=0;
      else if (node==head_neigbours[1])
        node_id=1;
      else if (node==head_neigbours[2])
        node_id=2;
      Node * new_root = new Node();
      new_root->setName("R");
      if (( node_id==2 && head->hasFather() ) || head->isLeaf())
        {
	      if (!head->hasDistanceToFather())
          head->setDistanceToFather(0.1);
	      double d=head->getDistanceToFather();
	      node->removeSon(head);
	      node->addSon(new_root);
	      new_root->addSon(head);
	      new_root->setDistanceToFather(d/2.);
	      head->setDistanceToFather(d/2.);
        }
      else
        {
	      if (!node->hasDistanceToFather())
          node->setDistanceToFather(0.1);
	      
	      double d=node->getDistanceToFather();
	      head->removeSon(node);
	      head->addSon(new_root);
	      new_root->addSon(node);
	      new_root->setDistanceToFather(d/2.);
	      node->setDistanceToFather(d/2.);
        }     
      gin_trees[tree_i]->rootAt(new_root);
      */
      
      return_pair.second=TreeTemplateTools::treeToParenthesis(*trees[tree_i],false);
     // return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    else //DTL likelihood summed over roots and reconciliations.
      {
      return_pair.first=log(psum);
      return_pair.second=TreeTemplateTools::treeToParenthesis(*trees[tree_i],false);

   //   return_pair.second=TreeTemplateTools::treeToParenthesis(*gin_trees[tree_i],false,"ID");
      }
    
    return_v.push_back(return_pair);
    tree_i++;  
    }
  
  
  // delocate
  cherries.clear();
  cherry_reg.clear();
  leaf_hidden_event_name.clear();
  leaf_hidden_event_x.clear();
  max_hidden_events.clear();
  max_event_x.clear();
  max_event_xc1.clear();
  max_event_xc2.clear();
  max_event_name.clear();
  
  for (std::vector< std::vector< std::pair<node_type*,node_type*> > >::iterator it=run_vec.begin();it!=run_vec.end();it++)
    (*it).clear();
  run_vec.clear();  
  
  
  //for (std::vector< std::pair<node_type*,node_type*> > ::iterator jt=it.begin();jt!=it.end();jt++)
  for (std::map < std::pair<node_type *,node_type *>,  long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)
    delete ur_a[(*it).first];  
  ur_a.clear();
  
  for (std::vector< std::vector< std::pair<node_type*,node_type*> > > ::iterator it=trees_roots.begin();it!=trees_roots.end();it++)
    (*it).clear();
  trees_roots.clear();
  
  
  for ( std::vector<tree_type *>::iterator it= gin_trees.begin();it!=gin_trees.end();it++)
    delete (*it);
  gin_trees.clear();
  in_trees.clear();
  return return_v;
}



  
