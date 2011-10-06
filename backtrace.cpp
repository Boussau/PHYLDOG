#include "DTL.h"
#include <float.h>
using namespace std;
using namespace bpp;

void Species_tree:: backstep_m_N2(pair <Node*,Node*> dnode ,pair <Node*,Node*>  left_dnode,pair <Node*,Node*>  right_dnode)
{

  pair <Node*,Node*> child_1, child_2;
  child_1 = left_dnode;
  child_2 = right_dnode;
  dnode_c1[dnode]=child_1;
  dnode_c2[dnode]=child_2;

  bool is_leaf_1,is_leaf_2;
  is_leaf_1=child_1.first->isLeaf();
  is_leaf_2=child_2.first->isLeaf();

  bool imaroot= (dnode.first==dnode.second) || (is_leaf_1 && is_leaf_2);
  if (imaroot)
    root_max[dnode]=0;
  if (imaroot)
    root_max_s[dnode]=0;

  scalar_type max_root=0;
  int max_x=-1;

  scalar_type max_D_c=0;
  scalar_type max_T_c=0;

  int leaf_1=-1,leaf_2=-1,xp,xpp,x,z,y;
  long_vector_type *ax;
  int_vector_type *ax_max_event_name,*ax_max_event_x,*ax_max_event_xc1,*ax_max_event_xc2;
  map<int,vector<int> > * ax_max_hidden_events;
  vector_type *fQef_xp,*fQef_xpp,*fQef_x;

  long_vector_type *a_c1,*a_c2;
  scalar_type sev,sew,sfv,sfw,sum_tmp,tv_sum_tmp,tw_sum_tmp,T_sum,D_sum,scalar_tmp,tmp_Qxy,tmp_ay;
  scalar_type max_sev,max_sew,max_sfv,max_sfw,tmp_1,tmp_2,max_T_sum,max_D_sum,max_tv_sum_tmp,max_tw_sum_tmp,max_sum_tmp;
  scalar_type nd2_tv_sum_tmp,nd2_tw_sum_tmp,nd2_T_sum,nd2_D_sum;
  int T_y,D_y,tv_x,tw_x,nd2_tw_x,nd2_tv_x;

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

	  int tmp=leaf_1;
	  leaf_1=leaf_2;
	  leaf_2=tmp;
	  bool tmb=is_leaf_1;
	  is_leaf_1=is_leaf_2;
	  is_leaf_2=tmb;
	  dnode_c1[dnode]=child_2;
	  dnode_c2[dnode]=child_1;
  
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
	  if (cherry_a.count(cherry)==0 ||  cherry_max_event_name.count(cherry)==0 ||  cherry_max_event_x.count(cherry)==0 ||  cherry_max_event_xc1.count(cherry)==0 ||  cherry_max_event_xc2.count(cherry)==0)
	    {
	      cherry_a[cherry] = new long_vector_type(lin_size);

	      cherry_max_event_name[cherry]=new int_vector_type(lin_size);
	      cherry_max_event_x[cherry]=new int_vector_type(lin_size);
	      cherry_max_event_xc1[cherry]=new int_vector_type(lin_size);
	      cherry_max_event_xc2[cherry]=new int_vector_type(lin_size);
	      for ( map <int , vector<int> > ::iterator it=cherry_max_hidden_events[cherry].begin();it!=cherry_max_hidden_events[cherry].end();it++)
		(*it).second.clear();
	      cherry_max_hidden_events[cherry].clear(); 
	    }
	  else
	    {
	      cherry_a[cherry]->resize(lin_size);

	      cherry_max_event_name[cherry]->resize(lin_size);
	      cherry_max_event_x[cherry]->resize(lin_size);
	      cherry_max_event_xc1[cherry]->resize(lin_size);
	      cherry_max_event_xc2[cherry]->resize(lin_size);
	      for ( map <int , vector<int> > ::iterator it=cherry_max_hidden_events[cherry].begin();it!=cherry_max_hidden_events[cherry].end();it++)
		(*it).second.clear();
	      cherry_max_hidden_events[cherry].clear(); 	      
	    }
	  ax = cherry_a[cherry];

	  ax_max_event_name = cherry_max_event_name[cherry];
	  ax_max_event_xc1 = cherry_max_event_xc1[cherry];
	  ax_max_event_xc2 = cherry_max_event_xc2[cherry];
	  ax_max_event_x = cherry_max_event_x[cherry];
	  ax_max_hidden_events = &cherry_max_hidden_events[cherry]; 

	  cherries[dnode]=cherry;
	  cherry_reg[cherry]=1;
	}
    }
  else
    {
      for ( map <int , vector<int> > ::iterator it=max_hidden_events[dnode].begin();it!=max_hidden_events[dnode].end();it++)
	(*it).second.clear();
      max_hidden_events[dnode].clear(); 
      //XXZZ
      if (ur_a.count(dnode)==0 || max_event_name.count(dnode)==0 || max_event_x.count(dnode)==0 || max_event_xc1.count(dnode)==0 || max_event_xc2.count(dnode)==0 )
	{
	  ur_a[dnode] = new long_vector_type(lin_size);

	  max_event_name[dnode]=new int_vector_type(lin_size);
	  max_event_x[dnode]=new int_vector_type(lin_size);
	  max_event_xc1[dnode]=new int_vector_type(lin_size);
	  max_event_xc2[dnode]=new int_vector_type(lin_size);	  
	}
      else
	ur_a[dnode]->resize(lin_size);

      ax = ur_a[dnode];

      ax_max_event_name = max_event_name[dnode];
      ax_max_event_xc1 = max_event_xc1[dnode];
      ax_max_event_xc2 = max_event_xc2[dnode];
      ax_max_event_x = max_event_x[dnode];
      ax_max_hidden_events = &max_hidden_events[dnode]; 
      
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
  int T_sum_x;
  int D_sum_x;

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
	max_tv_sum_tmp=tv_sum_tmp;
	nd2_tv_x=tv_x;
	tv_x=xpp;
      }
    //tv_sum_tmp +=scalar_tmp*sfv;
    tw_sum_tmp = scalar_tmp*sfw;
    if (max_tw_sum_tmp<tw_sum_tmp) 
      {
	nd2_tw_sum_tmp=max_tw_sum_tmp;
	max_tw_sum_tmp=tw_sum_tmp;		  
	nd2_tw_x=tw_x;
	tw_x=xpp;
      }
    //tw_sum_tmp +=scalar_tmp*sfw;
    // ############## T #############
  }

  double norm=lin_size-ur_S_N_leaves,max_S=0;
  for (x=lin_size-ur_S_N_leaves;x>=0;x--){
    //BBBB//
    (*ax)(x)=0;
    if (type_of_x[x] == 1){
      // e and f
      xp = below_son0[x];
      xpp = below_son1[x];		   
      if (is_leaf_1 || is_leaf_2){
	fQef_xp = & flat_Qef[xp];
	fQef_xpp = & flat_Qef[xpp];}
      // a speciation
      fQef_x = & flat_Qef[x];
      // ############## S #############		  
      // place v in e and f
      if (is_leaf_1){
	sev = (*fQef_xp)( leaf_1 - xp);
	sfv = (*fQef_xpp)( leaf_1 - xpp);}
      else{
	//sev=(*fQef_x)(xp-x)*(*a_c1)(xp);
	//sfv=(*fQef_x)(xpp-x)*(*a_c1)(xpp);
	sev=(*a_c1)(xp);
	sfv=(*a_c1)(xpp);
      }
      // place w in e and f
      if (is_leaf_2){
	sew = (*fQef_xp)( leaf_2 - xp);
	sfw = (*fQef_xpp)( leaf_2 - xpp);}
      else{
	//sew=(*fQef_x)(xp-x)*(*a_c2)(xp);
	//sfw=(*fQef_x)(xpp-x)*(*a_c2)(xpp);
	sew=(*a_c2)(xp);
	sfw=(*a_c2)(xpp);
	
      }
      //BBBB//
      //max_sum_tmp=max(sev*sfw , sfv*sew);
      max_sum_tmp=Qegp[xp]*Qegp[xpp]*max(sev*sfw , sfv*sew);

      if (max_S<max_sum_tmp || 1)
	{
	  (*ax_max_event_name)(x)=1;
	  (*ax_max_event_x)(x)=x;     
	  max_S=max_sum_tmp;
	  if (sev*sfw > sfv*sew)
	    {
	      (*ax_max_event_xc1)(x)=xp;
	      (*ax_max_event_xc2)(x)=xpp;
	    }	
	  else
	    {
	      (*ax_max_event_xc1)(x)=xpp;
	      (*ax_max_event_xc2)(x)=xp;
	    }
	}	

      // ############## S #############
      // ############## SL #############
      sum_tmp=(*fQef_x)(xp-x)*(*ax)(xp);	      
      //BBBB//
      if (max_sum_tmp<sum_tmp) 	{
	  max_sum_tmp=sum_tmp;
	  (*ax_max_event_xc1)(x)=(*ax_max_event_xc1)(xp);
	  (*ax_max_event_xc2)(x)=(*ax_max_event_xc2)(xp);
	  (*ax_max_event_name)(x)=2;
	  (*ax_max_hidden_events)[x].push_back(xp);     
	  (*ax_max_event_x)(x)=x;
	}
      //sum_tmp+=(*fQef_x)(xp-x)*(*ax)(xp);	   
      sum_tmp=(*fQef_x)(xpp-x)*(*ax)(xpp);	      
      //BBBB//
      if (max_sum_tmp<sum_tmp) 	
	{
	  max_sum_tmp=sum_tmp;
	  (*ax_max_event_xc1)(x)=(*ax_max_event_xc1)(xpp);
	  (*ax_max_event_xc2)(x)=(*ax_max_event_xc2)(xpp);
	  (*ax_max_event_name)(x)=2;
	  (*ax_max_hidden_events)[x].push_back(xpp);     
	  (*ax_max_event_x)(x)=x;
	}
      //sum_tmp+=(*fQef_x)(xpp-x)*(*ax)(xpp);	   
      // ############## SL #############
      //BBBB//
      (*ax)(x)=max_sum_tmp;

      if (imaroot)
	{
	  root_max_s[dnode]+=(*ax)(x)*p_omega(x);
	  //cout << ":"<<dnode.first<<dnode.second << " " << (*ax)(x)/norm <<endl;
	  if(max_root<=(*ax)(x)*p_omega(x) )//XX
	    {
	      max_root=(*ax)(x)*p_omega(x);       
	      max_x=x;
	    }
	}

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

	tv_sum_tmp = scalar_tmp*sfv;
	if (max_tv_sum_tmp<tv_sum_tmp) 
	  {
	    nd2_tv_sum_tmp=max_tv_sum_tmp;
	    max_tv_sum_tmp=tv_sum_tmp;
	    nd2_tv_x=tv_x;
	    tv_x=xpp;
	  }
	//tv_sum_tmp +=scalar_tmp*sfv;
	tw_sum_tmp = scalar_tmp*sfw;
	if (max_tw_sum_tmp<tw_sum_tmp) 
	  {
	    nd2_tw_sum_tmp=max_tw_sum_tmp;
	    max_tw_sum_tmp=tw_sum_tmp;		  
	    nd2_tw_x=tw_x;
	    tw_x=xpp;
	  }
	T_sum=(*fQef_xpp)(y-xpp) * scalar_tmp * tmp_ay;
	if (max_T_sum<T_sum) 
	  {
	    max_T_sum=T_sum;       
	    T_sum_x=y;
	  }
	D_sum= tau_dt[y] * tmp_ay;	
	if (max_D_sum<D_sum) 
	  {
	    max_D_sum=D_sum;       
	    D_sum_x=y;	   
	  }
	// ############## T #############
      }
    }
    else {
      sum_tmp = 0;
      //BBBB//
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
      //BBBB//
      if (max_sum_tmp<sum_tmp) 
	{
	  max_sum_tmp=sum_tmp;
	  (*ax_max_event_name)(x)=3;
	  (*ax_max_event_x)(x)=x;     
	  (*ax_max_event_xc1)(x)=x;
	  (*ax_max_event_xc2)(x)=x;
	}
      //sum_tmp +=  delta_dt[x] * sev*sew;	   
      // ############## D #############
	
      // ############## T #############	  
      int tmp_vx,tmp_wx;
      if (x==tv_x && x==tw_x)
	{
	  sum_tmp = max(nd2_tv_sum_tmp*sew , nd2_tw_sum_tmp*sev);
	  tmp_vx=nd2_tv_x;tmp_wx=nd2_tw_x; 
	  if (max_tv_sum_tmp*sew > max_tw_sum_tmp*sev)
	    {tmp_vx=tv_x;tmp_wx= x;} 
	  else
	    {tmp_vx= x;tmp_wx=tw_x;}
	}
      else if (x==tv_x)
	{
	  sum_tmp = max(nd2_tv_sum_tmp*sew , max_tw_sum_tmp*sev);
	  tmp_vx=nd2_tv_x;tmp_wx=tw_x; 
	  if (max_tv_sum_tmp*sew > max_tw_sum_tmp*sev)
	    {tmp_vx=tv_x;tmp_wx= x;}
	  else
	    {tmp_vx= x;tmp_wx=x;}
	}
      else if (x==tw_x )
	{
	  sum_tmp = max(max_tv_sum_tmp*sew , nd2_tw_sum_tmp*sev);
	  tmp_vx=tv_x;tmp_wx=nd2_tw_x; 
	  if (max_tv_sum_tmp*sew > max_tw_sum_tmp*sev)
	    {tmp_vx=tv_x;tmp_wx= x;} 
	  else
	    {tmp_vx= x;tmp_wx=tw_x;} 
	}
      else
	{
	  sum_tmp = max(max_tv_sum_tmp*sew , max_tw_sum_tmp*sev);
	  if (max_tv_sum_tmp*sew > max_tw_sum_tmp*sev)
	    {tmp_vx=tv_x;tmp_wx= x;} 
	  else
	    {tmp_vx= x;tmp_wx=tw_x;} 	   
	}

      //sum_tmp = max(max_tv_sum_tmp*sew , max_tw_sum_tmp*sev) - 1*(tau_dt[x]*sew*sev);	 
      
      //BBBB//
      if (max_sum_tmp<sum_tmp) 
	{
	  max_sum_tmp=sum_tmp;
	  (*ax_max_event_name)(x)=4;
	  (*ax_max_event_x)(x)=x;     

	  if (max_tv_sum_tmp*sew > max_tw_sum_tmp*sev)
	    {
	      (*ax_max_event_xc1)(x)=tv_x;
	      (*ax_max_event_xc2)(x)=x;
	    }
	  else
	    {
	      (*ax_max_event_xc1)(x)=x;
	      (*ax_max_event_xc2)(x)=tw_x;
	    }
	}

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
	//BBBB//

	if (max_sum_tmp<sum_tmp) 
	  {
	    max_sum_tmp=sum_tmp;
	    (*ax_max_event_xc1)(x)=(*ax_max_event_xc1)(y);
	    (*ax_max_event_xc2)(x)=(*ax_max_event_xc2)(y);
	    (*ax_max_event_name)(x)=(*ax_max_event_name)(y);
	    (*ax_max_event_x)(x)=(*ax_max_event_x)(y);
	  }
	//sum_tmp+=scalar_tmp;	      
	sum_tmp= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (max_T_sum);	      
	//BBBB//

	if (max_sum_tmp<sum_tmp) 
	  {
	    max_sum_tmp=sum_tmp;
	    (*ax_max_event_xc1)(x)=(*ax_max_event_xc1)(T_sum_x);
	    (*ax_max_event_xc2)(x)=(*ax_max_event_xc2)(T_sum_x);
	    (*ax_max_event_name)(x)=5;
	    (*ax_max_hidden_events)[x].push_back(T_sum_x);     
	    (*ax_max_event_x)(x)=x;
	  }

	//sum_tmp+= Qe[lin_slice[x]][dc_mod][lin_edge[x]] * (T_sum -  tau_dt[x] * scalar_tmp);	      
	sum_tmp= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * max_D_sum );     
	//BBBB//
	if (max_sum_tmp<sum_tmp)
	  {
	    max_sum_tmp=sum_tmp;
	    (*ax_max_event_xc1)(x)=(*ax_max_event_xc1)(D_sum_x);
	    (*ax_max_event_xc2)(x)=(*ax_max_event_xc2)(D_sum_x);
	    (*ax_max_event_name)(x)=5;
	    (*ax_max_hidden_events)[x].push_back(D_sum_x);     
	    (*ax_max_event_x)(x)=x;
	  }

	//sum_tmp+= Qe[lin_slice[y]][dc_mod][lin_edge[y]] * ( tmp_Qxy * D_sum - scalar_tmp * tau_dt[y]) ;
      }
      // ############## 0&TL #############	  
      //BBBB//
      (*ax)(x) = max_sum_tmp;
      if (imaroot)
	{
	  root_max_s[dnode]+=(*ax)(x)*p_omega(x);
	  //cout << ":"<<dnode.first<<dnode.second << " " << (*ax)(x)/norm <<endl;	  
	  if(max_root<=(*ax)(x)*p_omega(x))//XX
	    {
	      max_root=(*ax)(x)*p_omega(x);  
	      max_x=x;
	    }
	}
    }

  }
  if (imaroot)//XX
    {
      root_max[dnode]=max_root;
      root_max_x[dnode]=max_x;
    }

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
	  int tmp=leaf_1;
	  leaf_1=leaf_2;
	  leaf_2=tmp;
	  bool tmb=is_leaf_1;
	  is_leaf_1=is_leaf_2;
	  is_leaf_2=tmb;
	  dnode_c1[dnode]=child_2;
	  dnode_c2[dnode]=child_1;
  
	}      
      cherry_root_max[cherry]=root_max[dnode];
      cherry_root_max_s[cherry]=root_max_s[dnode];
      cherry_root_max_x[cherry]=root_max_x[dnode];
    }
}

int Species_tree::x_down_f(int x)
{
  while (type_of_x[x]!=1 && lin_slice[x]>0)
    x=below_son0[x];
  if (lin_slice[x]==0)
    return lin_edge[x]+N_slices;
  else
    return lin_slice[0]-lin_slice[x]+1;
}


int  Species_tree::node_down_f(int x)
{
  while (type_of_x[x]!=1 && lin_slice[x]>0)
    x=below_son0[x];
  return  x;
}

void  Species_tree::x_DFS()
{
  x_DFS(tree->getRootNode());
  return;
}
void  Species_tree::x_DFS(node_type * node)
{
  if (node->isLeaf())
    {
      (*event_stream)<<x_down(flat_node_map[node]) << ",";
      return;
    }
  vector <Node*> sons=node->getSons();
  if (sons.size()==1)
    {
      x_DFS(sons[0]);
      return;
    }
  x_DFS(sons[0]);
  (*event_stream)<<x_down(flat_node_map[node]) << ",";
  x_DFS(sons[1]);
  return;
}

string  Species_tree::legend()
{
  stringstream out;
  out << "# " << "node_id" << "  " << "branch_id" << "  " << "time_slice" << "  " << "name" << endl;

  for (int x=0;x<lin_size;x++)
    {
      out << "# " << x << " " << x_down(x) << " " << lin_slice[0]-lin_slice[x]+1 << " " << type_of_x[x] << " " <<inverse_flat_node_map[x]->getName() << endl;
    }
  return out.str();
}
void Species_tree::backtrace(pair <Node*,Node*> root_dnode, tree_type * tree,int in_o_x)
{  

  pair <Node*,Node*> dnode = root_dnodes[root_dnode];
  int e_type=0;
  stringstream out;
  if (cherries.count(root_dnode)==0)
    {
      int o_x = in_o_x;//root_max_x[root_dnode];
      string event_name="";
      string hidden_name="";

      e_type=(*max_event_name[root_dnode])(o_x);
      int nx;
      //####################################################################################################
      while (e_type==2 || e_type==5 )
	{
	  while  (max_hidden_events[root_dnode].count(o_x)==0)
	    o_x=(*max_event_x[root_dnode])(o_x);	  	  

	  nx=(max_hidden_events[root_dnode][o_x]).back();

	  if (e_type==2)
	    {
	      stringstream hout;
	      hout<<"&SL@";	     
	      hout<<lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x) <<"|"<< o_x;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	  else
	    {
	      stringstream hout;
	      hout<<"&TL@";
	      hout<<lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x) << "|"<< x_down(nx) <<"|"<< o_x;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	    
	  //max_hidden_events[root_dnode][nx].pop_back();
	  e_type=(*max_event_name[root_dnode])(nx);
	  o_x=(*max_event_x[root_dnode])(nx);	  	  
	}
      //####################################################################################################
      out<<  lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x);

      int at_x=o_x;

      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	{
	  event_name="T";
	  if ((*max_event_xc1[root_dnode])(o_x)==o_x)
	    {
	      at_x=(*max_event_xc2[root_dnode])(o_x);
	      out << "|" <<x_down((*max_event_xc2[root_dnode])(o_x));
	    }
	  else  
	    {
	      at_x=(*max_event_xc1[root_dnode])(o_x);
	      out << "|" <<x_down((*max_event_xc1[root_dnode])(o_x));
	    }
	}
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();

      Node* head=dnode.first;
      Node* node=dnode.second;
      vector <Node*> head_neigbours = head->getSons();
      if (head->hasFather())
	head_neigbours.push_back(head->getFather());
      vector <Node*> node_neigbours = node->getSons();
      //typo? BB
      //if (head->hasFather())
      if (head->hasFather())

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
      new_root->setBranchProperty("ID",BppString(event_name+"@"+node_name+hidden_name));

      backtrace(dnode_c1[root_dnode],(*max_event_xc1[root_dnode])(o_x));
      if (dnode_c1[root_dnode].first->isLeaf())
	name_leaf(dnode_c1[root_dnode].first,(*max_event_xc1[root_dnode])(o_x));
      backtrace(dnode_c2[root_dnode],(*max_event_xc2[root_dnode])(o_x));    
      if (dnode_c2[root_dnode].first->isLeaf())
	name_leaf(dnode_c2[root_dnode].first,(*max_event_xc2[root_dnode])(o_x));
  
      //(*event_stream)<< TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
      //cout<< TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
      //(*event_stream) << stream_xL(tree) << endl;
    }
  else
    {
      pair<int,int> cherry=cherries[root_dnode];
      //int o_x = cherry_root_max_x[cherry];
      int o_x = in_o_x;//root_max_x[root_dnode];

      string event_name="";
      string hidden_name="";

      e_type=(*cherry_max_event_name[cherry])(o_x);
      int nx;
      //####################################################################################################
      while (e_type==2 || e_type==5 )
	{
	  while  (cherry_max_hidden_events[cherry].count(o_x)==0)
	    o_x=(*cherry_max_event_x[cherry])(o_x);	  	  

	  nx=cherry_max_hidden_events[cherry][o_x].back();

	  if (e_type==2)
	    {
	      stringstream hout;
	      hout<<"&SL@";
	      hout<<lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x) <<"|"<< o_x;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	  else
	    {
	      stringstream hout;
	      hout<<"&TL@";	      
	      hout<<lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x) << "|"<< x_down(nx) <<"|"<< o_x;
	      hidden_name=hout.str()+hidden_name;	  
	    }	  
	  //cherry_max_hidden_events[cherry][nx].pop_back();
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  o_x=(*cherry_max_event_x[cherry])(nx);	  	  
	}
      //####################################################################################################
      out<<  lin_slice[0]-lin_slice[o_x]+1 << "|" << x_down(o_x);
      int at_x=o_x;

      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	{
	  event_name="T";
	  if ((*cherry_max_event_xc1[cherry])(o_x)==o_x)
	    {
	      at_x=(*cherry_max_event_xc2[cherry])(o_x);
	      out << "|" <<x_down((*cherry_max_event_xc2[cherry])(o_x));
	    }
	  else  
	    {
	      at_x=(*cherry_max_event_xc1[cherry])(o_x);
	      out << "|" <<x_down((*cherry_max_event_xc1[cherry])(o_x));
	    }
	}
      else if (e_type==5)
	event_name="0TL";
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();
      tree->getRootNode()->setBranchProperty("ID",BppString(event_name+"@"+node_name+hidden_name));
      name_leaf(dnode_c1[root_dnode].first,(*cherry_max_event_xc1[cherry])(o_x));
      name_leaf(dnode_c2[root_dnode].first,(*cherry_max_event_xc2[cherry])(o_x));

      //(*event_stream)<< TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
      //cout<<TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
    }
} 
void Species_tree::backtrace(pair <Node*,Node*> dnode,int x)
{
  
  if (cherries.count(dnode)==0)
    {
      if (max_event_name.count(dnode)==0)
	return;
      int e_type=(*max_event_name[dnode])(x);
      int ex=(*max_event_x[dnode])(x);
      string hidden_name;
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  nx=(max_hidden_events[dnode][ex]).back();

	  if (e_type==2)
	    {
	      stringstream hout;
	      hout<<"&SL@";
	      hout<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex) <<"|"<< ex;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	  else
	    {
	      stringstream hout;
	      hout<<"&TL@";
	      hout<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex) << "|"<< x_down(nx) <<"|"<< ex;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	    
	  //max_hidden_events[dnode][nx].pop_back();
	  e_type=(*max_event_name[dnode])(nx);
	  ex=(*max_event_x[dnode])(nx);	  	  
	}
      //####################################################################################################


      stringstream out;
      out<<lin_slice[0]-lin_slice[ex]+1;
      int at_x=ex;
      string event_name;      

      if (e_type==1)
	{
	  event_name="S";
	  out << "|" << x_down(ex);
	}
      else if (e_type==3)
	{
	  event_name="D";
	  out << "|" << x_down(ex);	  
	}
      else if (e_type==4)
	{
	  event_name="T";
	  if ((*max_event_xc1[dnode])(ex)==ex)
	    {
	      at_x=(*max_event_xc2[dnode])(ex);
	      out<< "|" << x_down(ex) << "|" <<x_down((*max_event_xc2[dnode])(ex)) ;
	    }
	  else  
	    {
	      at_x=(*max_event_xc1[dnode])(ex);
	      out << "|" << x_down(ex) << "|" <<x_down((*max_event_xc1[dnode])(ex));
	    }
	}
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();
      dnode.first->setBranchProperty("ID",BppString(event_name+"@"+node_name+hidden_name));
      backtrace(dnode_c1[dnode],(*max_event_xc1[dnode])(x));
      if (dnode_c1[dnode].first->isLeaf())
	name_leaf(dnode_c1[dnode].first,(*max_event_xc1[dnode])(ex));
      backtrace(dnode_c2[dnode],(*max_event_xc2[dnode])(x));
      if (dnode_c2[dnode].first->isLeaf())
	name_leaf(dnode_c2[dnode].first,(*max_event_xc2[dnode])(ex));
	
    }
  else 
    {

      pair<int,int> cherry=cherries[dnode];
      if (cherry_max_event_name.count(cherry)==0)
	return;

      string hidden_name;
      int e_type=(*cherry_max_event_name[cherry])(x);
      int ex=(*cherry_max_event_x[cherry])(x);
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  //oo
	  nx=cherry_max_hidden_events[cherry][ex].back();

	  if (e_type==2)
	    {
	      stringstream hout;
	      hout<<"&SL@";
	      hout<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex) <<"|"<< ex;
	      hidden_name=hout.str()+hidden_name;	  
	    }
	  else
	    {
	      stringstream hout;
	      hout<<"&TL@";
	      hout<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex) << "|"<< x_down(nx) <<"|"<< ex;
	      hidden_name=hout.str()+hidden_name;	  
	    }	  
	  //cherry_max_hidden_events[cherry][nx].pop_back();
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  ex=(*cherry_max_event_x[cherry])(nx);	  	  
	}
      //####################################################################################################

      stringstream out;
      out<<lin_slice[0]-lin_slice[ex]+1;
      int at_x=ex;
      
      string event_name;

      if (e_type==1)
	{
	  event_name="S";
	  out << "|" << x_down(ex);	  	  
	}
      else if (e_type==3)
	{
	  event_name="D";
	  out << "|" << x_down(ex);	  	 
	}
      else if (e_type==4)
	{
	  event_name="T";
	  if ((*cherry_max_event_xc1[cherry])(ex)==ex)
	    {
	      at_x=(*cherry_max_event_xc2[cherry])(ex);
	      out << "|" << x_down(ex)<< "|" <<x_down((*cherry_max_event_xc2[cherry])(ex));
	    }
	  else  
	    {
	      at_x=(*cherry_max_event_xc1[cherry])(ex);	      
	      out << "|" << x_down(ex)<< "|" <<x_down((*cherry_max_event_xc1[cherry])(ex));
	    }
	}
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();
      dnode.first->setBranchProperty("ID",BppString(event_name+"@"+node_name+hidden_name));
      name_leaf(dnode_c1[dnode].first,(*cherry_max_event_xc1[cherry])(ex));
      name_leaf(dnode_c2[dnode].first,(*cherry_max_event_xc2[cherry])(ex));
    }  

}


void Species_tree::name_leaf(Node * leaf,int x,int leaf_x)
{
  // cherry switch solution
  int other_leaf_x=leaf_x;
  if (leaf_x==-1)
    leaf_x=ur_sigma[leaf->getName()];
  
    
  int nx=x;
  while (nx!=leaf_x)
    {
      int e_type=leaf_hidden_event_name[leaf_x][nx];
      string hidden_name="";

      if (e_type==2)
	{
	  stringstream hout;
	  hout<<"&SL@";
	  hout<<lin_slice[0]-lin_slice[nx]+1<< "|" << x_down(nx) <<"|"<< nx;
	  hidden_name=hout.str()+hidden_name;	  
	}
      else if (e_type==5 || x_down(nx)!=x_down(leaf_x))
	{
	  stringstream hout;
	  hout<<"&TL@";
	  hout<<lin_slice[0]-lin_slice[nx]+1<< "|" << x_down(nx) << "|"<< x_down(leaf_hidden_event_x[leaf_x][nx]) <<"|"<< nx;
	  hidden_name=hout.str()+hidden_name;	  
	}
      nx=leaf_hidden_event_x[leaf_x][nx];
      leaf->setName(leaf->getName()+hidden_name);
    }
}



void Species_tree::tracestream(pair <Node*,Node*> root_dnode, tree_type * tree, int x)
{  
  pair <Node*,Node*> dnode = root_dnodes[root_dnode];
  int e_type=0;
  stringstream out;
  if (cherries.count(root_dnode)==0)
    {
      int o_x;
      if (x==-1)
	o_x = root_max_x[root_dnode];
      else
	o_x=x;
      stringstream closer;
      int nx;
      string event_name="";
      e_type=(*max_event_name[root_dnode])(o_x);
      
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  //<<??>
	  while  (max_hidden_events[root_dnode].count(o_x)==0)
	    o_x=(*max_event_x[root_dnode])(o_x);	  	  
	  nx=(max_hidden_events[root_dnode][o_x]).back();
	  if (e_type==2)
	    {
	      int at_branch=x_down(o_x);
	      closer<<">@"<<at_branch<<"|";
	      (*event_stream)<<"<@"<<at_branch<<"|"<<"S@"<<at_branch<<"|";
	    }
	  else
	    {
	      int at_branch=x_down(nx);
	      closer<<"T>@"<<at_branch<<"|";
	      (*event_stream)<<"<T@"<<at_branch<<"|"<<"T@"<<at_branch<<"|";
	    }
	  //max_hidden_events[root_dnode][nx].pop_back();
	  e_type=(*max_event_name[root_dnode])(nx);
	  o_x=(*max_event_x[root_dnode])(nx);	  	  
	}
      //####################################################################################################
      out<<  lin_slice[0]-lin_slice[o_x]+1<< "|" << x_down(o_x);
      int at_x=o_x;
      int at_branch=x_down(at_x);
      
      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	event_name="T";
      else
	event_name="###########";
      out<<"|"<<at_x;
      
      string node_name=out.str();
      if (event_name=="S")
	{
	  closer<<">@"<<at_branch<<"|";
	  (*event_stream)<<"<@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";
	}
      else if (event_name=="T")
	{
	  if(x_down((*max_event_xc1[root_dnode])(o_x))==at_branch)
	    at_branch=x_down((*max_event_xc2[root_dnode])(o_x));
	  else
	    at_branch=x_down((*max_event_xc1[root_dnode])(o_x));	 
	  closer<<"T>@"<<at_branch<<"|";
	  (*event_stream)<<"<T@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";	  
	}
      else
	(*event_stream)<<event_name<<"@"<<at_branch<<"|";

      tracestream(dnode_c1[root_dnode],(*max_event_xc1[root_dnode])(o_x));

      if (dnode_c1[root_dnode].first->isLeaf())
	tracestream(dnode_c1[root_dnode].first,(*max_event_xc1[root_dnode])(o_x));
      tracestream(dnode_c2[root_dnode],(*max_event_xc2[root_dnode])(o_x));    


      if (dnode_c2[root_dnode].first->isLeaf())
	tracestream(dnode_c2[root_dnode].first,(*max_event_xc2[root_dnode])(o_x));
      (*event_stream)<<closer.str();
    }
  else
    {
      pair<int,int> cherry=cherries[root_dnode];
      int o_x;	
      if (x==-1)
	o_x = cherry_root_max_x[cherry];	
      else
	o_x=x;
      stringstream closer;

      string event_name="";

      e_type=(*cherry_max_event_name[cherry])(o_x);
      int nx;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  while  (cherry_max_hidden_events[cherry].count(o_x)==0)
	    o_x=(*cherry_max_event_x[cherry])(o_x);	  	  
	  nx=cherry_max_hidden_events[cherry][o_x].back();
	  
	  if (e_type==2)
	    {
	      int at_branch=x_down(o_x);
	      closer<<">@"<<at_branch<<"|";
	      (*event_stream)<<"<@"<<at_branch<<"|"<<"S@"<<at_branch<<"|";
	    }
	  else
	    {
	      int at_branch=x_down(nx);
	      closer<<"T>@"<<at_branch<<"|";
	      (*event_stream)<<"<T@"<<at_branch<<"|"<<"T@"<<at_branch<<"|";
	    }	  
	  //cherry_max_hidden_events[cherry][nx].pop_back();
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  o_x=(*cherry_max_event_x[cherry])(nx);	  	  
	}
      //####################################################################################################
      out<<  lin_slice[0]-lin_slice[o_x]+1 << "|" << x_down(o_x);
      int at_x=o_x;
      int at_branch=x_down(at_x);

      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	event_name="T";
      else if (e_type==5)
	event_name="0TL";
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();
      if (event_name=="S")
	{
	  closer<<">@"<<at_branch<<"|";
	  (*event_stream)<<"<@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";
	}
      else if (event_name=="T")
	{
	  if((*cherry_max_event_xc1[cherry])(o_x)==at_branch)
	    at_branch=(*cherry_max_event_xc2[cherry])(o_x);
	  else
	    at_branch=x_down((*cherry_max_event_xc1[cherry])(o_x));	 
	  closer<<"T>@"<<at_branch<<"|";
	  (*event_stream)<<"<T@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";	  
	}
      else
	(*event_stream)<<event_name<<"@"<<at_branch<<"|";
      
      tracestream(dnode_c1[root_dnode].first,(*cherry_max_event_xc1[cherry])(o_x));
      tracestream(dnode_c2[root_dnode].first,(*cherry_max_event_xc2[cherry])(o_x));
      (*event_stream)<<closer.str();

    }
  (*event_stream)<< endl;

} 
void Species_tree::tracestream(pair <Node*,Node*> dnode,int x)
{
  
  if (cherries.count(dnode)==0)
    {
      if (max_event_name.count(dnode)==0)
	return;
      int e_type=(*max_event_name[dnode])(x);
      int ex=(*max_event_x[dnode])(x);
      string hidden_name;
      stringstream closer;
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{	  
	  nx=(max_hidden_events[dnode][ex]).back();

	  if (e_type==2)
	    {
	      int at_branch=x_down(ex);
	      closer<<">@"<<at_branch<<"|";
	      (*event_stream)<<"<@"<<at_branch<<"|"<<"S@"<<at_branch<<"|";
	    }
	  else
	    {
	      int at_branch=x_down(nx);
	      closer<<"T>@"<<at_branch<<"|";
	      (*event_stream)<<"<T@"<<at_branch<<"|"<<"T@"<<at_branch<<"|";
	    }	    
	  //max_hidden_events[dnode][nx].pop_back();	
	  e_type=(*max_event_name[dnode])(nx);
	  int tmp=ex;ex=(*max_event_x[dnode])(nx);	  	  
	  if (tmp==ex) break;//strange exception occured for # HBG460549 cyano trees 3098
	}
      //####################################################################################################

      stringstream out;
      out<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex);
      int at_x=ex;
      int at_branch=x_down(ex);
      
      string event_name;      

      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	event_name="T";
      else
	event_name="###########";
      out<<"|"<<at_x;
      string node_name=out.str();
      if (event_name=="S")
	{
	  closer<<">@"<<at_branch<<"|";
	  (*event_stream)<<"<@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";
	}
      else if (event_name=="T")
	{
	  if(x_down((*max_event_xc1[dnode])(ex))==at_branch)
	    at_branch=x_down((*max_event_xc2[dnode])(ex));
	  else
	    at_branch=x_down((*max_event_xc1[dnode])(ex));	 
	  closer<<"T>@"<<at_branch<<"|";
	  (*event_stream)<<"<T@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";	  
	}
      else
	(*event_stream)<<event_name<<"@"<<at_branch<<"|";
      tracestream(dnode_c1[dnode],(*max_event_xc1[dnode])(x));
      if (dnode_c1[dnode].first->isLeaf())
	tracestream(dnode_c1[dnode].first,(*max_event_xc1[dnode])(ex));
      tracestream(dnode_c2[dnode],(*max_event_xc2[dnode])(x));
      if (dnode_c2[dnode].first->isLeaf())
	tracestream(dnode_c2[dnode].first,(*max_event_xc2[dnode])(ex));
      (*event_stream)<<closer.str();
    }
  else 
    {

      pair<int,int> cherry=cherries[dnode];
      if (cherry_max_event_name.count(cherry)==0)
	return;

      string hidden_name;
      stringstream closer;
      int e_type=(*cherry_max_event_name[cherry])(x);
      int ex=(*cherry_max_event_x[cherry])(x);
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  nx=(cherry_max_hidden_events[cherry][ex]).back();
	  if (e_type==2)
	    {
	      int at_branch=x_down(ex);
	      closer<<">@"<<at_branch<<"|";
	      (*event_stream)<<"<@"<<at_branch<<"|"<<"S@"<<at_branch<<"|";
	    }
	  else
	    {
	      int at_branch=x_down(nx);
	      closer<<"T>@"<<at_branch<<"|";
	      (*event_stream)<<"<T@"<<at_branch<<"|"<<"T@"<<at_branch<<"|";
	    }	  
	  //cherry_max_hidden_events[cherry][nx].pop_back();
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  ex=(*cherry_max_event_x[cherry])(nx);	  	  	 
	}
      //####################################################################################################

      stringstream out;
      out<<lin_slice[0]-lin_slice[ex]+1<< "|" << x_down(ex);
      int at_x=ex;
      int at_branch=x_down(ex);
      
      string event_name;
      if (e_type==1)
	event_name="S";
      else if (e_type==3)
	event_name="D";
      else if (e_type==4)
	event_name="T";
      else
	event_name="###########";
      if (event_name=="S" )
	{
	  closer<<">@"<<at_branch<<"|";
	  (*event_stream)<<"<@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";
	}
      else if (event_name=="T")
	{
	  if(x_down((*cherry_max_event_xc1[cherry])(ex))==at_branch)
	    at_branch=x_down((*cherry_max_event_xc2[cherry])(ex));
	  else
	    at_branch=x_down((*cherry_max_event_xc1[cherry])(ex));	 
	  closer<<"T>@"<<at_branch<<"|";
	  (*event_stream)<<"<T@"<<at_branch<<"|"<<event_name<<"@"<<at_branch<<"|";	  
	}
      else
	(*event_stream)<<event_name<<"@"<<at_branch<<"|";

      string node_name=out.str();

      tracestream(dnode_c1[dnode].first,(*cherry_max_event_xc1[cherry])(ex));
      tracestream(dnode_c2[dnode].first,(*cherry_max_event_xc2[cherry])(ex));
      (*event_stream)<<closer.str();
    }
}


void Species_tree::tracestream(Node * leaf,int x,int leaf_x)
{
  // cherry switch solution
  int other_leaf_x=leaf_x;
  if (leaf_x==-1)
    leaf_x=ur_sigma[leaf->getName()];
  string closer="";
  int nx=x;
  while (nx!=leaf_x)
    {
      int e_type=leaf_hidden_event_name[leaf_x][nx];
      stringstream tmp_closer;

      if (e_type==2)
	{
	  int at_branch=x_down(nx);
	  tmp_closer<<">@"<<at_branch<<"|";
	  (*event_stream)<<"<@"<<at_branch<<"|"<<"S@"<<at_branch<<"|";

	}
      else if (e_type==5 || x_down(nx)!=x_down(leaf_x))
	{
	  int at_branch=x_down(leaf_hidden_event_x[leaf_x][nx]);
	  tmp_closer<<"T>@"<<at_branch<<"|";
	  (*event_stream)<<"<T@"<<at_branch<<"|"<<"T@"<<at_branch<<"|";
	}
      nx=leaf_hidden_event_x[leaf_x][nx];	  
      if (e_type==2||e_type==5)
	{
	  //(*event_stream)<<closer.str();
	  closer=tmp_closer.str()+closer;
	}
    }
  (*event_stream)<<"O@"<<x_down(ur_sigma[leaf->getName()])<<"|";      

  (*event_stream)<<closer;
}


void Species_tree::tracecount(pair <Node*,Node*> root_dnode, tree_type * tree,int in_o_x)
{  
  //we record originations here
  // O
  O_counts[x_down(in_o_x)]+=tmp_pea;
  branch_zero[x_down(in_o_x)]+=tmp_pea;
  branch_count[x_down(in_o_x)]+=tmp_pea;
  

  pair <Node*,Node*> dnode = root_dnodes[root_dnode];
  int e_type=0;
  if (cherries.count(root_dnode)==0)
    {
      int o_x = in_o_x;//root_max_x[root_dnode];
      e_type=(*max_event_name[root_dnode])(o_x);
      int nx;
      //####################################################################################################
      while (e_type==2 || e_type==5 )
	{
	  while  (max_hidden_events[root_dnode].count(o_x)==0)
	    o_x=(*max_event_x[root_dnode])(o_x);	  	  
	  nx=(max_hidden_events[root_dnode][o_x]).back();

	  record_event(e_type, o_x,nx,-1,0);	  
	  e_type=(*max_event_name[root_dnode])(nx);
	  o_x=(*max_event_x[root_dnode])(nx);	  	  
	}
      //####################################################################################################

      int at_x=o_x;
      record_event(e_type, o_x,(*max_event_xc1[root_dnode])(o_x),(*max_event_xc2[root_dnode])(o_x),0);      
      tracecount(dnode_c1[root_dnode],(*max_event_xc1[root_dnode])(o_x));
      if (dnode_c1[root_dnode].first->isLeaf())
	tracecount(dnode_c1[root_dnode].first,(*max_event_xc1[root_dnode])(o_x));
      tracecount(dnode_c2[root_dnode],(*max_event_xc2[root_dnode])(o_x));    
      if (dnode_c2[root_dnode].first->isLeaf())
	tracecount(dnode_c2[root_dnode].first,(*max_event_xc2[root_dnode])(o_x));  
    }
  else
    {
      pair<int,int> cherry=cherries[root_dnode];
      //int o_x = cherry_root_max_x[cherry];
      int o_x = in_o_x;//root_max_x[root_dnode];

      e_type=(*cherry_max_event_name[cherry])(o_x);
      int nx;
      //####################################################################################################
      while (e_type==2 || e_type==5 )
	{
	  while  (cherry_max_hidden_events[cherry].count(o_x)==0)
	    o_x=(*cherry_max_event_x[cherry])(o_x);	  	  

	  nx=cherry_max_hidden_events[cherry][o_x].back();
	  record_event(e_type,o_x,nx,-1,0);
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  o_x=(*cherry_max_event_x[cherry])(nx);	  	  
	}
      //####################################################################################################
      int at_x=o_x;

      record_event(e_type,o_x,(*cherry_max_event_xc1[cherry])(o_x),(*cherry_max_event_xc2[cherry])(o_x),0);
      tracecount(dnode_c1[root_dnode].first,(*cherry_max_event_xc1[cherry])(o_x));
      tracecount(dnode_c2[root_dnode].first,(*cherry_max_event_xc2[cherry])(o_x));

    }
} 
void Species_tree::tracecount(pair <Node*,Node*> dnode,int x)
{
  
  if (cherries.count(dnode)==0)
    {
      if (max_event_name.count(dnode)==0)
	return;
      int e_type=(*max_event_name[dnode])(x);
      int ex=(*max_event_x[dnode])(x);
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  
	  nx=(max_hidden_events[dnode][ex]).back();
	  record_event(e_type,ex,nx);
	  e_type=(*max_event_name[dnode])(nx);
	  ex=(*max_event_x[dnode])(nx);	  	  
	}
      //####################################################################################################


      int at_x=ex;
      record_event(e_type,ex,(*max_event_xc1[dnode])(ex),(*max_event_xc2[dnode])(ex));
      tracecount(dnode_c1[dnode],(*max_event_xc1[dnode])(x));
      if (dnode_c1[dnode].first->isLeaf())
	tracecount(dnode_c1[dnode].first,(*max_event_xc1[dnode])(ex));
      tracecount(dnode_c2[dnode],(*max_event_xc2[dnode])(x));
      if (dnode_c2[dnode].first->isLeaf())
	tracecount(dnode_c2[dnode].first,(*max_event_xc2[dnode])(ex));
	
    }
  else 
    {

      pair<int,int> cherry=cherries[dnode];
      if (cherry_max_event_name.count(cherry)==0)
	return;

      int e_type=(*cherry_max_event_name[cherry])(x);
      int ex=(*cherry_max_event_x[cherry])(x);
      int nx=x;
      //####################################################################################################
      while (e_type==2 || e_type==5)
	{
	  nx=cherry_max_hidden_events[cherry][ex].back();
	  record_event(e_type,ex,nx);
	  e_type=(*cherry_max_event_name[cherry])(nx);
	  ex=(*cherry_max_event_x[cherry])(nx);	  	  
	}
      //####################################################################################################

      int at_x=ex;
      record_event(e_type,ex,(*cherry_max_event_xc1[cherry])(ex),(*cherry_max_event_xc2[cherry])(ex));      
      tracecount(dnode_c1[dnode].first,(*cherry_max_event_xc1[cherry])(ex));
      tracecount(dnode_c2[dnode].first,(*cherry_max_event_xc2[cherry])(ex));
    }  

}


void Species_tree::tracecount(Node * leaf,int x,int leaf_x)
{
  // cherry switch solution
  int other_leaf_x=leaf_x;
  if (leaf_x==-1)
    leaf_x=ur_sigma[leaf->getName()];
  
    
  /*
  int nx=x;
  while (nx!=leaf_x)
    {
      int e_type=leaf_hidden_event_name[leaf_x][nx];
      string hidden_name="";

      if (e_type==2)
	{
	  stringstream hout;
	  hout<<"&SL@";
	  hout<<lin_slice[0]-lin_slice[nx]+1<< "|" << x_down(nx) <<"|"<< nx;
	  hidden_name=hout.str()+hidden_name;	  
	}
      else if (e_type==5 || x_down(nx)!=x_down(leaf_x))
	{
	  stringstream hout;
	  hout<<"&TL@";
	  hout<<lin_slice[0]-lin_slice[nx]+1<< "|" << x_down(nx) << "|"<< x_down(leaf_hidden_event_x[leaf_x][nx]) <<"|"<< nx;
	  hidden_name=hout.str()+hidden_name;	  
	}
      nx=leaf_hidden_event_x[leaf_x][nx];
      leaf->setName(leaf->getName()+hidden_name);
    }
  */

  int nx=x;  
  while (nx!=leaf_x)
    {
      int e_type=leaf_hidden_event_name[leaf_x][nx];
      if (e_type==2 )
	{
	  record_event(e_type,nx,nx);
	}
      else if (e_type==5 || x_down(nx)!=x_down(leaf_x))
	{
	  record_event(5,nx,leaf_hidden_event_x[leaf_x][nx]);
	}
      nx=leaf_hidden_event_x[leaf_x][nx];

    }
  int at_x=x_down(ur_sigma[leaf->getName()]);
  // L
  branch_sum[at_x]+=tmp_pea;
  branch_zero[at_x]-=tmp_pea;
  branch_genome_size[at_x]+=tmp_pea;
}


void Species_tree::record_event(int e_type, int x,int xc1,int xc2,int o)
{
  int at_x,from_x;
  if (e_type==1 || e_type==2)
    {//S & SL
      at_x=x_down(x);
	  
      int b0=below_branch0[at_x];
      int b1=below_branch1[at_x];
      branch_sum[at_x]+=tmp_pea;
      branch_zero[at_x]-=tmp_pea;
      branch_genome_size[at_x]+=tmp_pea;

      branch_zero[b0]+=tmp_pea;
      branch_zero[b1]+=tmp_pea;
      branch_count[b0]+=tmp_pea;
      branch_count[b1]+=tmp_pea;
    }
  else if (e_type==3)
    {//D
      at_x=x_down(xc1);
      D_counts[at_x]+=tmp_pea;
      branch_zero[at_x]+=tmp_pea;
    }

  else if (e_type==4)
    {//T
      if (xc1==x)
	{
	  at_x=x_down(xc2);
	  from_x=x_down(xc1);
	  T_counts[at_x]+=tmp_pea;
	  Ttf[from_x][at_x]+=tmp_pea;
	  branch_zero[at_x]+=tmp_pea;
	  branch_sum[at_x]-=tmp_pea;
	}
      else  
	{
	  at_x=x_down(xc1);
	  from_x=x_down(xc2);
	  T_counts[at_x]+=tmp_pea;
	  Ttf[from_x][at_x]+=tmp_pea;
	  branch_zero[at_x]+=tmp_pea;
	  branch_sum[at_x]-=tmp_pea;
	}
    }
  else if (e_type==5)
    {//TL      
      at_x=x_down(xc1);
      from_x=x_down(x);
      T_counts[at_x]+=tmp_pea;	  
      Ttf[from_x][at_x]+=tmp_pea;
      branch_zero[at_x]+=tmp_pea;
      branch_sum[at_x]-=tmp_pea;
    }

}

