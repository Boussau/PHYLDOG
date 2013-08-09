#include "DTL.h"
#include <float.h>
using namespace std;
using namespace bpp;
scalar_type scalar_one=1;

void Species_tree::count_labels(vector<string> trees)
{
  O_counts.clear();  D_counts.clear();  T_counts.clear(); from_count.clear();
  branch_sum.clear();  branch_count.clear();  branch_zero.clear();  branch_one.clear();  branch_two.clear();
  mem_branch_count.clear();  mem_branch_fam_count.clear(),from_count.clear();branch_genome_size.clear();
  //sim_Ds_count.clear();sim_Ts_count.clear();
  rec_Ds_count.clear();rec_Ts_count.clear();

  for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
    {
      from_count.push_back(0); O_counts.push_back(0); D_counts.push_back(0); T_counts.push_back(0);
      branch_sum.push_back(0); branch_count.push_back(0); branch_zero.push_back(0); branch_one.push_back(0); branch_two.push_back(0);branch_genome_size.push_back(0);
      mem_branch_count[i]=0;
    }
  scalar_type tmp_pea = 1;
  scalar_type Opg=0.,Dpg=0.,Tpg=0.,Lpg=0.;
  

  vector <scalar_type> mo,md,mt,ml;
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
    {mo.push_back(0.);md.push_back(0.);mt.push_back(0.);ml.push_back(0.);}
  for (vector<string>::iterator si=trees.begin();si!=trees.end();si++)
    {
      //cout << (*si) <<endl;
      for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
	{mo[i]=O_counts[i];md[i]=D_counts[i];mt[i]=T_counts[i];ml[i]=branch_zero[i];}
      vector <Node*> nodes=TreeTemplateTools::parenthesisToTree(*si,false,"ID")->getNodes();
      for (vector<Node*>::iterator ni=nodes.begin();ni!=nodes.end();ni++)
	{
	  string name;
	  if  ((*ni)->isLeaf())
	    name = (*ni)->getName();
	  else
	    name = (* (dynamic_cast<const BppString *>((*ni)->getBranchProperty("ID")))).toSTL();

	  vector<string> name_tokens;
	  Tokenize(name,name_tokens,"&");
	  for (vector<string>::iterator ti=name_tokens.begin();ti!=name_tokens.end();ti++)
	    if ((*ti).find("@")!=(*ti).npos)
	      {
		vector<string> event_tokens;
		Tokenize((*ti),event_tokens,"@");
		string event=event_tokens[0];
		vector<string> time_tokens;
		Tokenize(event_tokens[1],time_tokens,"|");
		int at_x = atoi(time_tokens[1].c_str());
		int from_x = at_x;
		if (event=="T")
		  from_x=atoi(time_tokens[2].c_str());
		//TO FROM
		if (!(*ni)->hasFather())
		  if (event=="T")
		    {
		      O_counts[from_x]+=tmp_pea;
		      branch_zero[from_x]+=1;
		    }
		  else
		    {
		      O_counts[at_x]+=tmp_pea;
		      branch_zero[at_x]+=1;
		    }
		if (event=="S" || event=="SL")
		  {
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
		else if (event=="D")
		  {
		    D_counts[at_x]+=tmp_pea;
		    branch_zero[at_x]+=tmp_pea;
		    }
		  else if (event=="T" || event=="TL")
		    {
		      T_counts[at_x]+=tmp_pea;
		      branch_zero[at_x]+=tmp_pea;
		      branch_sum[at_x]-=tmp_pea;		     
		    }
		}
	      else
		{
		  int at_x = x_down(ur_sigma[(*ti)]);		  
 		  branch_sum[at_x]+=tmp_pea;
		  branch_zero[at_x]-=tmp_pea;
		  branch_genome_size[at_x]+=tmp_pea;
		}	  	  
	}
      //for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
	//{Opg+=}
	//cout << i << " " << O_counts[i]-mo[i] << " " << D_counts[i]-md[i] << " " << T_counts[i]-mt[i] << " " << branch_zero[i] - ml[i]<<" "<< branch_genome_size[i] <<endl;
    }
}


scalar_type Species_tree::count_events(vector<string> events,bool homogenous,bool subT )
{
  scalar_type pea;
  vector <int> gene;  vector <int> state;

  O_counts.clear();  D_counts.clear();  T_counts.clear(); from_count.clear();
  branch_sum.clear();  branch_count.clear();  branch_zero.clear();  branch_one.clear();  branch_two.clear();
  mem_branch_count.clear();  mem_branch_fam_count.clear(),from_count.clear();gene.clear();state.clear();

  //sim_Ds_count.clear();sim_Ts_count.clear();
  rec_Ds_count.clear();rec_Ts_count.clear();


  for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
    {
      gene.push_back(0); state.push_back(0); 
      from_count.push_back(0); O_counts.push_back(0); D_counts.push_back(0); T_counts.push_back(0);
      branch_sum.push_back(0); branch_count.push_back(0); branch_zero.push_back(0); branch_one.push_back(0); branch_two.push_back(0);
      mem_branch_count[i]=0;
    }

  scalar_type pea_count=0; scalar_type tree_count=0.;
  for (vector<string>::iterator ei=events.begin();ei!=events.end();ei++)
    {
      if (peas.size()>0) pea=peas[tree_count]; else pea=1.;
      pea=1;//XX
      pea_count+=pea;
      map <int,map<int,scalar_type> > tmp_from_count;
      tree_count++;
      for (int s1=0;s1<N_slices;s1++)
	for (int s2=0;s2<N_slices;s2++)
	  tmp_from_count[s1][s2]=0;
      vector <string> tokens;      
      Tokenize((*ei),tokens,"@|");
     for (int i = 0;i<(tokens.size())/2;i++)
	{
	  string event= tokens[2*i];
	  int at_branch= atoi(tokens[2*i+1].c_str());
	  if (!i) O_counts[at_branch]+=1.*pea;
	  int b0=below_branch0[at_branch];
	  int b1=below_branch1[at_branch];
	  if (event=="<" )	   
	    {
	      if (state[b0]==0) {state[b0]=1;gene[b0]=0;}
	      if (state[b1]==0) {state[b1]=1;gene[b1]=0;}
	    }
	  else if (event==">")
	    {
	      if (state[b0]==1)
		{
		  state[b0]=0;
		  branch_sum[b0]+=gene[b0]*pea;mem_branch_count[b0]+=1*pea;//XX		  
		  for (map<int,int>::iterator st=branch_2_slices[b0].begin();st!=branch_2_slices[b0].end();st++)		    
		    for (map<int,int>::iterator it=slice_2_branches[(*st).first].begin();it!=slice_2_branches[(*st).first].end();it++)
		      if ((*it).first!=b0)
			tmp_from_count[(*it).first][(*st).first]=max(scalar_one,tmp_from_count[(*it).first][(*st).first]);
		  if (gene[b0]==0 ) branch_zero[b0]+=1*pea;
		  if (gene[b0]==1 ) branch_one[b0]+=1*pea;
		  if (gene[b0]==2 ) branch_two[b0]+=1*pea;
		  gene[b0]=0;
		}
	      if (state[b1]==1)
		{
		  state[b1]=0; branch_sum[b1]+=gene[b1]*pea; mem_branch_count[b1]+=1*pea;//XX
		  for (map<int,int>::iterator st=branch_2_slices[b1].begin();st!=branch_2_slices[b1].end();st++)		    
		    for (map<int,int>::iterator it=slice_2_branches[(*st).first].begin();it!=slice_2_branches[(*st).first].end();it++)
		      if ((*it).first!=b1)
			tmp_from_count[(*it).first][(*st).first]=max(scalar_one,tmp_from_count[(*it).first][(*st).first]);
		  if (gene[b1]==0 ) branch_zero[b1]+=1*pea;
		  if (gene[b1]==1 ) branch_one[b1]+=1*pea;
		  if (gene[b1]==2 ) branch_two[b1]+=1*pea;
		  gene[b1]=0;
		}
	    }
	  else if (event=="<T") {if (state[at_branch]==1) state[at_branch]=2;}
	  else if (event=="T>") {if (state[at_branch]==2) state[at_branch]=1;}
	  else if (event=="S" || event=="O")
	    {
	      branch_count[at_branch]+=pea;	      
	      if (state[at_branch]==1) gene[at_branch]+=1;
	      if ( event=="O")		
		for (map<int,int>::iterator st=branch_2_slices[at_branch].begin();st!=branch_2_slices[at_branch].end();st++)		    
		  for (map<int,int>::iterator it=slice_2_branches[(*st).first].begin();it!=slice_2_branches[(*st).first].end();it++)
		    if ((*it).first!=at_branch)
		      tmp_from_count[(*it).first][(*st).first]=max(scalar_one,tmp_from_count[(*it).first][(*st).first]);
	    }	
	  else if (event=="T") T_counts[at_branch]+=1*pea;
	  else if (event=="D") D_counts[at_branch]+=1*pea;

	  if (event=="T") if (in_trees.size()==0) sim_Ts_count[tree_count-1]++; else rec_Ts_count[tree_count-1]++;
	  if (event=="D") if (in_trees.size()==0) sim_Ds_count[tree_count-1]++; else rec_Ds_count[tree_count-1]++;
	    
	}
      for (int i = 0;i<N_slices;i++)
	for (map<int,int>::iterator jt=slice_2_branches[i].begin();jt!=slice_2_branches[i].end();jt++)
	  from_count[(*jt).first]+=tmp_from_count[(*jt).first][i]*pea/(scalar_type)branch_2_slices[(*jt).first].size();     
	  
      
      for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
	{

	  if (mem_branch_count[i]>0)
	    mem_branch_fam_count[i]+=1;
	}
    }

  mem_D_count.clear(); mem_T_count.clear(); mem_L_count.clear();mem_O_count.clear();
  //if (!MPI)
  //  apparent_rates(pea_count,homogenous,subT);
  return pea_count;
}
void Species_tree::apparent_rates(scalar_type pea_count,bool homogenous,bool subT )
{ 
  remember_rates();
  vector <scalar_type> broadcast_delta,broadcast_tau,broadcast_lambda;
  int server = 0;  
  int rank = world.rank();
  int size = world.size(); 

  if (rank==server)
    {

  map <int,int> lb;
  //######## subtree loss / convergence correction ###############
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      mem_branch_genome_size[i]= branch_genome_size[i];
      mem_D_count[i]=D_counts[i];
      mem_T_count[i]=T_counts[i];
      mem_L_count[i]=branch_zero[i];
      mem_O_count[i]=O_counts[i];
    }   


  for (int i = 2;subT&&(i<N_slices+tree->getNumberOfLeaves()-1);i++)
    {
      scalar_type tmp_branch_one=branch_sum[i]-D_counts[i];
      int x=branch_2_x[i];             
      scalar_type dt=branch_length[i];
      scalar_type t=branch_length[i];
      scalar_type pL=1-exp(-namewise_lambda[branchname[i]]/t*dt);
      scalar_type pT=1-exp(-namewise_tau[branchname[i]]/t*dt);
      scalar_type Qei=Qe[lin_slice[x]][0][lin_edge[x]]*(1-pL);

      int b0=below_branch0[i];
      int b1=below_branch1[i];
      scalar_type lost_below;
      //      if(in_trees.size()>1 )
      //{
      lost_below = min(branch_count[i]*Qei*(1-pL), branch_zero[i]-1);
      branch_zero[i]+=tmp_branch_one*pL*pT;	  
      branch_sum[i]-=tmp_branch_one*pL*pT;	  
      T_counts[i]+=tmp_branch_one*pL*pT;
	  //}
	  //else
	  //lost_below = min(branch_count[i]*Qei*(1-pL), branch_zero[i]-1);      
      D_counts[i]/=(1-Qei); T_counts[i]/=(1-Qei);
      branch_zero[i]-=lost_below;branch_sum[i]+=lost_below;tmp_branch_one+=lost_below;
      lb[b0]+=lost_below;lb[b1]+=lost_below;
      branch_count[b0]+=lost_below;branch_zero[b0]+=lost_below;     
      branch_count[b1]+=lost_below;branch_zero[b1]+=lost_below;
    }


  int r0=below_branch0[1];  int r1=below_branch1[1];  scalar_type D_sum=0.;  scalar_type T_sum=0.;  scalar_type L_sum=0.;  scalar_type delta_sum=0,delta_ssum=0;
  scalar_type tau_sum=0,tau_ssum=0;  scalar_type lambda_sum=0,lambda_ssum=0;  scalar_type sum_count=0;
  map <int,int> avg_rate;
  avg_rate[0]=1;  avg_rate[1]=1;
  
  for (int i = 2;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    if (avg_rate.count(i)==0)
      {
	from_count[i]=pea_count;
	int x=branch_2_x[i];
	scalar_type p= branch_zero[i]/branch_count[i];     
	scalar_type q= 1-(1-p)/(branch_sum[i]/branch_count[i]);
	if (p==q)
	  p-=1e-4;
	scalar_type t=branch_length[i];
	scalar_type lambda=branch_zero[i]/max(branch_count[i],scalar_one)/t;
	lambda=((p*log((-1 + q)/(-1 + p)))/((p - q)*t));
	if (lambda!=lambda)
	  lambda=branch_zero[i]/max(branch_count[i],scalar_one )/t;
	  
	//lambda = branch_zero[i]/max(branch_sum[i],1. )/t;
	scalar_type delta=D_counts[i]/max(branch_sum[i],scalar_one)/t;
	delta= (q*log((-1 + q)/(-1 + p)))/((p - q)*t);	
	if (delta!=delta)
	  delta=D_counts[i]/max(branch_sum[i],scalar_one)/t;

	//if (avg_rate.count(i)==0)
	//T_counts[i]+=O_counts[i];
	scalar_type tau=T_counts[i]/max(from_count[i],scalar_one)/t;
	//if (tau>0.05 && lambda > 1e-4) 
	tau=(T_counts[i]/max(from_count[i],scalar_one))*lambda/(1.-exp(-lambda*t) );	
	if (tau!=tau)
	   tau=T_counts[i]/max(from_count[i],scalar_one)/t;

	scalar_type new_delta=branchwise_delta[i]*w+(1-w)*min(max(delta,1e-5*scalar_one),scalar_one);
	scalar_type new_tau=branchwise_tau[i]*w+(1-w)*min(max(tau,1e-5*scalar_one),scalar_one);
	scalar_type new_lambda=branchwise_lambda[i]*w+(1-w)*min(max(lambda,1e-5*scalar_one),3*scalar_one);
	if (avg_rate.count(i)==0)
	{
	  sum_count+=t;
	  delta_sum+=new_delta*t;
	  tau_sum+=new_tau*t;
	  lambda_sum+=new_lambda*t;
	  delta_ssum+=(new_delta)*(new_delta)*t;
	  tau_ssum+=(new_tau)*(new_tau)*t;
	  lambda_ssum+=(new_lambda)*(new_lambda)*t;
	}
      branchwise_delta[i]=new_delta;
      branchwise_tau[i]=new_tau;
      branchwise_lambda[i]=new_lambda;
      //cout <<i << " "<< branch_zero[i] << " " << branch_count[i] <<" "<< branch_sum[i] << " " << from_count[i] << " : ";
      //cout <<D_counts[i]<<" "<<T_counts[i] << " " << endl;
      
      }
  delta_avg=delta_sum/sum_count;
  tau_avg=tau_sum/sum_count;
  lambda_avg=lambda_sum/sum_count;
  scalar_type delta_var=abs(delta_avg*delta_avg-delta_ssum/sum_count);
  scalar_type tau_var=abs(tau_avg*tau_avg-tau_ssum/sum_count);
  scalar_type lambda_var=abs(lambda_avg*lambda_avg-lambda_ssum/sum_count);
  scalar_type delta_k=delta_avg*delta_avg/delta_var;
  scalar_type delta_Omega=delta_var/delta_avg;
  scalar_type tau_k=tau_avg*tau_avg/tau_var;
  scalar_type tau_Omega=tau_var/tau_avg;
  scalar_type lambda_k=lambda_avg*lambda_avg/lambda_var;
  scalar_type lambda_Omega=lambda_var/lambda_avg;
  if (world.rank()==0)    
    {
      cout << delta_avg << ":D:" << sqrt(delta_var) << " ";
      cout << tau_avg << ":T:" << sqrt(tau_var) << " ";
      cout << lambda_avg << ":L:" << sqrt(lambda_var) << " ";
      cout <<endl;
    }
  if (nclasses<1)
    nclasses=5;
  if (rank==server)
    {cout << "#################NCLASSES: " << nclasses <<endl;}
  vector <scalar_type> delta_limits;
  vector <scalar_type> delta_xs;
  vector <scalar_type> tau_limits;
  vector <scalar_type> tau_xs;
  vector <scalar_type> lambda_limits;
  vector <scalar_type> lambda_xs;


  for (scalar_type d=1;d<=nclasses;d++)
    {
      scalar_type delta_limit,delta_x,tau_limit,tau_x,lambda_limit,lambda_x;
      scalar_type p_limit=d/nclasses;
      scalar_type p_value=d/nclasses-1./(2.*nclasses);
      
      if (p_limit<1)
	delta_limit=boost::math::gamma_p_inv( delta_k,p_limit)*delta_Omega;
      else
	delta_limit=1e20;
      delta_x=boost::math::gamma_p_inv( delta_k,p_value)*delta_Omega;

      if (p_limit<1)
	tau_limit=boost::math::gamma_p_inv( tau_k,p_limit)*tau_Omega;
      else
	tau_limit=1e20;
      tau_x=boost::math::gamma_p_inv( tau_k,p_value)*tau_Omega;

      if (p_limit<1)
	lambda_limit=boost::math::gamma_p_inv( lambda_k,p_limit)*lambda_Omega;
      else
	lambda_limit=1e20;
      lambda_x=boost::math::gamma_p_inv( lambda_k,p_value)*lambda_Omega;

      delta_limits.push_back(delta_limit);
      delta_xs.push_back(delta_x);
      tau_limits.push_back(tau_limit);
      tau_xs.push_back(tau_x);
      lambda_limits.push_back(lambda_limit);
      lambda_xs.push_back(lambda_x);
    }
  
  scalar_type dtltmp=delta_avg+tau_avg+lambda_avg;
  scalar_type DTL_sum=0.; 
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      D_sum+=D_counts[i];//-count_branch_D[i];
      T_sum+=T_counts[i];//-count_branch_T[i];
      L_sum+=branch_zero[i];//-count_branch_L[i];
      DTL_sum+=D_counts[i]+T_counts[i]+branch_zero[i];
    }
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      int interval=0;
      while (branchwise_delta[i]*t>delta_limits[interval])
	interval++;
      branchwise_delta[i]=delta_xs[interval];

      interval=0;
      while (branchwise_tau[i]*t>tau_limits[interval])
	interval++;
      branchwise_tau[i]=tau_xs[interval];

      interval=0;
      while (branchwise_lambda[i]*t>lambda_limits[interval])
	interval++;
      branchwise_lambda[i]=lambda_xs[interval];

    }
  scalar_type tmp=0;
  
 
     for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
       {
	 broadcast_delta.push_back(0);
	 broadcast_tau.push_back(0);
	 broadcast_lambda.push_back(0);
       }
     for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
       {
	 
	 string name=branchname[i];
	 if (avg_rate.count(i)!=0 )
	   {
	     branchwise_delta[i]=delta_avg;
	     branchwise_tau[i]=tau_avg;
	     branchwise_lambda[i]=lambda_avg;
	   }
	 else if (branchwise==0 || homogenous==1)
	   {
	     branchwise_delta[i]=branchwise_delta[i]*w+(1-w)*delta_avg;
	     branchwise_tau[i]=branchwise_tau[i]*w+(1-w)*tau_avg;
	     branchwise_lambda[i]=branchwise_lambda[i]*w+(1-w)*lambda_avg;
	   }      
	 else
	   {
	     branchwise_delta[i]=branchwise_delta[i];
	     branchwise_tau[i]=branchwise_tau[i];
	     branchwise_lambda[i]=branchwise_lambda[i];
	   }
	 scalar_type t=branch_length[i];
	 
	 namewise_delta[name]=branchwise_delta[i]*t;
	 namewise_tau[name]=branchwise_tau[i]*t;
	 namewise_lambda[name]=branchwise_lambda[i]*t;
	 namewise_omega[name]=branchwise_omega[i];
	 broadcast_delta[i]=branchwise_delta[i];
	 broadcast_tau[i]=branchwise_tau[i];
	 broadcast_lambda[i]=branchwise_lambda[i];
       }
     
 avg_rate.clear();
 delta_limits.clear();
 delta_xs.clear();
 tau_limits.clear();
 tau_xs.clear();
 lambda_limits.clear();
 lambda_xs.clear();

   }
 
 broadcast(world,broadcast_delta,server);
 broadcast(world,broadcast_tau,server);
 broadcast(world,broadcast_lambda,server);
 for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
   {
     branchwise_delta[i]=broadcast_delta[i];
     branchwise_tau[i]=broadcast_delta[i];
     branchwise_lambda[i]=broadcast_lambda[i];
     string name=branchname[i];
     scalar_type t=branch_length[i];
     namewise_delta[name]=branchwise_delta[i]*t;
     namewise_tau[name]=branchwise_tau[i]*t;
     namewise_lambda[name]=branchwise_lambda[i]*t;
     namewise_omega[name]=branchwise_omega[i];
   }
 branchwise=1;
 omega_avg=O_counts[1]/pea_count+O_counts[2]/pea_count;
 
}
void Species_tree::remember_rates()
{  
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name=branchname[i];
      namewise_omega[name]=branchwise_omega[i];
      mem_branchwise_delta[i]=branchwise_delta[i];
      mem_branchwise_tau[i]=branchwise_tau[i];
      mem_branchwise_lambda[i]=branchwise_lambda[i];
      mem_branchwise_omega[i]=branchwise_omega[i];
      mmem_D_count[i]=mem_D_count[i];
      mmem_T_count[i]=mem_T_count[i];
      mmem_L_count[i]=mem_L_count[i];
      mmem_O_count[i]=mem_O_count[i];
      mem_rec_O_profile[i]=rec_O_profile[i];
      mmem_branch_genome_size[i]=mem_branch_genome_size[i];
    }
}

void Species_tree::recall_rates()
{
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      branchwise_delta[i]=mem_branchwise_delta[i];
      branchwise_tau[i]=mem_branchwise_tau[i];
      branchwise_lambda[i]=mem_branchwise_lambda[i];
      branchwise_omega[i]=mem_branchwise_omega[i];
      //branch_count[i]=mem_branch_count[i];

      string name=branchname[i];
      scalar_type t=branch_length[i];
      namewise_delta[name]=branchwise_delta[i]*t;
      namewise_tau[name]=branchwise_tau[i]*t;
      namewise_lambda[name]=branchwise_lambda[i]*t;
      namewise_omega[name]=branchwise_omega[i];      

      mem_D_count[i]=mmem_D_count[i];
      mem_T_count[i]=mmem_T_count[i];
      mem_L_count[i]=mmem_L_count[i];
      mem_O_count[i]=mmem_O_count[i];
      rec_O_profile[i]=mem_rec_O_profile[i];
      mem_branch_genome_size[i]=mmem_branch_genome_size[i];

    }
}

void Species_tree::print_rates()
{
  delta_avg=0,tau_avg=0,lambda_avg=0;
  scalar_type delta_tavg=0,tau_tavg=0,lambda_tavg=0,tmpt=0,tmp=0;
  if (0){
  tree_type * rate_tree = S_tree->clone();
  node_type * rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  vector <node_type*> nodes=rate_tree->getNodes();

  rate_tree = S_tree->clone();
  rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  nodes=rate_tree->getNodes();
  node_type * it;
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name =branchname[i];
      for (vector <node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
	{
	  if ((*nt)->getName()==name)
	    {
	      it=(*nt);
	      if (it->hasFather())
		it->setDistanceToFather(mem_L_count[i]);
	      stringstream out;	      
	      out<<mem_branch_sum[i];	  
	      if (it->isLeaf())
		{
		  it->setName(it->getName()+"="+out.str());
		}
	      else
		it->setBranchProperty("ID",BppString(out.str())); 	      
	    }
	  else
	    it=rate_root;      
	}
    }
  rate_trees["lambda"]= TreeTemplateTools::treeToParenthesis(*rate_tree,false,"ID");
  cout << "L counts "<<endl;
  cout<<rate_trees["lambda"];

  rate_tree = S_tree->clone();
  rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  nodes=rate_tree->getNodes();
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name =branchname[i];
      for (vector <node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
	{
	  if ((*nt)->getName()==name)
	    {
	      it=(*nt);
	      if (it->hasFather())
		it->setDistanceToFather(mem_T_count[i]);
	      stringstream out;	      
	      out<<mem_branch_sum[i];	  
	      if (it->isLeaf())
		{
		  it->setName(it->getName()+"="+out.str());
		}
	      else
		it->setBranchProperty("ID",BppString(out.str())); 	      
	    }
	  else
	    it=rate_root;      
	}
    }
  rate_trees["tau"]= TreeTemplateTools::treeToParenthesis(*rate_tree,false,"ID");
  cout << "T counts "<<endl;
  cout<<rate_trees["tau"];

  rate_tree = S_tree->clone();
  rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  nodes=rate_tree->getNodes();
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name =branchname[i];
      for (vector <node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
	{
	  if ((*nt)->getName()==name)
	    {
	      it=(*nt);
	      if (it->hasFather())
		it->setDistanceToFather(mem_D_count[i]);
	      stringstream out;	      
	      out<<mem_branch_sum[i];	  
	      if (it->isLeaf())
		{
		  it->setName(it->getName()+"="+out.str());
		}
	      else
		it->setBranchProperty("ID",BppString(out.str())); 	      
	    }
	  else
	    it=rate_root;      
	}
    }
  rate_trees["delta"]= TreeTemplateTools::treeToParenthesis(*rate_tree,false,"ID");
  cout << "D counts "<<endl;
  cout<<rate_trees["delta"];

  rate_tree = S_tree->clone();
  rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  nodes=rate_tree->getNodes();
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name =branchname[i];
      for (vector <node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
	{
	  if ((*nt)->getName()==name)
	    {
	      it=(*nt);
	      if (it->hasFather())
		it->setDistanceToFather(mem_O_count[i]);
	      stringstream out;	      
	      out<<mem_branch_sum[i];	  
	      if (it->isLeaf())
		{
		  it->setName(it->getName()+"="+out.str());
		}
	      else
		it->setBranchProperty("ID",BppString(out.str())); 	      
	    }
	  else
	    it=rate_root;      
	}
    }
  rate_trees["omega"]= TreeTemplateTools::treeToParenthesis(*rate_tree,false,"ID");
  cout << "O counts "<<endl;
  cout<<rate_trees["omega"];


  rate_tree = S_tree->clone();
  rate_root = rate_tree->getRootNode();
  name_internal(rate_root);
  nodes=rate_tree->getNodes();
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      string name =branchname[i];
      for (vector <node_type*>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
	{
	  if ((*nt)->getName()==name)
	    {
	      it=(*nt);
	      if (it->hasFather())
		it->setDistanceToFather( mem_branch_genome_size[i]);
	      stringstream out;	      
	      out<<mem_branch_genome_size[i];	  
	      if (it->isLeaf())
		{
		  it->setName(it->getName()+"="+out.str());
		}
	      else
		it->setBranchProperty("ID",BppString(out.str())); 	      
	    }
	  else
	    it=rate_root;      
	}
    }
  rate_trees["genomesize"]= TreeTemplateTools::treeToParenthesis(*rate_tree,false,"ID");
  cout << "genome size  "<<endl;
  cout<<rate_trees["genomesize"];

  }
  scalar_type ds=0;  scalar_type ts=0;  scalar_type ls=0;
  scalar_type sds=0;  scalar_type sts=0;  scalar_type sls=0;
  scalar_type s_ds=0;  scalar_type s_ts=0;  scalar_type s_ls=0;
  scalar_type s_sds=0;  scalar_type s_sts=0;  scalar_type s_sls=0;

  scalar_type sbs=0;  scalar_type bs=0;  
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
    {
      scalar_type t=branch_length[i];
      delta_avg+=branchwise_delta[i];tau_avg+=branchwise_tau[i];lambda_avg+=branchwise_lambda[i],tmp+=1;
      delta_tavg+=branchwise_delta[i]*t;tau_tavg+=branchwise_tau[i]*t;lambda_tavg+=branchwise_lambda[i]*t;tmpt+=t;
      sbs+=sim_O_count[i];
      bs+=mem_O_count[i];            
      sds+= sim_D_count[i];ds+= mem_D_count[i];
      sts+= sim_T_count[i];ts+= mem_T_count[i];
      sls+= sim_L_count[i];ls+= mem_L_count[i];
      s_sds+= sim_D_count[i]*sim_D_count[i];s_ds+= mem_D_count[i]*mem_D_count[i];
      s_sts+= sim_T_count[i]*sim_T_count[i];s_ts+= mem_T_count[i]*mem_T_count[i];
      s_sls+= sim_L_count[i]*sim_L_count[i];s_ls+= mem_L_count[i]*mem_L_count[i];
      

    }
  for (int i = 0;i<(N_slices+tree->getNumberOfLeaves()-1);i++)
    {
      cout << i<< " ";
      //cout << "@" << " ";      
      scalar_type t=branch_length[i];
      //cout << (1.-TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*inverse_flat_node_map[branch_2_x[i]],*inverse_flat_node_map[branch_2_x[1]]))+t/2. << " " << t << " ";
      //cout << branchwise_delta[i] << " ";
      //cout << branchwise_tau[i] << " ";
      //cout << branchwise_lambda[i] << " ";
      cout << sim_branch_genome_size[i] << " ";
      cout << mem_branch_genome_size[i] << "    ";
      cout << sim_O_count[i] << " ";
      cout << O_counts[i] << " ";
      //cout << branchwise_omega[i] << " ";
      //cout << sim_O_profile[i] << " ";
      //cout << rec_O_profile[i] << "    ";
      cout << sim_D_count[i] << " ";sds+= sim_D_count[i];
      cout << mem_D_count[i] << "    ";ds+= mem_D_count[i];
      cout << sim_T_count[i] << " ";sts+= sim_T_count[i];
      cout << mem_T_count[i] << "    ";ts+= mem_T_count[i];
      cout << sim_L_count[i] << " ";sls+= sim_L_count[i];
      cout << mem_L_count[i] << "    ";ls+= mem_L_count[i];
      cout << branchname[i] << " ";
      cout << endl;
    }

  delta_avg/=tmp ; tau_avg/=tmp;lambda_avg/=tmp;omega_avg=branchwise_omega[1]+branchwise_omega[2];
  cout << "AVG: " << delta_avg << "\t" << tau_avg << "\t" << lambda_avg <<endl;
  cout << "RAH: " << (delta_avg+tau_avg)/(delta_avg+tau_avg+lambda_avg) << "\t" << (tau_avg)/(delta_avg+tau_avg) << "\t" << (delta_avg+tau_avg+lambda_avg) <<endl;

  cout << "tAVG: " << delta_tavg/tmpt << "\t" << tau_tavg/tmpt << "\t" << lambda_tavg/tmpt <<endl;
  cout << "sums:" << sds << ":" << ds << " " <<  sts << ":" << ts << " " <<  sls << ":" << ls << endl;
  cout << "ratios:" << sds/(sds+sts+sls) << ":" << ds/(ds+ts+ls) << " " <<  sts/(sds+sts+sls) << ":" << ts/(ds+ts+ls) << " " <<  sls/(sds+sts+sls) << ":" << ls/(ds+ts+ls) << endl;
  cout << "pertree:" << sds/(sbs) << ":" << ds/(bs) << " " <<  sts/(sbs) << ":" << ts/(bs) << " " <<  sls/(sbs) << ":" << ls/(bs) << endl;
  cout << "var:";
  cout << sqrt((s_sds/sbs)-((sds/(sbs))*(sds/(sbs)))) << ":" << sqrt((s_ds/bs)-((ds/(bs))*(ds/(bs)))) << " ";
  cout << sqrt((s_sts/sbs)-((sts/(sbs))*(sts/(sbs)))) << ":" << sqrt((s_ts/bs)-((ts/(bs))*(ts/(bs)))) << " "; 
  cout << sqrt((s_sls/sbs)-((sls/(sbs))*(sls/(sbs)))) << ":" << sqrt((s_ls/bs)-((ls/(bs))*(ls/(bs)))) << " "; 
  cout <<endl;


  cout << "Br " << (sim_O_count[1]+sim_O_count[2])/sbs << " " << (mem_O_count[1]+mem_O_count[2])/bs << endl;

  if (world.rank()==0)
    {
      ofstream cout_stream( outstream_file_name.c_str(),ios::app); 
      cout_stream << "AVG: " << delta_avg << "\t" << tau_avg << "\t" << lambda_avg <<endl;
      cout_stream << "RAH: " << (delta_avg+tau_avg)/(delta_avg+tau_avg+lambda_avg) << "\t" << (tau_avg)/(delta_avg+tau_avg) << "\t" << (delta_avg+tau_avg+lambda_avg) <<endl;
      cout_stream << "tAVG: " << delta_tavg/tmpt << "\t" << tau_tavg/tmpt << "\t" << lambda_tavg/tmpt <<endl;
      cout_stream << "sums:" << sds << ":" << ds << " " <<  sts << ":" << ts << " " <<  sls << ":" << ls << endl;
      cout_stream << "Br " << sim_O_count[1]/sbs << " " << mem_O_count[1]/bs << endl;      
      cout_stream.close();     
    }
}


scalar_type Species_tree::DTLclock_height(int i)
{
  if (i==0) i=1;
  if (inverse_flat_node_map[branch_2_x[i]]->isLeaf())
    return 0;
  else
    {
      int b0=below_branch0[i];
      int b1=below_branch1[i];     
      return max((D_counts[b0]+T_counts[b0]+branch_zero[b0]+DTLclock_height(b0)+D_counts[b1]+T_counts[b1]+branch_zero[b1]+DTLclock_height(b1))/2.,max(DTLclock_height(b0)+1e-4,DTLclock_height(b1))+1e-4);
    }
}  


