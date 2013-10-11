#define  isorisnot( slice, si, edge, N_edges) (si==0 && (edge==N_edges-1 || slice == 0)) || (si>0)

#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include "DTL.h"

using namespace std;
using namespace bpp;
namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

//*****************************************************
//************** PUBLIC ** PUBLIC ** PUBLIC ***********
//*****************************************************


void Species_tree::init_L_improved(scalar_type delta_in,scalar_type tau_in,scalar_type lambda_in,scalar_type stem_omega,scalar_type out_tau,string mode)
{
  cherries.clear();  
  advance_from.clear();
  advance_until.clear();
  below_from.clear();
  below_until.clear();
  flat_Qef.clear();
  p_flat_Qef.clear();
  delta_dt.clear();
  tau_dt.clear();
  flat_node_map.clear();
  lin_si.clear();
  delta_table.clear(); tau_table.clear(); lambda_table.clear();
  id_delta_table.clear(); id_tau_table.clear(); id_lambda_table.clear();
  id_slice_table.clear();
  id_table.clear();
  l_table.clear();
  N_table.clear();
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


  infer_mode=true;
  //delete tree;

  init();  
  
  real_root = root;  root = new Node();  root -> setName("*Rstem");
  real_root->addSon(root);  root->setDistanceToFather(stem_omega);
  tree->rootAt(root);  

  bool out_group= false;//( out_tau>0 );// && (mode!="max");
  if ( out_group )
    {
      Node * out_group = new Node();  
      root->addSon(out_group);
      out_group->setDistanceToFather(1.+stem_omega*1.);
      out_group -> setName("OUT_GROUP|");
    }
  construct_slices();

  delta=delta_in;  tau=tau_in;  lambda=lambda_in; N_slices =  time_slices.size();

  map <int, map <int, Node * > > node_map;
  for (int slice = 0; slice<N_slices;slice++)
    {
      //cout << "# "<<slice <<  " : " ; 
      for (int edge=0;edge<time_slices[N_slices-slice-1].size();edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  node_map[slice][edge] = node;
	  //cout << node->getName() << " " << node->getDistanceToFather() << " | ";
	}
      //cout << endl;
    }
 
  flat_node_map.clear();
  inverse_flat_node_map.clear();

  int slice_offset[N_slices];
  lin_slice.clear();
  lin_si.clear();
  lin_edge.clear();

  //  scalar_type edges_in_slice[N_slices];
  // resolution for time integration
  ddc=10;
  double ddt = 1./((scalar_type) ddc-1.); 
  //construct vectors for Qe integration
  //vector_type  Qe[N_slices][ddc]; 
  //rate tree coding

  //LL
  vector_type  deltas[N_slices],  taus[N_slices], lambdas[N_slices], phis[N_slices];

  //topology coding

  //the indicies in slice-1 of the descendants of S nodes (i.e. g -> g'',g'') 
  //in any given slice the speciation node always has index -1 (i.e edges_in_slice[slice]-1 )
  int S_below[N_slices][2];
  //the indicies in slice-1 of the descendants of Sp nodes (i.e. h) 
  map<int , int> Sp_below[N_slices];
  //the index of nodes father in slice+1 
  map<int , int> S_above[N_slices];
  // index of node in it's slice
  map<Node * , int> index_;

  for (int slice = 0; slice<N_slices;slice++)
    {
      edges_in_slice[slice]=time_slices[N_slices-slice-1].size();
      for (int edge=0;edge<edges_in_slice[slice];edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  index_[node] = edge;
	}
    }
  for (int slice = 0; slice<N_slices;slice++)
    {
      int tmp=0;
      for (int edge=0;edge<edges_in_slice[slice];edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  if(slice > 0)
	    {
	      if(node->getNumberOfSons()==2)
		{
		  S_below[slice][0] = index_[node->getSons()[0]];
		  S_below[slice][1] = index_[node->getSons()[1]];
		}	    
	      else
		{
		  Sp_below[slice][edge] = index_[node->getSons()[0]];
		}	    	         
	    }
	  if (slice<N_slices-1)
	    S_above[slice][edge] = index_[node->getFather()];
	}
    }

  //resolution for gene tree reconciliation
  lin_size=0; inverse_flat_node_map.clear();lin_slice.clear();  lin_edge.clear();
  dc=2;
  dc_mod = (ddc-1)/(dc-1);

  // construct index for flat iteration forward in time over
  // S'' nodes that are in S or are not contemporary to a node in S
  for (int from_slice = N_slices-1; from_slice>-1; from_slice--)
    {
      slice_offset[from_slice] = lin_size;
      int N_from_edges = edges_in_slice[from_slice];            
      for (int si = dc-1; si > -1; si--)
	for (int edge_from = 0; edge_from < N_from_edges; edge_from++)
	  // S'' nodes that are in S or are not contemporary to a node in S
	  if ( isorisnot(from_slice,si,edge_from,N_from_edges) )
	    {
	      flat_map[from_slice][si][edge_from]=lin_size;
	      flat_node_map[node_map[from_slice][edge_from]] = lin_size;
	      inverse_flat_node_map[lin_size] = node_map[from_slice][edge_from];
	      
	      lin_slice[lin_size] = from_slice;
	      lin_edge[lin_size] = edge_from;
	      lin_si[lin_size] = si; 
	      lin_size++;
	    }
    }
  //NEWLIN SIZE 


  // MEM: N^2 x ddc
  for (int slice = 0; slice<N_slices;slice++)
    {
      scalar_type C = edges_in_slice[slice]-1;
      deltas[slice].resize(edges_in_slice[slice]);
      taus[slice].resize(edges_in_slice[slice]);
      lambdas[slice].resize(edges_in_slice[slice]);
      phis[slice].resize(edges_in_slice[slice]);       
      for (int ci=0;ci<ddc;ci++)
	Qe[slice][ci].resize(edges_in_slice[slice]);
      for (int edge=0;edge<edges_in_slice[slice];edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  scalar_type d = node->getDistanceToFather(); 
	  //if (node==root)
	  //  d=0.;
	  deltas[slice](edge) = pow(d * delta,beta);
	  if (C>0)
	    taus[slice](edge) = pow(d * tau / C,beta);
	  else
	    taus[slice](edge) = 0;
	  if (node->getName()[node->getName().size()-1]=='|')
	    lambdas[slice](edge) = 1e-5;
	  else
	    lambdas[slice](edge) = pow(d * lambda,beta);
	  //phis[slice](edge) = d * (delta + tau + lambda) ;
	  phis[slice](edge) = deltas[slice](edge)+taus[slice](edge)*(C)+lambdas[slice](edge);
	  // <<<<<<-XX-V-XX->>>>>>>
	}
    }

  for (map<int,int>::iterator it=branch_2_x.begin();it!=branch_2_x.end();it++)
    if (lin_slice[branch_2_x[(*it).first]]>0)
      {

	
	below_branch0[(*it).first]=x_down(below_son0[branch_2_x[(*it).first]]);	
	below_branch1[(*it).first]=x_down(below_son1[branch_2_x[(*it).first]]);
	int x=below_son0[branch_2_x[(*it).first]];
	scalar_type d=0;
	while (type_of_x[x]!=1 && lin_slice[x]>0)
	  {	    
	    d+=inverse_flat_node_map[x]->getDistanceToFather();
	    x=below_son0[x];
	  }
	if (d==0)
	  d=inverse_flat_node_map[x]->getDistanceToFather();
	branch_length[x_down(below_son0[branch_2_x[(*it).first]])]=d;

	x=below_son1[branch_2_x[(*it).first]];
	d=0;
	while (type_of_x[x]!=1 && lin_slice[x]>0)
	  {	    
	    d+=inverse_flat_node_map[x]->getDistanceToFather();
	    x=below_son0[x];
	  }
	if (d==0)
	  d=inverse_flat_node_map[x]->getDistanceToFather();
	branch_length[x_down(below_son1[branch_2_x[(*it).first]])]=d;		
      }



  for (int x=0;x<lin_size;x++)
    {
      string name=branchname[x_down(x)];
      branchi[name]=x_down(x);
      scalar_type t=branch_length[x_down(x)];
      if (x_down(x)==1)
	t=stem_omega;      
      delta=namewise_delta[name]/t;
      tau=namewise_tau[name]/t;
      lambda=namewise_lambda[name]/t;
      branchwise_omega[x_down(x)]=namewise_omega[name];
      
      int edge=lin_edge[x];
      int slice=lin_slice[x];
      scalar_type C = edges_in_slice[slice]-1;
      Node * node=time_slices[N_slices-slice-1][edge];
      scalar_type d = node->getDistanceToFather();       

      deltas[slice](edge) = d*delta;
      if (C>0)
	taus[slice](edge) = d*tau/ C;
      else
	taus[slice](edge) = 0;
      if (node->getName()[node->getName().size()-1]=='|')
	lambdas[slice](edge) = 1e-5;
      else
	lambdas[slice](edge) = d*lambda;
      phis[slice](edge) = deltas[slice](edge)+taus[slice](edge)*(C)+lambdas[slice](edge);
      // <<<<<<-XX-V-XX->>>>>>>        
    }  
  //calculate Qe (integrating back in time)
  // PARAL: none, N2 x ddc
  //LL
  vector_type term1;  vector_type term2;  vector_type term3;  vector_type dQ;  
  vector_type k1; vector_type k2;  vector_type k3;  vector_type k4;  vector_type ktmp; 
  vector_type QefQe; vector_type Qsq;
  vector_type v_tmp; double s_tmp;

  for (int slice = 0; slice<N_slices;slice++)
    {
      if (slice==0)
	//init Qe at leaves
	for (int edge=0;edge<edges_in_slice[slice];edge++)
	  Qe[0][0](edge) = 0;
      else
	//init Qe from slice below
	{
	  //pass thru Sp nodes
	  for (int edge=0;edge<edges_in_slice[slice]-1;edge++)
	    Qe[slice][0](edge) = Qe[slice-1][ddc-1]( Sp_below[slice][edge] );
	  //and the single speciation
	  Qe[slice][0](edges_in_slice[slice]-1) 
	    = Qe[slice-1][ddc-1]( S_below[slice][0] ) * Qe[slice-1][ddc-1]( S_below[slice][1] );
	}

      term1.resize(edges_in_slice[slice]); term2.resize(edges_in_slice[slice]); term3.resize(edges_in_slice[slice]);
      dQ.resize(edges_in_slice[slice]);      
      k1.resize(edges_in_slice[slice]); k2.resize(edges_in_slice[slice]); k3.resize(edges_in_slice[slice]); k4.resize(edges_in_slice[slice]);
      ktmp.resize(edges_in_slice[slice]);      
      Qsq.resize(edges_in_slice[slice]);
      v_tmp.resize(edges_in_slice[slice]);
      // RK4: 4th order Runge-Kutta for y'=f(y)
      // using BLAS/UBLAS for "inner for loop"
      for (int ddi = 1; ddi<ddc; ddi++)
	{	    
	  //RK4 k1
	  // k1 = f(y[n]) 
	  k1*=0;ktmp*=0;
	  blas::axpy( 1. , Qe[slice][ddi-1] , ktmp );

	  term1*=0;term2*=0;term3*=0;	    
	  //  delta_e ( Q_e Q_e )  
	  Qsq = ublas::element_prod( ktmp , ktmp );
	  term1 = ublas::element_prod( deltas[slice] , Qsq );
	  //+ Sum_CEe Q_e (tau_f/C  Q_f)
	  // tau_f/C Q_f 
	  v_tmp = ublas::element_prod( ktmp , taus[slice] );
	  // Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  s_tmp = blas::asum( v_tmp );
	  // Q_e x Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  blas::axpy( s_tmp , ktmp , term2 );
	  // tau_e Q_e Q_e 
	  v_tmp = ublas::element_prod( Qsq , taus[slice] ); 
	  blas::axpy( -1. , v_tmp , term2);
	  //+ lambda_e - phi_e Q_e
	  term3 = - ublas::element_prod( ktmp , phis[slice] );
	  blas::axpy( 1. , lambdas[slice] , term3 );

	  blas::axpy( 1. , term1 , k1 );
	  blas::axpy( 1. , term2 , k1 );
	  blas::axpy( 1. , term3 , k1 );

	  //RK4 k2
	  // k2 = f(y[n]+h/2 k1)
	  k2*=0;ktmp*=0;
	  blas::axpy( 1. , Qe[slice][ddi-1] , ktmp );
	  blas::axpy( ddt/2. , k1 , ktmp );

	  term1*=0;term2*=0;term3*=0;	    
	  //  delta_e ( Q_e Q_e )  
	  Qsq = ublas::element_prod( ktmp , ktmp );
	  term1 = ublas::element_prod( deltas[slice] , Qsq );
	  //+ Sum_CEe Q_e (tau_f/C  Q_f)
	  // tau_f/C Q_f 
	  v_tmp = ublas::element_prod( ktmp , taus[slice] );
	  // Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  s_tmp = blas::asum( v_tmp );
	  // Q_e x Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  blas::axpy( s_tmp , ktmp , term2 );
	  // tau_e Q_e Q_e 
	  v_tmp = ublas::element_prod( Qsq , taus[slice] ); 
	  blas::axpy( -1. , v_tmp , term2);
	  //+ lambda_e - phi_e Q_e
	  term3 = - ublas::element_prod( ktmp , phis[slice] );
	  blas::axpy( 1. , lambdas[slice] , term3 );

	  blas::axpy( 1. , term1 , k2 );
	  blas::axpy( 1. , term2 , k2 );
	  blas::axpy( 1. , term3 , k2 );

	  //RK4 k3
	  // k3 = f(y[n]+h/2 k2)
	  k3*=0;ktmp*=0;
	  blas::axpy( 1. , Qe[slice][ddi-1] , ktmp );
	  blas::axpy( ddt/2. , k2 , ktmp );

	  term1*=0;term2*=0;term3*=0;	    
	  //  delta_e ( Q_e Q_e )  
	  Qsq = ublas::element_prod( ktmp , ktmp );
	  term1 = ublas::element_prod( deltas[slice] , Qsq );
	  //+ Sum_CEe Q_e (tau_f/C  Q_f)
	  // tau_f/C Q_f 
	  v_tmp = ublas::element_prod( ktmp , taus[slice] );
	  // Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  s_tmp = blas::asum( v_tmp );
	  blas::axpy( s_tmp , ktmp , term2 );
	  // tau_e Q_e Q_e 
	  v_tmp = ublas::element_prod( Qsq , taus[slice] );
	  // Q_e x Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e  
	  blas::axpy( -1. , v_tmp , term2);
	  //+ lambda_e - phi_e Q_e
	  term3 = - ublas::element_prod( ktmp , phis[slice] );
	  blas::axpy( 1. , lambdas[slice] , term3 );

	  blas::axpy( 1. , term1 , k3 );
	  blas::axpy( 1. , term2 , k3 );
	  blas::axpy( 1. , term3 , k3 );

	  //RK4 k4
	  // k4 = f(y[n]+h k3)
	  k4*=0;ktmp*=0;
	  blas::axpy( 1. , Qe[slice][ddi-1] , ktmp );
	  blas::axpy( ddt , k3 , ktmp );

	  term1*=0;term2*=0;term3*=0;	    
	  //  delta_e ( Q_e Q_e )  
	  Qsq = ublas::element_prod( ktmp , ktmp );
	  term1 = ublas::element_prod( deltas[slice] , Qsq );
	  //+ Sum_CEe Q_e (tau_f/C  Q_f)
	  // tau_f/C Q_f 
	  v_tmp = ublas::element_prod( ktmp , taus[slice] );
	  // Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e
	  s_tmp = blas::asum( v_tmp );
	  // Q_e x Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e 
	  blas::axpy( s_tmp , ktmp , term2 );
	  // tau_e Q_e Q_e 
	  v_tmp = ublas::element_prod( Qsq , taus[slice] ); 
	  blas::axpy( -1. , v_tmp , term2);
	  //+ lambda_e - phi_e Q_e
	  term3 = - ublas::element_prod( ktmp , phis[slice] );
	  blas::axpy( 1. , lambdas[slice] , term3 );

	  blas::axpy( 1. , term1 , k4 );
	  blas::axpy( 1. , term2 , k4 );
	  blas::axpy( 1. , term3 , k4 );

	  //RK4 step
	  // y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4) 
	  dQ*=0;
	  blas::axpy( 1. , k1 , dQ );
	  blas::axpy( 2. , k2 , dQ );
	  blas::axpy( 2. , k3 , dQ );
	  blas::axpy( 1. , k4 , dQ );
	  blas::axpy( 1. , Qe[slice][ddi-1] , Qe[slice][ddi] );
	  blas::axpy( ddt/6. , dQ , Qe[slice][ddi] );

	}	
    }
  //  cout << "# we have Qecalc" << endl;
  
  //construct vectors for Qef integration
  //LL

  vector_type deltas_s_t[N_slices]; vector_type taus_s_t[N_slices]; vector_type phis_s_t[N_slices]; vector_type Qe_s_t[N_slices][ddc];
  vector_type Qef_s_t[N_slices][ddc][ddc];
  
  for (int slice = 0; slice<N_slices;slice++)
    {
      int N_edges = edges_in_slice[slice];
      deltas_s_t[slice].resize(N_edges * N_edges);	    
      taus_s_t[slice].resize(N_edges * N_edges);	    
      phis_s_t[slice].resize(N_edges * N_edges);	    
      for (int si = 0; si<ddc; si++)
	{	    
	  Qe_s_t[slice][si].resize(N_edges * N_edges);	    
	  for (int ti = 0; ti<=si; ti++)
	    Qef_s_t[slice][ti][si].resize(N_edges * N_edges);
	}
      for (int edge_e = 0; edge_e<N_edges; edge_e++)
	for (int edge_f = 0; edge_f<N_edges; edge_f++)
	  {
	    deltas_s_t[slice](edge_e*N_edges + edge_f) = deltas[slice](edge_e);      
	    taus_s_t[slice](edge_e*N_edges + edge_f) = taus[slice](edge_e);      
	    phis_s_t[slice](edge_e*N_edges + edge_f) = phis[slice](edge_e);      
	    for (int si = 0; si<ddc; si++)	    	      
	      Qe_s_t[slice][si](edge_e*N_edges + edge_f) = Qe[slice][si](edge_e);      
	  }
    }

  //calculate Q_ef (integrating back in time)
  //PARAL: idependent time slices: N_slices x N^2 ddc^2 

  for (int slice = 0; slice<N_slices;slice++)
    {
      int N_edges = edges_in_slice[slice];
      //init Qef_es
      for (int edge_e = 0; edge_e<N_edges; edge_e++)
	for (int edge_f = 0; edge_f<N_edges; edge_f++)
	  for (int ti = 0; ti < ddc; ti++)
	    {
	      if (edge_e==edge_f)
		Qef_s_t[slice][ti][ti](edge_e*N_edges + edge_f) = 1.;
	      else
		Qef_s_t[slice][ti][ti](edge_e*N_edges + edge_f) = 0.;
	    }
      for(int ti = 0; ti<ddc; ti++ )	
	{
	  term1.resize(N_edges*N_edges);term2.resize(N_edges*N_edges);term3.resize(N_edges*N_edges);dQ.resize(N_edges*N_edges);	  
	  k1.resize(N_edges*N_edges);k2.resize(N_edges*N_edges);k3.resize(N_edges*N_edges);k4.resize(N_edges*N_edges);ktmp.resize(N_edges*N_edges);	  
	  QefQe.resize(N_edges*N_edges);
	  v_tmp.resize(N_edges*N_edges);

	  // RK4: 4th order Runge-Kutta for y'=f(y)
	  // using BLAS/UBLAS for "inner for loop"
	  for (int ddi = ti+1; ddi<ddc; ddi++)
	    {	    
	      //RK4 k1
	      // k1 = f(y[n]) 
	      k1*=0;ktmp*=0;
	      blas::axpy( 1. , Qef_s_t[slice][ti][ddi-1] , ktmp );

	      term1*=0;term2*=0;term3*=0;	    
	      //  delta_e ( Q_ef Q_e )  
	      QefQe = ublas::element_prod( ktmp , Qe_s_t[slice][ddi]);
	      term1 = ublas::element_prod( deltas_s_t[slice] , QefQe );

	      //+ Q_e Sum_CEg (tau_g/C  Q_gf) + Q_efSum_CEg tau_g/C Q_g

	      // Q_e Sum_CEg (tau_g/C  Q_gf)
	      // // (tau_g/C  Q_gf) 
	      v_tmp = ublas::element_prod( ktmp , taus_s_t[slice] );
	      // // Sum_CEg (tau_g/C  Q_gf) + (tau_e/C  Q_ef)
	      s_tmp = blas::asum( v_tmp );
	      // // Qe x Sum_CEg (tau_g/C  Q_gf) + (tau_e/C  Q_ef)
	      blas::axpy( s_tmp , Qe_s_t[slice][ddi] , term2 );
	      
	      // Q_ef Sum_CEg tau_g/C Q_g
	      // // (tau_g/C Q_g)
	      v_tmp = ublas::element_prod( Qe_s_t[slice][ddi] , taus_s_t[slice] );
	      // // Sum_CEe  Q_e (tau_f/C  Q_f) + tau_e Q_e
	      s_tmp = blas::asum( v_tmp );
	      // // Qef x Sum_CEg (tau_g/C  Q_g) + (tau_e/C  Q_e)
	      blas::axpy( s_tmp , Qef_s_t[slice][ti][ddi] , term2 );

	      // (tau_e/C  Q_efQ_e)
	      v_tmp = ublas::element_prod( QefQe , taus_s_t[slice] );
	      // Q_e Sum_CEg (tau_g/C  Q_gf) + Q_ef Sum_CEg tau_g/C Q_g
	      blas::axpy( -2. , v_tmp , term2 );
	      
	      //- phi_e Q_e
	      term3 = - ublas::element_prod( ktmp , phis_s_t[slice] );
	      
	      blas::axpy( 1. , term1 , k1 );
	      blas::axpy( 1. , term2 , k1 );
	      blas::axpy( 1. , term3 , k1 );

	      //RK4 k2
	      // k2 = f(y[n]+h/2 k1)
	      k2*=0;ktmp*=0;
	      blas::axpy( 1. , Qef_s_t[slice][ti][ddi-1] , ktmp );
	      blas::axpy( ddt/2. , k1 , ktmp );

	      blas::axpy( 1. , term1 , k2 );
	      blas::axpy( 1. , term2 , k2 );
	      blas::axpy( 1. , term3 , k2 );

	      //RK4 k3
	      // k3 = f(y[n]+h/2 k2)
	      k3*=0;ktmp*=0;
	      blas::axpy( 1. , Qef_s_t[slice][ti][ddi-1] , ktmp );
	      blas::axpy( ddt/2. , k2 , ktmp );

	      blas::axpy( 1. , term1 , k3 );
	      blas::axpy( 1. , term2 , k3 );
	      blas::axpy( 1. , term3 , k3 );

	      //RK4 k4
	      // k4 = f(y[n]+h k3)
	      k4*=0;ktmp*=0;
	      blas::axpy( 1. , Qef_s_t[slice][ti][ddi-1] , ktmp );
	      blas::axpy( ddt , k3 , ktmp );

	      blas::axpy( 1. , term1 , k4 );
	      blas::axpy( 1. , term2 , k4 );
	      blas::axpy( 1. , term3 , k4 );

	      //RK4 step
	      // y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4) 
	      dQ*=0;
	      blas::axpy( 1. , k1 , dQ );
	      blas::axpy( 2. , k2 , dQ );
	      blas::axpy( 2. , k3 , dQ );
	      blas::axpy( 1. , k4 , dQ );
	      blas::axpy( 1. , Qef_s_t[slice][ti][ddi-1] , Qef_s_t[slice][ti][ddi] );
	      blas::axpy( ddt/6. , dQ , Qef_s_t[slice][ti][ddi] );	  	  
	    }	  
	}      
    }
  scalar_type dt = 1./((scalar_type) dc -1.); 
  // construct Qef accros slices
  array_type Qef_hp[dc][N_slices],Qef_ep[dc][N_slices];
  for (int from_slice = 0; from_slice<N_slices; from_slice++)
    {

      int N_from_edges = edges_in_slice[from_slice];
      slice_2_slice[from_slice][from_slice].resize(N_from_edges,N_from_edges);

      for (int edge_e = 0; edge_e<N_from_edges; edge_e++)
	{
	  for (int edge_h = 0; edge_h<N_from_edges; edge_h++) 
	    {
	      // Qeh(-1,0)	     
	      slice_2_slice[from_slice][from_slice](edge_e,edge_h) = 
		Qef_s_t[from_slice][0][ddc-1](edge_e*N_from_edges + edge_h);	   
	      //XX	    
	    }
	}
    }  
  // calculate Qef accros neigbooring slices 
  // and construct Qeh(s,0) and Qh't(-1,t)   
  for (int from_slice = 1; from_slice<N_slices; from_slice++)
    {
      int to_slice = from_slice-1;
	
      int N_from_edges = edges_in_slice[from_slice];      
      int N_to_edges = edges_in_slice[to_slice];      

      /*
      */

      for (int ti=0; ti<dc; ti++)
	{
	  int tmp_i = ti * dc_mod;
	  // Qhf(-1,t) in Qh'f(-1,t) format
	  Qef_hp[ti][from_slice].resize(N_from_edges,N_to_edges );
	  for (int edge_f = 0; edge_f<N_to_edges; edge_f++)
	    {
	      for (int edge_h = 0; edge_h<N_from_edges -1; edge_h++)      
		{
		  //  * Qh'f(-1,t)
		  int hp =Sp_below[from_slice][edge_h];	     
		  Qef_hp[ti][from_slice](edge_h,edge_f) = 
		    Qef_s_t[to_slice][tmp_i][ddc-1](hp*N_to_edges + edge_f);	    
		}
	      int edge_h = N_from_edges -1;
	      int gp = S_below[from_slice][0];
	      int gpp = S_below[from_slice][1];
	      // * [ Qg''f(-1,t) Qg'(-1) + Qg'f(-1,t) Qg''(-1)]
	      Qef_hp[ti][from_slice](edge_h,edge_f) = 
		Qef_s_t[to_slice][tmp_i][ddc-1](gp*N_to_edges + edge_f) * Qe[to_slice][ddc-1](gpp) + 
		Qef_s_t[to_slice][tmp_i][ddc-1](gpp*N_to_edges + edge_f) * Qe[to_slice][ddc-1](gp) ;	    
	    }      
	  // Qeh(s,0) in Qeh'(-1,t) format
	  Qef_ep[ti][from_slice].resize(N_from_edges,N_to_edges );
	  for (int edge_e = 0; edge_e<N_from_edges; edge_e++)      	
	    {
	      for (int edge_f = 0; edge_f<N_from_edges-1; edge_f++)	    
		{
		  //  Qef'(s,0) *
		  int fp =Sp_below[from_slice][edge_f];	     
		  Qef_ep[ti][from_slice](edge_e,fp) = 
		    Qef_s_t[from_slice][0][tmp_i](edge_e*N_from_edges + edge_f);	    
		}	  
	      int edge_f = N_from_edges -1;
	      int gp = S_below[from_slice][0];
	      int gpp = S_below[from_slice][1];
	      // Qeg'(s,0) [Qg''(-1)]  *
	      Qef_ep[ti][from_slice](edge_e,gp) = 
		Qef_s_t[from_slice][0][tmp_i](edge_e*N_from_edges + edge_f) * Qe[to_slice][ddc-1](gpp) ;	    
	      // Qeg''(s,0) [Qg'(-1)]  *
	      Qef_ep[ti][from_slice](edge_e,gpp) = 
		Qef_s_t[from_slice][0][tmp_i](edge_e*N_from_edges + edge_f) * Qe[to_slice][ddc-1](gp) ;	    
	    }    
	}
      slice_2_slice[from_slice][to_slice].resize(N_from_edges,N_to_edges );
      // Qef(-1,0) = Qeh(-1,0) * Qh'f(-1,0) + Qeg(-1,0) * Qg'f(-1,0) * Qg''(-1) +  Qeg(-1,0) * Qg''f(-1,0) * Qg'(-1)    
      blas::gemm(slice_2_slice[from_slice][from_slice],Qef_hp[0][from_slice],slice_2_slice[from_slice][to_slice]);      
    }

  for (int from_slice = 0; from_slice<N_slices; from_slice++)
    {
      int to_slice = from_slice;
      //e2e_ft =& edge_2_edge[from_slice][to_slice];

      int N_from_edges = edges_in_slice[from_slice];            
      int N_to_edges = edges_in_slice[to_slice];      	  
      for (int si = 0; si<dc; si++)
	{
	  //e2e_fts =& (*e2e_ft)[si];
	  for (int ti = 0; ti<=si; ti++)
	    {
	      //e2e_ftst =& (*e2e_fts)[ti];		   	      
	      int d_si = si * dc_mod;
	      int d_ti = ti * dc_mod;
	      // (*e2e_ftst).resize(N_from_edges,N_to_edges);
	      edge_2_edge[from_slice][to_slice][si][ti].resize(N_from_edges,N_to_edges);
	      for (int edge_e = 0; edge_e<N_from_edges; edge_e++)
		for (int edge_f = 0; edge_f<N_to_edges; edge_f++) 			
		  {
		    edge_2_edge[from_slice][to_slice][si][ti](edge_e,edge_f) =  
		      Qef_s_t[from_slice][d_ti][d_si](edge_e*N_from_edges + edge_f);		   

		  }
	    }
	}
    }

  // calculate Qef between internal S'' nodes
  for (int from_slice = 1; from_slice<N_slices; from_slice++)
    {
      int N_from_edges = edges_in_slice[from_slice];            

      int to_slice=from_slice-1;
      int N_to_edges = edges_in_slice[to_slice];      	  
      array_type tmp;
      tmp.resize(N_from_edges,edges_in_slice[to_slice+1]);
      for (int si = 0; si<dc; si++)
	{
	  for (int ti = 0; ti<dc; ti++)	      
	    if ( !(si == dc-1 && ti == 0) )
	      {
		if (from_slice == to_slice+1)
		  {
		    edge_2_edge[from_slice][to_slice][si][ti].resize(N_from_edges,N_to_edges);
		    blas::gemm(edge_2_edge[from_slice][from_slice][si][0],Qef_hp[ti][to_slice+1],
		    	       edge_2_edge[from_slice][to_slice][si][ti]);   
		  }
		else
		  {
		    edge_2_edge[from_slice][to_slice][si][ti].resize(N_from_edges,N_to_edges);
		    blas::gemm(Qef_ep[si][from_slice],slice_2_slice[from_slice-1][to_slice+1],
			       tmp);
		    blas::gemm(tmp,Qef_hp[ti][to_slice+1],
			       edge_2_edge[from_slice][to_slice][si][ti]);
		  }
	      }
	}

      tmp.resize(0,0);      

    }
  //cout << "# across all slices " <<t.elapsed() << endl;

  
  //cout << "# lin_size " <<t.elapsed() << endl;
  //t.restart();

  // construct flatQef for efficiinet calculation
  flat_Qef.clear();  delta_dt.clear();  tau_dt.clear();
  // speciation 1 or internal edge 0
  //int type_of_x[lin_size];
   type_of_x.clear();  below_son0.clear();  below_son1.clear();  
  type_of_x.clear();
  // for type 1 xp and xpp
  // for type 0 .
  below_from.clear();
  below_until.clear();
  Qegp.clear();

  for (int from_slice = N_slices-1; from_slice>-1; from_slice--)
    {
      int lin_from = slice_offset[from_slice];
      int N_from_edges = edges_in_slice[from_slice];            
      int fromm1_slice_off=slice_offset[from_slice-1];
      int from_slice_off=slice_offset[from_slice];
      
      for (int si = dc-1; si > -1; si--)
	for (int edge_from = 0; edge_from < N_from_edges; edge_from++)
	  {
	    if ( isorisnot(from_slice,si,edge_from,N_from_edges) )
	      {
		//delta_dt[lin_from] =deltas[from_slice](edge_from) * dt;// 1.-exp(-deltas[from_slice](edge_from)*dt);//(*deltas_slice)(edge_from) * dt;
		delta_dt[lin_from] = 1.-exp(-deltas[from_slice](edge_from)*dt);
		//tau_dt[lin_from] = taus[from_slice](edge_from) * dt;//1.-exp(-taus[from_slice](edge_from)*dt);//(*taus_slice)(edge_from) * dt;
		tau_dt[lin_from] = 1.-exp(-taus[from_slice](edge_from)*dt);
		if ( si == 0 && from_slice>0 )
		  {
		    type_of_x[lin_from] = 1;
		    int xp_edge = S_below[from_slice][0];
		    int xpp_edge = S_below[from_slice][1];
		    below_from[lin_from] = fromm1_slice_off + xp_edge ;//flat_map[from_slice-1][1][xp_edge];//
		    below_until[lin_from] = fromm1_slice_off + xpp_edge ;//flat_map[from_slice-1][1][xpp_edge];//
		    
		    int this_slice=lin_slice[ fromm1_slice_off + xp_edge];
		    int N_edges = edges_in_slice[this_slice];    
		    Qegp[fromm1_slice_off + xp_edge]=Qef_s_t[this_slice][0][ddc-1](xp_edge* N_edges +  xp_edge);
		    Qegp[fromm1_slice_off + xpp_edge]=Qef_s_t[this_slice][0][ddc-1](xpp_edge* N_edges + xpp_edge);		    		    
		    
		    //cout << this_slice <<": "<< Qegp[fromm1_slice_off + xp_edge] << " " << Qegp[fromm1_slice_off + xpp_edge] <<" . "<< dc_mod << " " << ddc-1<< endl;

		    below_son0[lin_from] = fromm1_slice_off + xp_edge ;//flat_map[from_slice-1][1][xp_edge];//
		    below_son1[lin_from] = fromm1_slice_off + xpp_edge ;//flat_map[from_slice-1][1][xpp_edge];//

		    advance_from[lin_from]=flat_map[from_slice][1][0];
		    advance_until[lin_from]=flat_map[from_slice][1][edges_in_slice[from_slice]-1];
		    //cout << lin_from << " o " <<  advance_from[lin_from] << " " << advance_until[lin_from] <<endl;		    
		  }
		else	    
		  {
		    type_of_x[lin_from] = 0;
		    below_from[lin_from] = lin_from - edge_from ;
		    below_until[lin_from] = lin_from - edge_from + N_from_edges;			    
		    if (from_slice==0)
		      {
			advance_from[lin_from]=0;
			advance_until[lin_from]=-1;		       
		      }
		    else if (si==1 && from_slice>0)
		      {
			
			advance_from[lin_from]=flat_map[from_slice][0][N_from_edges-1];
			advance_until[lin_from]=flat_map[from_slice-1][dc-1][edges_in_slice[from_slice-1]-1];		       
			below_son0[lin_from]=fromm1_slice_off +Sp_below[from_slice][edge_from];
			if (edge_from==N_from_edges-1)
			  below_son0[lin_from]=flat_map[from_slice][0][N_from_edges-1];
		      }
		    else
		      {
			advance_from[lin_from]=flat_map[from_slice][si-1][0];
			advance_until[lin_from]=flat_map[from_slice][si-1][N_from_edges-1];
			below_son0[lin_from]=lin_from+ fromm1_slice_off;
		      }
		  }
		flat_Qef[lin_from].resize(lin_size-lin_from+1);
		int N_to_edges = N_from_edges;      	  		  
		int lin_to = from_slice_off;
		for (int ti = dc-1; ti > -1; ti--)
		  for (int edge_to = 0; edge_to < N_to_edges; edge_to++)		  
		    if ( isorisnot(from_slice,ti,edge_to,N_to_edges) )		    
		      {
			if (lin_to>=lin_from)//(ti<si)
			  {
			    if ( !(si == dc-1 && ti == 0) )
			      flat_Qef[lin_from](lin_to-lin_from) = 
				edge_2_edge[from_slice][from_slice][si][ti](edge_from,edge_to);
			    else
			      flat_Qef[lin_from](lin_to-lin_from) = 
				slice_2_slice[from_slice][from_slice](edge_from,edge_to);
			  }
			//XX
			lin_to++;			    
		      }
		lin_from++;
	      }
	  }
      int to_slice = from_slice-1;
      if (from_slice>0  && (to_slice==from_slice-1 || to_slice==0) )
	{
	  lin_from = slice_offset[from_slice];
	  
	  for (int si = dc-1; si > -1; si--)
	    for (int edge_from = 0; edge_from < N_from_edges; edge_from++)
	      if ( isorisnot(from_slice,si,edge_from,N_from_edges) )
		{
		  //e2e_fts=& (*e2e_ft)[si];
		  int lin_to = slice_offset[to_slice];
		  int N_to_edges = edges_in_slice[to_slice];      	  		  
		  for (int ti = dc-1; ti > -1; ti--)
		    for (int edge_to = 0; edge_to < N_to_edges; edge_to++)		  
		      if ( isorisnot(to_slice,ti,edge_to,N_to_edges) )			  
			{
			  //e2e_ftst=& (*e2e_fts)[ti];
			  if ( !(si == dc-1 && ti == 0) )			      
			    flat_Qef[lin_from](lin_to-lin_from) =			    
			      edge_2_edge[from_slice][to_slice][si][ti](edge_from,edge_to);
			  //(*e2e_ftst)(edge_from,edge_to);			      
			  else
			    flat_Qef[lin_from](lin_to-lin_from) =
			      slice_2_slice[from_slice][to_slice](edge_from,edge_to);						       
			  //(*s2s_ft)(edge_from,edge_to);
			  lin_to++;			    
			}
		  lin_from++;
		}
	}
    }
  x_down.resize(lin_size);
  node_down.resize(lin_size);

  for (int x=0;x<lin_size;x++)
    {
      x_down(x)=x_down_f(x);
      node_down(x)=node_down_f(x);
      branch_2_x[x_down(x)]=x;
      //branch_2_slices[x_down(x)][lin_slice[x]]=1;
      //slice_2_branches[lin_slice[x]][x_down(x)]=1;
    }


  //map <int, map <int, Node * > > node_map;
  for (map <int, map <int, Node * > > ::iterator it=node_map.begin();it!=node_map.end();it++)
    {
      (*it).second.clear();        
    }
  node_map.clear();

  //map<int , int> Sp_below[N_slices];
  //map<int , int> S_above[N_slices];
  for (int i = 0;i<N_slices;i++)
    {
      Sp_below[i].clear();
      S_above[i].clear();
      deltas[i].resize(0);
      taus[i].resize(0);
      lambdas[i].resize(0);
      phis[i].resize(0);
      //Qef_s_tail[i].resize(0);
    }
  //map<Node * , int> index_;
  index_.clear();

  //  vector_type  deltas[N_slices],  taus[N_slices], lambdas[N_slices], phis[N_slices];

  term1.resize(0);
  term2.resize(0);
  term3.resize(0);
  dQ.resize(0);
  
  k1.resize(0);
  k2.resize(0);
  k3.resize(0);
  k4.resize(0);
  ktmp.resize(0);
  
  Qsq.resize(0);
  QefQe.resize(0);
  v_tmp.resize(0);

  //array_type Qef_s_tail[N_slices];
  //array_type Qef_hp[dc][N_slices],Qef_ep[dc][N_slices];
  onetree=0;
  llonetree=0;
  zerotree= Qe[lin_slice[0]][0][lin_edge[0]];
  vector < Node * > S_leaves = tree->getLeaves();
  ur_S_N_leaves = tree->getLeaves().size();
  for (vector <Node*>::iterator jt = S_leaves.begin();jt != S_leaves.end(); jt++)
    {
      ur_sigma[(*jt)->getName()]=flat_node_map[(*jt)];

      if ((*jt)->getName()[0]!='*') 
	{
	  if (mode=="intbck" || mode =="bck" || mode=="tree")
	    Qef_bck_N2(*jt);
	  else
	    Qef_N2(*jt);
	    //Qef_m_N2(*jt);
	  onetree+=flat_Qef[0](flat_node_map[(*jt)]);
	  llonetree+=flat_Qef[0](flat_node_map[(*jt)]);
	}
    }
  
  // vector_type  deltas[N_slices],  taus[N_slices], lambdas[N_slices], phis[N_slices];
  
  //vector_type term1;  vector_type term2;  vector_type term3;  vector_type dQ;  
  //vector_type k1; vector_type k2;  vector_type k3;  vector_type k4;  vector_type ktmp; 
  //vector_type QefQe; vector_type Qsq;
  //vector_type v_tmp; 
  term1.resize(0);term2.resize(0);term3.resize(0);
  k1.resize(0); k2.resize(0); k3.resize(0); k4.resize(0); ktmp.resize(0);
  QefQe.resize(0); Qsq.resize(0);
  v_tmp.resize(0);

  term1.clear();term2.clear();term3.clear();
  k1.clear(); k2.clear(); k3.clear(); k4.clear(); ktmp.clear();
  QefQe.clear(); Qsq.clear();
  v_tmp.clear();
  for (int i=0;i< N_slices;i++)
    {
      deltas[i].resize(0);taus[i].resize(0);lambdas[i].resize(0);
      deltas[i].clear();taus[i].clear();lambdas[i].clear();
      
      deltas_s_t[i].resize(0);taus_s_t[i].resize(0); phis_s_t[i].resize(0);
      deltas_s_t[i].clear();taus_s_t[i].clear(); phis_s_t[i].clear();
      for (int j=0;j< ddc;j++)
	{Qe_s_t[i][j].resize(0);Qe_s_t[i][j].clear();}
      for (int j=0;j< ddc;j++)
	for (int k=0;k< ddc;k++)
	  {Qef_s_t[i][j][k].resize(0);Qef_s_t[i][j][k].clear();}      
    }
  S_leaves.clear();
  //vector_type deltas_s_t[N_slices]; vector_type taus_s_t[N_slices]; vector_type phis_s_t[N_slices]; vector_type Qe_s_t[N_slices][ddc];
  //vector_type Qef_s_t[N_slices][ddc][ddc];
  

}


void Species_tree::init_x()
{
  ur_S_N_leaves = tree->getLeaves().size();
  
  double stem_omega=0.1;
  cherries.clear();
  
  advance_from.clear();
  advance_until.clear();
  below_from.clear();
  below_until.clear();
  flat_node_map.clear();
  lin_si.clear();

  id_slice_table.clear();
  id_table.clear();
  l_table.clear();
  N_table.clear();

  infer_mode=true;
  //delete tree;

  init();    
  real_root = root;  root = new Node();  root -> setName("*Rstem1");
  real_root->addSon(root);  root->setDistanceToFather(stem_omega);
  tree->rootAt(root);  

  real_root = root;  root = new Node();  root -> setName("*Rstem2");
  real_root->addSon(root);  root->setDistanceToFather(stem_omega);  
  tree->rootAt(root);

  real_root = root;  root = new Node();  root -> setName("*Rstem3");  
  real_root->addSon(root);  root->setDistanceToFather(stem_omega);
  tree->rootAt(root);

  real_root = root;  root = new Node();  root -> setName("*Rstem4");  
  real_root->addSon(root); root->setDistanceToFather(stem_omega);
  tree->rootAt(root);

  real_root = root; root = new Node();  root -> setName("*Rstem5");
  real_root->addSon(root);  root->setDistanceToFather(stem_omega);
  tree->rootAt(root);

  bool out_group= false;//( out_tau>0 );// && (mode!="max");
  if ( out_group )
    {
      Node * out_group = new Node();  
      root->addSon(out_group);
      out_group->setDistanceToFather(1.+stem_omega*5.);
      out_group -> setName("OUT_GROUP");
    }
  
  construct_slices();
  N_slices =  time_slices.size();
  map <int, map <int, Node * > > node_map;
  for (int slice = 0; slice<N_slices;slice++)
    {
      //cout << "# "<<slice <<  " : " ; 
      for (int edge=0;edge<time_slices[N_slices-slice-1].size();edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  node_map[slice][edge] = node;
	  //cout << node->getName() << " " << node->getDistanceToFather() << " | ";
	}
      //cout << endl;
    }
 
  flat_node_map.clear();
  inverse_flat_node_map.clear();

  int slice_offset[N_slices];
  //topology coding

  //the indicies in slice-1 of the descendants of S nodes (i.e. g -> g'',g'') 
  //in any given slice the speciation node always has index -1 (i.e edges_in_slice[slice]-1 )
  int S_below[N_slices][2];
  //the indicies in slice-1 of the descendants of Sp nodes (i.e. h) 
  map<int , int> Sp_below[N_slices];
  //the index of nodes father in slice+1 
  map<int , int> S_above[N_slices];
  // index of node in it's slice
  map<Node * , int> index_;

  for (int slice = 0; slice<N_slices;slice++)
    {
      edges_in_slice[slice]=time_slices[N_slices-slice-1].size();
      for (int edge=0;edge<edges_in_slice[slice];edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  index_[node] = edge;
	}
    }
  for (int slice = 0; slice<N_slices;slice++)
    {
      int tmp=0;
      for (int edge=0;edge<edges_in_slice[slice];edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  if(slice > 0)
	    {
	      if(node->getNumberOfSons()==2)
		{
		  S_below[slice][0] = index_[node->getSons()[0]];
		  S_below[slice][1] = index_[node->getSons()[1]];
		}	    
	      else
		{
		  Sp_below[slice][edge] = index_[node->getSons()[0]];
		}	    	         
	    }
	  if (slice<N_slices-1)
	    S_above[slice][edge] = index_[node->getFather()];
	}
    }

  //resolution for gene tree reconciliation
  lin_size=0; inverse_flat_node_map.clear();lin_slice.clear();  lin_edge.clear();
  dc=2;
  dc_mod = (ddc-1)/(dc-1);

  // construct index for flat iteration forward in time over
  // S'' nodes that are in S or are not contemporary to a node in S
  for (int from_slice = N_slices-1; from_slice>-1; from_slice--)
    {
      slice_offset[from_slice] = lin_size;
      int N_from_edges = edges_in_slice[from_slice];  
      for (int si = dc-1; si > -1; si--)
	for (int edge_from = 0; edge_from < N_from_edges; edge_from++)
	  // S'' nodes that are in S or are not contemporary to a node in S
	  if ( isorisnot(from_slice,si,edge_from,N_from_edges) )
	    {
	      flat_map[from_slice][si][edge_from]=lin_size;
	      flat_node_map[node_map[from_slice][edge_from]] = lin_size;
	      inverse_flat_node_map[lin_size] = node_map[from_slice][edge_from];
	      
	      lin_slice[lin_size] = from_slice;
	      lin_edge[lin_size] = edge_from;
	      lin_si[lin_size] = si; 
	      lin_size++;
	    }
    }
  //NEWLIN SIZE 

  type_of_x.clear();
  // for type 1 xp and xpp
  // for type 0 .
  below_from.clear();
  below_until.clear();
  
  for (int from_slice = N_slices-1; from_slice>-1; from_slice--)
    {
      int lin_from = slice_offset[from_slice];
      int N_from_edges = edges_in_slice[from_slice];            
      int fromm1_slice_off=slice_offset[from_slice-1];
      int from_slice_off=slice_offset[from_slice];      
      for (int si = dc-1; si > -1; si--)
	for (int edge_from = 0; edge_from < N_from_edges; edge_from++)
	  {
	    if ( isorisnot(from_slice,si,edge_from,N_from_edges) )
	      {
		if ( si == 0 && from_slice>0 )
		  {
		    type_of_x[lin_from] = 1;
		    int xp_edge = S_below[from_slice][0];
		    int xpp_edge = S_below[from_slice][1];
		    below_from[lin_from] = fromm1_slice_off + xp_edge ;//flat_map[from_slice-1][1][xp_edge];//
		    below_until[lin_from] = fromm1_slice_off + xpp_edge ;//flat_map[from_slice-1][1][xpp_edge];//
		    

		    below_son0[lin_from] = fromm1_slice_off + xp_edge ;//flat_map[from_slice-1][1][xp_edge];//
		    below_son1[lin_from] = fromm1_slice_off + xpp_edge ;//flat_map[from_slice-1][1][xpp_edge];//

		    advance_from[lin_from]=flat_map[from_slice][1][0];
		    advance_until[lin_from]=flat_map[from_slice][1][edges_in_slice[from_slice]-1];
		    //cout << lin_from << " o " <<  advance_from[lin_from] << " " << advance_until[lin_from] <<endl;		    
		  }
		else	    
		  {
		    type_of_x[lin_from] = 0;
		    below_from[lin_from] = lin_from - edge_from ;
		    below_until[lin_from] = lin_from - edge_from + N_from_edges;			    
		    if (from_slice==0)
		      {
			advance_from[lin_from]=0;
			advance_until[lin_from]=-1;		       
		      }
		    else if (si==1 && from_slice>0)
		      {			
			advance_from[lin_from]=flat_map[from_slice][0][N_from_edges-1];
			advance_until[lin_from]=flat_map[from_slice-1][dc-1][edges_in_slice[from_slice-1]-1];		       
			below_son0[lin_from]=fromm1_slice_off +Sp_below[from_slice][edge_from];
			if (edge_from==N_from_edges-1)
			  below_son0[lin_from]=flat_map[from_slice][0][N_from_edges-1];
		      }
		    else
		      {
			advance_from[lin_from]=flat_map[from_slice][si-1][0];
			advance_until[lin_from]=flat_map[from_slice][si-1][N_from_edges-1];
			below_son0[lin_from]=lin_from+ fromm1_slice_off;
		      }
		  }
		lin_from++;
	      }
	  }
    }
  x_down.resize(lin_size);
  node_down.resize(lin_size);
  branch_2_slices.clear();
  slice_2_branches.clear();

  for (int x=0;x<lin_size;x++)
    {
      x_down(x)=x_down_f(x);
      branchname[x_down(x)]=inverse_flat_node_map[node_down(x)]->getName();
      branchi[inverse_flat_node_map[node_down(x)]->getName()]=x_down(x);
      node_down(x)=node_down_f(x);
      branch_2_x[x_down(x)]=x;
      branch_2_slices[x_down(x)][lin_slice[x]]=1;
      slice_2_branches[lin_slice[x]][x_down(x)]=1;
    }

  //NEWLIN SIZE 



}
void Species_tree::reset(scalar_type delta,scalar_type tau,scalar_type lambda,scalar_type stem_omega)
{
  init_x();
  init_L_improved(delta,tau,lambda,0.1,0.01,"estimate");   
  init_L_improved(delta,tau,lambda,0.1,0.01,"intbck");   

  O_counts.clear();  D_counts.clear();  T_counts.clear(); from_count.clear();
  branch_sum.clear();  branch_count.clear();  branch_zero.clear();branch_genome_size.clear();
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
    {
      from_count.push_back(0); O_counts.push_back(0); D_counts.push_back(0); T_counts.push_back(0);
      branch_sum.push_back(0); branch_count.push_back(0); branch_zero.push_back(0);branch_genome_size.push_back(0);
    }  
  /*
  cherries.clear();
  for (std::map < std::pair<node_type *,node_type *>,long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)  
    delete (*it).second;
  ur_a.clear();
  cherry_a.clear();
  cherry_reg.clear();
  leaf_hidden_event_name.clear();
  leaf_hidden_event_x.clear();
  max_hidden_events.clear();
  max_event_x.clear();
  max_event_xc1.clear();
  max_event_xc2.clear();
  max_event_name.clear();
  */
  for (int x=0;x<lin_size;x++)
    {
      scalar_type t=branch_length[x_down(x)];
      if (x_down(x)==1)
	t=stem_omega;
      string name=branchname[x_down(x)];
      branchi[name]=x_down(x);

      if (x_down(x)==1)
	{
	  branchwise_omega[x_down(x)]=1;
	  namewise_omega[name]=1;	  
	}
      else
	{
	  branchwise_omega[x_down(x)]=1;
	  namewise_omega[name]=1;
	}
      namewise_delta[name]=delta*t;
      namewise_tau[name]=tau*t;
      namewise_lambda[name]=lambda*t;
    }  
  p_omega.resize(lin_size);
  scalar_type norm=0;
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)    
    if (type_of_x[x]!=1)
      {
	p_omega(x)=branchwise_omega[x_down(x)];	  
	norm+=p_omega(x);
      }
    else
      p_omega(x)=0.;
  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
    {
      p_omega(x)/=norm;
      branchwise_omega[x_down(x)]=p_omega(x);
      string name=branchname[x_down(x)];
      namewise_omega[name]=p_omega(x);
      
    }
  
  //for (map<string,scalar_type>::iterator it=namewise_delta.begin();it!=namewise_delta.end();it++) (*it).second=0.01; 
  //for (map<string,scalar_type>::iterator it=namewise_tau.begin();it!=namewise_tau.end();it++) (*it).second=0.01; 
  //for (map<string,scalar_type>::iterator it=namewise_lambda.begin();it!=namewise_lambda.end();it++) (*it).second=0.01; 

  
}
