#include "DTL.h"
#include <float.h>
using namespace std;
using namespace bpp;


#define  isorisnot( slice, si, edge, N_edges) (si==0 && (edge==N_edges-1 || slice == 0)) || (si>0)

void Species_tree:: emit_G_trees(int N_trees , vector<string> * G_trees,vector<string> * codes,vector<string> * stream)
{
  tree_norm=0;
  for (int tree_i=0;tree_i<N_trees;)
    {
      //cout << tree_i << " **-**" <<endl;
      map<int,int> tmp_count_branch_omega,tmp_count_branch_D,tmp_count_branch_T,tmp_count_branch_L;
      tree_norm+=1;
  
      Node * G_root= new Node();
      TreeTemplate<Node>  G_tree(G_root);
      
      scalar_type expected_number_of_originations = 0;      
      for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	{
	  if (type_of_x[x]!=1 )	      
	    {
	      Node* node=inverse_flat_node_map[x];
	      if (node!=root)
		expected_number_of_originations+=node->getDistanceToFather()*branch_omega[node_down(x)];
	    }
	}
      expected_number_of_originations = omega_root;
      int root_x = 1;
      // ########## ORIGINATION ################
      //Poission(lambda,0) = exp(-lambda)
      //if one or more origination events occured 
      if( exp(-expected_number_of_originations)  < RandomTools::giveRandomNumberBetweenZeroAndEntry(1.))
	{
	  //we chose a time and place at random
	  scalar_type sum_tmp=0;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    if (type_of_x[x]!=1 )	      	    
	      sum_tmp+=branch_omega[node_down(x)];
	  scalar_type stop_tmp=RandomTools::giveRandomNumberBetweenZeroAndEntry(sum_tmp); 
	  sum_tmp=0.;
	  for (int x=lin_size-ur_S_N_leaves;x>=0;x--)
	    if (type_of_x[x]!=1 )	      
	      {
		sum_tmp+=branch_omega[node_down(x)];
		if (stop_tmp<sum_tmp)
		  {
		    root_x=x;
		    break;
		  }
	      }	 	  
	  
	}
      sim_x=root_x;

      tmp_count_branch_omega[x_down(root_x)]+=1;

      // ########## ORIGINATION ################
      G_root->setName(inverse_flat_node_map[root_x]->getName());
      pair <Node *,int> root_gene(G_root,root_x);
      vector < pair <Node *,int> > genes;      
      vector < pair <Node *,int> > new_genes;

      genes.push_back(root_gene);
      
      int slice=lin_slice[root_x];
      //cout<<lin_slice[1] << "# slices"<<endl;
      Node * next_node;
      scalar_type next_t;
      scalar_type time_from_root=0;
      bool saw_zero=false;
      while (slice>=0)
	{
	  //cout << slice << endl;
	  //cout << N_slices-slice-1 << " " << time_slices[N_slices-slice-1].size()-1 <<" " << inverse_flat_node_map[root_x]->getName() << " " << x_down(root_x) <<endl;
	  next_node=time_slices[N_slices-slice-1][time_slices[N_slices-slice-1].size()-1];
	  //cout << "current " << G_root->getName() << " " << type_of_x[root_x]<< endl;
	  //cout << "next " << next_node->getName() << " " << type_of_x[flat_node_map[next_node] ]<< endl;


	  next_t=next_node->getDistanceToFather();	  
	  //cout << next_t << endl;
	  scalar_type t=0,rate_sum;
	  while (1)
	    {
	      
	      // ########## RATE SUM ################
	      rate_sum=0;      
	      //cout << "GS: ";
	      for ( vector < pair <Node *,int> >::iterator gi=genes.begin();gi!=genes.end();gi++)
		{
		  int branch_x=node_down((*gi).second);
		  rate_sum+=branch_delta[branch_x];
		  rate_sum+=branch_lambda[branch_x];
		  //cout << x_down((*gi).second) << "," << lin_slice[(*gi).second] << " : " ;
		}
	      //cout << endl;
	      for (int edge=0;edge<time_slices[N_slices-slice-1].size();edge++)
		{
		  scalar_type C=0.;		  
		  Node * node=time_slices[N_slices-slice-1][edge];
		  int branch_x = flat_node_map[node];
		  if (type_of_x[branch_x]!=1 && slice>0)
		    branch_x=below_son0[branch_x];
		  for ( vector < pair <Node *,int> >::iterator di=genes.begin();di!=genes.end();di++)
		    if ((*di).second!=branch_x)
		      C+=1.;

		  if (C>0)
		    for ( vector < pair <Node *,int> >::iterator di=genes.begin();di!=genes.end();di++)
		      if ((*di).second!=branch_x)
			rate_sum+=branch_tau[node_down(branch_x)]*Tef_var[node_down(branch_x)][node_down((*di).second)]/C;
		}

	      // ########## RATE SUM ################
	      scalar_type time_to_next_event=RandomTools::randExponential(1/rate_sum);
	      if (t+time_to_next_event<next_t)
		{
		  t+=time_to_next_event;

		  time_from_root+=time_to_next_event;
		  // ########## RATE RE-SUM ################
		  scalar_type stop_sum=RandomTools::giveRandomNumberBetweenZeroAndEntry(rate_sum); 
		  rate_sum=0;      
		  for ( vector < pair <Node *,int> >::iterator gi=genes.begin();gi!=genes.end();gi++)
		    {
		      int branch_x=node_down((*gi).second);
		      rate_sum+=branch_delta[branch_x];
		      if (stop_sum<rate_sum)
			{
			  if ((*gi).first->hasFather())
			    {
			      scalar_type dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*gi).first->getFather());
			      (*gi).first->setDistanceToFather(time_from_root-dtmp);
			    }
			  Node * son0=new Node();
			  Node * son1=new Node();		  
			  stringstream out;
			  out<<  lin_slice[0]-lin_slice[(*gi).second]+1 << "|" << x_down((*gi).second)<< "|" << (*gi).second;
			  string name="D@"+out.str();
			  tmp_count_branch_D[x_down((*gi).second)]+=1;
			  (*gi).first->setBranchProperty("ID",BppString(name));
			  (*gi).first->setName( inverse_flat_node_map[x_down((*gi).second)]->getName() );
			  // XX
			  //cout << " E===>"<< name << endl;

			  (*gi).first->addSon(son0);
			  (*gi).first->addSon(son1);

			  pair <Node *,int> gene0(son0,(*gi).second);
			  pair <Node *,int> gene1(son1,(*gi).second);
			  son0->setDistanceToFather(0);
			  son1->setDistanceToFather(0);
			  son0->setName((*gi).first->getName());
			  son1->setName((*gi).first->getName());
			      
			  genes.erase(gi);
			  genes.push_back(gene0);
			  genes.push_back(gene1);
			  break;
			}
		      rate_sum+=branch_lambda[branch_x];
		      if (stop_sum<rate_sum)
			{
			  if ((*gi).first->hasFather())
			    {
			      double dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*gi).first->getFather());
			      (*gi).first->setDistanceToFather(time_from_root-dtmp);
			    }
			  stringstream out;
			  out<<  lin_slice[0]-lin_slice[(*gi).second]+1 << "|" << x_down((*gi).second) << "|" << (*gi).second;
			  
			  string name="*L@"+out.str();
			  tmp_count_branch_L[x_down((*gi).second)]+=1;

			  (*gi).first->setName( name );
			  (*gi).first->setBranchProperty("ID",BppString(name));
			  // XX
			  //cout<< " E===>" << name << endl;

			  genes.erase(gi);
			  break;
			}		     
		    }
		  if (stop_sum>rate_sum)
		    for (int edge=0;edge<time_slices[N_slices-slice-1].size();edge++)
		      {
			double trf_sum=0.;
			double C=0.;		  
			Node * node=time_slices[N_slices-slice-1][edge];
			int branch_x = flat_node_map[node];
			if (type_of_x[branch_x]!=1 && slice>0)
			  branch_x=below_son0[branch_x];
			for ( vector < pair <Node *,int> >::iterator di=genes.begin();di!=genes.end();di++)
			  if (x_down((*di).second)!=x_down(branch_x))
			    C+=1.;
			if (C>0)
			  for ( vector < pair <Node *,int> >::iterator di=genes.begin();di!=genes.end();di++)
			    if (x_down((*di).second)!=x_down(branch_x))
			      {
				rate_sum+=branch_tau[node_down(branch_x)]*Tef_var[node_down(branch_x)][node_down((*di).second)]/C;
				if (stop_sum<rate_sum)
				  {
				    
				    if ((*di).first->hasFather())
				      {
					double dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*di).first->getFather());
					(*di).first->setDistanceToFather(time_from_root-dtmp);
				      }
				    Node * son0=new Node();
				    Node * son1=new Node();		  
				    stringstream out;
				    out<<  lin_slice[0]-lin_slice[(*di).second]+1 << "|" << x_down(branch_x) << "|" << x_down((*di).second) << "|" << branch_x;
				    string name="T@"+out.str();
				    tmp_count_branch_T[x_down((*di).second)]+=1;

				    (*di).first->setBranchProperty("ID",BppString(name));
				    (*di).first->setName( inverse_flat_node_map[branch_x]->getName() );
				    
				    (*di).first->addSon(son0);
				    (*di).first->addSon(son1);
				    
				    pair <Node *,int> gene0(son0,(*di).second);
				    pair <Node *,int> gene1(son1,branch_x);
				    son0->setDistanceToFather(0);
				    son1->setDistanceToFather(0);
				    son0->setName((*di).first->getName());
				    son1->setName(inverse_flat_node_map[branch_x]->getName());

				    // XX
				    //cout << " E===>" << name << endl;
				    
				    genes.erase(di);
				    genes.push_back(gene0);
				    genes.push_back(gene1);

				    break;
				  }				
			      }
			if (stop_sum<rate_sum)
			  break;

		      }
		  // ########## RATE RE-SUM ################
		}
	      else
		{
		  time_from_root+=next_t-t;
		  break;		  
		}
	    }

	  if (slice==0)
	    for ( vector < pair <Node *,int> >::iterator gi=genes.begin();gi!=genes.end();gi++)	    
	      {
		new_genes.push_back((*gi));
		(*gi).first->setName( inverse_flat_node_map[(*gi).second]->getName());
		(*gi).first->setBranchProperty("ID",BppString(inverse_flat_node_map[(*gi).second]->getName()));

	      }
	  else
	    for ( vector < pair <Node *,int> >::iterator gi=genes.begin();gi!=genes.end();gi++)	    
	      {
		if (type_of_x[(*gi).second]!=1 )
		  {
		    if (slice>=1)
		      {
			if (lin_slice[(*gi).second]>0)
			  {
			  (*gi).second=below_son0[(*gi).second];
			  if ((*gi).first->hasFather())
			    {
			      double dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*gi).first->getFather());
			      (*gi).first->setDistanceToFather(time_from_root-dtmp);
			    }
			  }
		      }
		    new_genes.push_back((*gi));
		  }
		else
		  {
		    Node * son0=new Node();
		    Node * son1=new Node();		  
		    stringstream out;
		    out<<  lin_slice[0]-lin_slice[(*gi).second]+1 << "|" << x_down((*gi).second)<< "|" << (*gi).second;
		    string name="S@"+out.str();
		    //XX 
		    //cout<< " E===>" << name <<endl;

		    (*gi).first->setBranchProperty("ID",BppString(name));
		    (*gi).first->setName( inverse_flat_node_map[x_down((*gi).second)]->getName() );

		    (*gi).first->addSon(son0);
		    (*gi).first->addSon(son1);
		    if ((*gi).first->hasFather())
		      {
			double dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*gi).first->getFather());
			(*gi).first->setDistanceToFather(time_from_root-dtmp);
		      }
		    int x0 = below_son0[below_son0[(*gi).second]];
		    int x1 = below_son0[below_son1[(*gi).second]];
		    if (slice==1)
		      {
			x0=below_son0[(*gi).second];
			x1=below_son1[(*gi).second];
		      }
		    pair <Node *,int> gene0(son0,x0);
		    pair <Node *,int> gene1(son1,x1);
		  
		    gene0.first->setName( inverse_flat_node_map[gene0.second]->getName() );
		    gene1.first->setName( inverse_flat_node_map[gene1.second]->getName() );

		    new_genes.push_back(gene0);
		    new_genes.push_back(gene1);

		  }

	      }
	  genes.clear();
	  int max_slice=0;
	  //cout << "NGS : ";
	  for ( vector < pair <Node *,int> >::iterator gi=new_genes.begin();gi!=new_genes.end();gi++)	    
	    {
	      if (lin_slice[(*gi).second]>max_slice)
		max_slice=lin_slice[(*gi).second];
	      genes.push_back((*gi));
	      //cout << x_down((*gi).second) << "," << lin_slice[(*gi).second] << " : " ;

	    }
	  //cout << endl;
	  new_genes.clear();
	  //slice=max_slice;
	  //if (saw_zero && slice==0) slice=-1;
	  //if (slice==0) saw_zero=true;
	  //cout << "  ! " << max_slice << endl;
	  slice-=1;

	}
      for ( vector < pair <Node *,int> >::iterator gi=genes.begin();gi!=genes.end();gi++)
	{
	  (*gi).first->setName(inverse_flat_node_map[(*gi).second]->getName());
	  if ((*gi).first->hasFather())
	    {
	      double dtmp=TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*G_root,*(*gi).first->getFather());
	      (*gi).first->setDistanceToFather(time_from_root-dtmp);
	    }
	}
      if (G_tree.getNumberOfNodes()>1)      
	{
	  //cout<< " VN " <<TreeTemplateTools::treeToParenthesis(G_tree,false,"ID");
	  vector <Node*> leaves = G_tree.getLeaves();
	  /*
	  for (vector <Node*>::iterator lit=leaves.begin(); lit!=leaves.end(); lit++)
	    {
	      string name=(*lit)->getName();
	      if (name.find("*Sp.")!=name.npos)
		{
		  vector<string> tokens;
		  Tokenize(name,tokens,".");
		  name=tokens[1];	      
		  (*lit)->setName(name);
		}
	    }
	  */
	  G_strip_virtual_leaves(&G_tree,true);
	  //if (G_tree.getNumberOfNodes()>1)      
	    //cout<< " SVN " <<TreeTemplateTools::treeToParenthesis(G_tree,false,"ID");
	    //	  else
	    //cout<< " SVN " << "ORPHAN" << endl;
	}

      int real_leaves=0;
      vector <string> leaves = G_tree.getLeavesNames();
      for (vector<string>::iterator nt=leaves.begin();nt!=leaves.end();nt++)
	{
	  vector<string> tokens;
	  Tokenize((*nt),tokens,"&");
	  string name=tokens[0];	      
	  if (name!="OUT_GROUP")
	    {
	      real_leaves++;
	    }
	}
      if (real_leaves>1)//ooo
	{
	  for (int i = 0;i<N_slices+tree->getNumberOfLeaves()-1;i++)
	    {
	      count_branch_D[i]+=tmp_count_branch_D[i];
	      count_branch_T[i]+=tmp_count_branch_T[i];
	      count_branch_L[i]+=tmp_count_branch_L[i];
	      count_branch_omega[i]+=tmp_count_branch_omega[i];	      
	    }

	  tree_i++;
	  string xL_tree= TreeTemplateTools::treeToParenthesis(G_tree,false,"ID");
	  //cout << xL_tree;
	  //cout << stream_xL(&G_tree) << endl;
	  if (stream!=NULL)
	    stream->push_back(stream_xL(&G_tree));

	  /*
	  G_root=G_tree.getRootNode();
	  vector <Node * > leaves = TreeTemplateTools::getLeaves(*G_root);
	  vector<string> tokens;
	  for (vector<Node *>::iterator it = leaves.begin(); it!=leaves.end(); it++ )
	    {
	      tokens.clear();
	      string name = (*it)->getName();
	      Tokenize(name,tokens,"&");
	      (*it)->setName(tokens[0]);	      
	    }
	  */
	  stringstream out;
	  out << tree_i;
	  string G_tree_string = TreeTemplateTools::treeToParenthesis(G_tree,false,"ID");	  
	  if (G_trees!=NULL)
	    G_trees->push_back(G_tree_string);

	  if (codes!=NULL)
	    codes->push_back( "# SBG"+out.str());
	  (*event_stream) << "# SBG" <<out.str()<<endl; 
	  (*event_stream) << G_tree_string;
	  (*event_stream) << (*stream)[(*stream).size()-1] << endl;

	  //XX//
	}
      //destroy the G_tree!
    }
  //cout<<TreeTemplateTools::treeToParenthesis(*tree,false,"ID");
}


string Species_tree::stream_xL(TreeTemplate<Node>* G_tree)
{
  node_type * root= G_tree->getRootNode();
  string name=  (* (dynamic_cast<const BppString *>(root->getBranchProperty("ID")))).toSTL();
  vector<string> tokens;
  Tokenize(name,tokens,"&");
  name=tokens[0];
  int slice=lin_slice[sim_x];     
  sim_O_profile[N_slices-slice]+=1./edges_in_slice[slice];
  //root->setBranchProperty("ID",BppString(tokens[0]));
  string out=stream_xL(G_tree->getRootNode());
  return out.erase(0,1);
}
string Species_tree::stream_xL(Node* node)
{
  string open="";
  string close="";

  
  if (node->isLeaf())
    {
      string name=  node->getName();
      vector<string> tokens;
      Tokenize(name,tokens,"&");
      for(vector<string>::reverse_iterator st=tokens.rbegin();st<tokens.rend();++st)
	{
	  if (string::npos == (*st).find("@"))
	    {
	      stringstream out;
	      out<<x_down(ur_sigma[(*st)]);
	      open+="|O@"+out.str();
	    }
	  else
	    {
	      vector<string> event_token;
	      Tokenize((*st),event_token,"|");
	      string branch= event_token[1];
	      string tmp=event_token[0];
	      event_token.clear();
	      Tokenize(tmp,event_token,"@");
	      string e_type= event_token[0];
	      string e_time= event_token[1];
	      if (e_type=="SL")
		e_type="S";
	      if (e_type=="TL")
		e_type="T";	      
	      if (e_type=="S")
		{		  
		  open +="|<@"+branch+"|"+e_type+"@"+branch;
		  close="|>@"+branch+close;
		}
	      else if (e_type=="T")
		{
		  open +="|<T@"+branch+"|"+e_type+"@"+branch;
		  close="|T>@"+branch+close;
		}
	      else if (e_type=="D" )
		{
		  open +="|D@"+branch;
		}

	    }
	}
      return open+close;
    }
  else
    {
      string name=  (* (dynamic_cast<const BppString *>(node->getBranchProperty("ID")))).toSTL();
      vector<string> tokens;
      Tokenize(name,tokens,"&");
      for(vector<string>::reverse_iterator st=tokens.rbegin();st<tokens.rend();++st)
	{
	  vector<string> event_token;
	  Tokenize((*st),event_token,"|");
	  string branch= event_token[1];
	  string tmp=event_token[0];
	  event_token.clear();
	  Tokenize(tmp,event_token,"@");
	  string e_type= event_token[0];
	  string e_time= event_token[1];

	  if (e_type=="SL")
	    e_type="S";
	  if (e_type=="TL")
	    e_type="T";	      
	  if (e_type=="S")
	    {
	      open +="|<@"+branch+"|"+e_type+"@"+branch;
	      close="|>@"+branch+close;
	    }
	  else if (e_type=="T")
	    {
	      open +="|<T@"+branch+"|"+e_type+"@"+branch;
	      close="|T>@"+branch+close;
	    }
	  else if (e_type=="D" )
	    {
	      open +="|D@"+branch;
	    }
	  
	}
      
      return open+stream_xL(node->getSon(0))+stream_xL(node->getSon(1))+close;
    }
}


void Species_tree::init_sim(scalar_type stem_len, scalar_type delta_mean, scalar_type tau_mean, scalar_type lambda_mean, scalar_type delta_sd, scalar_type tau_sd, scalar_type lambda_sd, scalar_type omega_mean, scalar_type omega_sd, scalar_type pair_sd, scalar_type out_tau)
{
  ur_S_N_leaves = tree->getLeaves().size();

  /*
    ur_sigma.clear();
    for (std::map < std::pair<node_type *,node_type *>, long_vector_type *>::iterator it=ur_a.begin();it!=ur_a.end();it++)  
    delete (*it).second;
    ur_a.clear();
    for (  std::map < std::pair<int,int>, long_vector_type *>::iterator it= cherry_a.begin();it!= cherry_a.end();it++)  
    delete (*it).second;
    cherry_a.clear();
  */  
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
  //boost::progress_timer t;
  //t.restart();  
  init();    
  //t.restart();
  real_root = root;
  root = new Node();
  root -> setName("*Rstem1");
  real_root->addSon(root);
  root->setDistanceToFather(stem_len);
  tree->rootAt(root);
  bool out_group= ( out_tau>0 );// && (mode!="max");
  if ( out_group &&0)
    {
      Node * out_group = new Node();  
      root->addSon(out_group);
      out_group->setDistanceToFather(1.+stem_len);
      out_group -> setName("OUT_GROUP");
    }
  
  construct_slices();
  //real_root->setDistanceToFather(0);
  delta=delta_mean;
  tau=tau_mean;
  lambda=lambda_mean;
  N_slices =  time_slices.size();

  // tmp
  // tmp
  // tmp
  //map <int, map <int, Node * > > node_map;
  for (int slice = 0; slice<N_slices;slice++)
    {
      for (int edge=0;edge<time_slices[N_slices-slice-1].size();edge++)
	{
	  Node * node=time_slices[N_slices-slice-1][edge];
	  node_map[slice][edge] = node;
	}
    }
  // tmp
  // tmp
  // tmp
 
  //  scalar_type edges_in_slice[N_slices];
  // resolution for time integration
  ddc=10;
  scalar_type ddt = 1./((scalar_type) ddc-1.); 
  //construct vectors for Qe integration
  //vector_type  Qe[N_slices][ddc]; 
  //rate tree coding
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
  // MEM: N^2 x ddc


  
  
  //t.restart();

 
  //resolution for gene tree reconciliation
  dc=2;
  dc_mod = (ddc-1)/(dc-1);
  scalar_type dt = 1./((scalar_type) dc -1.); 



	
    

  flat_node_map.clear();
  inverse_flat_node_map.clear();
  //map < int , map <int , map <int , int > > > flat_map;

  int slice_offset[N_slices];
  lin_size=0;
  lin_slice.clear();
  lin_si.clear();
  lin_edge.clear();

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
  

  // construct flatQef for efficiinet calculation
  // vector_type flat_Qef[lin_size];
  // scalar_type delta_dt[lin_size];
  // scalar_type tau_dt[lin_size];
  flat_Qef.clear();
  delta_dt.clear();
  tau_dt.clear();
  // speciation 1 or internal edge 0
  //int type_of_x[lin_size];
  type_of_x.clear();
  // for type 1 xp and xpp
  // for type 0 ..
  //int below_from[lin_size];
  //int below_until[lin_size];
  below_from.clear();
  below_until.clear();

  //vector_type * deltas_slice;
  //vector_type * taus_slice;

  //array_type * s2s_ff;
  //array_type * s2s_ft;

  for (int from_slice = N_slices-1; from_slice>-1; from_slice--)
    {
      int lin_from = slice_offset[from_slice];
      int N_from_edges = edges_in_slice[from_slice];            
      //deltas_slice = &deltas[from_slice];
      //taus_slice = &taus[from_slice];
      int fromm1_slice_off=slice_offset[from_slice-1];
      int from_slice_off=slice_offset[from_slice];
      //e2e_f=& edge_2_edge[from_slice];
      //e2e_ff=& (*e2e_f)[from_slice];
      //s2s_ff=& slice_2_slice[from_slice][from_slice];
      
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
		//flat_Qef[lin_from].resize(lin_size-lin_from+1);
		int N_to_edges = N_from_edges;      	  		  
		int lin_to = from_slice_off;
		for (int ti = dc-1; ti > -1; ti--)
		  for (int edge_to = 0; edge_to < N_to_edges; edge_to++)		  
		    if ( isorisnot(from_slice,ti,edge_to,N_to_edges) )		    
		      {
			lin_to++;			    
		      }
		lin_from++;
	      }
	  }
      int to_slice = from_slice-1;
      //for (int to_slice = from_slice-1; to_slice >= 0; to_slice--)
      if (from_slice>0  && (to_slice==from_slice-1 || to_slice==0) )
	{
	  //e2e_ft=& (*e2e_f)[to_slice];
	  //s2s_ft=& slice_2_slice[from_slice][to_slice];	   
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
      branch_2_slices[x_down(x)][lin_slice[x]]=1;
      slice_2_branches[lin_slice[x]][x_down(x)]=1;
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


  //array_type Qef_s_tail[N_slices];
  //array_type Qef_hp[dc][N_slices],Qef_ep[dc][N_slices];
  vector < Node * > S_leaves = tree->getLeaves();
  
  ur_S_N_leaves = S_leaves.size();
  
  for (vector <Node*>::iterator jt = S_leaves.begin();jt != S_leaves.end(); jt++)
    {
      ur_sigma[(*jt)->getName()]=flat_node_map[(*jt)];
    }
  branch_delta.resize(lin_size);
  branch_tau.resize(lin_size);
  branch_lambda.resize(lin_size);
  branch_omega.resize(lin_size);
  (*event_stream) << "#Branch-wise rates:" <<endl; 
  (*event_stream) << "#branch_id"<<"\t"<<"D_rate"<<"\t"<<"\t"<<"T_rate"<<"\t"<<"L_rate"<<"\t"<<"Orig._rate"<<endl; 
  tmp_branch_delta.clear();tmp_branch_tau.clear();tmp_branch_lambda.clear();tmp_branch_omega.clear();
  for (int x=0;x<lin_size;x++)
    {
      if (type_of_x[x]==1 || inverse_flat_node_map[x]->isLeaf())
	{
	  int branch_x=x_down(x);
	  string name=inverse_flat_node_map[node_down(x)]->getName();

	  if (delta_sd>0)
	    tmp_branch_delta[branch_x]=RandomTools::randGamma(delta_mean*delta_mean/delta_sd/delta_sd , delta_mean/delta_sd/delta_sd);
	  else
	    tmp_branch_delta[branch_x]=delta_mean;
	  if (tau_sd>0)
	    tmp_branch_tau[branch_x]=RandomTools::randGamma(tau_mean*tau_mean/tau_sd/tau_sd,tau_mean/tau_sd/tau_sd);
	  else 
	    tmp_branch_tau[branch_x]=tau_mean;
	  if (lambda_sd>0)
	    tmp_branch_lambda[branch_x]=RandomTools::randGamma(lambda_mean*lambda_mean/lambda_sd/lambda_sd,lambda_mean/lambda_sd/lambda_sd);
	  else
	    tmp_branch_lambda[branch_x]=lambda_mean;
	  if (omega_sd>0)
	    tmp_branch_omega[branch_x]=RandomTools::randGamma(omega_mean*omega_mean/omega_sd/omega_sd,omega_mean/omega_sd/omega_sd);
	  else
	    tmp_branch_omega[branch_x]=omega_mean;
	  omega_root=-log(omega_mean);
	  
	  //if (name=="OUT_GROUP")
	  // tmp_branch_tau[branch_x]=out_tau;
	  //if (name=="OUT_GROUP")
	  // tmp_branch_omega[branch_x]=omega_mean;

	  //if (name=="OUT_GROUP")
	  // tmp_branch_lambda[branch_x]=;


	  (*event_stream) <<x_down(branch_x)<<"\t"<<tmp_branch_delta[branch_x]<<"\t"<<tmp_branch_tau[branch_x]<<"\t"<<tmp_branch_lambda[branch_x]<<"\t"<<tmp_branch_omega[branch_x]<< " " << name<<endl; 
	}
    }
  for (int x=0;x<lin_size;x++)
    {
      int branch_x=x_down(x);
      branch_delta[x]=tmp_branch_delta[branch_x];
      branch_tau[x]=tmp_branch_tau[branch_x];
      branch_lambda[x]=tmp_branch_lambda[branch_x];
      branch_omega[x]=tmp_branch_omega[branch_x];
    }

  Tef_var.clear();
  if (pair_sd>0 || 1)
    {
      (*event_stream)  << "#to_branch_id" << "\t" << "from_branch_id" << "\t" << "Tef_factor"<<endl; 
      vector <Node *> nodes=tree->getNodes();
      for (int from_x=0;from_x<lin_size;from_x++)
	for (int to_x=0;to_x<lin_size;to_x++)    
	  if (lin_slice[from_x]==lin_slice[to_x])
	    if (Tef_var.count(node_down(to_x))==0 || Tef_var[node_down(to_x)].count(node_down(from_x))==0)
	      {
		if  (pair_sd>0)
		  Tef_var[node_down(to_x)][node_down(from_x)]=RandomTools::randGamma(1./pair_sd/pair_sd,1./pair_sd/pair_sd);
		else
		  Tef_var[node_down(to_x)][node_down(from_x)]=1.;	  
		(*event_stream)  << x_down(to_x) << "\t" << x_down(from_x) << "\t" << Tef_var[node_down(to_x)][node_down(from_x)]<<endl; 
	      }
    }

  for (int x=0;x<lin_size;x++)
    branch_2_x[x_down(x)]=x;
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
  //XX
  for (int i = 0;i<N_slices+tree->getNumberOfLeaves();i++)
    {
      sim_O_profile.push_back(0.);
    }
}


void Species_tree::G_strip_virtual_leaves(tree_type * G_tree, bool addxL)
{
  bool striped=false;
  G_root=G_tree->getRootNode();
  while (!striped)
    {
      vector <Node * > leaves = TreeTemplateTools::getLeaves(*G_root);
      vector<Node *>::iterator it;	
      if (G_tree->getNumberOfNodes()==1)
	break;
      for ( it = leaves.begin(); it!=leaves.end(); it++ )
	{
	  if ((*it)!=G_root && (*it)->getName()[0]=='*' )
	    {
	      Node * father =  (*it)->getFather();
	      father->removeSon((*it));
	      //delete node;
	      //remove_leaf((*it));
	    }
	}
      striped=true;
      for ( it = leaves.begin(); it!=leaves.end(); it++ )
      	if ((*it)!=G_root && (*it)->getName()[0]=='*' )
	  striped=false;
    }
  while ( G_root->getNumberOfSons()==1 )
    {
      Node * old_root = G_root;
      //XX
      G_root=G_root->getSon(0);
      if (addxL)
	{
	  add_xL(old_root,G_root);     
	  string name = (* (dynamic_cast<const BppString *>(old_root->getBranchProperty("ID")))).toSTL();
	  string append="";
	  vector<string> tokens;
	  Tokenize(name,tokens,"&");
	  for (int i=1;i<tokens.size();i++)
	    {
	      append+="&"+tokens[i];
	    }
	  name = (* (dynamic_cast<const BppString *>(G_root->getBranchProperty("ID")))).toSTL();
	  G_root->setBranchProperty("ID",BppString(name+append));
      
	}
      string name = (* (dynamic_cast<const BppString *>(G_root->getBranchProperty("ID")))).toSTL();
      G_tree->rootAt(G_root);
      G_root->setBranchProperty("ID",BppString(name));
      //remove_leaf(old_root);
      Node * father =  old_root->getFather();
      father->removeSon(old_root);

    }
  strip(G_root,addxL);
}

void Species_tree::strip(Node * node,bool addxL)
{
  vector<Node * > sons = node->getSons();
  for( int i = 0; i < node->getNumberOfSons(); i++ )
    strip(sons[i],addxL);
  // removes Sp nodes
  if (root != node && !(node->isLeaf()) && isSpNode(node))
    remove_node(node,addxL);
  return;
}

// remove a S prime type node from tree
void Species_tree::remove_node(Node * node,bool addxL)
{
  vector<Node *> sons = node->getSons();
  for( int i = 0; i < node->getNumberOfSons(); i++ )
    if (sons[i]->getDistanceToFather()==0)
      node->removeSon(sons[i]);
  
  Node * father =  node->getFather();  
  Node * child = node->getSon(0);
  scalar_type d = TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*father,*child);
  father->removeSon(node);  
  father->addSon(child);
  if (addxL)
    add_xL(node,child);     
  child->setDistanceToFather(d);

  return;
}

void Species_tree::add_xL(Node* node,Node*child)
{      
  string node_name = (* (dynamic_cast<const BppString *>(node->getBranchProperty("ID")))).toSTL();      
  vector<string> tokens;
  Tokenize(node_name,tokens,"&");
  node_name=tokens[0];
  tokens.clear();
  Tokenize(node_name,tokens,"@");
  string event_type=tokens[0];
  string tail=tokens[1];
	    
  node_name=tokens[1];
  tokens.clear();
  Tokenize(node_name,tokens,"|");
  string event_time=tokens[0];
      
  if (event_type=="S")
    {
      if (child->isLeaf())
	{
	  string child_name =child->getName();
	  child->setName(child_name+"&"+event_type+"L"+"@"+tail);
	}
      else
	{
	  string child_name = (* (dynamic_cast<const BppString *>(child->getBranchProperty("ID")))).toSTL();      
	  child->setBranchProperty("ID",BppString(child_name+"&"+event_type+"L"+"@"+tail));
	}
    }
  if (event_type=="T")
    {	  
      int to=atoi(tokens[1].c_str());
      int from=atoi(tokens[2].c_str());
	  
      int child_branch;
      string child_name;
      if (child->isLeaf())
	child_name =child->getName();
      else
	child_name = (* (dynamic_cast<const BppString *>(child->getBranchProperty("ID")))).toSTL();      
      tokens.clear();
      Tokenize(child_name,tokens,"&");
      child_name=tokens.back();
      if (string::npos == child_name.find("@"))
	{
	  child_branch=x_down(ur_sigma[child_name]);
	}
      else
	{
	  tokens.clear();
	  Tokenize(child_name,tokens,"|");
	  child_name=tokens[0];
	  tokens.clear();
	  Tokenize(child_name,tokens,"@");
	  child_branch=atoi(tokens[1].c_str());
	}
      if (to==child_branch)
	{
	  if (child->isLeaf())
	    {
	      string child_name =child->getName();
	      child->setName(child_name+"&"+event_type+"L"+"@"+tail);
	    }
	  else
	    {
	      string child_name = (* (dynamic_cast<const BppString *>(child->getBranchProperty("ID")))).toSTL();      
	      child->setBranchProperty("ID",BppString(child_name+"&"+event_type+"L"+"@"+tail));
	    }
	}
    }
}

  // strip S prime type nodes from tree
  void Species_tree::G_strip(Node * node, tree_type * G_tree)
  {
    G_root=G_tree->getRootNode();  
    vector<Node * > sons = node->getSons();
    for( int i = 0; i < node->getNumberOfSons(); i++ )
      strip(sons[i]);
    // removes Sp nodes
    if (G_root != node && !(node->isLeaf()) && isSpNode(node))
      remove_node(node);
    return;
  }
