#include "DTL.h"
#include "DTL_mpi.h"
namespace mpi = boost::mpi;
using namespace std;
using namespace bpp;


pair<Species_tree *,string > sim_init_LL_mpi(string tree_file,string sim_file,int outgroup,
		 const mpi::communicator  world)
{
  //########################## GLOBAL #########################
  string S_tree_string;
  tree_type* S;
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  scalar_type delta,tau,lambda,omega,delta_var,tau_var,lambda_var,omega_var,noise;
  int N;
  long seed;
  vector < vector <scalar_type> > gather_sim_O_counts;
  vector < vector <scalar_type> > gather_sim_D_counts;
  vector < vector <scalar_type> > gather_sim_T_counts;
  vector < vector <scalar_type> > gather_sim_L_counts;
  vector < vector <scalar_type> > gather_sim_O_profile;
  //########################## GLOBAL #########################

  if (rank==server)
    {
      seed=good_seed();
      cout << "#SEED:" << seed <<endl;
    }
  broadcast(world,seed,server);
  RandomTools::setSeed(seed*rank+rank);
  if (rank == server) 
    {       
      string sim_file_name=sim_file;
      string tree_file_name=tree_file;
      vector <string> forest,codes;
      string line;
      int i=0;
      ifstream file_stream (sim_file_name.c_str());
      scalar_type alpha=0.5;
      scalar_type rho=0.5;
      scalar_type height=0.5;
      noise=0;
      omega=0.5;

      delta_var=0.;
      tau_var=0.;
      lambda_var=0.;
      omega_var=0.;
      N=20;
      if (file_stream.is_open())  //  ########## read sim params ############
	{
	  while (! file_stream.eof() )
	    {
	      getline (file_stream,line);
	      if (line.find("#N")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  N=atoi(tokens[1].c_str());		    
		}

	      if (line.find("#noise")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  noise=atof(tokens[1].c_str());		    
		}
	      
	      if (line.find("#ALPHA")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  alpha=atof(tokens[1].c_str());		    
		}
	      if (line.find("#RHO")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  rho=atof(tokens[1].c_str());		    
		}
	      if (line.find("#HEIGHT")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  height=atof(tokens[1].c_str());		    
		}
	      if (line.find("#OMEGA")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  omega=atof(tokens[1].c_str());		    
		}
	      if (line.find('#VAR')!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  delta_var=atof(tokens[1].c_str());		    
		  tau_var=atof(tokens[1].c_str());		    
		  lambda_var=atof(tokens[1].c_str());		    
		}
	      if (line.find("#OMEGAVAR")!=line.npos)
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  omega_var=atof(tokens[1].c_str());		    
		}

	    }
	  file_stream.close();
	}
      delta=(1-alpha)*(rho)*height;
      tau=(alpha)*(rho)*height;
      lambda=(1-rho)*height;
      if (rank==server) cout <<"D,T,L,O,N: "<< delta << " " << tau << " " << lambda << " " << omega << " " << N <<endl;
      if (rank==server) cout <<"Dv,Tv,Lv,Ov,NOISE: "<< delta_var << " " << tau_var << " " << lambda_var << " " << omega_var<<" " << noise<<endl;

      
      string S_string="";
      file_stream.open(tree_file_name.c_str());
      if (file_stream.is_open())  //  ########## read an input S tree ############
	{
	  while (! file_stream.eof() )
	    {
	      getline (file_stream,line);
	      if (line.find('(')!=line.npos && S_string=="")
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  for (int i=0;i<tokens.size();i++)
		    if (tokens[i].find('(')!=tokens[i].npos)		
		      S_string=tokens[i];
		}
	    }     
	  file_stream.close();
	}
      //########################## DISTRIBUTE ######################
      S_tree_string=S_string;
    }
  //########################## DISTRIBUTE ######################
  // SERVER CHECK POINT 1 - client trees sent -
  
  broadcast(world, N, server);   
  broadcast(world, delta, server);   
  broadcast(world, tau, server);   
  broadcast(world, lambda, server);   
  broadcast(world, omega, server);   
  broadcast(world, delta_var, server);   
  broadcast(world, tau_var, server);   
  broadcast(world, lambda_var, server);   
  broadcast(world, omega_var, server);   
  broadcast(world, noise, server);   
  
  broadcast(world, S_tree_string, server);   

  //########################## infer_tree init #######################
  
  S = TreeTemplateTools::parenthesisToTree(S_tree_string); 
  tree_type * So=S->clone(); 
  Node * root=So->getRootNode();
  double stem_len=0.1;
  if (outgroup || true)//  ########## add outgroup to root of S tree ############
    {
      Node * new_root = new Node();  
      Node * out_group = new Node();  
      root->addSon(new_root);
      new_root->addSon(out_group);      
      out_group->setDistanceToFather(1.+stem_len);
      out_group -> setName("OUT_GROUP");
      So->rootAt(new_root);
      Node * root=So->getRootNode();
      root->getSon(0)->setDistanceToFather(stem_len);
      root->getSon(1)->setDistanceToFather(stem_len);
    }
  Species_tree * infer_tree = new Species_tree(So);

  infer_tree->init_x();   

  //cout << "rank " << rank << " " << delta << " " << tau << " " << lambda << endl; 
  infer_tree->init_sim(0.1, delta, tau, lambda, delta_var*delta, tau_var*tau, lambda_var*lambda,omega,omega*omega_var,0.);      
  stringstream out;
  out<< rank << "trees.trees";
  ofstream trees_file(out.str().c_str());    
  infer_tree->event_stream=&trees_file;

  vector <string> client_trees;
  vector <string> client_codes;
  vector <string> client_stream;



  if (rank!=server)
    {
      infer_tree->emit_G_trees(N,&client_trees,&client_codes,&client_stream);         
      for (int i=0;i<client_trees.size();i++)
	{
	  if (noise>0)
	    {
	      scalar_type p0 = exp(-noise); 
	      scalar_type p=1;
	      int k=-1;
	      bool stop=false;
	      while(!stop)
		{
		  k++;
		  scalar_type u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
		  p*=u;
		  if (p<=p0)
		    stop=true;
		}
	      for (int nni=0;nni<k;nni++)
		{
		  vector <string> allnnis=all_NNIs(client_trees[i]);
		  client_trees[i]=allnnis[RandomTools::giveIntRandomNumberBetweenZeroAndEntry(allnnis.size())];
		}
	    }
	}
      infer_tree->unrooted_register(client_trees);  
      for (int i=0;i<client_trees.size();i++)
	infer_tree->codes[client_trees[i]]=client_codes[i];  

    }
  //########################## infer_tree init #######################
  if (noise==0)
    {
      infer_tree->world=world;
      infer_tree->branchwise=1;
      infer_tree->count_labels(client_trees);  
      scalar_type pea_count = sum_counts(infer_tree,world);
      for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves();i++)
	{
	  infer_tree->sim_O_count[i]=infer_tree->O_counts[i];
	  infer_tree->sim_D_count[i]=infer_tree->D_counts[i];
	  infer_tree->sim_T_count[i]=infer_tree->T_counts[i];
	  infer_tree->sim_L_count[i]=infer_tree->branch_zero[i];
	  infer_tree->sim_branch_sum[i]=infer_tree->branch_sum[i];
	  infer_tree->sim_branch_count[i]=infer_tree->branch_count[i];
	  infer_tree->sim_branch_genome_size[i]=infer_tree->branch_genome_size[i];
	  
	}
    }
  vector < vector <string> > gather_trees;
  gather(infer_tree->world, infer_tree->in_trees, gather_trees, server);    
  if (rank==server)
    {
      vector <string> forest;
      //GG
	for (int j=1;j<size;j++)
	  for (vector<string>::iterator tit=gather_trees[j].begin();tit!=gather_trees[j].end();tit++)
	    forest.push_back((*tit));
	infer_tree->panj_string= PANJ(S_tree_string,forest);
	cout << "#PANJ "<< infer_tree->panj_string <<endl;
	cout << "#PANJ "<< TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(S_tree_string),*TreeTemplateTools::parenthesisToTree(infer_tree->panj_string)) << endl;
      }

  pair<Species_tree *,string > return_pair;
  return_pair.first=infer_tree;      
  return_pair.second=S_tree_string;       
  return return_pair;
  //########################## CLIENT #########################

}





pair<Species_tree *,string > init_LL_mpi(string tree_file,string forest_file,int outgroup,
		 const mpi::communicator  world)
{
  //########################## GLOBAL #########################
  tree_type * S;
  string S_tree_string;
  string panj_string;
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  vector <string> client_trees;
  vector <string> client_codes;
  vector < vector <string> > scatter_trees;
  vector < vector <string> > scatter_codes;
  //########################## GLOBAL #########################

  if (rank == server) 
    { 
    //########################## SERVER #########################
      string forest_file_name=forest_file;
      string tree_file_name=tree_file;
      vector <string> forest,codes;
      string line;
      int i=0;
      ifstream file_stream (forest_file_name.c_str());
      if (file_stream.is_open())  //  ########## read a forest ############
	{
	  while (! file_stream.eof() )
	    {
	      getline (file_stream,line);
	      if (line.find('(')!=line.npos)
		{
		  i++;
		  forest.push_back(line);
		}
	      else if (line.find('#')!=line.npos)
		codes.push_back(line);
	    }
	  file_stream.close();
	}
      string S_string="";
      file_stream.open(tree_file_name.c_str());
      if (file_stream.is_open())  //  ########## read an input S tree ############
	{
	  while (! file_stream.eof() )
	    {
	      getline (file_stream,line);
	      if (line.find('(')!=line.npos && S_string=="")
		{
		  vector<string> tokens;
		  Tokenize(line,tokens," ");
		  for (int i=0;i<tokens.size();i++)
		    if (tokens[i].find('(')!=tokens[i].npos)		
		      S_string=tokens[i];
		}
	    }     
	  file_stream.close();
	}
      //########################## DISTRIBUTE ######################
      /// SAMPLE
      vector <int> lottery;
      for (int i=0;i<forest.size();i++)
	lottery.push_back(i);
      vector <int> winners;
      winners.resize(forest.size()/10.);
      RandomTools::getSample(lottery,winners);
      vector<string> sample_forest;
      vector<string> sample_codes;
      for (int i=0;i<winners.size();i++)
	{
	  int j=winners[i];
	  sample_forest.push_back(forest[j]);
	  sample_codes.push_back(codes[j]);
	}
      forest.clear();
      codes.clear();
      for (int i=0;i<winners.size();i++)
	{
	  forest.push_back(sample_forest[i]);
	  codes.push_back(sample_codes[i]);
	}
      /// SAMPLE

      S_tree_string=S_string;
      panj_string= PANJ(S_string,forest);
      cout << "#PANJ "<<panj_string <<endl;
      cout << "#PANJ "<< TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(S_string),*TreeTemplateTools::parenthesisToTree(panj_string)) << endl;

      map <scalar_type,string> sort_trees;
      map <scalar_type,string> sort_codes;
      
      for (int i=0;i<forest.size();i++)
	{
	  tree_type * S=TreeTemplateTools::parenthesisToTree(forest[i]);
	  scalar_type treesize=S->getNumberOfLeaves();
	  while ( sort_trees.count(treesize) )
	    treesize+=RandomTools::giveRandomNumberBetweenZeroAndEntry(1e-2);
	  sort_trees[treesize]=forest[i];
	  sort_codes[treesize]=codes[i];
	  delete S;
	  
	}
      int new_i=0;
      for ( map <scalar_type,string> ::iterator it=sort_trees.begin();it!=sort_trees.end();it++)
	{
	  forest[new_i]=it->second;
	  codes[new_i]=sort_codes[it->first];
	  new_i+=1;
	}
      map <int,int> large_trees;
      for (int i=0;i<forest.size();i++)
	{
	  tree_type * S=TreeTemplateTools::parenthesisToTree(forest[i]);
	  if  (S->getNumberOfLeaves()>150)
	    {
	      large_trees[i]=S->getNumberOfLeaves();
	    }
	  //large_trees[S->getNumberOfLeaves()]++;
	  delete S;
	}



      //int cum=0;
      //for ( map <int,int> ::iterator it=large_trees.begin();it!=large_trees.end();it++)
      //	{
      //	  cum+=(*it).second;
      //	  cout << (*it).first << " " <<  (*it).second << " " << cum << endl;
      //	}
      map <int,int> scatter_sizes;
      for (int j=0;j<size;j++)
	{ 
	  vector <string> tmp;
	  scatter_trees.push_back(tmp);
	  vector <string> tmp2;
	  scatter_codes.push_back(tmp2);
	}
      //cout << "distribute" <<endl;
      
      for (int i=0;i<forest.size();i++)
	{
	  int j = i%(size-1);
	  if (large_trees.count(i)==0)
	    {
	      scatter_trees[j+1].push_back(forest[i]);
	      scatter_codes[j+1].push_back(codes[i]);
	      tree_type * S=TreeTemplateTools::parenthesisToTree(forest[i]);
	      scatter_sizes[j]+=S->getNumberOfLeaves();
	      delete S;
	    }
	  
	}
      for ( map <int,int> ::iterator it=scatter_sizes.begin();it!=scatter_sizes.end();it++)
	cout << it->first << " " << it->second << endl;
      //cout << "distribute" <<endl;

      //.optimization could come here
      //########################## DISTRIBUTE ######################

      // SERVER CHECK POINT 1 - client trees sent -
      broadcast(world, S_tree_string, server);   
      broadcast(world, panj_string, server);   
      
      scatter(world, scatter_trees, client_trees, server);
      scatter(world, scatter_codes, client_codes, server);

      //########################## count_tree init #######################
      S = TreeTemplateTools::parenthesisToTree(S_tree_string); 
      tree_type * So=S->clone(); 
      Node * root=So->getRootNode();
      double stem_len=0.1;
      //!!
      if (outgroup || true)//  ########## add outgroup to root of S tree ############
	{
	  Node * new_root = new Node();  
	  Node * out_group = new Node();  
	  root->addSon(new_root);
	  new_root->addSon(out_group);      
	  out_group->setDistanceToFather(1.+stem_len);
	  out_group -> setName("OUT_GROUP");
	  So->rootAt(new_root);
	  Node * root=So->getRootNode();
	  if (root->getSon(0)->isLeaf())
	    root->getSon(1)->setDistanceToFather(stem_len);
	  else
	    root->getSon(0)->setDistanceToFather(stem_len);
	}
      Species_tree * count_tree = new Species_tree(So);
      count_tree->panj_string=panj_string;

      count_tree->init_x();   
      //########################## count_tree init #######################
      // SERVER CHECK POINT 1 - client trees sent -
      count_tree->world=world;
      pair<Species_tree *,string > return_pair;
      return_pair.first=count_tree;      
      return_pair.second=S_tree_string;       
      return return_pair;
  //########################## SERVER #########################
    }
  else
    {
  //########################## CLIENT #########################

      // CLIENT CHECK POINT 1 - client trees recived -
      broadcast(world, S_tree_string, server);   
      broadcast(world, panj_string, server);   

      scatter(world, scatter_trees, client_trees, server);
      scatter(world, scatter_codes, client_codes, server);

          //########################## infer_tree init #######################

      S = TreeTemplateTools::parenthesisToTree(S_tree_string); 
      tree_type * So=S->clone(); 
      Node * root=So->getRootNode();
      double stem_len=0.1;
      if (outgroup || true)//  ########## add outgroup to root of S tree ############
	{
	  Node * new_root = new Node();  
	  Node * out_group = new Node();  
	  root->addSon(new_root);
	  new_root->addSon(out_group);      
	  out_group->setDistanceToFather(1.+stem_len);
	  out_group -> setName("OUT_GROUP");
	  So->rootAt(new_root);
	  Node * root=So->getRootNode();
	  root->getSon(0)->setDistanceToFather(stem_len);
	  root->getSon(1)->setDistanceToFather(stem_len);
	}
      Species_tree * infer_tree = new Species_tree(So);
      infer_tree->panj_string=panj_string;     
      infer_tree->init_x();   
      infer_tree->unrooted_register(client_trees);  
      
      for (int i=0;i<client_trees.size();i++)
	infer_tree->codes[client_trees[i]]=client_codes[i];
  
      //########################## infer_tree init #######################
      // CLIENT CHECK POINT 1 - client trees recived -
      infer_tree->world=world;
      pair<Species_tree *,string > return_pair;
      return_pair.first=infer_tree;      
      return_pair.second=S_tree_string;       
      return return_pair;
      //########################## CLIENT #########################
    }
}


scalar_type LL_mpi(Species_tree * infer_tree,string Sstring,string mode,bool outgroup,scalar_type delta,scalar_type tau, scalar_type lambda)
{ 

  const mpi::communicator  world = infer_tree->world;
  int server = 0;  
  int rank = world.rank();
  int size = world.size();
  scalar_type ll;
  scalar_type return_ll; 

  
  if (rank==server) cout << mode << " " << Sstring <<endl; 

  broadcast(world, Sstring,server);
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);

  if (outgroup)
    {
      Node * root=S->getRootNode();
      
      vector <Node*> nodes=S->getNodes();
      for (vector<Node * > :: iterator it=nodes.begin();it!=nodes.end();it++)
	if ((*it)->hasFather())
	  if ((*it)->getDistanceToFather()==0)
	    (*it)->setDistanceToFather(0.001);
      
      infer_tree->node_makeclocklike_height((S->getRootNode()));    
      S->scaleTree(1./TreeTemplateTools::getDistanceBetweenAnyTwoNodes(*root,*(S->getLeaves()[0])));      
      double stem_len=0.1;
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
      //root->getSon(1)->setDistanceToFather(stem_len);
      infer_tree->reconstruct(S);      
    }
  vector <string> modes;
  Tokenize(mode,modes,"_");
  if (modes[0]=="clock" )
    {
      
      infer_tree->intp=0;
      infer_tree->w=0.;
      infer_tree->reset();
      infer_tree->branchwise=0;      
      string clock_string=Sstring;
      if (rank==server) cout << modes[0] << " " << modes[1] <<endl; 
      clock_string=DTL_clock(clock_string,infer_tree,modes[1]);
      //ll=LL_mpi(infer_tree, clock_string, "Hest",outgroup,delta,tau, lambda);
      ll=infer_tree->tmp_ll;
      infer_tree->tmp_ll=ll;
      delete S;
      if (rank==server) cout <<"$$  "<< clock_string <<endl; 
      infer_tree->tmp_Sstring=clock_string;
      infer_tree->tmp_ll=ll;
      return ll;
    }
  if (mode=="Nest" )
    {
      infer_tree->intp=0;
      infer_tree->w=0.;
      infer_tree->reset(delta,tau,lambda);
      infer_tree->branchwise=0;      
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");
      //ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
      infer_tree->omega_avg=infer_tree->branchwise_omega[0]+infer_tree->branchwise_omega[1]+infer_tree->branchwise_omega[2];
      delete S;
      infer_tree->tmp_ll=ll;
      return ll;
    }
  if (mode=="noestimate" )
    {
      ll= oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;
      delete S;
      return ll;
    }
  if (mode=="Pest" )
    {
      infer_tree->intp=0.9;
      infer_tree->w=0.;
      // !! no T
      infer_tree->reset();
      infer_tree->branchwise=0;      
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
      scalar_type pea_count = sum_counts(infer_tree,world);
      delete S;
      infer_tree->tmp_ll=ll;
      return ll;
    }
  if (mode=="Test" )
    {
      infer_tree->intp=0.9;
      infer_tree->w=0.;
      // !! no T
      infer_tree->reset();
      infer_tree->branchwise=0;      
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
      scalar_type pea_count = sum_counts(infer_tree,world);
      delete S;
      scalar_type Tll=0;
      for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	Tll+=infer_tree->T_counts[i];      
      infer_tree->tmp_ll=-Tll;
      return -Tll;
    }

  if (mode=="FASTest" )
    {
      infer_tree->intp=0.9;
      infer_tree->w=0.;
      // !! no T
      infer_tree->reset();
      infer_tree->branchwise=0;      
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
      scalar_type pea_count = sum_counts(infer_tree,world);
      infer_tree->apparent_rates(pea_count,true,true);

      //if (rank==server) cout<<"#p "<< ll  << " " << " "<< infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");
      infer_tree->tmp_ll=ll;
      infer_tree->omega_avg=infer_tree->branchwise_omega[0]+infer_tree->branchwise_omega[1]+infer_tree->branchwise_omega[2];
      delete S;
      return ll;
    }
  
  if (mode=="Hest" || mode=="IHest")
    {
      if (mode=="Hest")
	infer_tree->nclasses=5;
      
      // !! no T
      //tau=1e-30;
      //infer_tree->reset(0.01,1e-30,0.01);
      infer_tree->reset();
      infer_tree->intp=0.9;
      infer_tree->w=0.;
      infer_tree->branchwise=0;      	  
      ll=oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
      infer_tree->tmp_ll=ll;
      scalar_type pea_count = sum_counts(infer_tree,world);
      if (rank==server) cout<<"#p "<< ll  << " " << " "<< infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg <<" "<<pea_count << endl ;

      infer_tree->apparent_rates(pea_count,mode=="Hest",true);	
      
      //if (rank==server) infer_tree->print_rates();

      // ### sum events       
      scalar_type ll2=ll-2;
      while (ll>ll2 && abs(ll-ll2)>1 )
	{
	  ll2=ll;
	  // given rates estimate p_omega
	  ll=oLL(infer_tree,S,delta,tau,lambda,1,"estimate");
	  infer_tree->tmp_ll=ll;	 
	  if (ll2>ll)
	    {
	      infer_tree->recall_rates();	  
	      //scalar_type re_ll=oLL(infer_tree,S,delta,tau,lambda,1,"noestimate") ;	  
	      ll=ll2;
	      break;
	    }	  
	  infer_tree->omega_avg=infer_tree->branchwise_omega[0]+infer_tree->branchwise_omega[1]+infer_tree->branchwise_omega[2];
	  // given p_omega estimate rates
	  // formal position of remember 
	  oLL(infer_tree,S,delta,tau,lambda,1,"intbck");
	  scalar_type pea_count = sum_counts(infer_tree,world);
	  infer_tree->apparent_rates(pea_count,mode=="Hest",true);
	  
	  if (rank==server) cout<<"#e "<< ll  << " " << " "<< infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl ;
	  
	  // ### sum events      
	}
      return_ll=ll;
      infer_tree->tmp_ll=ll;
      delete S;
      return return_ll;
    }
}

void sum_ttf(Species_tree * infer_tree)
{
  int server = 0;  
  int rank = infer_tree->world.rank();
  int size = infer_tree->world.size();
  vector<  vector< vector<scalar_type> > > gather_Ttf;
  gather(infer_tree->world, infer_tree->Ttf, gather_Ttf, server);       
  //########################## SERVER #########################
  if (rank == server) 
    {      
      for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	for (int k = 0;k<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;k++)
	  infer_tree->Ttf[i][k]=0;
      for (int j=1;j<size;j++)
	for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	  for (int k = 0;k<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;k++)
	    infer_tree->Ttf[i][k]+=gather_Ttf[j][i][k];
    }  
  
}
scalar_type sum_counts(Species_tree * infer_tree,const mpi::communicator  world)
{
  //########################## GLOBAL #########################
  int server = 0;  
  int rank = world.rank();
  int size = world.size();

  scalar_type client_ll,broadcast_ll;
  vector<scalar_type> gather_ll;
  scalar_type client_pea_count;
  scalar_type broadcast_pea_count;
  vector<scalar_type> gather_pea_count;
  scatter_scalar gather_O_counts;
  scatter_scalar gather_D_counts;
  scatter_scalar gather_T_counts;
  scatter_scalar gather_L_counts;
  scatter_scalar gather_B_counts;
  scatter_scalar gather_S_counts;
  scatter_scalar gather_F_counts;
  scatter_scalar gather_G_counts;
  //vector<  vector< vector<scalar_type> > > gather_Ttf;
  
  vector<scalar_type> broadcast_O_counts;
  vector<scalar_type> broadcast_D_counts;
  vector<scalar_type> broadcast_T_counts;
  vector<scalar_type> broadcast_branch_zero;
  vector<scalar_type> broadcast_branch_count;
  vector<scalar_type> broadcast_branch_sum;
  vector<scalar_type> broadcast_from_count;
  vector<scalar_type> broadcast_genome_size;
  int start=1;
  //########################## GLOBAL #########################

  //########################## SERVER #########################
  //########################## SERVER #########################
  if (rank != server) 
    {
      //########################## CLIENT #########################
      client_pea_count=infer_tree->tmp_pea_sum;
      client_ll=infer_tree->tmp_ll;
      //########################## CLIENT #########################

    }
  //########################## COMMON #########################
  gather(world, infer_tree->O_counts, gather_O_counts, server);           
  gather(world, infer_tree->D_counts, gather_D_counts, server);           
  gather(world, infer_tree->T_counts, gather_T_counts, server);           
  gather(world, infer_tree->branch_zero, gather_L_counts, server);           
  gather(world, infer_tree->branch_count, gather_B_counts, server);                 
  gather(world, infer_tree->branch_sum, gather_S_counts, server);              
  gather(world, infer_tree->from_count, gather_F_counts, server);       
  gather(world, infer_tree->branch_genome_size, gather_G_counts, server);       

  //TTF
  // gather(world, infer_tree->Ttf, gather_Ttf, server);       
       
  gather(world, client_pea_count, gather_pea_count, server);                 
  gather(world, client_ll, gather_ll, server);      
  //########################## COMMON #########################
  if (rank == server) 
    {

      //########################## SERVER #########################
      for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	{
	  broadcast_O_counts.push_back(0);
	  broadcast_D_counts.push_back(0);
	  broadcast_T_counts.push_back(0);
	  broadcast_branch_zero.push_back(0);
	  broadcast_branch_count.push_back(0);
	  broadcast_branch_sum.push_back(0);
	  broadcast_from_count.push_back(0);      
	  broadcast_genome_size.push_back(0);      
	}
      scalar_type tmp_O=0;
      scalar_type tmp_D=0;
      scalar_type tmp_T=0;
      scalar_type tmp_z=0;
      scalar_type tmp_c=0;
      scalar_type tmp_s=0;
      scalar_type tmp_f=0;
      //cout << " start summing - sum " <<endl;     
      //for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
      //	for (int k = 0;k<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;k++)
      //	  infer_tree->Ttf[i][k]=0;
      for (int j=1;j<size;j++)
	{
	  for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	    {	      
	      //for (int k = 0;k<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;k++)
	      //	infer_tree->Ttf[i][k]+=gather_Ttf[j][i][k];
	      broadcast_O_counts[i]+=gather_O_counts[j][i];
	      broadcast_D_counts[i]+=gather_D_counts[j][i];
	      broadcast_T_counts[i]+=gather_T_counts[j][i];
	      broadcast_branch_zero[i]+=gather_L_counts[j][i];
	      broadcast_branch_count[i]+=gather_B_counts[j][i];
	      broadcast_branch_sum[i]+=gather_S_counts[j][i];	 
	      broadcast_from_count[i]+=gather_F_counts[j][i];	 
	      broadcast_genome_size[i]+=gather_G_counts[j][i];	 	      
	      tmp_O+=gather_O_counts[j][i];
	      tmp_D+=gather_D_counts[j][i];
	      tmp_T+=gather_T_counts[j][i];
	      tmp_z+=gather_L_counts[j][i];
	      tmp_c+=gather_B_counts[j][i];
	      tmp_s+=gather_S_counts[j][i];
	      tmp_f+=gather_F_counts[j][i];	     
	    }
	}

      /*
      for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
	{
	  cout << "tag "<<i;
	  cout << " " <<  broadcast_O_counts[i];
	  cout << " " <<  broadcast_D_counts[i];
	  cout << " " <<  broadcast_T_counts[i];
	  cout << " " <<  broadcast_branch_zero[i];
	  cout << " " <<  broadcast_branch_sum[i];
	  cout << " " <<  broadcast_branch_count[i];
	  cout << " " <<  broadcast_from_count[i];
	  cout << endl;
	}
      */
      //if (rank==server) cout << rank << " : "<<tmp_O << " " << tmp_D << " " << tmp_T << " " << tmp_z << " " << tmp_c << " " << tmp_s << " " << tmp_f << endl; 
      broadcast_pea_count=0.;
      broadcast_ll=0.;
      //cout << " start summing - llsum " << endl;      

      for (int j=1;j<size;j++)	
	{
	  broadcast_pea_count+=gather_pea_count[j];      
	  broadcast_ll+=gather_ll[j];
	}
      // ### sum events
      //########################## SERVER #########################
      //cout << " end summing"    << endl;   

    }
  //########################## CLIENT #########################
  //########################## CLIENT #########################

  //########################## COMMON #########################  
broadcast(world, broadcast_O_counts,server);
  broadcast(world, broadcast_D_counts,server);
  broadcast(world, broadcast_T_counts,server);
  broadcast(world, broadcast_branch_zero,server);
  broadcast(world, broadcast_branch_count,server);
  broadcast(world, broadcast_branch_sum,server);
  broadcast(world, broadcast_from_count,server);
  broadcast(world, broadcast_genome_size,server);
  broadcast(world, broadcast_pea_count,server);
  broadcast(world, broadcast_ll,server);
  //########################## COMMON #########################

  //########################## SEVER #########################
  //########################## CLIENT #########################
  //  if (rank==server)
  //  {
  for (int i = 0;i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1;i++)
    {
      infer_tree->O_counts[i]=broadcast_O_counts[i];
      infer_tree->D_counts[i]=broadcast_D_counts[i];
      infer_tree->T_counts[i]=broadcast_T_counts[i];
      infer_tree->branch_zero[i]=broadcast_branch_zero[i];
      infer_tree->branch_count[i]=broadcast_branch_count[i];
      infer_tree->branch_sum[i]=broadcast_branch_sum[i];
      infer_tree->from_count[i]=broadcast_from_count[i];
      infer_tree->branch_genome_size[i]=broadcast_genome_size[i];
    }
  //infer_tree->remember_rates();
  // }
  //infer_tree->apparent_rates(broadcast_pea_count,true,true);
  //########################## CLIENT #########################
  //########################## SERVER #########################
  
 
  gather_ll.clear();
  gather_pea_count.clear();

  for (scatter_scalar::iterator it=gather_O_counts.begin();it!=gather_O_counts.end();it++) (*it).clear();
  gather_O_counts.clear();
  for (scatter_scalar::iterator it=gather_D_counts.begin();it!=gather_D_counts.end();it++) (*it).clear();
  gather_D_counts.clear();
  for (scatter_scalar::iterator it=gather_T_counts.begin();it!=gather_T_counts.end();it++) (*it).clear();
  gather_T_counts.clear();
  for (scatter_scalar::iterator it=gather_L_counts.begin();it!=gather_L_counts.end();it++) (*it).clear();
  gather_L_counts.clear();
  for (scatter_scalar::iterator it=gather_B_counts.begin();it!=gather_B_counts.end();it++) (*it).clear();
  gather_B_counts.clear();
  for (scatter_scalar::iterator it=gather_S_counts.begin();it!=gather_S_counts.end();it++) (*it).clear();
  gather_S_counts.clear();
  for (scatter_scalar::iterator it=gather_F_counts.begin();it!=gather_F_counts.end();it++) (*it).clear();
  gather_F_counts.clear();
  for (scatter_scalar::iterator it=gather_G_counts.begin();it!=gather_G_counts.end();it++) (*it).clear();
  gather_G_counts.clear();

  broadcast_O_counts.clear();
  broadcast_D_counts.clear();
  broadcast_T_counts.clear();
  broadcast_branch_zero.clear();
  broadcast_branch_count.clear();
  broadcast_branch_sum.clear();
  broadcast_from_count.clear();
  broadcast_genome_size.clear();

  return broadcast_pea_count;

}

string strip_out(string Sstring)
{
  tree_type * S=TreeTemplateTools::parenthesisToTree(Sstring);   
  Node * root=S->getRootNode();
  vector <Node *> sons = root->getSons();
  while (sons.size()==1)
    {
      root=root->getSon(0);
      sons = root->getSons();
    }
  sons = root->getSons();
  Node * out;
  Node * other;
  if (sons[0]->isLeaf() && sons[0]->getName()=="OUT_GROUP" )
    other=sons[1];
  else
    other=sons[0];
  string tmp=TreeTemplateTools::nodeToParenthesis(*other)+";";
  delete S;
  return tmp;
}

pair<scalar_type,string> SPR_step(scalar_type pre_ll,Species_tree * infer_tree, string Sstring,string mode,bool greedy,bool root_step, bool skip_step)
{
  int server = 0;  
  int rank = infer_tree->world.rank();
  int size = infer_tree->world.size();  
  vector <string> spr_strings;
  vector <string> spr_names;
  map < string ,string> correlates;
  scalar_type oLL;
  if (pre_ll==1) oLL=LL_mpi(infer_tree,Sstring,mode);  
  else oLL=pre_ll;
  sum_ttf(infer_tree);

  pair<scalar_type,string> return_pair;
  return_pair.first=-1;
  return_pair.second=Sstring;
  if (root_step) skip_step=false;
  

  if (rank==server)
    {
      map < scalar_type, vector < int > > from;
      map < scalar_type, vector < int > > to;    

      for (int i = 0; i<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1; i++)
	for (int j = 0; j<infer_tree->N_slices+infer_tree->tree->getNumberOfLeaves()-1; j++)	  
	  {
	    to[infer_tree->Ttf[i][j]].push_back(i);
	    from[infer_tree->Ttf[i][j]].push_back(j);	      
	  }	  

      if ( root_step)
	{
      tree_type * S = TreeTemplateTools::parenthesisToTree(Sstring);

      Node * root = S->getRootNode();

      Node * left=root->getSon(0);
      Node * right=root->getSon(1);      
      scalar_type d_left=left->getDistanceToFather();
      scalar_type d_right=right->getDistanceToFather();

      Node * left_left;
      Node * left_right;     
      scalar_type d_left_left;
      scalar_type d_left_right;
      
      bool  left_is_leaf=left->isLeaf();
      bool  right_is_leaf=right->isLeaf();

      if (!left_is_leaf)
	{
	  left_left=left->getSon(0);
	  left_right=left->getSon(1);     
	  d_left_left=left_left->getDistanceToFather();
	  d_left_right=left_right->getDistanceToFather();
	}
      
      Node * right_right;
      Node * right_left;
      scalar_type d_right_left;
      scalar_type d_right_right;
      if (!right_is_leaf)
	{
	  right_right=right->getSon(1);     
	  right_left=right->getSon(0);
	  d_right_right=right_right->getDistanceToFather();
	  d_right_left=right_left->getDistanceToFather();
	}
      	

      if (!left_is_leaf)
	{
	  root->removeSon(left);
	  root->removeSon(right);

	  left->removeSon(left_left);
	  left->addSon(right);
	  
	  root->addSon(left);
	  root->addSon(left_left);

	  left_left->setDistanceToFather(d_left_left/2.);
	  left->setDistanceToFather(d_left_left/2.);
	  right->setDistanceToFather(d_left+d_right);
	  //cout << TreeTemplateTools::treeToParenthesis( * S);
	  string root_tmp=TreeTemplateTools::treeToParenthesis( * S);
	  spr_strings.push_back(root_tmp);	
	  correlates[root_tmp]="#root";

	  root->removeSon(left);
	  root->removeSon(left_left);	 
 	  left->addSon(left_left);

	  left->removeSon(left_right);
	 
	  root->addSon(left);
	  root->addSon(left_right);
 
	  left_right->setDistanceToFather(d_left_right/2.);
	  left->setDistanceToFather(d_left_right/2.);
	  root_tmp=TreeTemplateTools::treeToParenthesis( * S);
	  spr_strings.push_back(root_tmp);		    
	  correlates[root_tmp]="#root ";

	  left->removeSon(right);
	  root->addSon(right);

	  root->removeSon(left_right);	  	  
	  left->addSon(left_right);
	}

      if (!right_is_leaf)
	{
	  root->removeSon(left);
	  root->removeSon(right);

	  right->removeSon(right_left);
	  right->addSon(left);

	  root->addSon(right);
	  root->addSon(right_left);

	  right_left->setDistanceToFather(d_right_left/2.);
	  right->setDistanceToFather(d_right_left/2.);
	  left->setDistanceToFather(d_left+d_right);
	  string root_tmp=TreeTemplateTools::treeToParenthesis( * S);
	  spr_strings.push_back(root_tmp);		    
	  correlates[root_tmp]="#root ";

	  root->removeSon(right);
	  root->removeSon(right_left);	  
	  right->addSon(right_left);

	  right->removeSon(right_right);
	 
	  root->addSon(right);
	  root->addSon(right_right);
 
	  right_right->setDistanceToFather(d_right_right/2.);
	  right->setDistanceToFather(d_right_right/2.);
	  root_tmp=TreeTemplateTools::treeToParenthesis( * S);
	  spr_strings.push_back(root_tmp);		    
	  correlates[root_tmp]="#root ";

	  right->removeSon(left);
	  root->addSon(left);
	  root->removeSon(right_right);	  	  
	  right->addSon(right_right);
	}
      delete S;		
	}
 
      for ( map < scalar_type,vector < int > >::reverse_iterator it=from.rbegin(); it!=from.rend(); it++ )
	for (int i=0;i<(*it).second.size();i++)
	  {
	    int tmp_to,tmp_from;
	    scalar_type tmp_Ttf;
	    tmp_Ttf=(*it).first;
	    tmp_from=(*it).second[i];
	    tmp_to=to[(*it).first][i];	    
	    string to_name= infer_tree->branchname[tmp_to];
	    string from_name= infer_tree->branchname[tmp_from];
	      
	    if (tmp_Ttf>0 &&  from_name!="OUT_GROUP" && to_name!="OUT_GROUP")
	      {
		tree_type * S = TreeTemplateTools::parenthesisToTree(Sstring);
		Node * root=S->getRootNode();
		infer_tree->name_internal(root);
		vector < Node *> nodes=S->getNodes();
		Node * to_node;
		Node * from_node;
		for (vector < Node *>::iterator nt=nodes.begin();nt!=nodes.end();nt++)
		    if ((*nt)->getName()==to_name)
		      to_node=(*nt);
		    else if ((*nt)->getName()==from_name)
		      from_node=(*nt);
		  
		//cout << tmp_to << " " << tmp_from << " " << tmp_Ttf <<endl;
		//cout << to_name << " " << from_name << endl;	  
		if (to_node->hasFather() && to_node->getFather()!=from_node->getFather() )
		  {
		    Node * to_father=to_node->getFather();
		    Node * other_son=to_father->getSon(0);
		    if (other_son==to_node)
		      other_son=to_father->getSon(1);
		    bool re_root=false;
		    if (to_father->hasFather())
		      {
			to_father->removeSon(to_node);
			Node * to_grandfather = to_father->getFather();
			scalar_type d_1=to_father->getDistanceToFather();
			scalar_type d_2=other_son->getDistanceToFather();
			to_father->removeSon(other_son);
			to_grandfather->removeSon(to_father);
			to_grandfather->addSon(other_son);
			other_son->setDistanceToFather(d_1+d_2);
		      }
		    else
		      {
			re_root=true;
			to_father->removeSon(to_node);
		        S->rootAt(other_son);
			other_son->removeSon(to_father);
		      }
		    
		    Node * from_father=from_node->getFather();
		    from_father->removeSon(from_node);

		    Node * new_node=new Node();
		    from_father->addSon(new_node);
		    new_node->addSon(from_node);
		    new_node->addSon(to_node);
		    
		    to_node->setDistanceToFather(1.);
		    from_node->setDistanceToFather(1.);
		    new_node->setDistanceToFather(1.);
		    string tmp_spr=TreeTemplateTools::treeToParenthesis( * S);
		    //CHECK IF WE JUST TRIED THIS RECENTLY !!

		    string  spr_name=to_name+"->"+from_name;		    
		    if (skip_step) spr_names.push_back(spr_name);		    		      
		    spr_strings.push_back(tmp_spr);		    
		    stringstream out;
		    out << tmp_Ttf << " ";
		    out << infer_tree->T_counts[tmp_to] << " ";
		    out << infer_tree->T_counts[tmp_from] << " ";
		    out << infer_tree->branch_zero[tmp_from] << " ";
		    out << infer_tree->branch_zero[tmp_to] << " ";


		    correlates[tmp_spr]=out.str();
		  }
		delete S;		
	      }
	  }      
    }
  broadcast(infer_tree->world,spr_strings,server);
  if (skip_step) broadcast(infer_tree->world,spr_names,server);
  int c=0;
  for (vector<string>::iterator st=spr_strings.begin();st!=spr_strings.end();st++)
    {
      c++;
      if ( skip_step==false || infer_tree->tried_spr.count(spr_names[c-1])==0 )
	{
	  scalar_type tryLL=LL_mpi(infer_tree,(*st),mode);
	  scalar_type dLL=tryLL-oLL;	  
	  if (rank==server) 
	    {
	      ofstream cout_stream(infer_tree->outstream_file_name.c_str(),ios::app);  
	      cout_stream<< c << " " << correlates[(*st)] << " " << dLL << " " << oLL<< " " << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree((*st)),*TreeTemplateTools::parenthesisToTree(infer_tree->REF_Sstring))<< endl; 
	      cout_stream.close();
	    }
	  if (dLL>0)
	    {
	      return_pair.first=dLL;
	      return_pair.second=(*st);
	      infer_tree->tmp_ll=tryLL;
	      if (greedy) break;	  
	    }
	  else if (greedy && skip_step)
	    infer_tree->tried_spr[spr_names[c-1]]=1;
	}
      else
	{
	  if (rank==server) {ofstream cout_stream(infer_tree->outstream_file_name.c_str(),ios::app); cout_stream<< c <<" (*) "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree((*st)),*TreeTemplateTools::parenthesisToTree(infer_tree->REF_Sstring))<<endl;  cout_stream.close();}
	}
      if ((c>60 && mode=="clock_Test" )||(c>60+root_step*4 && !skip_step)) break;
    }  
  spr_strings.clear();
  spr_names.clear();
  if (return_pair.first<=0)
    {
      infer_tree->tmp_ll=oLL;
      infer_tree->tmp_Sstring=Sstring;
    }
  return return_pair; 
  
}


