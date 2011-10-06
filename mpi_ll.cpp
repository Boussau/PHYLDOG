// From the STL:
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "DTL.h"
#include "DTL_mpi.h"

//From Boost:
//#include <mpi.h> 
#include <boost/mpi.hpp>
#include <boost/progress.hpp>
#include <NumCalc/EigenValue.h>

//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/collectives.hpp>
//#include <boost/serialization/string.hpp>
//#include <boost/mpi/packed_iarchive.hpp>
using namespace std;
using namespace bpp;
namespace mpi = boost::mpi;


string   ML_SPR(string cout_name,Species_tree * infer_tree,string Sstring, mpi::communicator world,string restart_file_name,string explore_mode="noestimate",int random_nnis=0,bool greedy = false,string S0string="",bool clock=false,bool root_step=false,bool skip_step=false)
{
  int rank = world.rank();
  int size = world.size();
  int server =0;
  scalar_type ll,trialll;
  //tree_type * S0=S->clone(); 
  if (root_step) skip_step=false;
  infer_tree->tried_spr.clear();

  if (rank==server)  
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream <<"#---------------------"<<endl;
      cout_stream <<"#SPR  " << "  exp: " << explore_mode << " greedy: " << greedy << " clock: " << clock << endl;
      cout_stream <<"#---------------------"<<endl<<endl;
      cout_stream <<"#S0 "<< Sstring <<endl;
      cout_stream <<"#S1 "<< Sstring;
      cout_stream <<"#S0 - S1 "  <<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(Sstring),*TreeTemplateTools::parenthesisToTree(S0string)) <<endl<<endl;
      cout_stream.close();
    }
  //if (rank==server) Sstring=randomize_nni(Sstring,random_nnis);

  broadcast(world, Sstring,server);

  pair <double,string > step;
  string save;

  save=Sstring;
  
  //if (clock) Sstring=DTL_clock(Sstring,infer_tree);
  
  boost::timer * t;
  if (rank==server) 
    {t=new boost::timer();}
  
  step=SPR_step(1,infer_tree,Sstring,explore_mode,greedy,root_step,skip_step);     
  if (clock) step.second=infer_tree->tmp_Sstring;
  
  if (rank==server)  
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream << "#TIME "<<t->elapsed() <<endl;
      cout_stream << "#SPR "<< step.second << " " << step.first <<endl;
      cout_stream << "#SPR "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;
      cout_stream.close();

      ofstream restart_file(restart_file_name.c_str());    
      restart_file << "#SPR "<< step.second;
      restart_file << "#SPR "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;
      restart_file << "#rate "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg << " " << infer_tree->omega_avg << endl;
      restart_file <<"#mode " << explore_mode <<endl;
      restart_file.close();
    }


  while (step.first>0)
    {
      save=step.second;
      boost::timer * t; if (rank==server) {t=new boost::timer();}
      step=SPR_step(infer_tree->tmp_ll,infer_tree,step.second,explore_mode,greedy,root_step,skip_step);           
      if (clock) step.second=infer_tree->tmp_Sstring;
      
      if (rank==server)  
	{
	  ofstream cout_stream(cout_name.c_str(),ios::app);
	  cout_stream << "#TIME "<<t->elapsed() <<endl;
	  cout_stream << "#SPR "<< step.second << " " << step.first <<endl;
	  cout_stream << "#SPR "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;
	  cout_stream.close();
	  
	  ofstream restart_file(restart_file_name.c_str());    
	  restart_file << "#SPR "<< step.second;
	  restart_file << "#SPR "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;
	  restart_file << "#rate "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg << " " << infer_tree->omega_avg << endl;
	  restart_file <<"#mode " << explore_mode  <<endl;
	  restart_file.close();
	}      
      
      
    }
  step.second=DTL_clock( step.second,infer_tree);
  if (rank==server)  
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream << "#O_SPR "<< step.second << " " << step.first <<endl;
      cout_stream.close();
    }
  return step.second;
  
}


string   ML_NNI(string cout_name,Species_tree * infer_tree,string Sstring, mpi::communicator world,string mode,string explore_mode="noestimate",int random_nnis=0,bool greedy = false,string S0string="",bool clock=true)
{
  int rank = world.rank();
  int size = world.size();
  int server =0;
  scalar_type ll,trialll;
  //tree_type * S0=S->clone(); 
  
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream <<"#---------------------"<<endl;
      cout_stream <<"#NNI  " << mode << "  exp: " << explore_mode << " greedy: " << greedy << " clock: " << clock << endl;
      cout_stream <<"#---------------------"<<endl<<endl;
      cout_stream <<"#S0 "<< Sstring <<endl;
      cout_stream.close();
    }
  //if (rank==server && random_nnis) Sstring=randomize_nni(Sstring,random_nnis);

  broadcast(world, Sstring,server);
  pair <double,string > step,trialstep;
  string save;
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream <<"#S1 "<< Sstring;
      cout_stream <<"#S0 - S1 "  <<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(Sstring),*TreeTemplateTools::parenthesisToTree(S0string)) <<endl<<endl;      
      cout_stream.close();
    }

  save=Sstring;    
  ll=LL_mpi(infer_tree,Sstring,mode);
  trialstep=ML_nni_step(infer_tree,Sstring,explore_mode,greedy);     
  trialll=LL_mpi(infer_tree,trialstep.second,mode);
  if (trialll>ll)
    {
      step.first=trialstep.first;
      step.second=trialstep.second;
    }
  else
    {
      step.first=-1;
      step.second=save;
    }
  
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream << "#NNI "<< step.second << " " << step.first <<endl;
      cout_stream << "#NNI "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;
      cout_stream.close();
    }

  while (step.first>0)
    {

      save=step.second;
      ll=LL_mpi(infer_tree,step.second,mode);
      trialstep=ML_nni_step(infer_tree,step.second,explore_mode,greedy);     
      trialll=LL_mpi(infer_tree,trialstep.second,mode);
      if (trialll>ll)
	{
	  step.first=trialstep.first;
	  step.second=trialstep.second;
	}
      else
	{
	  step.first=-1;
	  step.second=save;
	  break;
	}
      
      if (rank==server)
	{
	  ofstream cout_stream(cout_name.c_str(),ios::app);
	  cout_stream<< "#NNI "<< step.second << " " << step.first <<endl;
	  cout_stream << "#NNI "<< infer_tree->tmp_ll << " "<< step.first <<" "<<TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(step.second),*TreeTemplateTools::parenthesisToTree(S0string))  <<endl;	  
	  
	}
    }
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);      
      cout_stream << "#O_NNI "<< step.second << " " << step.first <<endl;      
    }
  return step.second;
  
}


string   ML_run(string cout_name,Species_tree * infer_tree,string Sstring, mpi::communicator world,string mode,string restart_file_name,string explore_mode="noestimate",int random_steps=0,bool greedy = false)
{

  int rank = world.rank();
  int size = world.size();
  int server =0;
  string to;
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);      
      cout_stream <<"#---------------------"<<endl;
      cout_stream <<"#TO  " << mode << "  exp: " << explore_mode << " greedy: " << greedy<<endl;
      cout_stream <<"#---------------------"<<endl; 
      cout_stream <<"#S0 "<< Sstring;
      cout_stream <<"#S0 "<< print_time_order(Sstring) <<endl;  
      cout_stream <<"#S0 "<< Sstring;
      cout_stream.close();
    }
  //if (rank==server && random_steps>0) Sstring=randomize_time_order(Sstring,random_steps,"noroot");
  broadcast(world, Sstring,server);
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);      
      cout_stream <<"#S1 "<< Sstring;
      cout_stream<<"#S1 "<< print_time_order(Sstring) <<endl;  
      cout_stream.close();
    }

  double ll,lll,trialll; 
  int count =1,c=0;
  map <string,int> end;
  pair <double,string > step,trialstep;
  
  string save=Sstring;  
  if (mode=="Hest" )
    {
      ll=LL_mpi(infer_tree,Sstring,mode);
      save=Sstring;
      trialstep=ML_step(infer_tree,Sstring,explore_mode,"order",greedy);
      trialll=LL_mpi(infer_tree,trialstep.second, mode);      
      if (trialll>ll )
	{
	  step.first=trialstep.first;
	  step.second=trialstep.second;
	}
      else
	{
	  step.first=-1;
	  step.second=save;
	}

      if (rank==server)
	{
	  ofstream cout_stream(cout_name.c_str(),ios::app);
	  cout_stream << "#LTO "<< step.second << " " << step.first <<endl;
	  cout_stream << "#LTO "<< infer_tree->tmp_ll << " "<< step.first  <<endl;
	  cout_stream.close();
	}

    
      while (step.first>0 || c==1)
	{
	  count++;
	  //step=ML_step(infer_tree,step.second,mode,"root",greedy);
	  
	  to=print_time_order(step.second);
	  if (rank==server)
	    {
	      ofstream restart_file(restart_file_name.c_str());    
	      restart_file << "#M "<< step.second;
	      restart_file << "#rate "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg << " " << infer_tree->omega_avg << endl;
	      restart_file <<"#mode " << " " <<mode<< " "<< explore_mode <<" "<< to <<endl;
	      restart_file.close();
	    }

	  ll=LL_mpi(infer_tree,step.second,mode);
	  save=step.second;
	  trialstep=ML_step(infer_tree,step.second,explore_mode,"order",greedy);
	  trialll=LL_mpi(infer_tree,trialstep.second,mode);
	  if (trialll>ll )
	    {
	      step.first=trialstep.first;
	      step.second=trialstep.second;
	    }
	  else
	    {
	      step.first=-1;
	      step.second=save;
	    }

      if (rank==server)
	{
	  ofstream cout_stream(cout_name.c_str(),ios::app);
	  cout_stream << "#LTO "<< step.second << " " << step.first <<endl;
	  cout_stream << "#LTO "<< infer_tree->tmp_ll << " "<< step.first <<endl;
	  cout_stream.close();
	}

      to=print_time_order(step.second);	  
      if (end.count(to)==0)
	end[to]=1;
      else
	break;
	}
    }
  
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream << "#O_LTO "<< step.second << " " << step.first <<endl;
      cout_stream << "#O_LTO "<< infer_tree->tmp_ll << " "<< step.first  <<endl;
      cout_stream.close();
    }

  ll=LL_mpi(infer_tree,Sstring,mode);  
  save=Sstring;
  trialstep=last_step(infer_tree,Sstring,explore_mode,greedy);
  trialll=LL_mpi(infer_tree,trialstep.second,mode);
  if (trialll>ll )
    {
      step.first=trialstep.first;
      step.second=trialstep.second;
    }
  else
    {
      step.first=-1;
      step.second=save;
    }
  end.clear();


  while (step.first>0 || c==1)
    {
      count++;
      //step=ML_step(infer_tree,step.second,mode,"root",greedy);
      to=print_time_order(step.second);
      if (rank==server)
	{
	  ofstream restart_file(restart_file_name.c_str());    
	  restart_file << "#Ml "<< step.second;
	  restart_file << "#rate "<< ll  << " "  << infer_tree->delta_avg << " " << infer_tree->tau_avg << " " << infer_tree->lambda_avg<< " " << infer_tree->omega_avg << endl;
	  restart_file <<"#mode " << " " <<mode<< " "<< explore_mode <<" "<< to <<endl;
	  restart_file.close();
	}
      ll=LL_mpi(infer_tree,step.second,mode);
      save=step.second;
      trialstep=last_step(infer_tree,step.second,explore_mode,greedy);
      trialll=LL_mpi(infer_tree,trialstep.second, mode);
      if (trialll>ll )
	{
	  step.first=trialstep.first;
	  step.second=trialstep.second;
	}
      else
	{
	  step.first=-1;
	  step.second=save;
	}      

      if (rank==server)
	{
	  ofstream cout_stream(cout_name.c_str(),ios::app);
	  cout_stream << "#FTO "<< step.second << " " << step.first <<endl;
	  cout_stream << "#FTO "<< infer_tree->tmp_ll << " "<< step.first   <<endl;
	  cout_stream.close();
	}
      
      
      to=print_time_order(step.second);
      if (end.count(to)==0)
	end[to]=1;
      else
	break;
    }

  to=print_time_order(step.second);
  if (rank==server)
    {
      string tmp=restart_file_name+"."+mode+"_"+explore_mode;

      ofstream restart_file(tmp.c_str());    
      restart_file << "#O "<< step.second;
      restart_file << "#rate "<< ll  << " "  <<  infer_tree->delta_avg << " " << infer_tree->tau_avg<< " " << infer_tree->lambda_avg << " " << infer_tree->omega_avg << endl;      
      restart_file <<"#mode " << " " <<mode<< " "<< explore_mode <<" "<< to <<endl;
      //restart_file << infer_tree->rate_trees["D"] << endl << infer_tree->rate_trees["tau"] << endl <<infer_tree->rate_trees["lambda"] << endl << infer_tree->rate_trees["omega"] << endl;
      restart_file.close();
    }
  
  if (rank==server)
    {
      ofstream cout_stream(cout_name.c_str(),ios::app);
      cout_stream << "#O_TO "<< step.second << " " << step.first <<endl;
      cout_stream << "#O_TO "<< infer_tree->tmp_ll << " "<< step.first <<endl;
      cout_stream.close();
    }

  return step.second;
}



int main (int argc, char *  argv[]) 
{

  int rank, size;
  int server = 0;  
  //Using BOOST :
  mpi::environment env(argc, argv);
  mpi::communicator world;
  rank = world.rank();
  size = world.size();
  if (size==1) 
    {
      if (rank==server) cout <<"\n\n\n\t\tError: this program can only run if 2 or more processes are used."<<endl;
      if (rank==server) cout <<"\t\tUse 'mpirun -np k mpi_test ...', where k>=2"<<endl;
      exit(-1);
    }	
  if (argc<4)
    {
      if (rank==server) cout<< "mpirun -np k mpi_test forest S outgroup \n"; 
      return 1;
    } 
  //RowMatrix<double> mym(3,3);
  //mym(0,0)=1;
  //mym(1,1)=1;
  //mym(2,2)=1;
  //EigenValue<double> ev(mym);
  
  int outgroup=0;//atoi(argv[3]);
  pair<Species_tree *,string> tree_pair;
  string outstream_file_name=argv[3];      
  outstream_file_name+="cout";
  
  tree_pair=init_LL_mpi(argv[2],argv[1],outgroup,world);
  //tree_pair=sim_init_LL_mpi(argv[2],argv[1],1,world);
  if (rank==server &&0)
    {
      ofstream cout_stream(outstream_file_name.c_str(),ios::app);  
      cout << "#trees simulated" <<endl;     
      cout_stream << "#trees simulated" <<endl;
      cout_stream.close();     
    }
 
  Species_tree * infer_tree=tree_pair.first;
  infer_tree->MPI=true;
  string Sstring=tree_pair.second;
  string S0string=Sstring;

  infer_tree->outstream_file_name=outstream_file_name;
  string restart_file_name=argv[3];
  restart_file_name+="restart";

  boost::timer * t;

  if (rank==server)
    {
      ofstream cout_stream(outstream_file_name.c_str(),ios::app);  
      cout_stream << "#trees scattered" <<endl;
      cout << "#trees scattered" <<endl;
      cout_stream.close();     
      t=new boost::timer();
    }



  scalar_type ll=LL_mpi(infer_tree,Sstring,"Hest");
  if (rank==server)
    {
      ofstream cout_stream(outstream_file_name.c_str(),ios::app);  
      cout_stream << ll << " Hest " << t->elapsed()<<endl;
      cout << ll << " Hest " << t->elapsed()<<endl;
      
      cout_stream.close();
      infer_tree->print_rates();
    }

  ll=LL_mpi(infer_tree,Sstring,"IHest");
  if (rank==server)
    {
      ofstream cout_stream(outstream_file_name.c_str(),ios::app);  
      cout_stream << ll << " IHest " << t->elapsed()<<endl;
      cout << ll << " Hest " << t->elapsed()<<endl;
      
      cout_stream.close();
      infer_tree->print_rates();
    }

  infer_tree->REF_Sstring=S0string;
  // -- ref 1. --
  Sstring=infer_tree->panj_string;    
  Sstring=ML_SPR(outstream_file_name,infer_tree,Sstring,world,restart_file_name,"clock_Test",0,true,S0string,true,false,false);
  Sstring=ML_SPR(outstream_file_name,infer_tree,Sstring,world,restart_file_name,"clock_Hest",0,true,S0string,false,true);
  Sstring=ML_NNI(outstream_file_name,infer_tree,Sstring,world,"Hest","noestimate",0,true,S0string,false);
  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"Hest","tmp_restart","noestimate",0,true);
  // -- ref 1. --
  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"Hest","tmp_restart","Hest",0,true);
  Sstring=ML_NNI(outstream_file_name,infer_tree,Sstring,world,"Hest","Hest",0,true,S0string,false);
  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"Hest","tmp_restart","Hest",0,true);
  Sstring=ML_NNI(outstream_file_name,infer_tree,Sstring,world,"IHest","IHest",0,true,S0string,false);
  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"IHest","tmp_restart","IHest",0,true);

  return 1;
  Sstring=DTL_clock(Sstring,infer_tree,"Hest");


  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"Hest","tmp_restart","noestimate",0,true);
  Sstring=ML_NNI(outstream_file_name,infer_tree,Sstring,world,"Hest","noestimate",0,true,S0string,false);

  Sstring=ML_run(outstream_file_name, infer_tree,Sstring,world,"Hest","tmp_restart","Hest",0,true);

  return 1;
 
}
