/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a bioinformatics computer program whose purpose is to
simultaneously build gene and species trees when gene families have
undergone duplications and losses. It can analyze thousands of gene
families in dozens of genomes simultaneously, and was presented in
an article in Genome Research. Trees and parameters are estimated
in the maximum likelihood framework, by maximizing theprobability
of alignments given the species tree, the gene trees and the parameters
of duplication and loss.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;


#include "SpeciesTreeLikelihood.h"
#include "ClientComputingGeneLikelihoods.h"




/*
This program takes a species tree file, sequence alignments and files describing the relations between sequences and species.
To compile it, use the Makefile !
*/

/******************************************************************************/

void help()
{
    (*ApplicationTools::message << "__________________________________________________________________________").endLine();
    (*ApplicationTools::message << "phyldog parameter1_name=parameter1_value parameter2_name=parameter2_value"   ).endLine();
    (*ApplicationTools::message << "      ... param=option_file").endLine();
    (*ApplicationTools::message << "Example of some options: ").endLine();
    (*ApplicationTools::message).endLine();
    
    /*SequenceApplicationTools::printInputAlignmentHelp();
     PhylogeneticsApplicationTools::printInputTreeHelp();
     PhylogeneticsApplicationTools::printSubstitutionModelHelp();
     PhylogeneticsApplicationTools::printRateDistributionHelp();
     PhylogeneticsApplicationTools::printCovarionModelHelp();
     PhylogeneticsApplicationTools::printOptimizationHelp(true, false);
     PhylogeneticsApplicationTools::printOutputTreeHelp();*/
    //  (*ApplicationTools::message << "output.infos                      | file where to write site infos").endLine();
    //  (*ApplicationTools::message << "output.estimates                  | file where to write estimated parameter values").endLine();
    (*ApplicationTools::message << "PATH                                 | path to the directory containing input and output files").endLine();
    (*ApplicationTools::message << "init.species.tree                    | user or random or mrp").endLine();
    (*ApplicationTools::message << "species.tree.file                    | if the former is set to \"user\", path to a species tree" ).endLine();
    (*ApplicationTools::message << "species.names.file                   | if instead it was set to \"random\", or \"mrp\" path to a list of species names ").endLine();
    (*ApplicationTools::message << "starting.tree.file                   | file where to write the initial species tree").endLine();
    (*ApplicationTools::message << "output.tree.file                     | file where to write the end species tree").endLine();
    (*ApplicationTools::message << "output.temporary.tree.file           | file where to write the temporary species tree in case the job is ended because of time constraints.").endLine();

    //  (*ApplicationTools::message << "heuristics.level                  | 0, 1, 2, or 3; 0: DR exact algorithm (default, best); 1 fast heuristics; 2: exhaustive and fast; 3: exhaustive and slow").endLine();
    //(*ApplicationTools::message << "species.id.limit.for.root.position| Threshold for trying root positions").endLine();
    (*ApplicationTools::message << "genelist.file                        | file containing a list of gene option files to analyse").endLine();
    (*ApplicationTools::message << "branch.expected.numbers.optimization | average, branchwise, average_then_branchwise or no: how we optimize duplication and loss probabilities").endLine();
    (*ApplicationTools::message << "genome.coverage.file                 | file giving the percent coverage of the genomes used").endLine();
    (*ApplicationTools::message << "spr.limit                            | integer giving the breadth of SPR movements, in number of nodes. 0.1* number of nodes in the species tree might be OK.").endLine();
    (*ApplicationTools::message << "reconciliation.model                 | 'DL' or 'COAL' giving the type of model to reconcile gene trees against the species tree.").endLine();

    (*ApplicationTools::message << "  Refer to the README file or the Bio++ Program Suite Manual for a list of supplementary options.").endLine();
    (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}


/******************************************************************************/
/***************************Adding two std::vectors*********************************/


std::vector <int> operator + (std::vector <int> x, std::vector <int> y) {
  unsigned int temp=x.size();
  if (temp!=y.size()) {
    std::cout <<"problem : adding two std::vectors of unequal sizes."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
    exit (-1);
  }
 std::vector<int> result;
  for(unsigned int i = 0 ; i<temp ; i++){
    result.push_back(x[i]+y[i]);
  }
  return result;
}

std::vector <double> operator + (std::vector <double> x, std::vector <double> y) {
  unsigned int temp=x.size();
  if (temp!=y.size()) {
    std::cout <<"problem : adding two std::vectors of unequal sizes."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
    exit (-1);
  }
 std::vector<double> result;
  for(unsigned int i = 0 ; i<temp ; i++){
    result.push_back(x[i]+y[i]);
  }
  return result;
}





/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
//////////////////////////////////////////////////////MAIN/////////////////////////////////////////////////
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


int main(int args, char ** argv)
{
  if(args == 1)
  {
    help();
    exit(0);
  }

  unsigned int rank, size;
  unsigned int server = 0;

  //Using BOOST :
  mpi::environment env(args, argv);
  mpi::communicator world;
  rank = world.rank();
  size = world.size();
  std::string line;

  if (size==1) {
    std::cout <<"\n\n\n\t\tError: this program can only run if 2 or more processes are used."<<std::endl;
    std::cout <<"\t\tUse 'mpirun -np k phyldog ...', where k>=2"<<std::endl;
    MPI::COMM_WORLD.Abort(1);
    exit(-1);
  }
  
  
  try {
    ApplicationTools::startTimer();
        //All processors parse the main options
    std::map<std::string, std::string> params = AttributesTools::parseOptions(args, argv);

    
        //##################################################################################################################
        //##################################################################################################################
        //############################################# IF AT THE SERVER NODE ##############################################
        //##################################################################################################################
        //##################################################################################################################
        
    if (rank == server) { 
            
            /*
            int z = 0;
            //   char hostname[256];
            //gethostname(hostname, sizeof(hostname));
            std::cout <<"PID: "<<getpid()<<std::endl;
            std::cout <<"z: "<<z<<std::endl;
            //printf("PID %d on %s ready for attach\n", getpid(), hostname);
            // fflush(stdout);
            while (0 == z){
                std::cout << z <<std::endl;
                sleep(5);
            }*/
            
            
            std::cout << "******************************************************************" << std::endl;
            std::cout << "*                   PHYLDOG, version 1.1.0                       *" << std::endl;
            std::cout << "* Author: B. Boussau                            Created 16/07/07 *" << std::endl;
            std::cout << "******************************************************************" << std::endl;
            std::cout << std::endl;
            
            std::cout <<"Server of rank "<<rank <<" with PID "<< TextTools::toString((int)getpid())<<std::endl;

            SpeciesTreeLikelihood spTL = SpeciesTreeLikelihood(world, server, size, params);
            //initialization, first communications, first likelihood computation

            spTL.initialize();

            spTL.MLSearch();
                        
      std::cout << "PHYLDOG's done. Bye." << std::endl;
      ApplicationTools::displayTime("Total execution time:");
            MPI_Barrier(world);  
            MPI::Finalize( );
    }//End if at the server node
      
        
        //##################################################################################################################
        //##################################################################################################################
        //############################################# IF AT A CLIENT NODE ################################################
        //##################################################################################################################
        //##################################################################################################################
        if (rank >server)
	  {
               //  Code useful to use GDB on clients.
 	/*	int z = 0;
 	      //   char hostname[256];
 	      //gethostname(hostname, sizeof(hostname));
 	      std::cout <<"- PID: "<<getpid()<<std::endl;
 	      std::cout <<"-   z: "<<z<<std::endl;
 	      //printf("PID %d on %s ready for attach\n", getpid(), hostname);
 	      // fflush(stdout);
 	      while (0 == z){
                 std::cout << " z=" << z <<std::endl;
                 sleep(5);
 	      }*/
      
	      ApplicationTools::startTimer();
	      bool debug = ApplicationTools::getBooleanParameter("debug",params,false);
	      string path = ApplicationTools::getStringParameter("PATH", params, "", "", true, false);
	      string outputFile = path + "Client_"+TextTools::toString(rank)+".out";
	      streambuf *psbuf, *backup, *backupcerr;
	      ofstream filestr;
	      if (!debug) {
		  //Redirecting stdout to a specific file for this client
		  filestr.open (outputFile.c_str());
		  backup = cout.rdbuf();     // back up cout's streambuf 
		  backupcerr = cerr.rdbuf(); // back up cerr's streambuf
		  psbuf = filestr.rdbuf();   // get file's streambuf                                                                                                                                                               
		  cout.rdbuf(psbuf);         // assign streambuf to cout     
		  cerr.rdbuf(psbuf);
	      }
	      
	      //initialization, first communications, first likelihood computation
	      ClientComputingGeneLikelihoods client = ClientComputingGeneLikelihoods(world, server, rank, params);
	      
	      //Main loop, computation of the gene tree likelihoods given species trees sent by the server.
	      client.MLSearch();

			if (!debug) {
				cerr.rdbuf(backupcerr);    // restore cerr's original streambuf
				cout.rdbuf(backup);        // restore cout's original streambuf                                                                                                                                                  
				filestr.close();
			}
			MPI_Barrier(world);  
			MPI::Finalize( );
		}//end if a client node
  }
  catch(std::exception & e)
  {
    std::cout << e.what() << std::endl;
        MPI::COMM_WORLD.Abort(1);
    exit(-1);
  }
  return (0);
}










/*
 // This bit of code is useful to use GDB on clients, when put into the client's code:
 //launch the application, which will output the client pid
 //then launch gdb, attach to the given pid ("attach pid" or "gdb phyldog pid"), 
 //use "up" to go up the stacks, and set the variable z to !=0 to get out of the loop with "set var z = 8".
 int z = 0;
 //   char hostname[256];
 //gethostname(hostname, sizeof(hostname));
 std::cout <<"PID: "<<getpid()<<std::endl;
 std::cout <<"z: "<<z<<std::endl;
 //printf("PID %d on %s ready for attach\n", getpid(), hostname);
 // fflush(stdout);
 while (0 == z){
 std::cout <<z<<std::endl;
 sleep(5);
 }
 */





/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////RUBBISH////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


/*This code permits outputting trees with numbers of duplication or loss events at branches
 //set total numbers of loss and duplications on branches
 setLossesAndDuplications(*bestTree, bestLossNumbers, bestDuplicationNumbers);
 std::string dupTree = ApplicationTools::getStringParameter("output.duplications.tree.file", params, "AllDuplications.tree", "", false, false); 
 std::ofstream out (dupTree.c_str(), std::ios::out);
 out << treeToParenthesisWithIntNodeValues (*bestTree, false, DUPLICATIONS)<<std::endl;
 out.close();
 std::string lossTree = ApplicationTools::getStringParameter("output.losses.tree.file", params, "AllLosses.tree", "", false, false);
 out.open (lossTree.c_str(), std::ios::out);
 out << treeToParenthesisWithIntNodeValues (*bestTree, false, LOSSES)<<std::endl;
 out.close();
 */








