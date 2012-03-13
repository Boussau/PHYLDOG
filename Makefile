#Where the boost include files and directories can be found:
BOOST_INCLUDE = /usr/local/include/boost

#Where the include files for Bio++ can be found:
BIOPP_INCLUDE = /usr/local/include

#Where the Bio++ libraries can be found:
BIOPP_LIBRARIES = /usr/local/lib

#Where the Boost libraries can be found:
BOOST_LIBRARIES = /usr/local/lib

#Where the sources are:
SRC = ./

#Debugging compilation options
CC = mpic++ -m64 -g -Wall -I $(BOOST_INCLUDE) 

#Efficient compilation options
#CC = mpic++ -m64 -O3 -I $(BOOST_INCLUDE)


#LL =  -L/usr/local/lib -lboost_serialization -lboost_mpi
LL =  -L/usr/local/lib -L$(BIOPP_LIBRARIES) -L$(BOOST_LIBRARIES) -lboost_serialization -lboost_mpi

all : RearrangeGeneTreeDTL phyldog

###DTL model

LINK =  -L. -lDTL -lblas  


INCLUDE = -I$(BOOST_INCLUDE) -I$(BIOPP_INCLUDE)

RearrangeGeneTreeDTL : $(SRC)RearrangeGeneTreeDTL.cpp DTLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o libDTL.a mpi_time_order.o
	$(CC) $(INCLUDE) -o RearrangeGeneTreeDTL $(SRC)RearrangeGeneTreeDTL.cpp DTLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o libDTL.a mpi_time_order.o -lbpp-phyl -lbpp-seq -lbpp-core $(LL) $(LINK)

DTLGeneTreeLikelihood.o : $(SRC)DTLGeneTreeLikelihood.cpp $(SRC)DTLGeneTreeLikelihood.h ReconciliationTools.o GeneTreeAlgorithms.o libDTL.a mpi_time_order.o
	$(CC) $(INCLUDE) -c $(SRC)DTLGeneTreeLikelihood.cpp $(SRC)DTLGeneTreeLikelihood.h $(SRC)GeneTreeAlgorithms.h 

libDTL.a :  simDTL.o inferDTL_improved.o inferDTL_unrooted.o  time_order.o nni.o backtrace.o simulator.o branchwise.o $(SRC)DTL.h
	ar rcs libDTL.a simDTL.o nni.o time_order.o inferDTL_unrooted.o inferDTL_improved.o backtrace.o simulator.o branchwise.o
        
simDTL.o : $(SRC)simDTL.cpp $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)simDTL.cpp  

inferDTL_improved.o : $(SRC)inferDTL_improved.cpp $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)inferDTL_improved.cpp  

inferDTL_unrooted.o : $(SRC)inferDTL_unrooted.cpp  $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)inferDTL_unrooted.cpp  

backtrace.o : $(SRC)backtrace.cpp  $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)backtrace.cpp  

time_order.o : $(SRC)time_order.cpp  $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)time_order.cpp  

nni.o : $(SRC)nni.cpp  $(SRC)Makefile $(SRC)DTL.h
	$(CC) $(INCLUDE) -c $(SRC)nni.cpp  

simulator.o : $(SRC)simulator.cpp  $(SRC)DTL.h     
	$(CC) $(INCLUDE) -c $(SRC)simulator.cpp  

branchwise.o : $(SRC)branchwise.cpp  $(SRC)DTL.h   
	$(CC) $(INCLUDE) -c $(SRC)branchwise.cpp  

mpi_time_order.o : $(SRC)mpi_time_order.cpp  $(SRC)DTL.h $(SRC)DTL_mpi.h
	$(CC) $(INCLUDE) -c $(SRC)mpi_time_order.cpp  



####DL model only

#phyldog : ReconcileDuplications.o SpeciesTreeLikelihood.o ReconciliationTreeLikelihood.o DLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o
phyldog : ReconcileDuplications.o SpeciesTreeLikelihood.o DLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o                                                                                                                                                                                                              

#	$(CC) $(INCLUDE) -o phyldog ReconcileDuplications.o SpeciesTreeLikelihood.o ReconciliationTreeLikelihood.o DLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o -lbpp-phyl -lbpp-seq -lbpp-core $(LL)

	$(CC) $(INCLUDE) -o phyldog ReconcileDuplications.o SpeciesTreeLikelihood.o DLGeneTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o GeneTreeAlgorithms.o -lbpp-phyl -lbpp-seq -lbpp-core $(LL)

ReconcileDuplications.o : $(SRC)ReconcileDuplications.cpp SpeciesTreeLikelihood.o ReconciliationTools.o SpeciesTreeExploration.o GenericTreeExplorationAlgorithms.o
	$(CC) $(INCLUDE) -c $(SRC)ReconcileDuplications.cpp $(SRC)ReconciliationTools.h 
  
GeneTreeAlgorithms.o : $(SRC)GeneTreeAlgorithms.cpp $(SRC)GeneTreeAlgorithms.h ReconciliationTools.o
	$(CC) $(INCLUDE) -c $(SRC)GeneTreeAlgorithms.cpp $(SRC)GeneTreeAlgorithms.h $(SRC)ReconciliationTools.h
  
SpeciesTreeLikelihood.o : $(SRC)SpeciesTreeLikelihood.cpp $(SRC)SpeciesTreeLikelihood.h $(SRC)SpeciesTreeExploration.h $(SRC)ReconciliationTools.h
	$(CC) $(INCLUDE) -c $(SRC)SpeciesTreeLikelihood.cpp $(SRC)SpeciesTreeLikelihood.h $(SRC)SpeciesTreeExploration.h $(SRC)ReconciliationTools.h

DLGeneTreeLikelihood.o : $(SRC)DLGeneTreeLikelihood.cpp $(SRC)DLGeneTreeLikelihood.h ReconciliationTools.o $(SRC)GenericTreeExplorationAlgorithms.h
	$(CC) $(INCLUDE) -c $(SRC)DLGeneTreeLikelihood.cpp $(SRC)DLGeneTreeLikelihood.h $(SRC)ReconciliationTools.h 


#ReconciliationTreeLikelihood.o : $(SRC)ReconciliationTreeLikelihood.cpp $(SRC)ReconciliationTreeLikelihood.h ReconciliationTools.o $(SRC)GenericTreeExplorationAlgorithms.h
#	$(CC) $(INCLUDE) -c $(SRC)ReconciliationTreeLikelihood.cpp $(SRC)ReconciliationTreeLikelihood.h $(SRC)ReconciliationTools.h 

SpeciesTreeExploration.o : $(SRC)SpeciesTreeExploration.cpp $(SRC)SpeciesTreeExploration.h $(SRC)ReconciliationTools.h GenericTreeExplorationAlgorithms.o
	$(CC) $(INCLUDE) -c $(SRC)SpeciesTreeExploration.cpp $(SRC)SpeciesTreeExploration.h $(SRC)ReconciliationTools.h

GenericTreeExplorationAlgorithms.o : $(SRC)GenericTreeExplorationAlgorithms.cpp $(SRC)GenericTreeExplorationAlgorithms.h ReconciliationTools.o
	$(CC) $(INCLUDE) -c $(SRC)GenericTreeExplorationAlgorithms.cpp $(SRC)GenericTreeExplorationAlgorithms.h $(SRC)ReconciliationTools.h

ReconciliationTools.o : $(SRC)ReconciliationTools.cpp $(SRC)ReconciliationTools.h 
	$(CC) $(INCLUDE) -c $(SRC)ReconciliationTools.cpp $(SRC)ReconciliationTools.h 


### Other useful commands
clean :
	rm *.o *.gch
  
archive : 
	tar cfz archive.tgz *.cpp *.h Makefile Makefile_MAC INSTALL README

phyldogarchive :
	mkdir PHYLDOGArchive
	cp ReconciliationTools.cpp ReconciliationTools.h GenericTreeExplorationAlgorithms.cpp GenericTreeExplorationAlgorithms.h DLGeneTreeLikelihood.cpp DLGeneTreeLikelihood.h SpeciesTreeLikelihood.cpp SpeciesTreeLikelihood.h GeneTreeAlgorithms.cpp GeneTreeAlgorithms.h ReconcileDuplications.cpp  SpeciesTreeExploration.h SpeciesTreeExploration.cpp Makefile INSTALL README PHYLDOGArchive
	tar cfz PHYLDOGArchive.tgz PHYLDOGArchive
#ReconciliationTools.cpp ReconciliationTools.h GenericTreeExplorationAlgorithms.cpp GenericTreeExplorationAlgorithms.h DLGeneTreeLikelihood.cpp DLGeneTreeLikelihood.h SpeciesTreeLikelihood.cpp SpeciesTreeLikelihood.h GeneTreeAlgorithms.cpp GeneTreeAlgorithms.h ReconcileDuplications.cpp  SpeciesTreeExploration.h SpeciesTreeExploration.cpp Makefile

examplearchive :
	#rm example/*_5Mammals
	#rm example/*_starting.tree
	#rm example/*_starting.tree_sp_names
	#rm example/Client_*
	tar cfz example.tgz example

