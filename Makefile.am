AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = ReconcileDuplications
if STATIC
ReconcileDuplications_LDFLAGS = -static
else 
ReconcileDuplications_LDFLAGS = 
endif
ReconcileDuplications_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl $(BOOST_SERIALIZATION_LIB) $(BOOST_MPI_LIB)
AM_CXXFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)
H_FILES = \
ReconciliationTools.h\
SpeciesTreeExploration.h\
ReconciliationTreeLikelihood.h

CPP_FILES = \
ReconciliationTools.cpp\
SpeciesTreeExploration.cpp\
ReconciliationTreeLikelihood.cpp\
ReconcileDuplications.cpp



ReconcileDuplications_SOURCES = $(H_FILES) $(CPP_FILES)

