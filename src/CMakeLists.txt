SET(PHYLDOG_SRCS
  SpeciesTreeLikelihood.cpp
  SpeciesTreeLikelihood.h
  SpeciesTreeExploration.h
  SpeciesTreeExploration.cpp
  ReconciliationTools.h
  ReconciliationTools.cpp
  DLGeneTreeLikelihood.cpp
  DLGeneTreeLikelihood.h
  COALGeneTreeLikelihood.cpp
  COALGeneTreeLikelihood.h
  GeneTreeLikelihood.cpp
  GeneTreeLikelihood.h
  COALTools.cpp COALTools.h
  GenericTreeExplorationAlgorithms.h
  GenericTreeExplorationAlgorithms.cpp
  GeneTreeAlgorithms.h
  GeneTreeAlgorithms.cpp
  FastRHomogeneousTreeLikelihood.h
  FastRHomogeneousTreeLikelihood.cpp
  ClientComputingGeneLikelihoods.h
  ClientComputingGeneLikelihoods.cpp
  )
ADD_EXECUTABLE(phyldog ReconcileDuplications.cpp ${PHYLDOG_SRCS})
target_link_libraries(phyldog 
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${LIBS}
)

install(TARGETS phyldog DESTINATION bin)
