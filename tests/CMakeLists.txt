SET(PHYLDOG_SRCS
  ../src/SpeciesTreeLikelihood.cpp
  ../src/SpeciesTreeLikelihood.h
  ../src/SpeciesTreeExploration.h
  ../src/SpeciesTreeExploration.cpp
  ../src/ReconciliationTools.h
  ../src/ReconciliationTools.cpp
  ../src/DLGeneTreeLikelihood.cpp
  ../src/DLGeneTreeLikelihood.h
  ../src/COALGeneTreeLikelihood.cpp
  ../src/COALGeneTreeLikelihood.h
  ../src/GeneTreeLikelihood.cpp
  ../src/GeneTreeLikelihood.h
  ../src/COALTools.cpp
  ../src/COALTools.h
  ../src/GenericTreeExplorationAlgorithms.h
  ../src/GenericTreeExplorationAlgorithms.cpp
  ../src/GeneTreeAlgorithms.h
  ../src/GeneTreeAlgorithms.cpp
  ../src/FastRHomogeneousTreeLikelihood.h
  ../src/FastRHomogeneousTreeLikelihood.cpp
  ../src/ClientComputingGeneLikelihoods.h
  ../src/ClientComputingGeneLikelihoods.cpp
  ../src/LikelihoodEvaluator.h
  ../src/LikelihoodEvaluator.cpp
)

ADD_EXECUTABLE(test_likelihoodEvaluator test_likelihoodEvaluator.cpp ${PHYLDOG_SRCS})
target_link_libraries(test_likelihoodEvaluator 
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${PLL_LIBRARIES}
  ${BPP_LIBRARIES}
)

ADD_EXECUTABLE(test_SPRs test_SPRs.cpp ${PHYLDOG_SRCS})
target_link_libraries(test_SPRs
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${PLL_LIBRARIES}
  ${BPP_LIBRARIES}
)

install(TARGETS test_SPRs test_likelihoodEvaluator DESTINATION tests)
