SET(PHYLDOG_SRCS
  Constants.h
  MPI_SpeciesTreeLikelihood.cpp
  MPI_SpeciesTreeLikelihood.h
  MPI_SpeciesTreeExploration.h
  MPI_SpeciesTreeExploration.cpp
  ReconciliationTools.h
  ReconciliationTools.cpp
  DLGeneTreeLikelihood.cpp
  DLGeneTreeLikelihood.h
  COALGeneTreeLikelihood.cpp
  COALGeneTreeLikelihood.h
  GeneTreeLikelihood.cpp
  GeneTreeLikelihood.h
  COALTools.cpp
  COALTools.h
  GenericTreeExplorationAlgorithms.h
  GenericTreeExplorationAlgorithms.cpp
  GeneTreeAlgorithms.h
  GeneTreeAlgorithms.cpp
  FastRHomogeneousTreeLikelihood.h
  FastRHomogeneousTreeLikelihood.cpp
  ClientComputingGeneLikelihoods.h
  ClientComputingGeneLikelihoods.cpp
  LikelihoodEvaluator.h
  LikelihoodEvaluator.cpp
  )
  
ADD_EXECUTABLE(phyldog Phyldog.cpp ${PHYLDOG_SRCS})
target_link_libraries(phyldog 
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${PLL_LIBRARIES}
  ${BPP_LIBRARIES}
)

SET(PHYLDOGLIGHT_SRCS
  Constants.h
  ReconciliationTools.h
  ReconciliationTools.cpp
  GenericTreeExplorationAlgorithms.h
  GenericTreeExplorationAlgorithms.cpp
  GeneTreeAlgorithms.h
  GeneTreeAlgorithms.cpp
  FastRHomogeneousTreeLikelihood.h
  FastRHomogeneousTreeLikelihood.cpp
  LikelihoodEvaluator.h
  LikelihoodEvaluator.cpp
  )


ADD_EXECUTABLE(phyldog_light Phyldog_light.cpp ${PHYLDOGLIGHT_SRCS})
target_link_libraries(phyldog_light
  ${PLL_LIBRARIES}
  ${BPP_LIBRARIES}
)



# 
# ADD_EXECUTABLE(phyldog_static ReconcileDuplications.cpp ${PHYLDOG_SRCS})
# target_link_libraries(phyldog_static 
#   ${Boost_SERIALIZATION_LIBRARY_STATIC}
#   ${Boost_MPI_LIBRARY_STATIC}
#   ${MPI_LIBRARIES}
#   ${PLL_LIBRARIES_STATIC}
#   ${BPP_LIBRARIES_STATIC}
# )

install(TARGETS phyldog DESTINATION bin)
