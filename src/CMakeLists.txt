# CMake script for Phyldog
# Author: Bastien Boussau
# Created: 12/05/2010

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(phyldog CXX)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF(NOT NO_VIRTUAL_COV)

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF(NO_VIRTUAL_COV)

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()



IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE Distribution)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -w")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -w")



#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})

  IF(BUILD_STATIC)
    FIND_PATH(${lib_name}_LIBRARY_DIR lib${lib_name}.a HINTS ${CMAKE_PREFIX_PATH} ${CMAKE_PREFIX_PATH}/lib)
    IF(${lib_name}_LIBRARY_DIR)
      SET(${lib_name}_LIBRARY ${${lib_name}_LIBRARY_DIR}/lib${lib_name}.a)
    ELSE(${lib_name}_LIBRARY_DIR)
      MESSAGE(FATAL_ERROR "${lib_name} static library could not be found.")
    ENDIF(${lib_name}_LIBRARY_DIR)
  ELSE(BUILD_STATIC)
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  ENDIF(BUILD_STATIC)

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   libraries: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ELSE(
  SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" ))
ENDIF(CMAKE_INSTALL_PREFIX)

#Finidng the appropriate openmp flags
# ENABLE_LANGUAGE(C)
# FIND_PACKAGE( OpenMP)
# if(OPENMP_FOUND)
# message("OPENMP FOUND")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()


#Find the libraries. The order is very important for static linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

#Boost libraries:
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
include(FindBoostAlt)
find_package(BoostAlt COMPONENTS serialization mpi REQUIRED)

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH}) 
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})


#Compiling executables
SET(PHYLDOG_SRCS SpeciesTreeLikelihood.cpp SpeciesTreeLikelihood.h SpeciesTreeExploration.h SpeciesTreeExploration.cpp ReconciliationTools.h ReconciliationTools.cpp DLGeneTreeLikelihood.cpp DLGeneTreeLikelihood.h COALGeneTreeLikelihood.cpp COALGeneTreeLikelihood.h GeneTreeLikelihood.cpp GeneTreeLikelihood.h COALTools.cpp COALTools.h GenericTreeExplorationAlgorithms.h GenericTreeExplorationAlgorithms.cpp GeneTreeAlgorithms.h GeneTreeAlgorithms.cpp FastRHomogeneousTreeLikelihood.h FastRHomogeneousTreeLikelihood.cpp)
ADD_EXECUTABLE(phyldog ReconcileDuplications.cpp ${PHYLDOG_SRCS})


target_link_libraries(phyldog 
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${LIBS}
)
# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(RECDUP_VERSION_CURRENT "1")
SET(RECDUP_VERSION_REVISION "0")
SET(RECDUP_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR RECDUP_VERSION_MAJOR "${RECDUP_VERSION_CURRENT} - ${RECDUP_VERSION_AGE}")
SET(RECDUP_VERSION_MINOR ${RECDUP_VERSION_AGE})
SET(RECDUP_VERSION_PATCH ${RECDUP_VERSION_REVISION})
SET(RECDUP_VERSION "${RECDUP_VERSION_MAJOR}.${RECDUP_VERSION_MINOR}.${RECDUP_VERSION_PATCH}")

# Subdirectories
#ADD_SUBDIRECTORY(Utils) No subdirectory, the code is at the root.

# Doxygen
#FIND_PACKAGE(Doxygen)
#IF (DOXYGEN_FOUND)
#  ADD_CUSTOM_TARGET (apidoc ${DOXYGEN_EXECUTABLE} Doxyfile)
#ENDIF (DOXYGEN_FOUND)

# Packager
SET(CPACK_PACKAGE_NAME "phyldog")
SET(CPACK_PACKAGE_VENDOR "Bastien Boussau")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Joint species and gene trees inference.")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 "/CVS/"
 ".*\\\\.gz"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 "\\\\..*\\\\.swp"
 "phyldog"
 ${CPACK_SOURCE_IGNORE_FILES}
)
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(deb fakeroot debian/rules binary)
ADD_DEPENDENCIES(deb dist)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist)

#Compiling executables
#ADD_EXECUTABLE(phyldog phyldog)
#TARGET_LINK_LIBRARIES(phyldog ${LIBS})


ENDIF()

