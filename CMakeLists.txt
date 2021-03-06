CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(phyldog CXX)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF(NOT NO_VIRTUAL_COV)

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF(NO_VIRTUAL_COV)




IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Distribution)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -w -DNDEBUG")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -w -DNDEBUG")



#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  IF(${BUILD_STATIC})
    SET(${lib_name}_NAMES lib${lib_name}.a)
  ELSE(${BUILD_STATIC})
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  ENDIF(${BUILD_STATIC})
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} ${lib_type} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   library: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)


#Find the libraries. The order is very important for static linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(BPP_LIBRARIES bpp-core Bpp/Clonable.h)

IF(NOT PLL_LIBNAME)
  SET(PLL_LIBNAME "pll-sse3")
ENDIF(NOT PLL_LIBNAME)


IMPROVED_FIND_LIBRARY(PLL_LIBRARIES ${PLL_LIBNAME} pll/pll.h)

# #Boost libraries:
# SET(Boost_USE_STATIC_LIBS ON)
# find_package(Boost COMPONENTS mpi serialization REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})
# set(Boost_LIBRARIES_STATIC ${Boost_LIBRARIES})
# set(Boost_SERIALIZATION_LIBRARY_STATIC ${Boost_SERIALIZATION_LIBRARY})
# set(Boost_MPI_LIBRARY_STATIC ${Boost_MPI_LIBRARY})
# UNSET(Boost_LIBRARIES)
# UNSET(Boost_SERIALIZATION_LIBRARY_STATIC)
# UNSET(Boost_MPI_LIBRARY_STATIC)
# 


SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})

find_package(Boost COMPONENTS mpi serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH}) 
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Source Directory
add_subdirectory(src ${PHYLDOG_BINARY_DIR})

# Tests
add_subdirectory(tests)


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
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_BUILD_TOOL} package_source)
