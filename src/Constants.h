#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <string>

#include <Bpp/Text/TextTools.h>


//Defining a DEBUG macro to print debug messages.
#ifndef NDEBUG
    #define WHEREAMI(x,y) (  std::cout <<"DEBUG: "<< x << " : "<< TextTools::toString<int>(y) << std::endl )
    //std::cerr << "PASSED: "<< x << " : "<< TextTools::toString<int>(y)<< std::endl; ; std::cerr.flush() ; std::cout.flush(); 
  //  #define D(x) ( std::cerr << "PASSED: "<< x << std::endl;  std::cout << "PASSED: " << x << std::endl; std::cerr.flush() ; std::cout::flush(); )
#else 
    #define WHEREAMI(x,y)  
#endif



const std::string SPECIESID="SPECIESID";
const std::string EVENT="EVENT";
const std::string LOSSES="L";
const std::string DUPLICATIONS="D";
const std::string EVENTSPROBA="EVENTSPROBA";
const std::string LOWLIK="LOWLIK";
const std::string NUMGENES="NUMGENES";
const std::string NUMLINEAGES="NUMLINEAGES";

const double UNLIKELY=-100000000000000000000.0;
const double SMALLPROBA=0.0000000001;
const double BIGPROBA=0.9999999999;
const int MAXFILENAMESIZE = 500;
const int MAXSPECIESTREESIZE = 10000; //size of the species tree, in number of CHARs, as it is written in Newick format
const double DIST = 0.1;


#endif  //_CONSTANTS_H_