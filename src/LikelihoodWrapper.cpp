


extern "C" {
#include <pll/pll.h>
}

#include "LikelihoodWrapper.h"


using namespace std;

LikelihoodWrapper::loadPLLtree(){
  // using pllNewickParseString (""); 
  
  alignmentData_PLL() = pllParseAlignmentFile(PLL_FORMAT_FASTA, fastaForPLL);
}