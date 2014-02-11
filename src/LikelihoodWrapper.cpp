


extern "C" {
#include <pll/pll.h>
}

#include "LikelihoodWrapper.h"


using namespace std;
using namespace bpp;

loadPLLtree(){
  alignmentData_PLL = pllParseAlignmentFile(PLL_FORMAT_FASTA, fastaForPLL);
}