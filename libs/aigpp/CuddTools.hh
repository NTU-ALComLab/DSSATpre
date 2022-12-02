#ifndef AIGPP_CUDDTOOLS_HH
#define AIGPP_CUDDTOOLS_HH

#include <iostream>

struct DdManager;

void printBDDStatistics(DdManager* bddman, std::ostream& os = std::cout);

#endif /* AIGPP_CUDDTOOLS_HH */
