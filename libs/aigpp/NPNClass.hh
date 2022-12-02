/**************************************************************
 *
 *       AIGPP // NPNClass.hh
 *
 *       Copyright (C) 2008 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 452 $
 *         $Author: pigorsch $
 *
 ***************************************************************/
#ifndef AIGPP_NPNCLASS_HH
#define AIGPP_NPNCLASS_HH

#include <list>
#include <vector>

#include "DynamicRewritingTable.hh"
#include "FourInputFunction.hh"

namespace aigpp {
enum
{
    NumberOfNPNClasses = 222,
    InvalidNPNClass    = NumberOfNPNClasses + 1
};

struct NPNClassAndTransformations
{
    unsigned short int npnClass;
    NPNTransformation  t1, t2;

    bool operator==(const NPNClassAndTransformations& ntt) const
    {
        return (npnClass == ntt.npnClass && t1 == ntt.t1 && t2 == ntt.t2);
    }
};

struct NPNClassImplementation
{
    NPNClassImplementation() : size(-1), usedInputsMask(0), usedInputs(0), implementations() {}

    int                                     size;
    unsigned int                            usedInputsMask;
    unsigned int                            usedInputs;
    std::list<DynamicRewritingTable::Index> implementations;
};

std::vector<NPNClassAndTransformations> computeNPNClasses();

}  // namespace aigpp

#endif /* AIGPP_NPNCLASS_HH */
