/**************************************************************
 *
 *       AIGPP // FourInputFunction.hh
 *
 *       Copyright (C) 2008 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 535 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#ifndef AIGPP_FOURINPUTFUNCTION_HH
#define AIGPP_FOURINPUTFUNCTION_HH

#include <iostream>

namespace aigpp {
using FourInputFunction = unsigned int;
void printFunction(const FourInputFunction& f, std::ostream& os = std::cout);

inline bool
dep0(const FourInputFunction& f)
{
    return ((f & 0x5555) != ((f >> 1) & 0x5555));
}
inline bool
dep1(const FourInputFunction& f)
{
    return ((f & 0x3333) != ((f >> 2) & 0x3333));
}
inline bool
dep2(const FourInputFunction& f)
{
    return ((f & 0x0f0f) != ((f >> 4) & 0x0f0f));
}
inline bool
dep3(const FourInputFunction& f)
{
    return ((f & 0x00ff) != ((f >> 8) & 0x00ff));
}

/*
 * lower 13 bits:
 * bit 0: output negation
 * bit 1-4: input negations (bit 1: first input, bit 2: second input, ...)
 * bit 5-12: input permutation indices (bits 5+6: first input, bits 7+8: second
 * input)
 */
typedef unsigned int NPNTransformation;
void                 printTransformation(const NPNTransformation& t, std::ostream& os = std::cout);

NPNTransformation initTransformation(bool nout, bool n0, bool n1, bool n2, bool n3, unsigned int p0, unsigned int p1,
                                     unsigned int p2, unsigned int p3);

inline bool
isOutputNegated(const NPNTransformation& t)
{
    return (t & 1) != 0;
}

inline bool
isInputNegated(const NPNTransformation& t, unsigned int i)
{
    return ((t >> 1) & (1 << i)) != 0;
}

inline unsigned int
getInputPermutation(const NPNTransformation& t, unsigned int i)
{
    return ((t >> 5) >> (2 * i)) & 3;
}

NPNTransformation inverseTranformation(const NPNTransformation& t);
NPNTransformation concatTranformations(const NPNTransformation& t1, const NPNTransformation& t2);

FourInputFunction applyTransformation(const FourInputFunction& f, const NPNTransformation& t);
}  // namespace aigpp

#endif /* AIGPP_FOURINPUTFUNCTION_HH */
