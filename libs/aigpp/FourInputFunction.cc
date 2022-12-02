/**************************************************************
 *
 *       AIGPP // FourInputFunction.cc
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

#include "FourInputFunction.hh"

void
aigpp::printFunction(const aigpp::FourInputFunction& f, std::ostream& os)
{
    for (unsigned int i = 0; i != 16; ++i) {
        os << ((f >> i) & 1);
    }
}

void
aigpp::printTransformation(const aigpp::NPNTransformation& t, std::ostream& os)
{
    os << "nout=" << isOutputNegated(t) << " ninput=" << isInputNegated(t, 0) << ":" << isInputNegated(t, 1) << ":"
       << isInputNegated(t, 2) << ":" << isInputNegated(t, 3) << " pinput=" << getInputPermutation(t, 0) << ":"
       << getInputPermutation(t, 1) << ":" << getInputPermutation(t, 2) << ":" << getInputPermutation(t, 3);
}

aigpp::NPNTransformation
aigpp::initTransformation(bool nout, bool n0, bool n1, bool n2, bool n3, unsigned int p0, unsigned int p1,
                          unsigned int p2, unsigned int p3)
{
    NPNTransformation t = 0;
    /*
    t |= ( nout ? 1 : 0 );
    t |= ( n0 ? ( 1 << 0 ) : 0 ) << 1;
    t |= ( n1 ? ( 1 << 1 ) : 0 ) << 1;
    t |= ( n2 ? ( 1 << 2 ) : 0 ) << 1;
    t |= ( n3 ? ( 1 << 3 ) : 0 ) << 1;
    t |= ( p0 << 0 ) << 5;
    t |= ( p1 << 2 ) << 5;
    t |= ( p2 << 4 ) << 5;
    t |= ( p3 << 6 ) << 5;
    */

    t |= (nout ? 1 : 0);
    t |= (n0 ? (1 << 1) : 0);
    t |= (n1 ? (1 << 2) : 0);
    t |= (n2 ? (1 << 3) : 0);
    t |= (n3 ? (1 << 4) : 0);
    t |= (p0 << 5);
    t |= (p1 << 7);
    t |= (p2 << 9);
    t |= (p3 << 11);

    return t;
}

aigpp::NPNTransformation
aigpp::inverseTranformation(const aigpp::NPNTransformation& t)
{
    unsigned int p0 = getInputPermutation(t, 0);
    unsigned int p1 = getInputPermutation(t, 1);
    unsigned int p2 = getInputPermutation(t, 2);
    unsigned int p3 = getInputPermutation(t, 3);

    unsigned int ip[4] = {5, 5, 5, 5};
    bool         in[4] = {false, false, false, false};

    ip[p0] = 0;
    ip[p1] = 1;
    ip[p2] = 2;
    ip[p3] = 3;

    in[p0] = isInputNegated(t, 0);
    in[p1] = isInputNegated(t, 1);
    in[p2] = isInputNegated(t, 2);
    in[p3] = isInputNegated(t, 3);

    return initTransformation(isOutputNegated(t), in[0], in[1], in[2], in[3], ip[0], ip[1], ip[2], ip[3]);
}

aigpp::NPNTransformation
aigpp::concatTranformations(const aigpp::NPNTransformation& t1, const aigpp::NPNTransformation& t2)
{
    int pt1[4] = {5, 5, 5, 5};
    for (int i = 0; i != 4; ++i) {
        pt1[i] = getInputPermutation(t1, i);
    }

    return initTransformation(
        isOutputNegated(t1) != isOutputNegated(t2), (isInputNegated(t1, 0) != isInputNegated(t2, pt1[0])),
        (isInputNegated(t1, 1) != isInputNegated(t2, pt1[1])), (isInputNegated(t1, 2) != isInputNegated(t2, pt1[2])),
        (isInputNegated(t1, 3) != isInputNegated(t2, pt1[3])), getInputPermutation(t2, pt1[0]),
        getInputPermutation(t2, pt1[1]), getInputPermutation(t2, pt1[2]), getInputPermutation(t2, pt1[3]));
}

aigpp::FourInputFunction
aigpp::applyTransformation(const aigpp::FourInputFunction& f, const aigpp::NPNTransformation& t)
{
    FourInputFunction ft = f;
    if (isOutputNegated(t)) {
        ft ^= 0xFFFF;
    }

    if (isInputNegated(t, 0)) {
        ft = ((ft & 0xAAAA) >> 1) | ((ft & 0x5555) << 1);
    }

    if (isInputNegated(t, 1)) {
        ft = ((ft & 0xCCCC) >> 2) | ((ft & 0x3333) << 2);
    }

    if (isInputNegated(t, 2)) {
        ft = ((ft & 0xF0F0) >> 4) | ((ft & 0x0F0F) << 4);
    }

    if (isInputNegated(t, 3)) {
        ft = ((ft & 0xFF00) >> 8) | ((ft & 0x00FF) << 8);
    }

    FourInputFunction ftp = 0;

    unsigned int p0 = getInputPermutation(t, 0);
    unsigned int p1 = getInputPermutation(t, 1);
    unsigned int p2 = getInputPermutation(t, 2);
    unsigned int p3 = getInputPermutation(t, 3);

    for (int i3 = 0; i3 != 2; ++i3) {
        for (int i2 = 0; i2 != 2; ++i2) {
            for (int i1 = 0; i1 != 2; ++i1) {
                for (int i0 = 0; i0 != 2; ++i0) {
                    int oldindex = (i0 << 0) | (i1 << 1) | (i2 << 2) | (i3 << 3);
                    int newindex = (i0 << p0) | (i1 << p1) | (i2 << p2) | (i3 << p3);

                    ftp |= ((ft >> oldindex) & 1) << newindex;
                }
            }
        }
    }

    return ftp & 0xFFFF;
}
