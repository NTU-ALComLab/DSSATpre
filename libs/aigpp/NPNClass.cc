/**************************************************************
 *
 *       AIGPP // NPNClass.cc
 *
 *       Copyright (C) 2008 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include <iostream>

#include "NPNClass.hh"

std::vector<aigpp::NPNClassAndTransformations>
aigpp::computeNPNClasses()
{
    NPNClassAndTransformations invalidNPNClassAndTransformations;
    invalidNPNClassAndTransformations.npnClass = InvalidNPNClass;

    std::vector<NPNClassAndTransformations> f2n(1 << 16, invalidNPNClassAndTransformations);

    unsigned int nextNPNClass = 0;

    for (FourInputFunction f = 0; f != (1 << 16); ++f) {
        /* NPN class already computed */
        if (f2n[f].npnClass != InvalidNPNClass) {
            continue;
        }

        /* apply all possible transformations to f */
        for (unsigned int negations = 0; negations != (1 << 5); ++negations) {
            bool n0   = (negations & (1 << 0)) != 0;
            bool n1   = (negations & (1 << 1)) != 0;
            bool n2   = (negations & (1 << 2)) != 0;
            bool n3   = (negations & (1 << 3)) != 0;
            bool nout = (negations & (1 << 4)) != 0;

            for (int p0 = 0; p0 != 4; ++p0) {
                for (int p1 = 0; p1 != 4; ++p1) {
                    if (p1 == p0) {
                        continue;
                    }
                    for (int p2 = 0; p2 != 4; ++p2) {
                        if (p2 == p0 || p2 == p1) {
                            continue;
                        }
                        for (int p3 = 0; p3 != 4; ++p3) {
                            if (p3 == p0 || p3 == p1 || p3 == p2) {
                                continue;
                            }

                            NPNTransformation t1 = initTransformation(nout, n0, n1, n2, n3, p0, p1, p2, p3);
                            FourInputFunction ft = applyTransformation(f, t1);

                            NPNClassAndTransformations& ntt = f2n[ft];

                            if (ntt.npnClass == nextNPNClass) {
                                continue;
                            }
                            assert(ntt.npnClass == InvalidNPNClass);

                            NPNTransformation t2 = inverseTranformation(t1);

                            assert(nextNPNClass <= 222);
                            ntt.npnClass = static_cast<unsigned short int>(nextNPNClass);
                            ntt.t1       = t2;
                            ntt.t2       = t1;
                        }
                    }
                }
            }
        }

        ++nextNPNClass;
    }

    /* QUIET version */
#if 0
    std::cout << "found " << nextNPNClass << " NPN classes" << std::endl;
#endif

    assert(nextNPNClass == NumberOfNPNClasses);

    return f2n;
}
