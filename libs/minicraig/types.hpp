/**********************************************************************
 * types.hpp
 *
 * Stefan Kupferschmid <skupfers@informatik.uni-freiburg.de>
 *
 **********************************************************************/
#ifndef MINICRAIG_TYPES_HPP
#define MINICRAIG_TYPES_HPP

namespace minicraig {

    enum ClauseType{A_CLAUSE = 0,
                    B_CLAUSE = 1,
                    L_CLAUSE = 2
    };

    enum VarType{UNDEF = 0,
                 A_LOCAL = 1,
                 B_LOCAL = 2,
                 GLOBAL = 3
    };

    enum ClauseInfo{INITIAL_CLAUSE = 0,
                    TARGET_CLAUSE = 1,
                    TRANS_CLAUSE = 2,
                    CRAIG_CLAUSE = 3,
                    LEARNT_CLAUSE = 4
    };

    enum BMCResult{UNSAT = 0,
                   SAT = 1,
                   OVER_UNSAT = 2,
                   OVER_SAT = 3,
                   UNDEF_RESULT = 4,
                   SAFE = 5
    };

}

#endif /* TYPES_HPP */
