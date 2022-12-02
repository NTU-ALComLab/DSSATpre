#include "CuddTools.hh"

#include <cudd.h>

void
printBDDStatistics(DdManager* bddman, std::ostream& os)
{
    if (bddman == nullptr) {
        os << "no bdd manager present\n";
        return;
    }

    os << "BDD: current nodes = " << Cudd_ReadNodeCount(bddman) << '\n';
    os << "BDD: current dead nodes = " << Cudd_ReadDead(bddman) << '\n';
    os << "BDD: peak live nodes = " << Cudd_ReadPeakLiveNodeCount(bddman) << '\n';
    os << "BDD: live node limit = " << Cudd_ReadMaxLive(bddman) << '\n';

    Cudd_ReorderingType rmethod;
    if (Cudd_ReorderingStatus(bddman, &rmethod) != 0) {
        os << "BDD: reordering enabled with method ";

        switch (rmethod) {
            case CUDD_REORDER_SIFT:
                os << "SIFT";
                break;
            case CUDD_REORDER_LAZY_SIFT:
                os << "LAZY_SIFT";
                break;
            default:
                os << "UNKNOWN(" << rmethod << ")";
        }

        os << '\n';

        os << "BDD: next reordering = " << Cudd_ReadNextReordering(bddman) << '\n';
    } else {
        os << "BDD: reordering disabled\n";
    }

    os << "BDD: garbage collection " << (Cudd_GarbageCollectionEnabled(bddman) != 0 ? "enabled" : "disabled") << '\n';

    os << std::flush;
}
