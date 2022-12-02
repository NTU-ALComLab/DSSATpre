#ifndef AIGPP_EQUIVALENCECHECKER_HH
#define AIGPP_EQUIVALENCECHECKER_HH

#include <vector>

#include "EdgeRef.hh"
#include "Manager.hh"

namespace aigpp {
/*
 *
 * Class for checking the equivalence of two (or more) vectors of EdgeRefs.
 *
 * Usage:
 *
 * ...
 * EquivalenceChecker eq;
 * eq.setOriginalRoots( someVectorOfEdgeRefs );
 * ... do something, e.g. rewriting, etc.
 * bool equivalent = eq.checkNewRoots( someOtherVectorOfEdgeRefs );
 * if( !equivalent ) ...
 */
class EquivalenceChecker
{
   public:
    /* initialize the equivalence checker by some vector of EdgeRefs */
    void setOriginalRoots(const std::vector<EdgeRef>& roots);
    void setOriginalRoots(const std::vector<InternalEdgeRef>& roots);

    /* perform pairwise equivalence checks of the elements of "roots" with the
     * elements of the vector set by "setOriginalRoots" return "true" if all pairs
     * are equivalent
     */
    bool checkNewRoots(const std::vector<EdgeRef>& roots);
    bool checkNewRoots(const std::vector<InternalEdgeRef>& roots);

   private:
    InternalEdgeRef copyCone(const InternalEdgeRef& root, std::map<Node*, InternalEdgeRef>& mapping);
    void            copyVariables(Manager* manager);
    Manager*        getManager(const std::vector<EdgeRef>& roots) const;
    Manager*        getManager(const std::vector<InternalEdgeRef>& roots) const;

   private:
    Manager              man_;
    std::vector<EdgeRef> roots_;
};
}  // namespace aigpp

#endif /* AIGPP_EQUIVALENCECHECKER_HH */
