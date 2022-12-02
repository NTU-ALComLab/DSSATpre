/**************************************************************
 *
 *       AIGPP Package // ExtRefTable.hh
 *
 *       Copyright (C) 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: $
 *         $Author: $
 *
 ***************************************************************/

#ifndef AIGPP_EXTREFTABLE_HH
#define AIGPP_EXTREFTABLE_HH

/* stdlib + stl */
#include <map>
#include <stack>
#include <vector>

/* lrabsutil */
#include <lrabsutil/Assert.hh>

/* aigpp */
#include "InternalEdgeRef.hh"

namespace aigpp {
class ExtRefTable
{
   public:
    ExtRefTable()                   = default;
    ExtRefTable(const ExtRefTable&) = delete;
    ExtRefTable(ExtRefTable&&)      = delete;
    ExtRefTable& operator=(const ExtRefTable&) = delete;
    ExtRefTable& operator=(ExtRefTable&&) = delete;
    ~ExtRefTable();

    int                    insert(const InternalEdgeRef& target);
    const InternalEdgeRef& lookup(int index) const;

    void clear();
    void ref(int index);
    void deref(int index);

    void rebuildMapping();

   private:
    class Item
    {
       public:
        Item() : _refCount(0), _target() {}

        Item(const InternalEdgeRef& target, int refCount) : _refCount(refCount), _target(target) {}

        int             _refCount;
        InternalEdgeRef _target;
    };

    std::vector<Item>             _targets;
    std::map<Node*, unsigned int> _indexMap;
    std::stack<unsigned int>      _freeIndexes;

    friend class Manager;
    friend class RewritingManager;
    friend class FRAIGManager;
};
}  // namespace aigpp

#endif /* AIGPP_EXTREFTABLE_HH */
