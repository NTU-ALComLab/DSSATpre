/**************************************************************
 *
 *       AIGPP // DynamicRewritingTable.hh
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

#ifndef AIGPP_DYNAMICREWRITINGTABLE_HH
#define AIGPP_DYNAMICREWRITINGTABLE_HH

#include "FourInputFunction.hh"

#include <cassert>
#include <map>
#include <stack>
#include <vector>

namespace aigpp {
class DynamicRewritingTable
{
   public:
    DynamicRewritingTable();
    DynamicRewritingTable(DynamicRewritingTable&& other)      = default;
    DynamicRewritingTable(const DynamicRewritingTable& other) = delete;
    ~DynamicRewritingTable();
    DynamicRewritingTable& operator=(DynamicRewritingTable&& other) = default;
    DynamicRewritingTable& operator=(const DynamicRewritingTable& other) = delete;

    typedef int Index;
    enum
    {
        InvalidIndex = -1
    };

    typedef std::pair<Index, Index> IndexPair;

    static bool  isInverted(const Index& i);
    static Index getNonInvertedIndex(const Index& i);
    static Index invert(const Index& i);

    static bool         isInput(const Index& i);
    static unsigned int getInput(const Index& i);
    Index               lookupInput(unsigned int input) const;

    const IndexPair& getAnd(const Index& i) const;
    Index            lookupAnd(const Index& i1, const Index& i2) const;
    /* no lookup is performed! */
    Index insertAnd(const Index& i1, const Index& i2);

    void ref(const Index& i);
    void deref(const Index& i);

    std::size_t size() const { return _andMap.size(); }

    void              print(const Index& i) const;
    FourInputFunction getFunction(const Index& i) const;
    std::size_t       getSize(const Index& i) const;

    int usedInputsMask(const Index& i) const;

    void clear();

   private:
    typedef std::vector<IndexPair>     AndMap;
    typedef std::map<IndexPair, Index> AndLookupMap;
    typedef std::vector<int>           RefMap;

    Index allocateIndex();

    Index             _nextIndex;
    std::stack<Index> _freeIndexes;

    AndMap       _andMap;
    AndLookupMap _andLookupMap;
    RefMap       _refMap;

    int _refs, _derefs;

    friend class Manager;
};
}  // namespace aigpp

#endif /* AIGPP_DYNAMICREWRITINGTABLE_HH */
