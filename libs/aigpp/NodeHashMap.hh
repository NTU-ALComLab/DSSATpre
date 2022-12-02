/**************************************************************
 *
 *       AIGPP Package // NodeHashMap.hh
 *
 *       Copyright (C) 2007 Florian Pigorsch
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

#ifndef AIGPP_NODEHASHMAP_HH
#define AIGPP_NODEHASHMAP_HH

/* std */
#include <stack>
#include <utility>

/* LRABS Utilities */
#include <lrabsutil/Assert.hh>

/* local */
#include "InternalEdgeRef.hh"
#include "Node.hh"
#include "Primes.hh"

namespace aigpp {
template <typename VALUE>
struct NodeHashMapEntry
{
    NodeHashMapEntry(Node* k, const VALUE& v, NodeHashMapEntry* n) : key(k), value(v), next(n){};

    void clearValue();

    Node* key;
    VALUE value;

    NodeHashMapEntry* next;

   private:
    /* copying is not allowed */
    NodeHashMapEntry(const NodeHashMapEntry&);
    NodeHashMapEntry& operator=(const NodeHashMapEntry&);
};

template <typename VALUE, int FLAG, int CLEARFLAG = 0>
class NodeHashMap
{
   public:
    NodeHashMap();
    NodeHashMap(const NodeHashMap&) = delete;
    NodeHashMap(NodeHashMap&&)      = delete;
    NodeHashMap& operator=(const NodeHashMap&) = delete;
    NodeHashMap& operator=(NodeHashMap&&) = delete;
    ~NodeHashMap();

    void         insert(Node* key, const VALUE& v);
    const VALUE& lookup(Node* key) const;
    void         remove(Node* key);
    void         clear();
    void         resize(int newPrimesIndex);
    std::size_t  index(Node* key) const;

   private:
    typedef NodeHashMapEntry<VALUE> Entry;

    int         _primesIndex;
    std::size_t _tableSize;
    Entry**     _table;

    std::size_t _entries;

    std::stack<Entry*> _freeEntries;
};
}  // namespace aigpp

#include "NodeHashMap.icc"

#endif /* AIGPP_NODEHASHMAP_HH */
