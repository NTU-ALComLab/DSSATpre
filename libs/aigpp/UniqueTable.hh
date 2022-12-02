/**************************************************************
 *
 *       AIGPP Package // UniqueTable.hh
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#ifndef AIGPP_UNIQUETABLE_HH
#define AIGPP_UNIQUETABLE_HH

#include <stack>

#include "Edge.hh"
#include "Node.hh"

namespace aigpp {
class UniqueTable
{
   public:
    UniqueTable();
    ~UniqueTable();

    /* copying is not allowed */
    UniqueTable(const UniqueTable&) = delete;
    UniqueTable& operator=(const UniqueTable&) = delete;

    /* moving is ok */
    UniqueTable(UniqueTable&&) noexcept;
    UniqueTable& operator=(UniqueTable&&) noexcept;

    bool lookup(const Edge& a, const Edge& b, Edge& result) const;

    void insert(Node* n);

    void remove(Node* n);

    void garbageCollect();
    void clear();

    void checkIntegrity(Node* nodeslist) const;

   private:
    void grow();

    unsigned long index(const Edge& a, const Edge& b) const { return ((a.hash() + b.hash()) % _tableSize); };

   private:
    std::size_t _primesIndex;
    std::size_t _tableSize;
    Node**      _table;

    std::size_t _entries;
};
}  // namespace aigpp

#endif /* AIGPP_UNIQUETABLE_HH */
