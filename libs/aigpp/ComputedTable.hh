/**************************************************************
 *
 *       AIGPP Package // ComputedTable.hh
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 379 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#ifndef AIGPP_COMPUTEDTABLE_HH
#define AIGPP_COMPUTEDTABLE_HH

//#define USE_COMPUTEDTABLE

#ifdef USE_COMPUTEDTABLE
#    include <stack>

#    include "Edge.hh"
#    include "Node.hh"

namespace aigpp {

class ComputedTable
{
   public:
    ComputedTable();
    ComputedTable(const ComputedTable&)     = delete;
    ComputedTable(ComputedTable&&) noexcept = delete;
    ~ComputedTable();
    ComputedTable& operator=(const ComputedTable&) = delete;
    ComputedTable& operator=(ComputedTable&&) = delete;

    enum Operation
    {
        OP_AND,
        OP_COFACTOR
    };

    void setMaxSize(long m) { _maxSize = m; };

    void insert(Operation op, const Edge& op1, const Edge& op2, const Edge& result);

    void insertAnd(const Edge& op1, const Edge& op2, const Edge& result);

    void insertCofactor(const Edge& n, const Edge& var, const Edge& result);

    bool lookup(Operation op, const Edge& op1, const Edge& op2, Edge& result) const;

    bool lookupAnd(const Edge& op1, const Edge& op2, Edge& result) const;

    bool lookupCofactor(const Edge& n, const Edge& var, Edge& result) const;

    void garbageCollect();

    void clear();

   private:
    void grow();

    long indexAnd(const Edge& op1, const Edge& op2) const { return ((op1.hash() + 3 * op2.hash()) % _tableSize); };

    long indexCofactor(Node* n, const Edge& var) const
    {
        assert(n != 0);

        return (n->index() + var.hash() * 3) % _tableSize;
    };

    long index(Operation op, const Edge& op1, const Edge& op2) const
    {
        switch (op) {
            case OP_AND:
                return indexAnd(op1, op2);
                break;
            case OP_COFACTOR:
                assert(!(op1.isInverted()));
                return indexCofactor(op1.node(), op2);
                break;
            default:
                /* can not get here */
                return 0;
        }
    };

    class Entry
    {
       public:
        Entry(Operation op_, const Edge& op1_, const Edge& op2_, const Edge& result_) :
            op(op_),
            op1(op1_),
            op2(op2_),
            result(result_){};

        void set(Operation op_, const Edge& op1_, const Edge& op2_, const Edge& result_)
        {
            op     = op_;
            op1    = op1_;
            op2    = op2_;
            result = result_;
        };

        bool matches(Operation op_, const Edge& op1_, const Edge& op2_) const
        {
            return (op == op_ && op1.structurallyEquivalent(op1_) && op2.structurallyEquivalent(op2_));
        };

        Operation op;
        Edge      op1;
        Edge      op2;
        Edge      result;
    };

    long    _primesIndex;
    long    _tableSize;
    Entry** _table;

    long _entries;

    long _maxSize;

    std::stack<Entry*> _freeEntries;
};

}  // namespace aigpp

#endif /* USE_COMPUTEDTABLE */

#endif /* AIGPP_COMPUTEDTABLE_HH */
