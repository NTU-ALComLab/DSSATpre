/**************************************************************
 *
 *       AIGPP Package // ComputedTable.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *	Author:
 *         Florian Pigorsch
 *	  University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *	Last revision:
 *         $Revision: 509 $
 *	  $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "ComputedTable.hh"

#ifdef USE_COMPUTEDTABLE

#    include <lrabsutil/Assert.hh>

#    include "Primes.hh"
#    include "StatManager.hh"

aigpp::ComputedTable::ComputedTable() :
    _primesIndex(0),
    _tableSize(PowerPrimes[_primesIndex]),
    _table(new Entry*[_tableSize]),
    _entries(0),
    _maxSize(10000),
    _freeEntries()
{
    std::fill(_table, _table + _tableSize, (Entry*)0);
}

aigpp::ComputedTable::~ComputedTable()
{
    Entry** e = _table;
    for (long i = 0; i != _tableSize; ++i, ++e) {
        if (*e != nullptr) {
            delete *e;
        }
    }

    while (!_freeEntries.empty()) {
        delete _freeEntries.top();
        _freeEntries.pop();
    }

    delete[] _table;
}

void
aigpp::ComputedTable::insert(aigpp::ComputedTable::Operation op, const aigpp::Edge& op1, const aigpp::Edge& op2,
                             const aigpp::Edge& result)
{
    switch (op) {
        case OP_AND:
            insertAnd(op1, op2, result);
            return;

        case OP_COFACTOR:
            insertCofactor(op1, op2, result);
            return;

        default:
            NEVER_GET_HERE;
    }
}

void
aigpp::ComputedTable::insertAnd(const aigpp::Edge& op1, const aigpp::Edge& op2, const aigpp::Edge& result)
{
    StatManager::instance().incComputedInsertions();

    grow();

    Entry** tableentry = _table + indexAnd(op1, op2);
    if (*tableentry != nullptr) {
        StatManager::instance().incComputedCollisions();
        (*tableentry)->set(OP_AND, op1, op2, result);
    } else {
        if (_freeEntries.empty()) {
            *tableentry = new Entry(OP_AND, op1, op2, result);
        } else {
            *tableentry = _freeEntries.top();
            _freeEntries.pop();
            (*tableentry)->set(OP_AND, op1, op2, result);
        }

        ++_entries;
    }
}

void
aigpp::ComputedTable::insertCofactor(const aigpp::Edge& n, const aigpp::Edge& var, const aigpp::Edge& result)
{
    assert(!(n.isConstant()));
    assert(!(var.isConstant()));
    assert(var.node()->isVar());

    StatManager::instance().incComputedInsertions();

    grow();

    Entry** tableentry = _table + indexCofactor(n.node(), var);
    if (*tableentry != nullptr) {
        StatManager::instance().incComputedCollisions();

        if (n.isInverted()) {
            (*tableentry)->set(OP_COFACTOR, !n, var, !result);
        } else {
            (*tableentry)->set(OP_COFACTOR, n, var, result);
        }
    } else {
        if (_freeEntries.empty()) {
            if (n.isInverted()) {
                *tableentry = new Entry(OP_COFACTOR, !n, var, !result);
            } else {
                *tableentry = new Entry(OP_COFACTOR, n, var, result);
            }
        } else {
            *tableentry = _freeEntries.top();
            _freeEntries.pop();

            if (n.isInverted()) {
                (*tableentry)->set(OP_COFACTOR, !n, var, !result);
            } else {
                (*tableentry)->set(OP_COFACTOR, n, var, result);
            }
        }

        ++_entries;
    }
}

bool
aigpp::ComputedTable::lookup(aigpp::ComputedTable::Operation op, const aigpp::Edge& op1, const aigpp::Edge& op2,
                             aigpp::Edge& result) const
{
    switch (op) {
        case OP_AND:
            return lookupAnd(op1, op2, result);

        case OP_COFACTOR:
            return lookupCofactor(op1, op2, result);

        default:
            NEVER_GET_HERE;
            return false;
    }
}

bool
aigpp::ComputedTable::lookupAnd(const aigpp::Edge& op1, const aigpp::Edge& op2, aigpp::Edge& result) const
{
    Entry** tableentry = _table + indexAnd(op1, op2);
    if (*tableentry != nullptr) {
        if ((*tableentry)->matches(OP_AND, op1, op2)) {
            result = (*tableentry)->result;

            StatManager::instance().incComputedHitsAnd();
            return true;
        } else {
            StatManager::instance().incComputedMissAnd();
            return false;
        }
    } else {
        StatManager::instance().incComputedMissAnd();
        return false;
    }
}

bool
aigpp::ComputedTable::lookupCofactor(const aigpp::Edge& n, const aigpp::Edge& var, aigpp::Edge& result) const
{
    Entry** tableentry = _table + indexCofactor(n.node(), var);
    if (*tableentry != 0) {
        if ((*tableentry)->matches(OP_COFACTOR, n.notIf(n.isInverted()), var)) {
            result = (*tableentry)->result.notIf(n.isInverted());

            StatManager::instance().incComputedHitsCofactor();
            return true;
        } else {
            StatManager::instance().incComputedMissCofactor();
            return false;
        }
    } else {
        StatManager::instance().incComputedMissCofactor();
        return false;
    }
}

void
aigpp::ComputedTable::garbageCollect()
{
    Entry** e = _table;

    for (long i = 0; i != _tableSize; ++i, ++e) {
        if (*e == 0) {
            continue;
        }

        if ((!(*e)->result.isConstant() && (*e)->result.node()->refCount() == 0) || (*e)->op1.node()->refCount() == 0
            || (*e)->op2.node()->refCount() == 0) {
            _freeEntries.push(*e);
            *e = 0;
            --_entries;
        }
    }
}

void
aigpp::ComputedTable::clear()
{
    Entry** e = _table;

    for (long i = 0; i != _tableSize; ++i, ++e) {
        if (*e == 0) {
            continue;
        }

        _freeEntries.push(*e);
        *e = 0;
    }

    _entries = 0;
}

void
aigpp::ComputedTable::grow()
{
    if (static_cast<double>(_entries) < 0.75 * _tableSize || _primesIndex + 1 >= NumPowerPrimes) {
        return;
    }

    if (_tableSize > _maxSize) return;

    long oldTableSize = _tableSize;

    ++_primesIndex;
    _tableSize = PowerPrimes[_primesIndex];

    Entry** newtable = new Entry*[_tableSize];
    std::fill(newtable, newtable + _tableSize, (Entry*)0);

    _entries = 0;

    Entry** tableentry = _table;
    for (long i = 0; i != oldTableSize; ++i, ++tableentry) {
        if (*tableentry == 0) {
            continue;
        }

        StatManager::instance().incComputedInsertions();

        Entry*  entry         = *tableentry;
        Entry** newtableentry = newtable + index(entry->op, entry->op1, entry->op2);

        if (*newtableentry == 0) {
            *newtableentry = entry;
            ++_entries;
        } else {
            StatManager::instance().incComputedCollisions();

            (*newtableentry)->set(entry->op, entry->op1, entry->op2, entry->result);
            _freeEntries.push(entry);
        }
    }

    delete[] _table;
    _table = newtable;
}

#endif /* USE_COMPUTEDTABLE */
