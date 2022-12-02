/**************************************************************
 *
 *       AIGPP Package // UniqueTable.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *         Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "UniqueTable.hh"

/* LRABS Utilities */
#include <lrabsutil/Assert.hh>

/* local */
#include "Manager.hh"
#include "Primes.hh"
#include "StatManager.hh"

aigpp::UniqueTable::UniqueTable() :
    _primesIndex(0),
    _tableSize(PowerPrimes[_primesIndex]),
    _table(new Node*[_tableSize]),
    _entries(0)
{
    std::fill(_table, _table + _tableSize, (Node*)nullptr);
}

aigpp::UniqueTable::UniqueTable(aigpp::UniqueTable&& other) noexcept :
    _primesIndex(other._primesIndex),
    _tableSize(other._tableSize),
    _table(other._table),
    _entries(other._entries)
{
    other._table = nullptr;
}

aigpp::UniqueTable&
aigpp::UniqueTable::operator=(aigpp::UniqueTable&& other) noexcept
{
    if (&other != this) {
        _primesIndex = other._primesIndex;
        _tableSize   = other._tableSize;
        _table       = other._table;
        _entries     = other._entries;
        other._table = nullptr;
    }
    return *this;
}

aigpp::UniqueTable::~UniqueTable()
{
    delete[] _table;
}

bool
aigpp::UniqueTable::lookup(const aigpp::Edge& a, const aigpp::Edge& b, aigpp::Edge& result) const
{
    /* search collision chain of hash value (a, b) */
    /* iterate through chain members and find match */
    for (Node* entry = _table[index(a, b)]; entry != nullptr; entry = entry->_nextUniqueTableEntry) {
        if ((entry->parent1().structurallyEquivalent(a) && entry->parent2().structurallyEquivalent(b))
            || (entry->parent1().structurallyEquivalent(b) && entry->parent2().structurallyEquivalent(a))) {
            StatManager::instance().incUniqueHits();

            result = Edge(entry, entry->isNAND());
            return true;
        }
    }

    /* key pair (a, b) not found */
    return false;
}

void
aigpp::UniqueTable::insert(aigpp::Node* n)
{
    assert(n != 0);
    assert(!(n->isVar()));

    grow();

    ++_entries;

    Node** tableentry = _table + index(n->parent1(), n->parent2());

    n->_prevUniqueTableEntry = nullptr;
    n->_nextUniqueTableEntry = *tableentry;
    if (*tableentry != nullptr) {
        (*tableentry)->_prevUniqueTableEntry = n;
    }
    *tableentry = n;
}

void
aigpp::UniqueTable::remove(Node* n)
{
    assert(n != 0);
    assert(!(n->isVar()));

    Node** tableentry = _table + index(n->parent1(), n->parent2());

    assert(*tableentry != 0);

    /* first element of collision chain matches */
    if (*tableentry == n) {
        --_entries;

        *tableentry = n->_nextUniqueTableEntry;
        if (*tableentry != nullptr) {
            (*tableentry)->_prevUniqueTableEntry = nullptr;
        }
        n->_nextUniqueTableEntry = nullptr;
    }
    /* search chain for match */
    else {
        assert(n->_prevUniqueTableEntry != 0);

        n->_prevUniqueTableEntry->_nextUniqueTableEntry = n->_nextUniqueTableEntry;
        if (n->_nextUniqueTableEntry != nullptr) {
            n->_nextUniqueTableEntry->_prevUniqueTableEntry = n->_prevUniqueTableEntry;
        }

        --_entries;
    }
}

void
aigpp::UniqueTable::checkIntegrity(aigpp::Node* nodeslist) const
{
    /* check that every node is properly doubly linked */
    Node** tableentry = _table;
    for (std::size_t i = 0; i < _tableSize; ++i, ++tableentry) {
        for (Node* entry = *tableentry; entry != nullptr; entry = entry->_nextUniqueTableEntry) {
            if (entry != *tableentry) {
                assert(entry->_prevUniqueTableEntry != 0);
                assert(entry->_prevUniqueTableEntry->_nextUniqueTableEntry == entry);
            } else {
                assert(entry->_prevUniqueTableEntry == 0);
            }
        }
    }

    /* check that all nodes from the nodeslist are in the unique table */
    for (Node* n = nodeslist; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            continue;
        }

        Edge result;
        if (lookup(n->parent1(), n->parent2(), result)) {
            if (result.node() != n) {
                std::cout << "uniquetable: lookup of node " << n->index() << " returned node " << result.node()->index()
                          << std::endl;

                std::cout << *n << std::endl << *(result.node()) << std::endl;

                assert(false);
            }
        } else {
            std::cout << "node " << n->index() << " not in unique table" << std::endl;
            std::cout << "parent1 " << (n->parent1().isInverted() ? "!" : "") << n->parent1().node()->index()
                      << std::endl
                      << "parent2 " << (n->parent2().isInverted() ? "!" : "") << n->parent2().node()->index()
                      << std::endl;

            assert(false);
        }
    }

    /* check that all nodes in the unique table are actually contained in the
     * nodeslist */
    std::set<Node*> nodesset;
    for (Node* n = nodeslist; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            continue;
        }
        nodesset.insert(n);
    }

    tableentry = _table;
    for (std::size_t i = 0; i < _tableSize; ++i, ++tableentry) {
        for (Node* entry = *tableentry; entry != nullptr; entry = entry->_nextUniqueTableEntry) {
            if (nodesset.find(entry) == nodesset.end()) {
                std::cout << "node " << entry->index() << " found in the uniquetable but not in the aig!";
                assert(false);
            }
        }
    }
}

void
aigpp::UniqueTable::garbageCollect()
{
    Node** tableentry = _table;
    for (std::size_t i = 0; i < _tableSize; ++i, ++tableentry) {
        Node* entry = *tableentry;

        while (entry != nullptr) {
            if (entry->refCount() == 0) {
                --_entries;

                if (entry->_prevUniqueTableEntry == nullptr) {
                    *tableentry = entry->_nextUniqueTableEntry;
                    if (entry->_nextUniqueTableEntry != nullptr) {
                        entry->_nextUniqueTableEntry->_prevUniqueTableEntry = nullptr;
                    }
                } else {
                    entry->_prevUniqueTableEntry->_nextUniqueTableEntry = entry->_nextUniqueTableEntry;
                    if (entry->_nextUniqueTableEntry != nullptr) {
                        entry->_nextUniqueTableEntry->_prevUniqueTableEntry = entry->_prevUniqueTableEntry;
                    }
                }
            }

            entry = entry->_nextUniqueTableEntry;
        }
    }
}

void
aigpp::UniqueTable::clear()
{
    std::fill(_table, _table + _tableSize, (Node*)nullptr);
    _entries = 0;
}

void
aigpp::UniqueTable::grow()
{
    if (static_cast<double>(_entries) < 0.75 * static_cast<double>(_tableSize) || _primesIndex + 1 >= NumPowerPrimes) {
        return;
    }

    std::size_t oldTableSize = _tableSize;

    ++_primesIndex;
    _tableSize = PowerPrimes[_primesIndex];

    auto** newtable = new Node*[_tableSize];
    std::fill(newtable, newtable + _tableSize, (Node*)nullptr);

    Node** tableentry = _table;
    for (std::size_t i = 0; i < oldTableSize; ++i, ++tableentry) {
        Node* entry = *tableentry;

        while (entry != nullptr) {
            Node* next = entry->_nextUniqueTableEntry;

            Node** newtableentry = newtable + index(entry->parent1(), entry->parent2());

            entry->_prevUniqueTableEntry = nullptr;
            entry->_nextUniqueTableEntry = *newtableentry;
            if (*newtableentry != nullptr) {
                (*newtableentry)->_prevUniqueTableEntry = entry;
            }
            *newtableentry = entry;

            entry = next;
        }
    }

    delete[] _table;
    _table = newtable;
}
