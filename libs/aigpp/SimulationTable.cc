/**************************************************************
 *
 *       AIGPP Package // SimulationTable.cc
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

#include "SimulationTable.hh"

#include "Node.hh"
#include "Primes.hh"
#include "SimVector.hh"

#define _unused(x) ((void)(x))

aigpp::SimulationTable::SimulationTable() : _capacity(PowerPrimes[0]), _entries(0), _primesIndex(0)
{
    _table = new Node*[_capacity];
    std::fill(_table, _table + _capacity, (Node*)nullptr);
}

aigpp::SimulationTable::~SimulationTable()
{
    assert(_entries == 0);
    delete[] _table;
}

aigpp::Node*
aigpp::SimulationTable::getClass(aigpp::Node* n) const
{
    assert(n->_nEqualSim != 0 && n->_pEqualSim != 0);

    Node* simclass = n;

    while (simclass->_nEqualSimHash == nullptr) {
        simclass = simclass->_nEqualSim;
        if (simclass == n) {
            break;
        }
    }

    assert(simclass->_pEqualSimHash != 0);

    return simclass;
}

aigpp::Node*
aigpp::SimulationTable::getInvertedClass(aigpp::Node* nonInvertedClass) const
{
    return getClass(nonInvertedClass)->_invertedSimClass;
}

aigpp::Node*
aigpp::SimulationTable::lookup(const aigpp::SimVector& sim) const
{
    Node* simtableentry = _table[sim.hash() % _capacity];
    if (simtableentry == nullptr) {
        return nullptr;
    }

    Node* simclass = simtableentry;
    do {
        if (simclass->sim() == sim) {
            return simclass;
        }

        simclass = simclass->_nEqualSimHash;
    } while (simclass != simtableentry);

    return nullptr;
}

aigpp::Node*
aigpp::SimulationTable::lookupInverted(const aigpp::SimVector& sim, aigpp::Node* nonInvertedClass) const
{
    if (nonInvertedClass != nullptr) {
        /* assure "nonInvertedClass" is class head */
        assert(nonInvertedClass->_nEqualSimHash != 0);
        return nonInvertedClass->_invertedSimClass;
    }

    Node* simtableentry = _table[sim.hashInverted() % _capacity];
    if (simtableentry == nullptr) {
        return nullptr;
    }

    Node* simclass = simtableentry;
    do {
        if (simclass->sim().equalInverted(sim)) {
            return simclass;
        }

        simclass = simclass->_nEqualSimHash;
    } while (simclass != simtableentry);

    return nullptr;
}

void
aigpp::SimulationTable::insert(aigpp::Node* n)
{
    n->_pEqualSimHash    = nullptr;
    n->_nEqualSimHash    = nullptr;
    n->_pEqualSim        = nullptr;
    n->_nEqualSim        = nullptr;
    n->_invertedSimClass = nullptr;

    // cast ok
    if (static_cast<double>(_entries) > 0.75 * static_cast<double>(_capacity)) {
        resize(_primesIndex + 1);
    }

    Node** simtableentry = _table + n->sim().hash() % _capacity;

    /*
     * no entry with an appropriate hash value
     */
    if (*simtableentry == nullptr) {
        ++_entries;

        *simtableentry    = n;
        n->_pEqualSimHash = n;
        n->_nEqualSimHash = n;
        n->_pEqualSim     = n;
        n->_nEqualSim     = n;

        Node* inverted       = lookupInverted(n->sim());
        n->_invertedSimClass = inverted;
        if (inverted != nullptr) {
            assert(inverted->_invertedSimClass == 0);
            inverted->_invertedSimClass = n;
        }

        return;
    }

    /*
     * entry with matching hash value found
     */

    /*
     * search for simulation class
     */
    Node* simclass = *simtableentry;
    do {
        if (simclass->sim() == n->sim()) {
            /*
             * simulation class found. insert as second element of simclass.
             */
            simclass->_nEqualSim->_pEqualSim = n;
            n->_pEqualSim                    = simclass;
            n->_nEqualSim                    = simclass->_nEqualSim;
            simclass->_nEqualSim             = n;

            return;
        }

        simclass = simclass->_nEqualSimHash;
    } while (simclass != *simtableentry);

    /*
     * no matching simulation class found.
     * create a new one and insert as second element of simtableentry.
     */
    ++_entries;

    n->_pEqualSim = n;
    n->_nEqualSim = n;

    (*simtableentry)->_nEqualSimHash->_pEqualSimHash = n;
    n->_pEqualSimHash                                = (*simtableentry);
    n->_nEqualSimHash                                = (*simtableentry)->_nEqualSimHash;
    (*simtableentry)->_nEqualSimHash                 = n;

    Node* inverted       = lookupInverted(n->sim());
    n->_invertedSimClass = inverted;
    if (inverted != nullptr) {
        assert(inverted->_invertedSimClass == 0);
        inverted->_invertedSimClass = n;
    }
}

void
aigpp::SimulationTable::clear()
{
    std::fill(_table, _table + _capacity, (Node*)nullptr);
    _entries = 0;
}

void
aigpp::SimulationTable::clearRememberClasses()
{
    Node** simtableentry = _table;

    for (std::size_t i = 0; i < _capacity; ++i, ++simtableentry) {
        if (*simtableentry == nullptr) {
            continue;
        }

        Node* simclass = *simtableentry;

        do {
            Node* nextClass = simclass->_nEqualSimHash;

            Node* n = simclass;
            do {
                Node* next = n->_nEqualSim;

                n->_pEqualSim        = nullptr;
                n->_nEqualSim        = nullptr;
                n->_pEqualSimHash    = nullptr;
                n->_nEqualSimHash    = nullptr;
                n->_invertedSimClass = nullptr;

                n->_nEqualSim = simclass;

                n = next;
            } while (n != simclass);

            simclass = nextClass;
        } while (simclass != *simtableentry);

        *simtableentry = nullptr;
    }

    clear();
}

void
aigpp::SimulationTable::rebuild(aigpp::Node* nodes, bool useRememberedClasses)
{
    if (useRememberedClasses) {
        /* insert remembered class heads first */
        for (Node* n = nodes; n != nullptr; n = n->next()) {
            assert(n->_nEqualSim != 0);
            assert(n->_pEqualSim == 0);

            if (n->_nEqualSim == n) {
                n->_nEqualSim = nullptr;
                insert(n);

                assert(n->_pEqualSim != 0);
            }
        }

        /* insert remaining nodes */
        for (Node* n = nodes; n != nullptr; n = n->next()) {
            /* node already inserted by previous loop */
            if (n->_pEqualSim != nullptr) {
                continue;
            }

            /* remembered class still correct */
            Node* rememberedClass = n->_nEqualSim;
            assert(rememberedClass->_nEqualSim != 0);
            assert(rememberedClass->_pEqualSim != 0);

            n->_nEqualSim = nullptr;
            if (rememberedClass->_nEqualSim->sim() == n->sim()) {
                n->_pEqualSim                           = rememberedClass;
                n->_nEqualSim                           = rememberedClass->_nEqualSim;
                rememberedClass->_nEqualSim->_pEqualSim = n;
                rememberedClass->_nEqualSim             = n;

                n->_invertedSimClass = nullptr;
            } else {
                insert(n);
            }
        }
    } else {
        for (Node* n = nodes; n != nullptr; n = n->next()) {
            n->_nEqualSim        = nullptr;
            n->_pEqualSim        = nullptr;
            n->_nEqualSimHash    = nullptr;
            n->_pEqualSimHash    = nullptr;
            n->_invertedSimClass = nullptr;

            insert(n);
        }
    }
}

void
aigpp::SimulationTable::resize(int newPrimesIndex)
{
    if (newPrimesIndex < 0 || newPrimesIndex >= NumPowerPrimes) {
        return;
    }

    std::size_t oldcapacity = _capacity;

    _primesIndex = newPrimesIndex;
    _capacity    = PowerPrimes[_primesIndex];

    auto** newtable = new Node*[_capacity];
    std::fill(newtable, newtable + _capacity, (Node*)nullptr);

    Node** simtableentry = _table;

    for (std::size_t i = 0; i < oldcapacity; ++i, ++simtableentry) {
        if (*simtableentry == nullptr) {
            continue;
        }

        std::vector<Node*> classes;

        Node* simclass = *simtableentry;
        do {
            classes.push_back(simclass);
            simclass = simclass->_nEqualSimHash;
        } while (simclass != *simtableentry);

        for (Node* n : classes) {
            Node* nn = n;

            Node** newsimtableentry = newtable + (nn->sim().hash() % _capacity);
            if (*newsimtableentry == nullptr) {
                nn->_nEqualSimHash = nn;
                nn->_pEqualSimHash = nn;
                *newsimtableentry  = nn;
            } else {
                nn->_nEqualSimHash                                  = (*newsimtableentry)->_nEqualSimHash;
                (*newsimtableentry)->_nEqualSimHash->_pEqualSimHash = nn;
                nn->_pEqualSimHash                                  = (*newsimtableentry);
                (*newsimtableentry)->_nEqualSimHash                 = nn;
            }
        }
    }

    delete[] _table;
    _table = newtable;
}

void
aigpp::SimulationTable::garbageCollect()
{
    /*
     * go through the whole table
     */
    std::vector<Node*> classes;
    std::vector<Node*> classnodes;

    std::size_t newEntries = 0;
    std::size_t oldNodes = 0, newNodes = 0;

    for (std::size_t i = 0; i != _capacity; ++i) {
        Node* simtableentry = _table[i];
        if (simtableentry == nullptr) {
            continue;
        }

        /*
         * go through all simclasses of this table entry
         */

        /* collect simclasses */
        classes.clear();
        {
            Node* n = simtableentry;
            do {
                classes.push_back(n);
                n = n->_nEqualSimHash;
            } while (n != simtableentry);
        }

        _table[i] = nullptr;

        for (Node* simclass : classes) {
            /* collect classnodes */
            classnodes.clear();
            {
                Node* n = simclass;
                do {
                    classnodes.push_back(n);
                    n = n->_nEqualSim;
                } while (n != simclass);
            }

            /* filter out dead nodes and fill "newclass" */
            Node* inverted = nullptr;
            Node* newclass = nullptr;
            for (Node* nn : classnodes) {
                ++oldNodes;

                Node* n = nn;

                if (n->_invertedSimClass != nullptr) {
                    assert(n == simclass);
                    inverted             = n->_invertedSimClass;
                    n->_invertedSimClass = nullptr;
                }

                /* dead node */
                if (n->refCount() == 0) {
                    n->_nEqualSim     = nullptr;
                    n->_pEqualSim     = nullptr;
                    n->_nEqualSimHash = nullptr;
                    n->_pEqualSimHash = nullptr;
                }
                /* alive node */
                else {
                    ++newNodes;

                    /* first node of class -> create newclass */
                    if (newclass == nullptr) {
                        newclass = n;

                        n->_nEqualSimHash = nullptr;
                        n->_pEqualSimHash = nullptr;

                        n->_nEqualSim = n;
                        n->_pEqualSim = n;
                    }
                    /* other nodes already in class -> insert after head */
                    else {
                        n->_nEqualSimHash = nullptr;
                        n->_pEqualSimHash = nullptr;

                        n->_pEqualSim                    = newclass;
                        n->_nEqualSim                    = newclass->_nEqualSim;
                        newclass->_nEqualSim->_pEqualSim = n;
                        newclass->_nEqualSim             = n;
                    }
                }
            }

            if (newclass != nullptr) {
                ++newEntries;

                /* first class in hash table entry -> set as head */
                if (_table[i] == nullptr) {
                    newclass->_nEqualSimHash = newclass;
                    newclass->_pEqualSimHash = newclass;

                    _table[i] = newclass;
                }
                /* other classes already in hash table entry -> insert after head */
                else {
                    newclass->_pEqualSimHash                  = _table[i];
                    newclass->_nEqualSimHash                  = _table[i]->_nEqualSimHash;
                    _table[i]->_nEqualSimHash->_pEqualSimHash = newclass;
                    _table[i]->_nEqualSimHash                 = newclass;
                }

                if (inverted != nullptr) {
                    inverted->_invertedSimClass = newclass;
                    newclass->_invertedSimClass = inverted;
                } else {
                    newclass->_invertedSimClass = nullptr;
                }
            } else {
                if (inverted != nullptr) {
                    inverted->_invertedSimClass = nullptr;
                }
            }
        }
    }

    _entries = newEntries;
}

std::size_t
aigpp::SimulationTable::simClassSize(aigpp::Node* simclass)
{
    std::size_t size = 0;

    Node* n = simclass;
    do {
        ++size;
        n = n->_nEqualSim;
    } while (n != simclass);

    return size;
}

std::size_t
aigpp::SimulationTable::maxSimClassSize() const
{
    std::size_t maxsize = 0;

    for (std::size_t i = 0; i != _capacity; ++i) {
        Node* simtableentry = _table[i];

        if (simtableentry == nullptr) {
            continue;
        }

        Node* simclass = simtableentry;
        do {
            std::size_t size = simClassSize(simclass);
            if (size > maxsize) {
                maxsize = size;
            }

            simclass = simclass->_nEqualSimHash;
        } while (simclass != simtableentry);
    }

    return maxsize;
}

void
aigpp::SimulationTable::checkIntegrity(aigpp::Node* nodeslist) const
{
    _unused(nodeslist);  // Is only used in debug mode.

    std::size_t simClassesFound = 0;
    std::size_t nodesFound      = 0;

    Node** simtableentry = _table;

    for (std::size_t i = 0; i < _capacity; ++i, ++simtableentry) {
        if (*simtableentry == nullptr) {
            continue;
        }

        Node* simclass = *simtableentry;
        assert((simclass->sim().hash() % _capacity) == i);

        do {
            ++simClassesFound;

            if (simclass->_invertedSimClass != nullptr) {
                assert(simclass->_invertedSimClass->_invertedSimClass == simclass);
            }

            Node* n = simclass;
            do {
                ++nodesFound;

                if (n != simclass) {
                    assert(n->_nEqualSimHash == 0);
                    assert(n->_pEqualSimHash == 0);
                    assert(n->_invertedSimClass == 0);
                }

                assert(n->sim() == simclass->sim());

                n = n->_nEqualSim;
            } while (n != simclass);

            simclass = simclass->_nEqualSimHash;
        } while (simclass != *simtableentry);
    }

    assert(_entries == simClassesFound);

#ifndef NDEBUG
    for (Node* n = nodeslist; n != 0; n = n->next()) {
        Node* lookupResult   = lookup(n->sim());
        Node* getclassResult = getClass(n);

        assert(lookupResult != 0);
        assert(lookupResult == getclassResult);
    }
#endif
}
