/**************************************************************
 *
 *       AIGPP // DynamicRewritingTable.cc
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

#include "DynamicRewritingTable.hh"

#include <iostream>
#include <set>
#include <stack>

#include <lrabsutil/SimpleStack.hh>

aigpp::DynamicRewritingTable::DynamicRewritingTable() :
    _nextIndex(4), /* reserve space for four inputs */
    _refMap(4, 1), /* add four non deletable inputs */
    _refs(0),
    _derefs(0)
{}

aigpp::DynamicRewritingTable::~DynamicRewritingTable()
{
    /* QUIET version */
#if 0
    if( _refs != _derefs )
    {
        std::cout << "refs " << _refs << " " << _derefs << std::endl;
    }
#endif

    /* no and nodes must be left */
    assert(_andLookupMap.empty());
    /* all allocated indexes must be freed (except the four inputs) */
    assert((int)_freeIndexes.size() + 4 == _nextIndex);
}

bool
aigpp::DynamicRewritingTable::isInverted(const aigpp::DynamicRewritingTable::Index& i)
{
    return ((i & 1u) == 1u);
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::getNonInvertedIndex(const aigpp::DynamicRewritingTable::Index& i)
{
    return i & ~1u;
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::invert(const aigpp::DynamicRewritingTable::Index& i)
{
    return i ^ 1u;
}

bool
aigpp::DynamicRewritingTable::isInput(const aigpp::DynamicRewritingTable::Index& i)
{
    return i < 8;
}

unsigned int
aigpp::DynamicRewritingTable::getInput(const aigpp::DynamicRewritingTable::Index& i)
{
    assert(isInput(i));
    return i / 2;
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::lookupInput(unsigned int i) const
{
    assert(i < 4);
    return Index(i * 2);
}

const aigpp::DynamicRewritingTable::IndexPair&
aigpp::DynamicRewritingTable::getAnd(const aigpp::DynamicRewritingTable::Index& i) const
{
    assert(i / 2 >= 4 && i / 2 - 4 < (int)_andMap.size());

    return (_andMap[i / 2 - 4]);
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::lookupAnd(const aigpp::DynamicRewritingTable::Index& i1,
                                        const aigpp::DynamicRewritingTable::Index& i2) const
{
    IndexPair ip(i1, i2);

    /* canonize order */
    if (ip.first > ip.second) {
        std::swap(ip.first, ip.second);
    }

    /* search for pair */
    const auto p = _andLookupMap.find(ip);

    /* not found -> return InvalidIndex */
    if (p == _andLookupMap.cend()) {
        return InvalidIndex;
    }
    /* found -> return matching index */
    else {
        return p->second;
    }
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::insertAnd(const aigpp::DynamicRewritingTable::Index& i1,
                                        const aigpp::DynamicRewritingTable::Index& i2)
{
    ref(i1);
    ref(i2);

    IndexPair ip(i1, i2);

    /* canonize order */
    if (ip.first > ip.second) {
        std::swap(ip.first, ip.second);
    }

    Index i = allocateIndex();

    if (i / 2 - 4 >= (int)_andMap.size()) {
        assert(i / 2 - 4 == (int)_andMap.size());
        _andMap.push_back(ip);
    } else {
        _andMap[i / 2 - 4] = ip;
    }

    _andLookupMap[ip] = i;

    return i;
}

void
aigpp::DynamicRewritingTable::ref(const aigpp::DynamicRewritingTable::Index& i)
{
    ++_refs;

    assert(i / 2 >= 0 && i / 2 < (int)_refMap.size());
    ++_refMap[i / 2];
}

void
aigpp::DynamicRewritingTable::deref(const aigpp::DynamicRewritingTable::Index& i)
{
    ++_derefs;

    assert(i / 2 >= 0 && i / 2 < (int)_refMap.size());
    int r = --_refMap[i / 2];
    assert(r >= 0);

    if (r == 0) {
        assert(!isInput(i));

        const IndexPair& ip = getAnd(i);
        _andLookupMap.erase(ip);

        _freeIndexes.push(i / 2);

        deref(ip.first);
        deref(ip.second);
    }
}

void
aigpp::DynamicRewritingTable::print(const aigpp::DynamicRewritingTable::Index& i) const
{
    std::stack<int> pending;
    std::set<int>   processed;

    std::cout << "impl: " << (isInverted(i) ? "!" : "") << i / 2 << std::endl;

    pending.push(i / 2);
    while (!pending.empty()) {
        if (processed.find(pending.top()) != processed.end()) {
            pending.pop();
            continue;
        }

        if (pending.top() < 4) {
            std::cout << pending.top() << ": variable" << std::endl;
            processed.insert(pending.top());
            pending.pop();
            continue;
        }

        const IndexPair& ip = getAnd(pending.top() * 2);

        if (processed.find(ip.first / 2) == processed.end()) {
            pending.push(ip.first / 2);
            continue;
        }

        if (processed.find(ip.second / 2) == processed.end()) {
            pending.push(ip.second / 2);
            continue;
        }

        std::cout << pending.top() << ": and " << (isInverted(ip.first) ? "!" : "") << ip.first / 2 << " "
                  << (isInverted(ip.second) ? "!" : "") << ip.second / 2 << std::endl;

        processed.insert(pending.top());
        pending.pop();
    }
}

aigpp::FourInputFunction
aigpp::DynamicRewritingTable::getFunction(const aigpp::DynamicRewritingTable::Index& i) const
{
    std::stack<int>             pending;
    std::map<int, unsigned int> mapping;

    mapping[0] = 0xAAAA;
    mapping[1] = 0xCCCC;
    mapping[2] = 0xF0F0;
    mapping[3] = 0xFF00;

    pending.push(i / 2);
    while (!pending.empty()) {
        if (mapping.find(pending.top()) != mapping.end()) {
            pending.pop();
            continue;
        }

        const IndexPair& ip = getAnd(pending.top() * 2);

        if (mapping.find(ip.first / 2) == mapping.end()) {
            pending.push(ip.first / 2);
            continue;
        }

        if (mapping.find(ip.second / 2) == mapping.end()) {
            pending.push(ip.second / 2);
            continue;
        }

        unsigned int s1 = mapping[ip.first / 2];
        if (isInverted(ip.first)) {
            s1 = ~s1;
        }
        unsigned int s2 = mapping[ip.second / 2];
        if (isInverted(ip.second)) {
            s2 = ~s2;
        }

        unsigned int s         = (s1 & s2);
        mapping[pending.top()] = s;

        pending.pop();
    }

    unsigned int result = mapping[i / 2];
    if (isInverted(i)) {
        result = ~result;
    }

    return result & 0xFFFF;
}

std::size_t
aigpp::DynamicRewritingTable::getSize(const aigpp::DynamicRewritingTable::Index& i) const
{
    std::set<int>   cone;
    std::stack<int> pending;

    pending.push(i / 2);
    while (!pending.empty()) {
        if (cone.find(pending.top()) != cone.end()) {
            pending.pop();
            continue;
        }

        if (pending.top() < 4) {
            cone.insert(pending.top());
            pending.pop();
            continue;
        }

        const IndexPair& ip = getAnd(pending.top() * 2);

        cone.insert(pending.top());
        pending.pop();

        if (cone.find(ip.first / 2) == cone.end()) {
            pending.push(ip.first / 2);
        }

        if (cone.find(ip.second / 2) == cone.end()) {
            pending.push(ip.second / 2);
        }
    }

    return cone.size();
}

int
aigpp::DynamicRewritingTable::usedInputsMask(const Index& i) const
{
    static lrabs::SimpleStack<int> pending;

    int mask = 0;

    pending.push(i / 2);
    while (!pending.empty()) {
        if (pending.top() < 4) {
            mask |= 1 << pending.top();
            pending.pop();
        } else {
            const IndexPair& ip = _andMap[pending.top() - 4];
            pending.pop();

            pending.push(ip.first / 2);
            pending.push(ip.second / 2);
        }
    }

    return mask;
}

void
aigpp::DynamicRewritingTable::clear()
{
    _nextIndex = 4;
    while (!_freeIndexes.empty()) {
        _freeIndexes.pop();
    }

    _andMap.clear();
    _andLookupMap.clear();
    _refMap.assign(4, 1);
    _refs   = 0;
    _derefs = 0;
}

aigpp::DynamicRewritingTable::Index
aigpp::DynamicRewritingTable::allocateIndex()
{
    /* allocate new index if no free is available */
    if (_freeIndexes.empty()) {
        assert((int)_refMap.size() == _nextIndex);

        _refMap.push_back(0);

        Index i = 2 * _nextIndex;
        ++_nextIndex;

        return i;
    }
    /* return free index */
    else {
        Index i = 2 * _freeIndexes.top();
        assert(_refMap[_freeIndexes.top()] == 0);

        _freeIndexes.pop();

        return i;
    }
}
