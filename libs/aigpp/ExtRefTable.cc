/**************************************************************
 *
 *       AIGPP Package // ExtRefTable.cc
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

#include "ExtRefTable.hh"

aigpp::ExtRefTable::~ExtRefTable()
{
    clear();
}

int
aigpp::ExtRefTable::insert(const aigpp::InternalEdgeRef& target)
{
    /* special handling for constant edges */
    if (target.isConstant()) {
        /* false -> index=-1 */
        if (!(target.isInverted())) {
            return -1;
        }
        /* true -> index = -2 */
        else {
            return -2;
        }
    }

    /* look for existing target */
    std::map<Node*, unsigned int>::const_iterator p = _indexMap.find(target._node);
    if (p != _indexMap.end()) {
        Item& t = _targets[p->second];
        ++(t._refCount);

        return p->second;
    }

    /* insert new target if not found */
    if (!_freeIndexes.empty()) {
        unsigned int f = _freeIndexes.top();
        _freeIndexes.pop();

        Item& t = _targets[f];
        assert(t._refCount == 0);

        t._target = target;
        ++(t._refCount);

        _indexMap[target._node] = f;

        return f;
    } else {
        _targets.emplace_back(target, 1);
        _indexMap[target._node] = static_cast<unsigned int>(_targets.size()) - 1;
        return static_cast<unsigned int>(_targets.size()) - 1;
    }
}

const aigpp::InternalEdgeRef&
aigpp::ExtRefTable::lookup(int index) const
{
    /* special handling for constant edges (index -1 or -2 ) */
    if (index == -1) {
        return const0;
    } else if (index == -2) {
        return const1;
    } else {
        const Item& t = _targets[index];
        assert(t._refCount > 0);

        return t._target;
    }
}

void
aigpp::ExtRefTable::clear()
{
    for (std::size_t i = 0; i != _targets.size(); ++i) {
        Item& t = _targets[i];
        if (t._refCount > 0) {
            t._refCount = 0;
            t._target.clear();

            _freeIndexes.push(static_cast<unsigned int>(i));
        }
    }

    _indexMap.clear();
}

void
aigpp::ExtRefTable::ref(int index)
{
    /* ignore constant nodes ( negative index ) */
    if (index >= 0) {
        Item& t = _targets[index];

        assert(t._refCount >= 0);

        /* QUIET version */
#if 0
        if( t._refCount == 0 )
        {
            std::cout << "WARNING: ExtRefTable: referencing dead node" << std::endl;
        }
#endif

        ++(t._refCount);
    }
}

void
aigpp::ExtRefTable::deref(int index)
{
    /* ignore constant nodes ( negative index ) */
    if (index >= 0) {
        Item& t = _targets[index];
        assert(t._refCount >= 1);

        --(t._refCount);

        if (t._refCount == 0) {
            _freeIndexes.push(index);
            _indexMap.erase(t._target._node);
            t._target.clear();
        }
    }
}

void
aigpp::ExtRefTable::rebuildMapping()
{
    _indexMap.clear();
    for (std::size_t i = 0; i != _targets.size(); ++i) {
        if (_targets[i]._refCount == 0) {
            continue;
        }
        if (_targets[i]._target.isConstant()) {
            continue;
        }

        _indexMap[_targets[i]._target._node] = static_cast<unsigned int>(i);
    }
}
