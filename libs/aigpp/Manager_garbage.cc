/**************************************************************
 *
 *       AIGPP Package // Manager_garbage.cc
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
 *
 ***************************************************************/

#include "Manager.hh"

/* lrabs utilities */
#include <lrabsutil/Resources.hh>
#include <lrabsutil/SimpleStack.hh>

/* local */
#include "StatManager.hh"

void
aigpp::Manager::reclaim(aigpp::Node* n)
{
    assert(n != nullptr);
    assert(n->refCount() == 0);

    std::size_t initialDead = _dead;

    static lrabs::SimpleStack<Node*> stk;
    stk.push(n);

    do {
        if (n->refCount() == 0) {
            assert(!(n->isVar()));

            n->ref();

            StatManager::instance().incReclaimed();
            --_dead;
            ++_nodeCount;

            stk.push(n->parent1().node());
            n = n->parent2().node();
        } else {
            n->ref();
            n = stk.top();
            stk.pop();
        }
    } while (!(stk.empty()));

    _reclaimed += (initialDead - _dead);

    assert(n->refCount() == 1);
    n->deref();
    --_nodeCount;
}

void
aigpp::Manager::reclaimIfDead(aigpp::Node* n)
{
    if (n != nullptr) {
        if (n->refCount() == 0) {
            reclaim(n);
        }
    }
}

void
aigpp::Manager::recursiveDeref(aigpp::Node* n)
{
    assert(n != nullptr);
    assert(n->refCount() != 0);

    if (n->refCount() > 1) {
        n->deref();
        return;
    }

    static lrabs::SimpleStack<Node*> stk;
    /* make sure we don't get to this point twice! (may happen, wenn recursively
     * deref caches) */
    assert(stk.size() == 0);
    stk.push(0);

    do {
        assert(n->refCount() != 0);

        if (n->refCount() == 1) {
            assert(!(n->isVar()));

            if (n->isCacheValid()) {
                /* schedule cache entry for dereferencing */
                Node* c = _cache.lookup(n).node();
                if (c != nullptr) {
                    c->ref();
                    stk.push(c);
                }

                /* delete cache entry (deref if called recursively!) */
                _cache.remove(n);
            }
            if (n->isZCacheValid()) {
                /* schedule cache entry for dereferencing */
                const std::pair<InternalEdgeRef, InternalEdgeRef>& p = _cacheZ.lookup(n);

                if (p.first.node() != nullptr) {
                    p.first.node();
                    stk.push(p.first.node());
                }

                if (p.second.node() != nullptr) {
                    p.second.node();
                    stk.push(p.second.node());
                }

                /* delete cache entry (deref if called recursively!) */
                _cacheZ.remove(n);
            }

            n->_refCount = 0;
            ++_dead;
            --_nodeCount;

            StatManager::instance().updateMaxDead(_dead);

            stk.push(n->parent1().node());
            n = n->parent2().node();
        } else {
            n->deref();
            n = stk.top();
            stk.pop();
        }
    } while (!(stk.empty()));
}

void
aigpp::Manager::garbageCollect()
{
    StatManager::instance().incGarbageCollections();

    const double startTime = lrabs::cpuTime();

    _unique.garbageCollect();
#ifdef USE_COMPUTEDTABLE
    _computed.garbageCollect();
#endif
    _simTable.garbageCollect();

    _lastNode = nullptr;

    Node* prev = nullptr;
    Node* n    = _nodes;
    while (n != nullptr) {
        if (n->refCount() == 0) {
            if (!(n->isVar())) {
                --_dead;
            }

            if (!(n->flag<Node::FLAG_ISREDUCED>())) {
                --_unreducedNodes;
            }

            StatManager::instance().incDeleted();

            assert(!(n->isCacheValid()));
            assert(!(n->isZCacheValid()));

            if (n->isSatVarValid()) {
                ++_deletedNodesInCNF;
                _nodesWithSATVars.erase(n);
            }

            Node* next = n->next();
            if (prev != nullptr) {
                prev->_nextNode = next;
            } else {
                _nodes = next;
            }

            releaseNode(n);

            n = next;
        } else {
            _lastNode = n;
            prev      = n;
            n         = n->next();
        }
    }

    StatManager::instance().incGarbageCollectionTime(lrabs::cpuTime() - startTime);
}
