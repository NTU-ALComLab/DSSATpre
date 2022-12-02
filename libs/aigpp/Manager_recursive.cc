/**************************************************************
 *
 *       AIGPP Package // Manager_recuresive.cc
 *
 *       Copyright (C) 2006, 2007 Florian Pigorsch
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 373 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "Manager.hh"

#include "StatManager.hh"

std::vector<aigpp::EdgeRef>
aigpp::Manager::vectorCompose(const std::vector<aigpp::EdgeRef>& functions,
                              const std::vector<aigpp::EdgeRef>& variables,
                              const std::vector<aigpp::EdgeRef>& replacements)
{
    std::vector<aigpp::InternalEdgeRef> functions2    = EdgeRef::toInternal(functions);
    std::vector<aigpp::InternalEdgeRef> variables2    = EdgeRef::toInternal(variables);
    std::vector<aigpp::InternalEdgeRef> replacements2 = EdgeRef::toInternal(replacements);

    std::vector<InternalEdgeRef> results2 = vectorCompose(functions2, variables2, replacements2);

    return EdgeRef::fromInternal(getExtRefTable(), results2);
}

std::vector<aigpp::InternalEdgeRef>
aigpp::Manager::vectorCompose(const std::vector<aigpp::InternalEdgeRef>& functions,
                              const std::vector<aigpp::InternalEdgeRef>& variables,
                              const std::vector<aigpp::InternalEdgeRef>& replacements)
{
    assert(variables.size() == replacements.size());

    for (auto v = variables.begin(); v != variables.end(); ++v) {
        assert(!(v->isConstant()));
        assert(v->node()->isVar());
    }

    if (functions.empty()) {
        return functions;
    }

    /* build cache for given variables */
    auto r = replacements.begin();
    for (const auto& v : variables) {
        assert(!(v.node()->isCacheValid()));
        _cache.insert(v.node(), r->notIf(v.isInverted()));
    }

    /* build cache for remaining variables */
    for (Node* n : _variables) {
        if (!(n->isCacheValid())) {
            _cache.insert(n, InternalEdgeRef(n));
        }
    }

    /* perform substitutions */
    std::vector<InternalEdgeRef> results;
    results.reserve(functions.size());

    for (const auto& f : functions) {
        if (f.isConstant()) {
            results.push_back(f);
        } else {
            buildRecursiveUsingCache(f.node());
            results.push_back(_cache.lookup(f.node()).notIf(f.isInverted()));
        }
    }

    /* clear cache */
    clearNodeCaches();

    return results;
}

void
aigpp::Manager::buildRecursiveUsingCache(aigpp::Node* n)
{
    std::stack<Node*>& pending = _recursiveNodeStack;

    pending.push(n);
    ref(n);

    while (!pending.empty()) {
        /* skip finished nodes */
        if (pending.top()->isCacheValid()) {
            deref(pending.top());
            pending.pop();
        }

        /* node must be an "and" node after the above check */

        /* check if first parent is finished */
        else if (!(pending.top()->parent1().node()->isCacheValid())) {
            ref(pending.top()->parent1().node());
            pending.push(pending.top()->parent1().node());
        }

        /* check if second parent is finished */
        else if (!(pending.top()->parent2().node()->isCacheValid())) {
            ref(pending.top()->parent2().node());
            pending.push(pending.top()->parent2().node());
        }

        /* process node */
        else {
            /* get first parent's result */
            InternalEdgeRef c1 = _cache.lookup(pending.top()->parent1().node());
            c1.toggleInverted(pending.top()->parent1().isInverted());

            /* get second parent's result */
            InternalEdgeRef c2 = _cache.lookup(pending.top()->parent2().node());
            c2.toggleInverted(pending.top()->parent2().isInverted());

            const bool      nand   = pending.top()->isNAND();
            InternalEdgeRef result = c1 * c2;
            result.toggleInverted(nand);

            /* set the node's cache */
            if (!(pending.top()->isCacheValid())) {
                _cache.insert(pending.top(), result);
            } else {
                std::cout << "buildRecursiveUsingCache: triggering sat check" << std::endl;

                assert(functionallyEquivalent(_cache.lookup(pending.top()), result));
            }

            deref(pending.top());
            pending.pop();
        }
    }
}

void
aigpp::Manager::buildCofactorRecursiveUsingCache(aigpp::Node* n, const Edge& /* var */)
{
    std::stack<Node*>& pending = _recursiveNodeStack;

    pending.push(n);
    ref(n);

    while (!pending.empty()) {
        /* skip finished nodes */
        if (pending.top()->isCacheValid()) {
            deref(pending.top());
            pending.pop();
        }
        /* node must be an "and" node after the above check */

#if 1
        /* one parent's cache is const0 (modulo inversion) */
        else if ((pending.top()->parent1().node()->flag<Node::FLAG_CACHECONST1>()
                  && pending.top()->parent1().isInverted())
                 || (pending.top()->parent1().node()->flag<Node::FLAG_CACHECONST0>()
                     && !(pending.top()->parent1().isInverted()))
                 || (pending.top()->parent2().node()->flag<Node::FLAG_CACHECONST1>()
                     && pending.top()->parent2().isInverted())
                 || (pending.top()->parent2().node()->flag<Node::FLAG_CACHECONST0>()
                     && !(pending.top()->parent2().isInverted()))) {
            /* set the node's cache */
            assert(!(pending.top()->isCacheValid()));

            if (!(pending.top()->isNAND())) {
                _cache.insert(pending.top(), const0);
                pending.top()->setFlag<Node::FLAG_CACHECONST0>();
            } else {
                _cache.insert(pending.top(), const1);
                pending.top()->setFlag<Node::FLAG_CACHECONST1>();
            }

            deref(pending.top());
            pending.pop();
        }
#endif
        /* check if first parent is finished */
        else if (!(pending.top()->parent1().node()->isCacheValid())) {
            ref(pending.top()->parent1().node());
            pending.push(pending.top()->parent1().node());
        }

        /* check if second parent is finished */
        else if (!(pending.top()->parent2().node()->isCacheValid())) {
            ref(pending.top()->parent2().node());
            pending.push(pending.top()->parent2().node());
        }

        /* process node */
        else {
            /* get first parent's result */
            InternalEdgeRef c1 = _cache.lookup(pending.top()->parent1().node());
            c1.toggleInverted(pending.top()->parent1().isInverted());

            /* get second parent's result */
            InternalEdgeRef c2 = _cache.lookup(pending.top()->parent2().node());
            c2.toggleInverted(pending.top()->parent2().isInverted());

            const bool      nand   = pending.top()->isNAND();
            InternalEdgeRef result = c1 * c2;
            result.toggleInverted(nand);

            /* set the node's cache */
            if (!(pending.top()->isCacheValid())) {
                _cache.insert(pending.top(), result);
            } else {
                std::cout << "buildCofactorRecursiveUsingCache: triggering sat check" << std::endl;
                assert(functionallyEquivalent(_cache.lookup(pending.top()), result));
            }

            if (result.isConstant()) {
                /* 0 */
                if (!(result.isInverted())) {
                    pending.top()->setFlag<Node::FLAG_CACHECONST0>();
                }
                /* 1 */
                else {
                    pending.top()->setFlag<Node::FLAG_CACHECONST1>();
                }
            }

            deref(pending.top());
            pending.pop();
        }
    }
}

void
aigpp::Manager::buildRecursiveZ(aigpp::Node* n)
{
    std::stack<Node*>& pending = _recursiveNodeStack;

    pending.push(n);
    ref(n);

    while (!pending.empty()) {
        if (pending.top()->isZCacheValid()) {
            deref(pending.top());
            pending.pop();
        } else if (!pending.top()->parent1().node()->isZCacheValid()) {
            ref(pending.top()->parent1().node());
            pending.push(pending.top()->parent1().node());
        } else if (!pending.top()->parent2().node()->isZCacheValid()) {
            ref(pending.top()->parent2().node());
            pending.push(pending.top()->parent2().node());
        } else {
            std::pair<InternalEdgeRef, InternalEdgeRef> p1 = _cacheZ.lookup(pending.top()->parent1().node());
            std::pair<InternalEdgeRef, InternalEdgeRef> p2 = _cacheZ.lookup(pending.top()->parent2().node());

            if (pending.top()->parent1().isInverted()) {
                p1.first.toggleInverted();
                p1.second.toggleInverted();
                std::swap(p1.first, p1.second);
            }
            if (pending.top()->parent2().isInverted()) {
                p2.first.toggleInverted();
                p2.second.toggleInverted();
                std::swap(p2.first, p2.second);
            }

            const bool                                  nand = pending.top()->isNAND();
            std::pair<InternalEdgeRef, InternalEdgeRef> result(p1.first * p2.first, p1.second * p2.second);

            if (nand) {
                result.first.toggleInverted();
                result.second.toggleInverted();
                std::swap(result.first, result.second);
            }

            if (!(pending.top()->isZCacheValid())) {
                _cacheZ.insert(pending.top(), result);
            }

            deref(pending.top());
            pending.pop();
        }
    }
}
