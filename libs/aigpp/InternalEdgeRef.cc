/**************************************************************
 *
 *       AIGPP Package // InternalEdgeRef.cc
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

#include "InternalEdgeRef.hh"

#include <cassert>
#include <fstream>
#include <set>
#include <utility>

#include <lrabsutil/Resources.hh>

#include "EdgeRef.hh"
#include "Manager.hh"
#include "Node.hh"
#include "StatManager.hh"
#include "VarSupport.hh"

aigpp::InternalEdgeRef::InternalEdgeRef() : Edge()
{}

aigpp::InternalEdgeRef::InternalEdgeRef(const aigpp::InternalEdgeRef& n) : Edge(n)
{
    if (!isConstant()) {
        node()->manager()->ref(node());
    }
}

aigpp::InternalEdgeRef::InternalEdgeRef(const aigpp::Edge& ndp) : Edge(ndp)
{
    if (!isConstant()) {
        node()->manager()->ref(node());
    }
}

aigpp::InternalEdgeRef::InternalEdgeRef(aigpp::Node* n, bool i) : Edge(n, i)
{
    if (!isConstant()) {
        node()->manager()->ref(node());
    }
}

aigpp::InternalEdgeRef::~InternalEdgeRef()
{
    if (!isConstant()) {
        node()->manager()->deref(node());
    }
}

void
aigpp::InternalEdgeRef::clear()
{
    if (!isConstant()) {
        node()->manager()->deref(node());
    }
    Edge::clear();
}

aigpp::InternalEdgeRef&
aigpp::InternalEdgeRef::operator=(const aigpp::InternalEdgeRef& n)
{
    setNode(n.node());
    setInverted(n.isInverted());

    return *this;
}

void
aigpp::InternalEdgeRef::setNode(aigpp::Node* n)
{
    assert(node() == nullptr || n == nullptr || (node()->manager() == n->manager()));

    if (node() != n) {
        if (!isConstant()) {
            node()->manager()->deref(node());
        }

        Edge::setNode(n);

        if (!isConstant()) {
            node()->manager()->ref(node());
        }
    }
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::existentialQ(const aigpp::InternalEdgeRef& var, bool internalCall) const
{
    assert(!var.isConstant());
    assert(var.node()->isVar());
    if (isConstant()) {
        return *this;
    }

    StatManager::instance().incQuantification();

    Manager* m = node()->manager();

    if (m->verbosity() >= 1) {
        std::cout << "quantifying variable " << var << std::endl;
    }

    if (!internalCall) {
        if (!(node()->varInSupport(var.node()->varIndex()))) {
            return *this;
        }
    }

    bool saveSimCreation = m->simCreation();
    if (!m->automaticFRAIGing()) {
        m->toggleSimCreation(false);
    }

    InternalEdgeRef result;

    InternalEdgeRef cofPos = cofactor(var, false);
    InternalEdgeRef cofNeg = cofactor(!var, false);
    result                 = cofPos + cofNeg;

    m->toggleSimCreation(saveSimCreation);

    if (!internalCall) {
        /* QUIET version */
#if 0
        std::cout << "size after: " << result.nodeCount() << std::endl;
#endif

        if (m->_bddSweepingDuringQuantification) {
            m->BDDSweepInternal(result);
        }

        /* QUIET version */
#if 0
        std::cout << "size after (BDD): " << result.nodeCount() << std::endl;
#endif

        if (!(result.isConstant())
            && m->settings().getRedundancyRemoval() == Settings::REDUNDANCYREMOVAL_QUANTIFICATION) {
            result.node()->findRedundantVars(true);
        }
    }

    return result;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::universalQ(const aigpp::InternalEdgeRef& var) const
{
    return !((!*this).existentialQ(var));
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::existentialQ(const std::vector<aigpp::InternalEdgeRef>& vars) const
{
    if (isConstant()) {
        return *this;
    }

    /* QUIET version */
#if 0
    std::cout << "quantified node: " << node()->refCount() << std::endl;
#endif

    InternalEdgeRef result = *this;

    Manager*  mgr       = node()->manager();
    const int verbosity = mgr->verbosity();

    if (mgr->settings().getQuantifierScheduling() == Settings::QUANTIFIERSCHEDULING_NONE) {
        for (auto var = vars.cbegin(); var != vars.cend(); ++var) {
            if (verbosity >= 1) {
                std::cout << "quantification: remaining variables " << vars.size() - (var - vars.cbegin()) << std::endl;
            }

            if (result.isConstant()) {
                break;
            }
            if (!(result.node()->varInSupport(var->node()->varIndex()))) {
                continue;
            }

            // result = result.existentialQ( *var, /*internalCall=*/true );
            mgr->existentialQInplace(result, *var, /*internalCall=*/true);

            /* QUIET version */
#if 0
            std::cout << "size after: " << result.nodeCount() << std::endl;
#endif
            if (mgr->_bddSweepingDuringQuantification) {
                mgr->BDDSweepInternal(result);
            }

            /* QUIET version */
#if 0
            std::cout << "size after (BDD): " << result.nodeCount() << std::endl;
#endif

            if (!(result.isConstant())
                && mgr->settings().getRedundancyRemoval() == Settings::REDUNDANCYREMOVAL_QUANTIFICATION) {
                result.node()->findRedundantVars(true);
            }
        }
    } else {
        std::set<int> v;
        for (const auto& var : vars) {
            v.insert(var.node()->varIndex());
        }

        while (!(v.empty())) {
            if (result.isConstant()) {
                break;
            }

            if (verbosity >= 1) {
                std::cout << "quantification: remaining variables " << v.size() << std::endl;
            }

            /* delete variables from the schedule that are not in the support */
            VarSupport sup = result.node()->support();

            for (auto p = v.begin(); p != v.end(); /**/) {
                if (!(sup.hasVar(*p))) {
                    v.erase(p++);
                } else {
                    ++p;
                }
            }

            /* quit if no variables are left */
            if (v.empty()) {
                break;
            }

            /* calculate heuristic values for all remaing variables */
            std::map<int, std::size_t> values;
            if (mgr->settings().getQuantifierScheduling() == Settings::QUANTIFIERSCHEDULING_REMAINING_NODES) {
                for (int p : v) {
                    values[p] = result.node()->remainingNodesAfterQuantification(p);
                }
            } else if (mgr->settings().getQuantifierScheduling() == Settings::QUANTIFIERSCHEDULING_MIN_DEPTH) {
                values = result.node()->minVarDepth();
            } else {
                NEVER_GET_HERE;
            }

            /* search for "best" variable according to heuristic value */
            int         bestVar   = -1;
            std::size_t bestValue = 0;

            for (int p : v) {
                assert(values.find(p) != values.end());
                std::size_t value = values[p];

                if (bestVar == -1 || value < bestValue) {
                    bestVar   = p;
                    bestValue = value;
                }
            }

            v.erase(bestVar);

            /* perform quantification */
            result = result.existentialQ(result.node()->manager()->variableInternal(bestVar),
                                         /*internalCall=*/true);

#if 0
            std::cout << "size after: " << result.nodeCount() << std::endl;
#endif

            if (mgr->_bddSweepingDuringQuantification) {
                mgr->BDDSweepInternal(result);
            }
#if 0
            std::cout << "size after (BDD): " << result.nodeCount() << std::endl;
#endif
            if (!(result.isConstant())
                && mgr->settings().getRedundancyRemoval() == Settings::REDUNDANCYREMOVAL_QUANTIFICATION) {
                result.node()->findRedundantVars(true);
            }
        }
    }

    return result;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::universalQ(const std::vector<aigpp::InternalEdgeRef>& vars) const
{
    return !((!*this).existentialQ(vars));
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::swapVariables(const std::vector<aigpp::InternalEdgeRef>& vars1,
                                      const std::vector<aigpp::InternalEdgeRef>& vars2) const
{
    assert(vars1.size() == vars2.size());

#ifndef NDEBUG
    for (const auto& var : vars1) {
        assert(!(var.isConstant()));
        assert(var.node()->isVar());
    }
    for (const auto& var : vars2) {
        assert(!(var.isConstant()));
        assert(var.node()->isVar());
    }
#endif

    if (isConstant()) {
        return *this;
    }

    Manager* mgr = node()->manager();

    /* build cache for given variables */
    auto var2 = vars2.begin();
    for (auto var1 = vars1.begin(); var1 != vars1.end(); ++var1, ++var2) {
        assert(!(var1->node()->isCacheValid()));
        mgr->_cache.insert(var1->node(), var2->notIf(var1->isInverted()));

        assert(!(var2->node()->isCacheValid()));
        mgr->_cache.insert(var2->node(), var1->notIf(var2->isInverted()));
    }

    /* build cache for remaining variables */
    for (Node* var : mgr->_variables) {
        if (!(var->isCacheValid())) {
            mgr->_cache.insert(var, InternalEdgeRef(var));
        }
    }

    /* perform substitution */
    mgr->buildRecursiveUsingCache(node());

    /* get result */
    InternalEdgeRef result = mgr->_cache.lookup(node());
    result.toggleInverted(isInverted());

    /* clear cache */
    mgr->clearNodeCaches();

    return result;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::cofactor(const aigpp::InternalEdgeRef& var, bool checkSupport) const
{
    assert(!var.isConstant());
    assert(var.node()->isVar());

    StatManager::instance().incCofactor();

    if (isConstant()) {
        return *this;
    }

    Manager* mgr = node()->manager();
    if (mgr->verbosity() >= 1) {
        std::cout << "cofactor " << var << std::endl;
        std::cout << "variables in cone: " << node()->support().vars() << std::endl;
        mgr->printStats();
    }

    if (checkSupport) {
        if (!(node()->varInSupport(var.node()->varIndex()))) {
            return *this;
        }
    }

    /* build cache for all variables */
    for (Node* v : mgr->_variables) {
        if (v == var.node()) {
            mgr->_cache.insert(v, const1.notIf(var.isInverted()));
            if (var.isInverted()) {
                v->setFlag<Node::FLAG_CACHECONST0>();
            } else {
                v->setFlag<Node::FLAG_CACHECONST1>();
            }
        } else {
            mgr->_cache.insert(v, InternalEdgeRef(v));
        }
    }

    /* perform cofactoring */
    mgr->_inCofactor  = true;
    mgr->_toBeRemoved = var.node()->varIndex();
    mgr->buildCofactorRecursiveUsingCache(node(), var);
    mgr->_toBeRemoved = -1;
    mgr->_inCofactor  = false;

    /* get result */
    InternalEdgeRef result = mgr->_cache.lookup(node());

    result.toggleInverted(isInverted());

    /* clear cache */
    mgr->clearNodeCaches();

    /* perform bdd sweeping on the result */
    if (mgr->settings().getBDDSweeping() == Settings::BDDSWEEPING_COFACTOR) {
        mgr->BDDSweepInternal(result);
    }

    return result;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::cofactor(const std::vector<aigpp::InternalEdgeRef>& vars) const
{
#ifndef NDEBUG
    for (const auto& var : vars) {
        assert(var.isConstant());
        assert(var.node()->isVar());
    }
#endif

    StatManager::instance().incCofactor();

    if (isConstant()) {
        return *this;
    }

    Manager* mgr = node()->manager();

    /* build cache for given variables */
    for (const auto& var : vars) {
        assert(!(var.node()->isCacheValid()));
        mgr->_cache.insert(var.node(), const1.notIf(var.isInverted()));
    }

    /* build cache for remaining variables */
    for (Node* v : mgr->_variables) {
        if (!(v->isCacheValid())) {
            mgr->_cache.insert(v, InternalEdgeRef(v));
        }
    }

    /* perform cofactoring */
    mgr->_inCofactor = true;
    mgr->buildRecursiveUsingCache(node());
    mgr->_inCofactor = false;

    /* get result */
    InternalEdgeRef result = mgr->_cache.lookup(node());
    result.toggleInverted(isInverted());

    /* clear cache */
    mgr->clearNodeCaches();

    /* perform bdd sweeping on the result */
    if (mgr->settings().getBDDSweeping() == Settings::BDDSWEEPING_COFACTOR) {
        mgr->BDDSweepInternal(result);
    }

    return result;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::compose(const std::vector<InternalEdgeRef>& vars,
                                const std::vector<InternalEdgeRef>& replacements) const
{
    assert(vars.size() == replacements.size());

#ifndef NDEBUG
    for (const auto& var : vars) {
        assert(!(var.isConstant()));
        assert(var.node()->isVar());
    }
#endif

    StatManager::instance().incCompose();

    if (isConstant()) {
        return *this;
    }

    Manager* mgr = node()->manager();

    /* build cache for given variables */
    auto repl = replacements.begin();
    for (auto var = vars.begin(); var != vars.end(); ++var, ++repl) {
        assert(!(var->node()->isCacheValid()));
        mgr->_cache.insert(var->node(), repl->notIf(var->isInverted()));
    }

    /* build cache for remaining variables */
    for (Node* v : mgr->_variables) {
        if (!(v->isCacheValid())) {
            mgr->_cache.insert(v, InternalEdgeRef(v));
        }
    }

    /* perform substitution */
    mgr->buildRecursiveUsingCache(node());

    /* get result */
    InternalEdgeRef result = mgr->_cache.lookup(node());
    result.toggleInverted(isInverted());

    /* clear cache */
    mgr->clearNodeCaches();

    return result;
}

std::vector<aigpp::InternalEdgeRef>
aigpp::InternalEdgeRef::cofactors(const std::vector<aigpp::InternalEdgeRef>& roots, const aigpp::InternalEdgeRef& var)
{
    assert(!var.isConstant());
    assert(var.node()->isVar());

    Manager* mgr = var.node()->manager();

    /* build cache for all variables */
    for (Node* v : mgr->_variables) {
        if (v == var.node()) {
            mgr->_cache.insert(v, const1.notIf(var.isInverted()));
        } else {
            mgr->_cache.insert(v, InternalEdgeRef(v));
        }
    }

    /* perform cofactoring */
    mgr->_inCofactor = true;
    for (const auto& r : roots) {
        if (r.isConstant()) {
            continue;
        }
        mgr->buildCofactorRecursiveUsingCache(r.node(), var);
    }
    mgr->_inCofactor = false;

    /* get result */
    std::vector<InternalEdgeRef> results;
    results.reserve(roots.size());

    for (const auto& r : roots) {
        if (r.isConstant()) {
            results.push_back(r);
        } else {
            InternalEdgeRef result = mgr->_cache.lookup(r.node());
            result.toggleInverted(r.isInverted());
            results.push_back(result);
        }
    }

    /* clear cache */
    mgr->clearNodeCaches();

    return results;
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::composeZ(const std::vector<InternalEdgeRef>& vars,
                                 const std::vector<InternalEdgeRef>& replacements,
                                 const aigpp::InternalEdgeRef&       Z) const
{
    assert(vars.size() == replacements.size());
    assert(!(Z.isConstant()) && !(Z.isInverted()) && (Z.node()->isVar()));

#ifndef NDEBUG
    for (const auto& var : vars) {
        assert(!(var.isConstant()));
        assert(var.node()->isVar());
    }
#endif

    StatManager::instance().incCompose();

    if (isConstant()) {
        return *this;
    }

    Manager* mgr = node()->manager();

    if (node()->isVar()) {
        for (std::size_t i = 0; i < vars.size(); ++i) {
            if (node()->varIndex() == vars[i].node()->varIndex()) {
                if (isInverted() == vars[i].isInverted()) {
                    return replacements[i];
                } else {
                    std::vector<InternalEdgeRef> vars2, repl2;
                    vars2.push_back(Z);
                    repl2.push_back(!Z);

                    return !(replacements[i].compose(vars2, repl2));
                }
            }
        }

        return *this;
    }

    /* build cache for given variables */
    auto repl = replacements.begin();
    for (auto var = vars.begin(); var != vars.end(); ++var, ++repl) {
        assert(!(var->node()->isZCacheValid()));

        mgr->_cacheZ.insert(var->node(),
                            std::pair<InternalEdgeRef, InternalEdgeRef>(repl->cofactor(!Z).notIf(var->isInverted()),
                                                                        repl->cofactor(Z).notIf(var->isInverted())));
    }

    /* build cache for remaining variables */
    for (Node* v : mgr->_variables) {
        if (!(v->isZCacheValid())) {
            mgr->_cacheZ.insert(v, std::pair<InternalEdgeRef, InternalEdgeRef>(InternalEdgeRef(v), InternalEdgeRef(v)));
        }
    }

    /* perform substitution */
    mgr->buildRecursiveZ(node());

    /* get result */
    std::pair<InternalEdgeRef, InternalEdgeRef> resultPair = mgr->_cacheZ.lookup(node());

    if (isInverted()) {
        resultPair.first.toggleInverted();
        resultPair.second.toggleInverted();
        std::swap(resultPair.first, resultPair.second);
    }
    InternalEdgeRef result = resultPair.first + Z * resultPair.second;

    /* clear cache */
    mgr->clearNodeCachesZ();

    return result;
}

aigpp::InternalEdgeRef aigpp::InternalEdgeRef::operator&(const aigpp::InternalEdgeRef& n2) const
{
    StatManager::instance().incAnd();

    if (isConstant()) {
        StatManager::instance().incAndConstant();

        /* 1 & x = x */
        if (isInverted()) {
            return n2;
        }
        /* 0 & x = 0 */
        else {
            return const0;
        }
    } else if (n2.isConstant()) {
        StatManager::instance().incAndConstant();

        /* x & 1 = x */
        if (n2.isInverted()) {
            return *this;
        }
        /* x & 0 = 0 */
        else {
            return const0;
        }
    }
    /* x & y */
    else {
        return node()->manager()->And(*this, n2);
    }
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::multiAnd(std::vector<aigpp::InternalEdgeRef> nodes)
{
    if (nodes.empty()) {
        return const1;
    }

    std::size_t size = nodes.size();

    while (size > 1) {
        std::size_t half = size / 2;
        if (size % 2 == 1) {
            for (std::size_t i = 0; i < half; ++i) {
                nodes[i] = nodes[2 * i] * nodes[2 * i + 1];
            }
            nodes[half] = nodes[size - 1];
            for (std::size_t i = half + 1; i < size; ++i) {
                nodes[i].clear();
            }
            size = half + 1;
        } else {
            for (std::size_t i = 0; i < half; ++i) {
                nodes[i] = nodes[2 * i] * nodes[2 * i + 1];
            }
            for (std::size_t i = half; i < size; ++i) {
                nodes[i].clear();
            }
            size = half;
        }
    }

    return nodes[0];
}

aigpp::InternalEdgeRef
aigpp::InternalEdgeRef::multiOr(std::vector<aigpp::InternalEdgeRef> nodes)
{
    if (nodes.empty()) {
        return const0;
    }

    std::size_t size = nodes.size();

    while (size > 1) {
        const std::size_t half = size / 2;
        if (size % 2 == 1) {
            for (std::size_t i = 0; i < half; ++i) {
                nodes[i] = nodes[2 * i] + nodes[2 * i + 1];
            }
            nodes[half] = nodes[size - 1];
            for (std::size_t i = half + 1; i < size; ++i) {
                nodes[i].clear();
            }
            size = half + 1;
        } else {
            for (std::size_t i = 0; i < half; ++i) {
                nodes[i] = nodes[2 * i] + nodes[2 * i + 1];
            }
            for (std::size_t i = half; i < size; ++i) {
                nodes[i].clear();
            }
            size = half;
        }
    }

    return nodes[0];
}
