/**************************************************************
 *
 *       AIGPP Package // Manager_BDDSweeping.cc
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

/* std */
#include <algorithm>
#include <cmath>
#include <stack>
#include <vector>

/* LRABS utilities */
#include <lrabsutil/Resources.hh>

#include <cuddInt.h>

/* local */
#include "CuddTools.hh"
#include "StatManager.hh"

#define USE_BDD_TIMEOUT

std::size_t
estimateAIGSize(DdManager* man, DdNode* bdd)
{
    if (Cudd_IsConstant(bdd)) {
        return 0;
    }

    DdNode* bdd1 = Cudd_ReadOne(man);
    DdNode* bdd0 = Cudd_Not(bdd1);

    std::set<DdNode*>   processed;
    std::stack<DdNode*> pending;

    std::size_t aigSize = 0;
    pending.push(Cudd_Regular(bdd));

    while (!pending.empty()) {
        if (processed.find(pending.top()) != processed.end()) {
            pending.pop();
        } else {
            DdNode* b = pending.top();
            pending.pop();
            processed.insert(b);

            /* special bdd node handling */

            /* (x->t)*(!x->0) ==  t*x */
            if (Cudd_E(b) == bdd0) {
                /* subcase (x->1)*(!x->0) == x */
                if (Cudd_IsConstant(Cudd_T(b))) {
                    /* t can not be complemented */
                    assert(Cudd_T(b) == bdd1);

                    /* nothing to do since we count variables separately */
                } else {
                    /* case (x->0)*(!x->0) can not occur since t can not be complemented
                     */
                    assert(!Cudd_IsConstant(Cudd_T(b)));

                    aigSize += 1;

                    pending.push(Cudd_T(b));
                }
            }
            /* (x->t)*(!x->1) ==  t+!x */
            else if (Cudd_E(b) == bdd1) {
                /* case (x->1)*(!x->1) = 1 can not occur (IMHO) */
                /* case (x->0)*(!x->1) can not occur since t can not be complemented */
                assert(!Cudd_IsConstant(Cudd_T(b)));

                aigSize += 1;

                pending.push(Cudd_T(b));
            }
            /* (x->1)*(!x->e) == e+x */
            else if (Cudd_T(b) == bdd1) {
                /* cases e==const handeled above */
                assert(!Cudd_IsConstant(Cudd_E(b)));

                aigSize += 1;

                pending.push(Cudd_Regular(Cudd_E(b)));
            }
            /* normal case */
            else {
                assert(!Cudd_IsConstant(Cudd_T(b)));
                assert(!Cudd_IsConstant(Cudd_E(b)));

                aigSize += 3;

                pending.push(Cudd_T(b));
                pending.push(Cudd_Regular(Cudd_E(b)));
            }
        }
    }

    aigSize += static_cast<std::size_t>(Cudd_SupportSize(man, bdd));

    return aigSize;
}

DdNode*
aigpp::Manager::buildBDDfromAIG(aigpp::Node* aig)
{
#ifdef ENABLE_BDD
    /* handle constant false node */
    if (aig == nullptr) {
        DdNode* bdd = Cudd_Not(Cudd_ReadOne(_bddManager));
        Cudd_Ref(bdd);
        return bdd;
    }

    _bddReachedNodeLimit = false;
    _bddReachedTimeLimit = false;

    double timeBDD = lrabs::cpuTime();

    lockFlags<1>();

    std::stack<Node*>& pending = _universalNodeStack;
    assert(pending.empty());
    std::vector<Node*>& processed = _universalNodeVector;
    assert(processed.empty());

    Node* n;

    /* count references/occurences of each node in the cone */
    using IntBDD  = std::pair<int, DdNode*>;
    using A2CBMap = std::map<Node*, IntBDD>;
    A2CBMap amap;

    pending.push(aig);

    while (!pending.empty()) {
        n = pending.top();

        assert(!(n->isVar()));

        /* already processed */
        if (n->flag<1>()) {
            pending.pop();
        }
        /* parent 1 not processed yet */
        else if (!(n->parent1().node()->isVar()) && !(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
        }
        /* parent 2 not processed yet */
        else if (!(n->parent2().node()->isVar()) && !(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
        }
        /* process node */
        else {
            /* increase parents' counts */
            if (!(n->parent1().node()->isVar())) {
                ++amap[n->parent1().node()].first;
            }
            if (!(n->parent2().node()->isVar())) {
                ++amap[n->parent2().node()].first;
            }

            amap[n] = IntBDD(0, 0);

            processed.push_back(n);
            n->setFlag<1>();
            pending.pop();
        }
    }

    for (Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    /* build bdd from aig, store intermediate bdd nodes in amap */

    bool reachedLimit = false;

#    ifdef USE_BDD_TIMEOUT
    /* set bdd timeout */
    const double maxBDDTime       = 20.0;
    double       bddTimeout       = lrabs::cpuTime() + maxBDDTime;
    bool         reachedTimeLimit = false;
    bool         reachedNodeLimit = false;
#    endif

    pending.push(aig);
    while (!pending.empty()) {
        n = pending.top();

        assert(!(n->isVar()));

        /* skip processed nodes */
        if (n->flag<1>()) {
            pending.pop();
        }
        /* process parent1 */
        else if (!(n->parent1().node()->isVar()) && !(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
        }
        /* process parent2 */
        else if (!(n->parent2().node()->isVar()) && !(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
        }
        /* process node */
        else {
            /* get parent1's bdd */
            DdNode* p1;
            if (n->parent1().node()->isVar()) {
                p1 = Cudd_bddIthVar(_bddManager, n->parent1().node()->varIndex());

                if (p1 == nullptr) {
                    reachedLimit     = true;
                    reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                    double t = bddTimeout - lrabs::cpuTime();
                    if (t < 0) {
                        reachedTimeLimit = true;
                    }
#    endif
                    break;
                }

                Cudd_Ref(p1);
            } else {
                auto cb = amap.find(n->parent1().node());
                assert(cb != amap.end());

                p1 = cb->second.second;
                assert(p1 != 0);

                --(cb->second.first);
                if (cb->second.first == 0) {
                    amap.erase(cb);
                } else {
                    Cudd_Ref(p1);
                }
            }
            if (n->parent1().isInverted()) {
                p1 = Cudd_Not(p1);
            }

            /* get parent2's bdd */
            DdNode* p2;
            if (n->parent2().node()->isVar()) {
                p2 = Cudd_bddIthVar(_bddManager, n->parent2().node()->varIndex());

                if (p2 == nullptr) {
                    Cudd_RecursiveDeref(_bddManager, p1);
                    reachedLimit     = true;
                    reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                    double t = bddTimeout - lrabs::cpuTime();
                    if (t < 0) {
                        reachedTimeLimit = true;
                    }
#    endif
                    break;
                }

                Cudd_Ref(p2);
            } else {
                auto cb = amap.find(n->parent2().node());
                assert(cb != amap.end());

                p2 = cb->second.second;
                assert(p2 != 0);

                --(cb->second.first);
                if (cb->second.first == 0) {
                    amap.erase(cb);
                } else {
                    Cudd_Ref(p2);
                }
            }
            if (n->parent2().isInverted()) {
                p2 = Cudd_Not(p2);
            }

            DdNode* b = Cudd_bddAnd(_bddManager, p1, p2);

            /* no bdd node was created -> node limit reached! */
            if (b == nullptr) {
                Cudd_RecursiveDeref(_bddManager, p1);
                Cudd_RecursiveDeref(_bddManager, p2);
                reachedLimit     = true;
                reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                double t = bddTimeout - lrabs::cpuTime();
                if (t < 0) {
                    reachedTimeLimit = true;
                }
#    endif
                break;
            }

            Cudd_Ref(b);
            Cudd_RecursiveDeref(_bddManager, p1);
            Cudd_RecursiveDeref(_bddManager, p2);

#    ifdef USE_BDD_TIMEOUT
            double t = bddTimeout - lrabs::cpuTime();
            if (t < 0) {
                Cudd_RecursiveDeref(_bddManager, b);
                reachedTimeLimit = true;
                reachedLimit     = true;
                break;
            }
#    endif

            if (n->isNAND()) {
                b = Cudd_Not(b);
            }

            auto cb = amap.find(n);
            assert(cb != amap.end());
            cb->second.second = b;

            processed.push_back(n);
            n->setFlag<1>();
        }
    }

    /* clear stack (stack may be non-empty if the bdd limit was reached) */
    assert(reachedLimit || pending.empty());
    while (!pending.empty()) {
        pending.pop();
    }

    /* get result */
    DdNode* reducedBDD = nullptr;
    if (!reachedLimit) {
        reducedBDD = amap[aig].second;
        assert(reducedBDD != 0);
        Cudd_Ref(reducedBDD);
    } else {
        _bddReachedNodeLimit = reachedNodeLimit;
        _bddReachedTimeLimit = reachedTimeLimit;
    }

    /* clear amap */
    for (const auto& p : amap) {
        if (p.second.second != nullptr) {
            Cudd_RecursiveDeref(_bddManager, p.second.second);
        }
    }
    amap.clear();

    for (Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    unlockFlags<1>();

    StatManager::instance().incBDDTime(lrabs::cpuTime() - timeBDD);

    return reducedBDD;
#else
    abort();
    return 0;
#endif
}

DdNode*
aigpp::Manager::buildBDDfromAIGlimited(aigpp::Node* aig, double growthLimit)
{
#ifdef ENABLE_BDD
    /* handle constant false node */
    if (aig == nullptr) {
        DdNode* bdd = Cudd_Not(Cudd_ReadOne(_bddManager));
        Cudd_Ref(bdd);
        return bdd;
    }

    _bddReachedGrowthLimit = false;
    _bddReachedNodeLimit   = false;
    _bddReachedTimeLimit   = false;

    const double timeBDD = lrabs::cpuTime();

    lockFlags<1>();

    std::stack<Node*>& pending = _universalNodeStack;
    assert(pending.empty());
    std::vector<Node*>& processed = _universalNodeVector;
    assert(processed.empty());

    Node* n = nullptr;

    /* count references/occurences of each node in the cone */
    using IntBDD  = std::pair<int, DdNode*>;
    using A2CBMap = std::map<Node*, IntBDD>;
    A2CBMap amap;

    pending.push(aig);

    while (!pending.empty()) {
        n = pending.top();

        assert(!(n->isVar()));

        /* already processed */
        if (n->flag<1>()) {
            pending.pop();
        }
        /* parent 1 not processed yet */
        else if (!(n->parent1().node()->isVar()) && !(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
        }
        /* parent 2 not processed yet */
        else if (!(n->parent2().node()->isVar()) && !(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
        }
        /* process node */
        else {
            /* increase parents' counts */
            if (!(n->parent1().node()->isVar())) {
                ++amap[n->parent1().node()].first;
            }
            if (!(n->parent2().node()->isVar())) {
                ++amap[n->parent2().node()].first;
            }

            amap[n] = IntBDD(0, nullptr);

            processed.push_back(n);
            n->setFlag<1>();
            pending.pop();
        }
    }

    for (Node* n : processed) {
        n->unsetFlag<1>();
    }
    processed.clear();

    /* build bdd from aig, store intermediate bdd nodes in amap */

    bool reachedLimit = false;

#    ifdef USE_BDD_TIMEOUT
    /* set bdd timeout */
    const double maxBDDTime         = 20.0;
    double       bddTimeout         = lrabs::cpuTime() + maxBDDTime;
    bool         reachedGrowthLimit = false;
    bool         reachedTimeLimit   = false;
    bool         reachedNodeLimit   = false;
#    endif

    std::size_t processedNodes = 0;

    pending.push(aig);
    while (!pending.empty()) {
        n = pending.top();

        assert(!(n->isVar()));

        /* skip processed nodes */
        if (n->flag<1>()) {
            pending.pop();
        }
        /* process parent1 */
        else if (!(n->parent1().node()->isVar()) && !(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
        }
        /* process parent2 */
        else if (!(n->parent2().node()->isVar()) && !(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
        }
        /* process node */
        else {
            /* get parent1's bdd */
            DdNode* p1 = nullptr;
            if (n->parent1().node()->isVar()) {
                p1 = Cudd_bddIthVar(_bddManager, n->parent1().node()->varIndex());

                if (p1 == nullptr) {
                    reachedLimit     = true;
                    reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                    const double t = bddTimeout - lrabs::cpuTime();
                    if (t < 0) {
                        reachedTimeLimit = true;
                    }
#    endif
                    break;
                }

                Cudd_Ref(p1);
            } else {
                auto cb = amap.find(n->parent1().node());
                assert(cb != amap.end());

                p1 = cb->second.second;
                assert(p1 != nullptr);

                --(cb->second.first);
                if (cb->second.first == 0) {
                    amap.erase(cb);
                } else {
                    Cudd_Ref(p1);
                }
            }
            if (n->parent1().isInverted()) {
                p1 = Cudd_Not(p1);
            }

            /* get parent2's bdd */
            DdNode* p2 = nullptr;
            if (n->parent2().node()->isVar()) {
                p2 = Cudd_bddIthVar(_bddManager, n->parent2().node()->varIndex());

                if (p2 == nullptr) {
                    Cudd_RecursiveDeref(_bddManager, p1);
                    reachedLimit     = true;
                    reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                    const double t = bddTimeout - lrabs::cpuTime();
                    if (t < 0) {
                        reachedTimeLimit = true;
                    }
#    endif
                    break;
                }

                Cudd_Ref(p2);
            } else {
                auto cb = amap.find(n->parent2().node());
                assert(cb != amap.end());

                p2 = cb->second.second;
                assert(p2 != nullptr);

                --(cb->second.first);
                if (cb->second.first == 0) {
                    amap.erase(cb);
                } else {
                    Cudd_Ref(p2);
                }
            }
            if (n->parent2().isInverted()) {
                p2 = Cudd_Not(p2);
            }

            DdNode* b = Cudd_bddAnd(_bddManager, p1, p2);

            /* no bdd node was created -> node limit reached! */
            if (b == nullptr) {
                Cudd_RecursiveDeref(_bddManager, p1);
                Cudd_RecursiveDeref(_bddManager, p2);
                reachedLimit     = true;
                reachedNodeLimit = true;

#    ifdef USE_BDD_TIMEOUT
                const double t = bddTimeout - lrabs::cpuTime();
                if (t < 0) {
                    reachedTimeLimit = true;
                }
#    endif
                break;
            }

            Cudd_Ref(b);
            Cudd_RecursiveDeref(_bddManager, p1);
            Cudd_RecursiveDeref(_bddManager, p2);

#    ifdef USE_BDD_TIMEOUT
            const double t = bddTimeout - lrabs::cpuTime();
            if (t < 0) {
                Cudd_RecursiveDeref(_bddManager, b);
                reachedTimeLimit = true;
                reachedLimit     = true;
                break;
            }
#    endif

            ++processedNodes;
            const int bddSize = static_cast<int>(Cudd_ReadNodeCount(_bddManager)) - Cudd_ReadSize(_bddManager);
            if (static_cast<double>(bddSize) > static_cast<double>(processedNodes) * growthLimit) {
                // std::cout << "bdd " << bddSize << " aig " << processedNodes <<
                // std::endl;
                Cudd_RecursiveDeref(_bddManager, b);
                reachedGrowthLimit = true;
                reachedLimit       = true;
                break;
            }

            if (n->isNAND()) {
                b = Cudd_Not(b);
            }

            auto cb = amap.find(n);
            assert(cb != amap.end());
            cb->second.second = b;

            processed.push_back(n);
            n->setFlag<1>();
        }
    }

    /* clear stack (stack may be non-empty if the bdd limit was reached) */
    assert(reachedLimit || pending.empty());
    while (!pending.empty()) {
        pending.pop();
    }

    /* get result */
    DdNode* reducedBDD = nullptr;
    if (!reachedLimit) {
        reducedBDD = amap[aig].second;
        assert(reducedBDD != nullptr);
        Cudd_Ref(reducedBDD);
    } else {
        _bddReachedGrowthLimit = reachedGrowthLimit;
        _bddReachedNodeLimit   = reachedNodeLimit;
        _bddReachedTimeLimit   = reachedTimeLimit;
    }

    /* clear amap */
    for (const auto& p : amap) {
        if (p.second.second != nullptr) {
            Cudd_RecursiveDeref(_bddManager, p.second.second);
        }
    }
    amap.clear();

    for (Node* n : processed) {
        n->unsetFlag<1>();
    }
    processed.clear();

    unlockFlags<1>();

    StatManager::instance().incBDDTime(lrabs::cpuTime() - timeBDD);

    return reducedBDD;
#else
    abort();
    return 0;
#endif
}

void
aigpp::Manager::buildAIGfromBDD(DdNode* bdd, aigpp::Node* originalAIG)
{
#ifdef ENABLE_BDD
    DdNode* bdd1 = Cudd_ReadOne(_bddManager);
    DdNode* bdd0 = Cudd_Not(bdd1);

    std::stack<DdNode*> pending;

    /* count references of bdd nodes */
    std::map<DdNode*, int> BddCount;

    pending.push(Cudd_Regular(bdd));

    while (!pending.empty()) {
        assert(!Cudd_IsComplement(pending.top()));

        /* already processed */
        if (BddCount.find(pending.top()) != BddCount.end()) {
            pending.pop();
            continue;
        }

        /* constant bdd */
        if (Cudd_IsConstant(pending.top())) {
            BddCount[pending.top()] = 0;
            pending.pop();
            continue;
        }

        assert(!Cudd_IsComplement(Cudd_T(pending.top())));

        /* variable encoding bdd "x*1 + !x*0" */
        if (Cudd_T(pending.top()) == bdd1 && Cudd_E(pending.top()) == bdd0) {
            BddCount[pending.top()] = 0;
            pending.pop();
            continue;
        }

        auto t = BddCount.find(Cudd_T(pending.top()));
        if (t == BddCount.end()) {
            pending.push(Cudd_T(pending.top()));
            continue;
        }

        auto e = BddCount.find(Cudd_Regular(Cudd_E(pending.top())));
        if (e == BddCount.end()) {
            pending.push(Cudd_Regular(Cudd_E(pending.top())));
            continue;
        }

        BddCount[pending.top()] = 0;
        if (!Cudd_IsConstant(Cudd_T(pending.top()))) {
            ++(t->second);
        }
        if (!Cudd_IsConstant(Cudd_Regular(Cudd_E(pending.top())))) {
            ++(e->second);
        }

        pending.pop();
    }

    /*
     * build the aig using the BDD structure
     */

    double timeAIG = lrabs::cpuTime();

    using BAMap = std::map<DdNode*, InternalEdgeRef>;
    BAMap bdd2aig;

    pending.push(Cudd_Regular(bdd));

    while (!pending.empty()) {
        assert(!Cudd_IsComplement(pending.top()));

        /* already processed */
        if (bdd2aig.find(pending.top()) != bdd2aig.end()) {
            pending.pop();
            continue;
        }

        /* constant bdd */
        if (Cudd_IsConstant(pending.top())) {
            assert(pending.top() == bdd1);
            bdd2aig[pending.top()] = aigpp::const1;
            pending.pop();
            continue;
        }

        assert(!Cudd_IsComplement(Cudd_T(pending.top())));

        /* variable encoding bdd "x->1 * !x->0 == x" */
        if (Cudd_T(pending.top()) == bdd1 && Cudd_E(pending.top()) == bdd0) {
            bdd2aig[pending.top()] = variableInternal(pending.top()->index);
            pending.pop();
            continue;
        }

        auto t = bdd2aig.find(Cudd_T(pending.top()));
        if (t == bdd2aig.end()) {
            pending.push(Cudd_T(pending.top()));
            continue;
        }

        auto e = bdd2aig.find(Cudd_Regular(Cudd_E(pending.top())));
        if (e == bdd2aig.end()) {
            pending.push(Cudd_Regular(Cudd_E(pending.top())));
            continue;
        }

        assert((pending.size() == 1) == (pending.top() == Cudd_Regular(bdd)));

        /* this is not the top bdd node! */
        if (pending.size() != 1) {
            /* special bdd node handling */

            /* (x->t)*(!x->0) ==  t*x */
            if (Cudd_E(pending.top()) == bdd0) {
                /* case (x->1)*(!x->0) handled above */
                /* case (x->0)*(!x->0) can not occur since t can not be complemented */
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));

                bdd2aig[pending.top()] = t->second * variableInternal(pending.top()->index);
            }
            /* (x->t)*(!x->1) ==  t+!x */
            else if (Cudd_E(pending.top()) == bdd1) {
                /* case (x->1)*(!x->1) = 1 can not occur (IMHO) */
                /* case (x->0)*(!x->1) can not occur since t can not be complemented */
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));
                bdd2aig[pending.top()] = t->second + !variableInternal(pending.top()->index);
            }
            /* (x->1)*(!x->e) == e+x */
            else if (Cudd_T(pending.top()) == bdd1) {
                /* cases e==const handeled above */
                assert(!Cudd_IsConstant(Cudd_E(pending.top())));

                bdd2aig[pending.top()] = e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top())))
                                         + variableInternal(pending.top()->index);
            }
            /* normal case */
            else {
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));
                assert(!Cudd_IsConstant(Cudd_E(pending.top())));

                bdd2aig[pending.top()] = variableInternal(pending.top()->index)
                                             .Mux(t->second, e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top()))));
            }
        }

        /* this is the top bdd node */
        else {
            /* special bdd node handling */
            /* (x->t)*(!x->0) ==  t*x */
            if (Cudd_E(pending.top()) == bdd0) {
                /* case (x->1)*(!x->0) handled above */
                /* case (x->0)*(!x->0) can not occur since t can not be complemented */
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));

                _nextAndSatResult           = originalAIG;
                _nextAndSatResultIsInverted = Cudd_IsComplement(bdd);
                _nextAndSkipSat             = true;
                InternalEdgeRef mux         = t->second * variableInternal(pending.top()->index);
                _nextAndSkipSat             = false;
                bdd2aig[pending.top()]      = mux;

            }
            /* (x->t)*(!x->1) ==  t+!x */
            else if (Cudd_E(pending.top()) == bdd1) {
                /* case (x->1)*(!x->1) = 1 can not occur (IMHO) */
                /* case (x->0)*(!x->1) can not occur since t can not be complemented */
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));

                _nextAndSatResult = originalAIG;
                /* complement because we calculate an "or"! */
                _nextAndSatResultIsInverted = !Cudd_IsComplement(bdd);
                _nextAndSkipSat             = true;
                InternalEdgeRef mux         = t->second + !variableInternal(pending.top()->index);
                _nextAndSkipSat             = false;
                bdd2aig[pending.top()]      = mux;
            }
            /* (x->1)*(!x->e) == e+x */
            else if (Cudd_T(pending.top()) == bdd1) {
                /* cases e==const handeled above */
                assert(!Cudd_IsConstant(Cudd_E(pending.top())));

                _nextAndSatResult = originalAIG;
                /* complement because we calculate an "or"! */
                _nextAndSatResultIsInverted = !Cudd_IsComplement(bdd);
                _nextAndSkipSat             = true;
                InternalEdgeRef mux         = e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top())))
                                      + variableInternal(pending.top()->index);
                _nextAndSkipSat        = false;
                bdd2aig[pending.top()] = mux;
            }
            /* normal case */
            else {
                assert(!Cudd_IsConstant(Cudd_T(pending.top())));
                assert(!Cudd_IsConstant(Cudd_E(pending.top())));

                InternalEdgeRef aigThen = !variableInternal(pending.top()->index) + t->second;

                InternalEdgeRef aigElse = variableInternal(pending.top()->index)
                                          + e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top())));

                _nextAndSatResult           = originalAIG;
                _nextAndSatResultIsInverted = Cudd_IsComplement(bdd);
                _nextAndSkipSat             = true;
                InternalEdgeRef mux         = aigThen * aigElse;
                _nextAndSkipSat             = false;
                bdd2aig[pending.top()]      = mux;
            }
        }

        /* decrease counts of bdd nodes and delete unused nodes */
        if (!Cudd_IsConstant(Cudd_T(pending.top()))) {
            auto c = BddCount.find(Cudd_T(pending.top()));

            --(c->second);
            if (c->second == 0) {
                bdd2aig.erase(t);
                BddCount.erase(c);
            }
        }
        if (!Cudd_IsConstant(Cudd_Regular(Cudd_E(pending.top())))) {
            auto c = BddCount.find(Cudd_Regular(Cudd_E(pending.top())));
            --(c->second);
            if (c->second == 0) {
                bdd2aig.erase(e);
                BddCount.erase(c);
            }
        }

        pending.pop();
    }

    /* clear cache */
    bdd2aig.clear();

    Cudd_RecursiveDeref(_bddManager, bdd);

    StatManager::instance().incBDDAIGTime(lrabs::cpuTime() - timeAIG);
#else
    abort();
#endif
}

aigpp::InternalEdgeRef
aigpp::Manager::createAIG(DdNode* bdd)
{
#ifdef ENABLE_BDD
    DdNode* bdd1 = Cudd_ReadOne(_bddManager);
    DdNode* bdd0 = Cudd_Not(bdd1);

    std::stack<DdNode*> pending;

    /* count references of bdd nodes */
    std::map<DdNode*, int> BddCount;

    pending.push(Cudd_Regular(bdd));

    while (!pending.empty()) {
        assert(!Cudd_IsComplement(pending.top()));

        /* already processed */
        if (BddCount.find(pending.top()) != BddCount.end()) {
            pending.pop();
            continue;
        }

        /* constant bdd */
        if (Cudd_IsConstant(pending.top())) {
            BddCount[pending.top()] = 0;
            pending.pop();
            continue;
        }

        assert(!Cudd_IsComplement(Cudd_T(pending.top())));

        /* variable encoding bdd "x*1 + !x*0" */
        if (Cudd_T(pending.top()) == bdd1 && Cudd_E(pending.top()) == bdd0) {
            BddCount[pending.top()] = 0;
            pending.pop();
            continue;
        }

        auto t = BddCount.find(Cudd_T(pending.top()));
        if (t == BddCount.end()) {
            pending.push(Cudd_T(pending.top()));
            continue;
        }

        auto e = BddCount.find(Cudd_Regular(Cudd_E(pending.top())));
        if (e == BddCount.end()) {
            pending.push(Cudd_Regular(Cudd_E(pending.top())));
            continue;
        }

        BddCount[pending.top()] = 0;
        if (!Cudd_IsConstant(Cudd_T(pending.top()))) {
            ++(t->second);
        }
        if (!Cudd_IsConstant(Cudd_Regular(Cudd_E(pending.top())))) {
            ++(e->second);
        }

        pending.pop();
    }

    /*
     * build the aig using the BDD structure
     */

    double timeAIG = lrabs::cpuTime();

    using BAMap = std::map<DdNode*, InternalEdgeRef>;
    BAMap bdd2aig;

    pending.push(Cudd_Regular(bdd));

    while (!pending.empty()) {
        assert(!Cudd_IsComplement(pending.top()));

        /* already processed */
        if (bdd2aig.find(pending.top()) != bdd2aig.end()) {
            pending.pop();
            continue;
        }

        /* constant bdd */
        if (Cudd_IsConstant(pending.top())) {
            assert(pending.top() == bdd1);
            bdd2aig[pending.top()] = aigpp::const1;
            pending.pop();
            continue;
        }

        assert(!Cudd_IsComplement(Cudd_T(pending.top())));

        /* variable encoding bdd "x->1 * !x->0 == x" */
        if (Cudd_T(pending.top()) == bdd1 && Cudd_E(pending.top()) == bdd0) {
            bdd2aig[pending.top()] = variableInternal(pending.top()->index);
            pending.pop();
            continue;
        }

        auto t = bdd2aig.find(Cudd_T(pending.top()));
        if (t == bdd2aig.end()) {
            pending.push(Cudd_T(pending.top()));
            continue;
        }

        auto e = bdd2aig.find(Cudd_Regular(Cudd_E(pending.top())));
        if (e == bdd2aig.end()) {
            pending.push(Cudd_Regular(Cudd_E(pending.top())));
            continue;
        }

        assert((pending.size() == 1) == (pending.top() == Cudd_Regular(bdd)));

        /* special bdd node handling */

        /* (x->t)*(!x->0) ==  t*x */
        if (Cudd_E(pending.top()) == bdd0) {
            /* case (x->1)*(!x->0) handled above */
            /* case (x->0)*(!x->0) can not occur since t can not be complemented */
            assert(!Cudd_IsConstant(Cudd_T(pending.top())));

            bdd2aig[pending.top()] = t->second * variableInternal(pending.top()->index);
        }
        /* (x->t)*(!x->1) ==  t+!x */
        else if (Cudd_E(pending.top()) == bdd1) {
            /* case (x->1)*(!x->1) = 1 can not occur (IMHO) */
            /* case (x->0)*(!x->1) can not occur since t can not be complemented */
            assert(!Cudd_IsConstant(Cudd_T(pending.top())));
            bdd2aig[pending.top()] = t->second + !variableInternal(pending.top()->index);
        }
        /* (x->1)*(!x->e) == e+x */
        else if (Cudd_T(pending.top()) == bdd1) {
            /* cases e==const handeled above */
            assert(!Cudd_IsConstant(Cudd_E(pending.top())));

            bdd2aig[pending.top()]
                = e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top()))) + variableInternal(pending.top()->index);
        }
        /* normal case */
        else {
            assert(!Cudd_IsConstant(Cudd_T(pending.top())));
            assert(!Cudd_IsConstant(Cudd_E(pending.top())));

            bdd2aig[pending.top()] = variableInternal(pending.top()->index)
                                         .Mux(t->second, e->second.notIf(Cudd_IsComplement(Cudd_E(pending.top()))));
        }

        /* decrease counts of bdd nodes and delete unused nodes */
        if (!Cudd_IsConstant(Cudd_T(pending.top()))) {
            auto c = BddCount.find(Cudd_T(pending.top()));

            --(c->second);
            if (c->second == 0) {
                bdd2aig.erase(t);
                BddCount.erase(c);
            }
        }
        if (!Cudd_IsConstant(Cudd_Regular(Cudd_E(pending.top())))) {
            auto c = BddCount.find(Cudd_Regular(Cudd_E(pending.top())));
            --(c->second);
            if (c->second == 0) {
                bdd2aig.erase(e);
                BddCount.erase(c);
            }
        }

        pending.pop();
    }

    InternalEdgeRef result = bdd2aig[Cudd_Regular(bdd)].notIf(Cudd_IsComplement(bdd));

    /* clear cache */
    bdd2aig.clear();

    StatManager::instance().incBDDAIGTime(lrabs::cpuTime() - timeAIG);

    return result;
#else
    abort();
    return const0;
#endif
}

aigpp::Manager::BDDSweepingResult
aigpp::Manager::BDDSweep(const aigpp::EdgeRef& bigNode, bool force)
{
#ifdef ENABLE_BDD
    return BDDSweepInternal(bigNode.getInternal(), force);
#else
    abort();
    return BDDSweeping_Skipped;
#endif
}

aigpp::Manager::BDDSweepingResult
aigpp::Manager::BDDSweepInternal(const aigpp::InternalEdgeRef& bigNode, bool force)
{
#ifdef ENABLE_BDD

    /* setup bdd */
    cleanBDD();

    if (_bddManager == nullptr) {
        return BDDSweeping_Skipped;
    }

    //    assert( static_cast<std::size_t>( Cudd_ReadNodeCount( _bddManager ) ) ==
    //    variableCount() + 1 );

    StatManager::instance().incBDDCalls();

    /* constants and variables cannot be reduced any further! */
    if (bigNode.isConstant() || bigNode.node()->isVar()) {
        StatManager::instance().incBDDSkipped();
        return BDDSweeping_Skipped;
    }

    std::size_t originalSize = bigNode.node()->nodeCount();
    if (!force) {
        /* skip if input aig is too small */
        if (originalSize < 1000) {
            StatManager::instance().incBDDSkipped();
            return BDDSweeping_Skipped;
        }

        /* skip if input aig is smaller than _nextReductionSize */
        if (originalSize < _nextReductionSize) {
            StatManager::instance().incBDDSkipped();

            _nextReductionSize
                = static_cast<std::size_t>(static_cast<double>(_nextReductionSize) * (1.0 - _BDDSimDecay));
            return BDDSweeping_Skipped;
        }
    }

    /* reset reordering next level */
    Cudd_SetNextReordering(_bddManager, static_cast<int>(originalSize));
    Cudd_AutodynEnable(_bddManager, CUDD_REORDER_SIFT);

    if (settings().getBDDSweepingInitializeOrder()) {
        /* calculate the depths of all variables occuring in bigNode */
        std::map<int, std::size_t> depthMap = bigNode.node()->minVarDepth();
        std::set<std::size_t>      depths;
        for (std::map<int, std::size_t>::const_iterator p = depthMap.begin(); p != depthMap.end(); ++p) {
            depths.insert(p->second);
        }

        auto*       newVarOrder = new int[variableCount()];
        std::size_t index       = 0;
        for (const std::size_t d : depths) {
            for (const auto& p : depthMap) {
                if (p.second == d) {
                    newVarOrder[index] = static_cast<int>(p.first);
                    ++index;
                }
            }
        }
        for (std::size_t v = 0; v < variableCount(); ++v) {
            if (depthMap.find(static_cast<int>(v)) == depthMap.end()) {
                newVarOrder[index] = static_cast<int>(v);
                ++index;
            }
        }

        /* set initial variable ordering */
        Cudd_SetMaxLive(_bddManager, 10000000);
        Cudd_ShuffleHeap(_bddManager, newVarOrder);
        delete[] newVarOrder;
    }

    /* limit number of intermediate bdd nodes */
    double BDDLimit = std::min(100.0 * static_cast<double>(originalSize),
                               10000000.0);  // conversion ok
    Cudd_SetMaxLive(_bddManager, static_cast<int>(BDDLimit));

    /* build bdd */
    DdNode* reducedBDD = buildBDDfromAIG(bigNode.node());

    /*
     * the bdd calculation was cancelled because the specified node limit was
     * reached
     */
    if (reducedBDD == nullptr) {
        cuddGarbageCollect(_bddManager, true);

        if (!force) {
            StatManager::instance().incBDDLimit();

            _nextReductionSize = 3 * std::max(_nextReductionSize, originalSize);

            _BDDSimDecay *= 0.5;
        }

        return BDDSweeping_BDDCouldNotBeConstructed;
    }

    int  bddSize         = Cudd_DagSize(reducedBDD);
    auto requiredBddSize = static_cast<int>((originalSize / 2));

    /*
     * optimize the resulting BDD by a final reordering step
     */
    if (settings().getBDDSweepingFinalReorder() && bddSize < 5 * requiredBddSize) {
        Cudd_ReduceHeap(_bddManager, CUDD_REORDER_LAZY_SIFT, 0);
        bddSize = Cudd_DagSize(reducedBDD);
    }

    /* the size of the resulting bdd is small enough */
    if (bddSize < requiredBddSize) {
        StatManager::instance().incBDDSuccess();

        /* create aig */
        buildAIGfromBDD(reducedBDD, bigNode.node());
        std::size_t reducedSize = bigNode.nodeCount();

        if (!force) {
            _bddreducedratio
                = 0.5 * _bddreducedratio + 0.5 * static_cast<double>(reducedSize) / static_cast<double>(bddSize);
            _nextReductionSize = reducedSize;
            _BDDSimDecay       = 0.01;
        }

        return BDDSweeping_Success;
    } else {
        Cudd_RecursiveDeref(_bddManager, reducedBDD);

        if (!force) {
            StatManager::instance().incBDDFailed();

            _nextReductionSize = static_cast<std::size_t>(1.5 * static_cast<double>(_nextReductionSize));
        }

        return BDDSweeping_BDDTooBig;
    }
#else
    abort();
    return BDDSweeping_Skipped;
#endif
}

void
aigpp::Manager::cleanBDD()
{
    if (_bddManager != nullptr) {
        /* check if there are any remaining bdd nodes */
        if (Cudd_ReadNodeCount(_bddManager) != 1 + Cudd_ReadSize(_bddManager)) {
            Cudd_Quit(_bddManager);
            _bddManager = nullptr;
        } else {
            return;
        }
    }

    if (_bddManager == nullptr) {
        _bddManager = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
        Cudd_AutodynEnable(_bddManager, CUDD_REORDER_SIFT);

        for (std::size_t v = 0; v != variableCount(); ++v) {
            DdNode* bddNode = Cudd_bddNewVar(_bddManager);
            if (bddNode == nullptr) {
                std::cout << "c cannot create bdd variable node -> disabling bdd manager" << std::endl;
                Cudd_Quit(_bddManager);
                return;
            } else {
                Cudd_Ref(bddNode);
            }
        }
    }
}

std::vector<DdNode*>
aigpp::Manager::computeBDD(const std::vector<aigpp::EdgeRef>& roots, DdManager* bddMan,
                           const std::map<int, DdNode*>& varMapping) const
{
    using IntBDD  = std::pair<int, DdNode*>;
    using A2CBMap = std::map<Node*, IntBDD>;

    A2CBMap           mapping;
    Node*             n;
    A2CBMap::iterator p1, p2;

    std::vector<Node*> processed;
    processed.reserve(10000);
    lockFlags<1>();

    for (const auto& p : varMapping) {
        assert(isVarIndexValid(p.first));

        n = variable(p.first).getInternal().node();
        assert(mapping.find(n) == mapping.end());

        mapping[n] = IntBDD(1, p.second);
        Cudd_Ref(p.second);

        processed.push_back(n);
        n->setFlag<1>();
    }

    std::stack<Node*> pending;

    for (const auto& p : roots) {
        if (!(p.isConstant())) {
            pending.push(p.getInternal().node());
        }
    }

    int remainingNodes = 0;

    /* count references/occurences of each node in the cone */
    while (!pending.empty()) {
        n = pending.top();

        /* already processed */
        if (mapping.find(n) != mapping.end()) {
            pending.pop();
        }
        /* parent 1 not processed yet */
        else if ((p1 = mapping.find(n->parent1().node())) == mapping.end()) {
            pending.push(n->parent1().node());
        }
        /* parent 2 not processed yet */
        else if ((p2 = mapping.find(n->parent2().node())) == mapping.end()) {
            pending.push(n->parent2().node());
        }
        /* process node */
        else {
            /* increase parents' counts */
            p1->second.first += 1;
            p2->second.first += 1;

            mapping[n] = IntBDD(0, 0);
            pending.pop();

            ++remainingNodes;
        }
    }

    for (const auto& p : roots) {
        if (!(p.isConstant())) {
            /* increase root nodes "used count" to avoid premature removal from
             * mapping */
            mapping[p.getInternal().node()].first++;
            pending.push(p.getInternal().node());
        }
    }

    while (!pending.empty()) {
        n = pending.top();

        if (n->flag<1>()) {
            pending.pop();
            continue;
        }

        assert(!(n->isVar()));

        if (!(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
        } else if (!(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
        } else {
            p1 = mapping.find(n->parent1().node());
            assert(p1 != mapping.end());

            DdNode* b1 = p1->second.second;
            assert(b1 != 0);

            if (n->parent1().isInverted()) {
                b1 = Cudd_Not(b1);
            }

            if ((--(p1->second.first)) == 0) {
                mapping.erase(p1);
            } else {
                Cudd_Ref(b1);
            }

            p2 = mapping.find(n->parent2().node());
            assert(p2 != mapping.end());

            DdNode* b2 = p2->second.second;
            assert(b2 != 0);

            if (n->parent2().isInverted()) {
                b2 = Cudd_Not(b2);
            }

            if ((--(p2->second.first)) == 0) {
                mapping.erase(p2);
            } else {
                Cudd_Ref(b2);
            }

            DdNode* b = Cudd_bddAnd(bddMan, b1, b2);
            assert(b != 0);

            if (n->isNAND()) {
                b = Cudd_Not(b);
            }
            Cudd_Ref(b);

            Cudd_RecursiveDeref(bddMan, b1);
            Cudd_RecursiveDeref(bddMan, b2);

            mapping[n].second = b;

            n->setFlag<1>();
            processed.push_back(n);

            pending.pop();

            --remainingNodes;
        }
    }

    std::vector<DdNode*> result;

    for (const auto& p : roots) {
        DdNode* b;

        if (!(p.isConstant())) {
            b = mapping[p.getInternal().node()].second;
            assert(b != 0);
        } else {
            b = Cudd_Not(Cudd_ReadOne(bddMan));
            assert(b != 0);
        }

        if (p.getInternal().isInverted()) {
            b = Cudd_Not(b);
        }
        Cudd_Ref(b);

        result.push_back(b);
    }

    for (const auto& p : mapping) {
        Cudd_RecursiveDeref(bddMan, p.second.second);
    }

    for (Node* p : processed) {
        p->unsetFlag<1>();
    }
    unlockFlags<1>();

    return result;
}
