/**************************************************************
 *
 *       AIGPP Package // Manager_and.cc
 *
 *       Copyright (C) 2006, 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "Manager.hh"

/* LRABS utilities */
#include <lrabsutil/Resources.hh>

/* local */
#include "StatManager.hh"

#define _unused(x) ((void)(x))

#define REWRITE

bool
aigpp::Manager::AndTrivialCase(const aigpp::InternalEdgeRef& n1, const aigpp::InternalEdgeRef& n2,
                               aigpp::InternalEdgeRef& result)
{
    /* checked before in Manager::And */
    assert(!(n1.isConstant() || n2.isConstant()));
    assert(n1.node() != n2.node());

    /*
     * level 1.5 trivial cases:
     * n1 is one of n2's parents or vice versa.
     */

    /* (!)a*(!)((!)a*b) or (!)a*(!)(b*(!)a) */
    if ((!(n2.node()->isVar())
         && (n1.node() == n2.node()->parent1().node() || n1.node() == n2.node()->parent2().node()))
        || (!(n1.node()->isVar())
            && (n2.node() == n1.node()->parent1().node() || n2.node() == n1.node()->parent2().node()))) {
        StatManager::instance().incAndTrivial15();

        /*
         * naming convention used here:
         * a1 := the level1 parent (n1 or n2) which is equal to the others second
         * level parent a2 := the level1 parent (n1 or n2) which is NOT equal to the
         * others second level parent b1 := the matching level 2 parent of a2 b2 :=
         * the NOT matching level 2 parent of a2
         */
        const Edge& a1 = (!(n2.node()->isVar())
                          && (n1.node() == n2.node()->parent1().node() || n1.node() == n2.node()->parent2().node()))
                             ? n1
                             : n2;
        const Edge& a2 = (structurallyEquivalent(a1, n1) ? n2 : n1);

        const Edge& b1 = (a1.node() == a2.node()->parent1().node()) ? a2.node()->parent1() : a2.node()->parent2();
        const Edge& b2 = (a1.node() == a2.node()->parent1().node()) ? a2.node()->parent2() : a2.node()->parent1();

        /* a1*(!)(b1[==a1]*b2) or !a1*(!)(!b1[==a1]*b2) */
        if (a1.isInverted() == b1.isInverted()) {
            /* a1*(b1[==a1]*b2) [== b1*b2 == a2] or !a1*(!b1[==a1]*b2) [== !b1*b2 ==
             * a2] */
            if (a2.isInverted() == a2.node()->isNAND()) {
                result = InternalEdgeRef(a2);
                return true;
            }
            /* a1*!(b1[==a1]*b2) [== a1*!b2] or !a1*!(!b1[==a1]*b2) [== !a1*!b2] */
            else {
#ifdef REWRITE
                result = InternalEdgeRef(a1) * InternalEdgeRef(!b2);
                return true;
#endif
            }
        }
        /* a1*(!)(!b1[==a1]*b2) or !a1*(!)(b1[==a1]*b2) */
        else {
            /* a1*(!b1[==a1]*b2) [== 0] or !a1*(b1[==a1]*b2) [== 0] */
            if (a2.isInverted() == a2.node()->isNAND()) {
                result = const0;
                return true;
            }
            /* a1*!(!b1[==a1]*b2) [== a1] or !a1*!(b1[==a1]*b2) [== !a1] */
            else {
                result = InternalEdgeRef(a1);
                return true;
            }
        }

        /* if we end up here, some of the above cases did not return */
        StatManager::instance().decAndTrivial15();
        return false;
    }

    /*
     * level 2.a trivial cases:
     * nodes share both parents.
     */
    /* (!)((!)a*(!)b)*(!)((!)a*(!)b) */
    if (!(n1.node()->isVar() || n2.node()->isVar())
        && ((n1.node()->parent1().node() == n2.node()->parent1().node()
             && n1.node()->parent2().node() == n2.node()->parent2().node())
            || (n1.node()->parent1().node() == n2.node()->parent2().node()
                && n1.node()->parent2().node() == n2.node()->parent1().node()))) {
        StatManager::instance().incAndTrivial2a();

        /* ((!)a*(!)b)*(!)((!)a*(!)b) */
        if (n1.isInverted() == n1.node()->isNAND()) {
            /* ((!)a*(!)b)*((!)a*(!)b) */
            if (n2.isInverted() == n2.node()->isNAND()) {
                /* since n1 and n2 cannot be equal (triv 1 check!), the result must be 0
                 */
                /* ^ wrong assumption (due to rewriting and non-fraiging) */
                // result = const0;

                /* nodes actually equivalent */
                if ((structurallyEquivalent(n1.node()->parent1(), n2.node()->parent1())
                     && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent2()))
                    || (structurallyEquivalent(n1.node()->parent1(), n2.node()->parent2())
                        && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent1()))) {
                    // std::cout << "WARNING: and-triv 2.a (1)" << std::endl;
                    result = n1;
                }
                /* at least one parent negation is inverse */
                else {
                    assert(structurallyAntivalent(n1.node()->parent1(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent1(), n2.node()->parent2())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent2()));
                    result = const0;
                }

                return true;
            }
            /* ((!)a*(!)b)*!((!)a*(!)b) */
            else {
                /* since n1 and n2 cannot be equal (triv 1 check!), the result must be
                 * n1 */
                /* ^ wrong assumption (due to rewriting and non-fraiging) */
                // result = n1;

                /* nodes actually equivalent */
                if ((structurallyEquivalent(n1.node()->parent1(), n2.node()->parent1())
                     && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent2()))
                    || (structurallyEquivalent(n1.node()->parent1(), n2.node()->parent2())
                        && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent1()))) {
                    // std::cout << "WARNING: and-triv 2.a (2)" << std::endl;
                    result = const0;
                }
                /* at least one parent negation is inverse */
                else {
                    assert(structurallyAntivalent(n1.node()->parent1(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent1(), n2.node()->parent2())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent2()));
                    result = n1;
                }

                return true;
            }
        }
        /* !((!)a*(!)b)*(!)((!)a*(!)b) */
        else {
            /* !((!)a*(!)b)*((!)a*(!)b) */
            if (n2.isInverted() == n2.node()->isNAND()) {
                /* since n1 and n2 cannot be equal (triv 1 check!), the result must be
                 * n2 */
                /* ^ wrong assumption (due to rewriting and non-fraiging) */
                // result = n2;

                /* nodes actually equivalent */
                if ((structurallyEquivalent(n1.node()->parent1(), n2.node()->parent1())
                     && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent2()))
                    || (structurallyEquivalent(n1.node()->parent1(), n2.node()->parent2())
                        && structurallyEquivalent(n1.node()->parent2(), n2.node()->parent1()))) {
                    // std::cout << "WARNING: and-triv 2.a (3)" << std::endl;
                    result = const0;
                }
                /* at least one parent negation is inverse */
                else {
                    assert(structurallyAntivalent(n1.node()->parent1(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent1(), n2.node()->parent2())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent1())
                           || structurallyAntivalent(n1.node()->parent2(), n2.node()->parent2()));
                    result = n2;
                }

                return true;
            }
            /* !((!)a*(!)b)*!((!)a*(!)b) */
            else {
                /*
                 * naming convention used here:
                 * a1 := n1 p1
                 * a2 := n1 p2
                 * b1 := n2's parent matching a1
                 * b2 := n2's parent matching a2
                 */

                const Edge& a1 = n1.node()->parent1();
                const Edge& a2 = n1.node()->parent2();
                const Edge& b1
                    = (a1.node() == n2.node()->parent1().node()) ? n2.node()->parent1() : n2.node()->parent2();
                const Edge& b2
                    = (a1.node() == n2.node()->parent1().node()) ? n2.node()->parent2() : n2.node()->parent1();

                /* !(a*b)*!(!a*b) = (!a+!b)*(a+!b) = 0+!a!b+a!b+!b = !b */
                if (a1.isInverted() != b1.isInverted() && a2.isInverted() == b2.isInverted()) {
                    result = InternalEdgeRef(!a2);
                    return true;
                }
                /* !(a*b)*!(a*!b) = (!a+!b)*(!a+b) = !a+!ab+!a!b+0 = !a */
                else if (a1.isInverted() == b1.isInverted() && a2.isInverted() != b2.isInverted()) {
                    result = InternalEdgeRef(!a1);
                    return true;
                }
                /* !(a*b)*!(a*b) can not occur (triv1!) */
                /* ^ wrong assumption (due to rewriting and non-fraiging) */
                else if (a1.isInverted() == b1.isInverted() && a2.isInverted() == b2.isInverted()) {
                    // std::cout << "WARNING: and-triv 2.a (4)" << std::endl;
                    result = n1;
                    return true;
                }

                /* remaining case here: !(a*b)*!(!a*!b) == xor( a, b ) -> not
                 * rewritable! */
            }
        }

        /* if we end up here, some of the above cases did not return */
        StatManager::instance().decAndTrivial2a();
        return false;
    }

    /*
     * level 2.b trivial cases:
     * nodes share one parent.
     */
    /* (!)((!)a*b)*(!)((!)a*c) and permutaions/invertions */
    if (!(n1.node()->isVar() || n2.node()->isVar())
        && ((n1.node()->parent1().node() == n2.node()->parent1().node())
            || (n1.node()->parent1().node() == n2.node()->parent2().node())
            || (n1.node()->parent2().node() == n2.node()->parent1().node())
            || (n1.node()->parent2().node() == n2.node()->parent2().node()))) {
        StatManager::instance().incAndTrivial2b();

        /*
         * naming convention used here:
         * a1 := n1's parent matching n2p1 or n2p2
         * a2 := n1' parent NOT matching -||-
         * b1 := n2's parent matching a1 or a2
         * b2 := n2's parent NOT matching a1 or 2
         */
        const Edge& a1 = ((n1.node()->parent1().node() == n2.node()->parent1().node())
                          || (n1.node()->parent1().node() == n2.node()->parent2().node()))
                             ? n1.node()->parent1()
                             : n1.node()->parent2();
        const Edge& a2
            = (structurallyEquivalent(a1, n1.node()->parent1())) ? n1.node()->parent2() : n1.node()->parent1();
        assert(a1.node() == n2.node()->parent1().node() || a1.node() == n2.node()->parent2().node());

        const Edge& b1 = (a1.node() == n2.node()->parent1().node()) ? n2.node()->parent1() : n2.node()->parent2();
        const Edge& b2 = (a1.node() == n2.node()->parent1().node()) ? n2.node()->parent2() : n2.node()->parent1();

        /* ((!)a*b)*(!)((!)a*c) */
        if (n1.isInverted() == n1.node()->isNAND()) {
            /* ((!)a*b)*((!)a*c) */
            if (n2.isInverted() == n2.node()->isNAND()) {
                /* (a*b)*(a*c) */
                if (a1.isInverted() == b1.isInverted()) {
#ifdef REWRITE
                    result = n1 * InternalEdgeRef(b2);
                    // alternative: result = InternalEdgeRef( a2 ) * n2;
                    return true;
#endif
                }
                /* (a*b)*(!a*c) */
                else {
                    result = const0;
                    return true;
                }
            }
            /* ((!)a*b)*!((!)a*c) */
            else {
                /* (a*b)*!(a*c) = (a*b)*!c */
                if (a1.isInverted() == b1.isInverted()) {
#ifdef REWRITE
                    result = n1 * InternalEdgeRef(!b2);
                    return true;
#endif
                }
                /* (a*b)*!(!a*c) = (a*b) */
                else {
                    result = n1;
                    return true;
                }
            }
        }
        /* !((!)a*b)*(!)((!)a*c) */
        else {
            /* !((!)a*b)*((!)a*c) */
            if (n2.isInverted() == n2.node()->isNAND()) {
                /* !(a*b)*(a*c) = !b*(a*c) */
                if (a1.isInverted() == b1.isInverted()) {
#ifdef REWRITE
                    result = n2 * InternalEdgeRef(!a2);
                    return true;
#endif
                }
                /* !(a*b)*(!a*c) = (!a*c) */
                else {
                    result = n2;
                    return true;
                }
            }
            /* !((!)a*b)*!((!)a*c) */
            else {
                /* !(a*b)*!(a*c) = !a + !b!c */
                if (a1.isInverted() == b1.isInverted()) {
                    /* this rewriting will create two new nodes! */
                    /*
          #ifdef REWRITE
                    result = InternalEdgeRef( !a1 ) + ( InternalEdgeRef( !a2 ) *
          InternalEdgeRef( !b2 ) ); return true; #endif
                    */
                }
                /* !(a*b)*!(!a*c) */
                else {
                    /* not rewritable! */
                }
            }
        }

        StatManager::instance().decAndTrivial2b();
        return false;
    }

    /* no trivial case */
    return false;
}

#ifdef USE_TEMPSIM
void
aigpp::Manager::removeUnequalTempSim(std::set<aigpp::Node*>& list, unsigned long tempSim) const
{
    for (auto c = list.begin(); c != list.end(); /**/) {
        if (tempSim != (*c)->tempSim()) {
            list.erase(c++);
        } else {
            ++c;
        }
    }
}
#endif

// TODO: adapt to cyclic simclasses
void
aigpp::Manager::initializeCandidates(std::set<aigpp::Node*>& list, aigpp::Node* simclass
#ifdef USE_TEMPSIM
                                     ,
                                     unsigned long tempSim
#endif
                                     ) const
{
    list.clear();

    if (simclass != nullptr) {
        Node* n = simclass;
        do {
            if (_skipTestWith != nullptr && n == _skipTestWith) {
                /*continue;*/
            }
#ifdef USE_TEMPSIM
            else if (tempSim != n->tempSim()) {
                /*continue;*/
            }
#endif
            else {
                list.insert(n);
            }

            n = n->_nEqualSim;
        } while (n != simclass);
    }
}

// TODO: adapt to cyclic simclasses
void
aigpp::Manager::updateCandidates(std::set<aigpp::Node*>& list, aigpp::Node* newsimclass) const
{
    if (newsimclass == nullptr) {
        list.clear();
    } else {
        std::set<Node*> newlist;

        Node* n = newsimclass;
        do {
            if (list.find(n) != list.end()) {
                newlist.insert(n);
            }

            n = n->_nEqualSim;
        } while (n != newsimclass);

        if (newlist.size() != list.size()) {
            list.swap(newlist);
        }
    }
}

aigpp::InternalEdgeRef
aigpp::Manager::And(const aigpp::InternalEdgeRef& n1, const aigpp::InternalEdgeRef& n2)
{
    if (verbosity() >= 2) {
        std::cout << "and: " << (n1.isInverted() ? "!" : "") << n1.node()->index() << "   "
                  << (n2.isInverted() ? "!" : "") << n2.node()->index() << std::endl;
    }

    /*
     * level 0 trivial cases:
     * n1 or n2 is constant
     * -> already handled
     */
    assert(!(n1.isConstant() || n2.isConstant()));

    /*
     * level 1 trivial cases:
     * n1 and n2 are equal or inverse.
     */
    if (n1.node() == n2.node()) {
        StatManager::instance().incAndTrivial1();
        return ((n1.isInverted() == n2.isInverted()) ? n1 : const0);
    }

    /*
     * canonize the edge order based on the nodes' indices
     */
    Edge n1p = n1;
    Edge n2p = n2;
    if (n1p.node()->index() > n2p.node()->index()) {
        n1p.swap(n2p);
    }

    /*
     * unique and computed table lookups
     */
    {
        Edge result;

        /*
         * unique table lookup
         */
        if (_unique.lookup(n1p, n2p, result)) {
            reclaimIfDead(result.node());
            if (verbosity() >= 3) {
                std::cout << "unique hit" << std::endl;
            }

            return InternalEdgeRef(result);
        }

#ifdef USE_COMPUTEDTABLE
        /*
         * computed table lookup
         */
        if (_computed.lookupAnd(n1p, n2p, result)) {
            reclaimIfDead(result.node());
            return InternalEdgeRef(result);
        }
#endif
    }

    /*
     * Apply garbage collection if necessary
     */
    if ((double)_dead > settings().getGCRatio() * (double)(_dead + _nodeCount)) {
        garbageCollect();
    }

    /* check trivial cases */
    if (!_nextAndSkipSat) {
        InternalEdgeRef result;
        if (AndTrivialCase(n1, n2, result)) {
            return result;
        }
    }

    /*
     * create a new temporary node. Use an abandoned Node object from _freeNodes
     * if possible.
     */
    Node* newnode = createNode(n1p, n2p);
    if (simCreation()) {
        assert(newnode->_sim != 0);
    }
    newnode->_index = _nextIndex;

    /*
     * the result of this and operation is known to be equivalent/antivalent to
     * _nextAndSatResult
     */
    if (_nextAndSkipSat) {
        StatManager::instance().incAndKnown();

        /*
         * is the the new node a better representation than the old node?
         */
        if (compareQuality(newnode, _nextAndSatResult) < 0) {
            if (_nextAndSatResultIsInverted) {
                resetSolver(true);
            }

            if (verbosity() >= 3) {
                std::cout << "replacement by BDD sweeping" << std::endl;
            }

            replace(_nextAndSatResult, newnode, _nextAndSatResultIsInverted);
        }

        releaseNode(newnode);

        return InternalEdgeRef(_nextAndSatResult, _nextAndSatResultIsInverted);
    }

    if (_automaticFRAIGing) {
        assert(simCreation() == true);

        /*
         * reset solver if needed
         */
        resetSolver();

        /*
         * Rebuild the simtable/simulation of all nodes if necessary.
         * And update the simulation of the temporary node since it is not
         * contained in the nodes list!!!
         */
        updateSimTable();
        newnode->updateSim();

#ifdef USE_TEMPSIM
        /* update tempsim if needed */
        propagateTempSim(newnode);
#endif

        /* candidates for equivalence */

        /* equivalence candidates before zero checking */
        std::set<Node*> candidates;
        bool            haveCandidatesAfterZero = false;

        /*
         * if the new node's simulation vector is zero,
         * it is probably equal to zero
         */

        if (
#ifdef USE_TEMPSIM
            newnode->tempSim() == 0ul &&
#endif
            newnode->sim().isZero()) {
            StatManager::instance().incSimZero();

            /* node IS zero */
            if (satEqualZero(newnode)) {
                StatManager::instance().incAndSATZero();
                StatManager::instance().incAndSAT();

#ifdef USE_COMPUTEDTABLE
                _computed.insertAnd(n1p, n2p, const0);
#endif

                _nodesWithSATVars.erase(newnode);
                newnode->invalidateSatVar();

                releaseNode(newnode);

                return const0;
            }

            /* node is not zero */

            haveCandidatesAfterZero = true;

#ifdef USE_TEMPSIM
            initializeCandidates(candidates, _simTable.lookup(newnode->sim()), newnode->tempSim());
#else
            initializeCandidates(candidates, _simTable.lookup(newnode->sim()));
#endif

#ifdef USE_TEMPSIM
            /*
             * immediately insert the new counterexample into the simtable
             * -> newnode is removed from the class of nodes with sim=0
             */
            if (!candidates.empty()) {
                updateTempSim(_newCounterExamples.back(), newnode);
                removeUnequalTempSim(candidates, newnode->tempSim());
            }
#endif
        }

        /*
         * check all nodes with the same simvector for equivalence
         */
        // TODO: adapt to cyclic simclasses
        Node* simclass = _simTable.lookup(newnode->sim());
        if (simclass != nullptr) {
            /*
             * collect initial set of equivalence candidates,
             * use candidates list from zero checking if it is available!
             */
            if (!haveCandidatesAfterZero) {
#ifdef USE_TEMPSIM
                initializeCandidates(candidates, simclass, newnode->tempSim());
#else
                initializeCandidates(candidates, simclass);
#endif
            } else {
                updateCandidates(candidates, simclass);
            }

            /*
             * check all candidates
             */
            while (!candidates.empty()) {
                /*
                 * If we collected more than SimVector::BinSize counterexamples
                 * integrate them and update the candidate list
                 */
                if (_newCounterExamples.size() >= SimVector::BinSize) {
                    updateSimTable();
                    newnode->updateSim();

                    /* update candidates set */
                    updateCandidates(candidates, _simTable.lookup(newnode->sim()));
                    if (candidates.empty()) {
                        break;
                    }
                }

                /*
                 * get and erase the first candidate
                 */
                Node* simnode = *(candidates.begin());
                candidates.erase(candidates.begin());

                /*
                 * check equivalence
                 */
                if (satEqual(newnode, simnode)) {
                    if (verbosity() >= 3) {
                        std::cout << "equivalent: " << simnode->index() << std::endl;
                    }

                    StatManager::instance().incAndSAT();

                    reclaimIfDead(simnode);

                    if (compareQuality(newnode, simnode) < 0) {
#ifdef USE_COMPUTEDTABLE
                        _computed.insertAnd(simnode->parent1(), simnode->parent2(), Edge(simnode));
#endif
                        replace(simnode, newnode, false);

                        simnode->unsetFlag<Node::FLAG_UNMODIFIED>();
                    }
#ifdef USE_COMPUTEDTABLE
                    else {
                        _computed.insertAnd(n1p, n2p, Edge(simnode));
                    }
#endif

                    if (newnode->isSatVarValid()) {
                        _nodesWithSATVars.erase(newnode);
                        newnode->invalidateSatVar();
                    }

                    releaseNode(newnode);

                    return InternalEdgeRef(simnode);
                }

                if (!candidates.empty()) {
                    if (_newCounterExamples.size() >= SimVector::BinSize) {
                        updateSimTable();
                        newnode->updateSim();

                        /* update candidates set */
                        updateCandidates(candidates, _simTable.lookup(newnode->sim()));
                    }
#ifdef USE_TEMPSIM
                    /* immediately propagate the latest counterexample to remove some
                       candidates */
                    else {
                        updateTempSim(_newCounterExamples.back(), newnode);
                        removeUnequalTempSim(candidates, newnode->tempSim());
                    }
#endif
                }
            }
        }

        /*
         * no equivalent node was found.
         * check nodes with complementary simvector for antivalence.
         */
        simclass = _simTable.lookupInverted(newnode->sim(), simclass);

        if (simclass != nullptr) {
            /* get initial set of antivalence candidates */
#ifdef USE_TEMPSIM
            initializeCandidates(candidates, simclass, ~(newnode->tempSim()));
#else
            initializeCandidates(candidates, simclass);
#endif

            while (!candidates.empty()) {
                /*
                 * If we collected more than SimVector::BinSize counterexamples
                 * integrate them and update the candidate list
                 */
                if (_newCounterExamples.size() >= SimVector::BinSize) {
                    updateSimTable();
                    newnode->updateSim();

                    /* update candidates set */
                    updateCandidates(candidates, _simTable.lookupInverted(newnode->sim()));
                    if (candidates.empty()) {
                        break;
                    }
                }

                /*
                 * get and erase the first candidate
                 */
                Node* simnode = *(candidates.begin());
                candidates.erase(candidates.begin());

                if (satEqualNegated(newnode, simnode)) {
                    if (verbosity() >= 3) {
                        std::cout << "antivalent: " << simnode->index() << std::endl;
                    }

                    StatManager::instance().incAndSAT();

                    reclaimIfDead(simnode);

                    if (compareQuality(newnode, simnode) < 0) {
                        resetSolver(true);

#ifdef USE_COMPUTEDTABLE
                        _computed.insertAnd(simnode->parent1(), simnode->parent2(), Edge(simnode));
#endif

                        replace(simnode, newnode, true);
                        simnode->unsetFlag<Node::FLAG_UNMODIFIED>();
                    }
#ifdef USE_COMPUTEDTABLE
                    else {
                        _computed.insertAnd(n1p, n2p, Edge(simnode, true));
                    }
#endif

                    if (newnode->isSatVarValid()) {
                        _nodesWithSATVars.erase(newnode);
                        newnode->invalidateSatVar();
                    }

                    releaseNode(newnode);

                    return InternalEdgeRef(simnode, true);
                }

                if (!candidates.empty()) {
                    if (_newCounterExamples.size() >= SimVector::BinSize) {
                        updateSimTable();
                        newnode->updateSim();

                        /* update candidates set */
                        updateCandidates(candidates, _simTable.lookupInverted(newnode->sim()));
                    }
#ifdef USE_TEMPSIM
                    /* immediately propagate the latest counterexample to remove some
                       candidates */
                    else {
                        updateTempSim(_newCounterExamples.back(), newnode);
                        removeUnequalTempSim(candidates, ~(newnode->tempSim()));
                    }
#endif
                }
            }
        }

        /*
         * no equivalent or antivalent node found.
         * node is reduced!
         */
        newnode->setFlag<Node::FLAG_ISREDUCED>();
    }

    /*
     * no equivalent or antivalent node found.
     * insert new node into simtable ans strashing structure.
     */

    if (verbosity() >= 3) {
        std::cout << "new: " << newnode->index() << std::endl;
    }

    newnode->unsetFlag<Node::FLAG_UNMODIFIED>();

    if (!(newnode->flag<Node::FLAG_ISREDUCED>())) {
        ++_unreducedNodes;
    }

    StatManager::instance().incAndNew();

    ref(n1.node());
    ref(n2.node());

#ifndef NDEBUG
    int nexti = nextIndex();
    assert(newnode->_index == nexti);
#endif

    if (simCreation()) {
        _simTable.insert(newnode);
    }

    addToNodesList(newnode);
    _unique.insert(newnode);

    return InternalEdgeRef(newnode);
}

aigpp::InternalEdgeRef
aigpp::Manager::SimpleAnd(const aigpp::InternalEdgeRef& n1, const aigpp::InternalEdgeRef& n2)
{
    if (verbosity() >= 2) {
        std::cout << "and: " << (n1.isInverted() ? "!" : "") << n1.node()->index() << "   "
                  << (n2.isInverted() ? "!" : "") << n2.node()->index() << std::endl;
    }

    /*
     * level 0 trivial cases:
     * n1 or n2 is constant
     */
    if (n1.isConstant()) {
        if (n1.isInverted()) /* constant "true" */
        {
            return n2;
        } else /* constant "false" */
        {
            return n1;
        }
    } else if (n2.isConstant()) {
        if (n2.isInverted()) /* constant "true" */
        {
            return n1;
        } else /* constant "false" */
        {
            return n2;
        }
    }

    /*
     * level 1 trivial cases:
     * n1 and n2 are equal or inverse.
     */
    if (n1.node() == n2.node()) {
        StatManager::instance().incAndTrivial1();
        return ((n1.isInverted() == n2.isInverted()) ? n1 : const0);
    }

    /*
     * canonize the edge order based on the nodes' indices
     */
    Edge n1p = n1;
    Edge n2p = n2;
    if (n1p.node()->index() > n2p.node()->index()) {
        n1p.swap(n2p);
    }

    /*
     * unique table lookup
     */
    {
        Edge result;

        /*
         * unique table lookup
         */
        if (_unique.lookup(n1p, n2p, result)) {
            reclaimIfDead(result.node());
            return InternalEdgeRef(result);
        }
    }

    /*
     * create a new temporary node
     */
    Node* newnode   = createNode(n1p, n2p);
    newnode->_index = nextIndex();

    newnode->unsetFlag<Node::FLAG_UNMODIFIED>();
    newnode->unsetFlag<Node::FLAG_ISREDUCED>();

    /*
     * insert new node into simtable ans strashing structure.
     */
    if (verbosity() >= 3) {
        std::cout << "new: " << newnode->index() << std::endl;
    }

    ++_unreducedNodes;

    StatManager::instance().incAndNew();

    ref(n1.node());
    ref(n2.node());

    if (simCreation()) {
        _simTable.insert(newnode);
    }

    addToNodesList(newnode);
    _unique.insert(newnode);

    return InternalEdgeRef(newnode);
}

aigpp::InternalEdgeRef
aigpp::Manager::DebugAnd(const aigpp::InternalEdgeRef& n1, const aigpp::InternalEdgeRef& n2, bool createNAND,
                         bool doTrivial, bool doUnique, bool doSAT)
{
    /* constant checks */
    if (n1.isConstant()) {
        return (n1.isInverted() ? n2 : const0).notIf(createNAND);
    } else if (n2.isConstant()) {
        return (n2.isInverted() ? n1 : const0).notIf(createNAND);
    }

    /*
     * level 1 trivial cases:
     * n1 and n2 are equal or inverse.
     */
    if (n1.node() == n2.node()) {
        return ((n1.isInverted() == n2.isInverted()) ? n1 : const0).notIf(createNAND);
    }

    /*
     * canonize the edge order based on the nodes' indices
     */
    Edge n1p = n1;
    Edge n2p = n2;
    if (n1p.node()->index() > n2p.node()->index()) {
        n1p.swap(n2p);
    }

    /*
     * unique table lookup
     */
    if (doUnique) {
        Edge result;

        /*
         * unique table lookup
         */
        if (_unique.lookup(n1p, n2p, result)) {
            reclaimIfDead(result.node());
            return InternalEdgeRef(result).notIf(createNAND);
        }
    }

    /*
     * Apply garbage collection if necessary
     */
    if ((double)_dead > settings().getGCRatio() * (double)(_dead + _nodeCount)) {
        garbageCollect();
    }

    /* check trivial cases */
    if (doTrivial) {
        InternalEdgeRef result;
        if (AndTrivialCase(n1, n2, result)) {
            return result;
        }
    }

    /*
     * create a new temporary node. Use an abandoned Node object from _freeNodes
     * if possible.
     */
    Node* newnode = createNode(n1p, n2p);
    if (simCreation()) {
        assert(newnode->_sim != 0);
    }
    newnode->_index = _nextIndex;

    if (createNAND) {
        newnode->setFlag<Node::FLAG_ISNAND>();
    }

    if (doSAT) {
        assert(simCreation() == true);

        /*
         * reset solver if needed
         */
        resetSolver();

        /*
         * Rebuild the simtable/simulation of all nodes if necessary.
         * And update the simulation of the temporary node since it is not
         * contained in the nodes list!!!
         */
        updateSimTable();
        newnode->updateSim();

#ifdef USE_TEMPSIM
        /* update tempsim if needed */
        propagateTempSim(newnode);
#endif

        /* candidates for equivalence */

        /* equivalence candidates before zero checking */
        std::set<Node*> candidates;
        bool            haveCandidatesAfterZero = false;

        /*
         * if the new node's simulation vector is zero,
         * it is probably equal to zero
         */

        if (
#ifdef USE_TEMPSIM
            newnode->tempSim() == 0ul &&
#endif
            newnode->sim().isZero()) {
            StatManager::instance().incSimZero();

            /* node IS zero */
            if (satEqualZero(newnode)) {
                StatManager::instance().incAndSATZero();
                StatManager::instance().incAndSAT();

                _nodesWithSATVars.erase(newnode);
                newnode->invalidateSatVar();

                releaseNode(newnode);

                return const0;
            }

            /* node is not zero */

            haveCandidatesAfterZero = true;

#ifdef USE_TEMPSIM
            // TODO: adapt to cyclic simclasses
            initializeCandidates(candidates, _simTable.lookup(newnode->sim()), newnode->tempSim());
#else
            initializeCandidates(candidates, _simTable.lookup(newnode->sim()));
#endif

#ifdef USE_TEMPSIM
            /*
             * immediately insert the new counterexample into the simtable
             * -> newnode is removed from the class of nodes with sim=0
             */
            if (!candidates.empty()) {
                updateTempSim(_newCounterExamples.back(), newnode);
                removeUnequalTempSim(candidates, newnode->tempSim());
            }
#endif
        }

        /*
         * check all nodes with the same simvector for equivalence
         */
        // TODO: adapt to cyclic simclasses
        Node* simclass = _simTable.lookup(newnode->sim());
        if (simclass != nullptr) {
            /*
             * collect initial set of equivalence candidates,
             * use candidates list from zero checking if it is available!
             */
            if (!haveCandidatesAfterZero) {
#ifdef USE_TEMPSIM
                initializeCandidates(candidates, simclass, newnode->tempSim());
#else
                initializeCandidates(candidates, simclass);
#endif
            } else {
                updateCandidates(candidates, simclass);
            }

            /*
             * check all candidates
             */
            while (!candidates.empty()) {
                /*
                 * If we collected more than SimVector::BinSize counterexamples
                 * integrate them and update the candidate list
                 */
                if (_newCounterExamples.size() >= SimVector::BinSize) {
                    updateSimTable();
                    newnode->updateSim();

                    /* update candidates set */
                    updateCandidates(candidates, _simTable.lookup(newnode->sim()));
                    if (candidates.empty()) {
                        break;
                    }
                }

                /*
                 * get and erase the first candidate
                 */
                Node* simnode = *(candidates.begin());
                candidates.erase(candidates.begin());

                /*
                 * check equivalence
                 */
                if (satEqual(newnode, simnode)) {
                    if (verbosity() >= 3) {
                        std::cout << "equivalent: " << simnode->index() << std::endl;
                    }

                    StatManager::instance().incAndSAT();

                    reclaimIfDead(simnode);

                    if (compareQuality(newnode, simnode) < 0) {
                        replace(simnode, newnode, false);
                        simnode->unsetFlag<Node::FLAG_UNMODIFIED>();
                    }

                    if (newnode->isSatVarValid()) {
                        _nodesWithSATVars.erase(newnode);
                        newnode->invalidateSatVar();
                    }

                    releaseNode(newnode);

                    return InternalEdgeRef(simnode);
                }

                if (!candidates.empty()) {
                    if (_newCounterExamples.size() >= SimVector::BinSize) {
                        updateSimTable();
                        newnode->updateSim();

                        /* update candidates set */
                        updateCandidates(candidates, _simTable.lookup(newnode->sim()));
                    }
#ifdef USE_TEMPSIM
                    /* immediately propagate the latest counterexample to remove some
                       candidates */
                    else {
                        updateTempSim(_newCounterExamples.back(), newnode);
                        removeUnequalTempSim(candidates, newnode->tempSim());
                    }
#endif
                }
            }
        }

        /*
         * no equivalent node was found.
         * check nodes with complementary simvector for antivalence.
         */
        // TODO: adapt to cyclic simclasses
        simclass = _simTable.lookupInverted(newnode->sim());

        if (simclass != nullptr) {
            /* get initial set of antivalence candidates */
#ifdef USE_TEMPSIM
            initializeCandidates(candidates, simclass, ~(newnode->tempSim()));
#else
            initializeCandidates(candidates, simclass);
#endif

            while (!candidates.empty()) {
                /*
                 * If we collected more than SimVector::BinSize counterexamples
                 * integrate them and update the candidate list
                 */
                if (_newCounterExamples.size() >= SimVector::BinSize) {
                    updateSimTable();
                    newnode->updateSim();

                    /* update candidates set */
                    updateCandidates(candidates, _simTable.lookupInverted(newnode->sim()));
                    if (candidates.empty()) {
                        break;
                    }
                }

                /*
                 * get and erase the first candidate
                 */
                Node* simnode = *(candidates.begin());
                candidates.erase(candidates.begin());

                if (satEqualNegated(newnode, simnode)) {
                    if (verbosity() >= 3) {
                        std::cout << "antivalent: " << simnode->index() << std::endl;
                    }

                    StatManager::instance().incAndSAT();

                    reclaimIfDead(simnode);

                    if (compareQuality(newnode, simnode) < 0) {
                        resetSolver(true);
                        replace(simnode, newnode, true);
                        simnode->unsetFlag<Node::FLAG_UNMODIFIED>();
                    }

                    if (newnode->isSatVarValid()) {
                        _nodesWithSATVars.erase(newnode);
                        newnode->invalidateSatVar();
                    }

                    releaseNode(newnode);

                    return InternalEdgeRef(simnode, true);
                }

                if (!candidates.empty()) {
                    if (_newCounterExamples.size() >= SimVector::BinSize) {
                        updateSimTable();
                        newnode->updateSim();

                        /* update candidates set */
                        updateCandidates(candidates, _simTable.lookupInverted(newnode->sim()));
                    }
#ifdef USE_TEMPSIM
                    /* immediately propagate the latest counterexample to remove some
                       candidates */
                    else {
                        updateTempSim(_newCounterExamples.back(), newnode);
                        removeUnequalTempSim(candidates, ~(newnode->tempSim()));
                    }
#endif
                }
            }
        }

        /*
         * no equivalent or antivalent node found.
         * node is reduced!
         */
        newnode->setFlag<Node::FLAG_ISREDUCED>();
    }

    /*
     * no equivalent or antivalent node found.
     * insert new node into simtable ans strashing structure.
     */

    newnode->unsetFlag<Node::FLAG_UNMODIFIED>();

    if (!(newnode->flag<Node::FLAG_ISREDUCED>())) {
        ++_unreducedNodes;
    }

    StatManager::instance().incAndNew();

    ref(n1.node());
    ref(n2.node());

#ifndef NDEBUG
    const int nexti = nextIndex();
    assert(newnode->_index == nexti);
#endif

    if (simCreation()) {
        _simTable.insert(newnode);
    }

    addToNodesList(newnode);
    _unique.insert(newnode);

    return InternalEdgeRef(newnode);
}
