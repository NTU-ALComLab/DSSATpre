#include "Manager.hh"

#include <lrabsutil/Resources.hh>
#include <lrabsutil/String.hh>

#include "EquivalenceChecker.hh"

#define _unused(x) ((void)(x))

struct EstimateRemainingNodes_Node
{
    int          _var  = 0;
    bool         _nand = false;
    std::size_t  _p1   = 0;
    bool         _p1n  = false;
    std::size_t  _p2   = 0;
    bool         _p2n  = false;
    unsigned int _sim  = 0;
};

std::vector<long>
aigpp::Manager::estimateRemainingNodesAfterQuantification(aigpp::Node*                     root,
                                                          const std::vector<aigpp::Node*>& variables) const
{
    // double startTime = lrabs::cpuTime();

    std::vector<long> result(variables.size(), -1);

    std::vector<EstimateRemainingNodes_Node> nodes;
    std::map<const Node*, std::size_t>       nodeindex;

    {
        lockFlags<1>();

        std::stack<const Node*>& pending = _universalConstNodeStack;
        assert(pending.empty());

        std::vector<const Node*>& processed = _universalConstNodeVector;
        assert(processed.empty());

        /* collect cone of root -> processed contains a topologically sorted list of
         * root's cone nodes */
        pending.push(root);
        while (!(pending.empty())) {
            const Node* n = pending.top();

            if (n->flag<1>()) {
                pending.pop();
                continue;
            } else if (n->isVar()) {
                n->setFlag<1>();
                processed.push_back(n);

                EstimateRemainingNodes_Node d;
                d._var       = n->varIndex();
                nodeindex[n] = nodes.size();
                nodes.push_back(d);

                pending.pop();
            } else {
                if (!(n->parent1().node()->flag<1>())) {
                    pending.push(n->parent1().node());
                }

                if (!(n->parent2().node()->flag<1>())) {
                    pending.push(n->parent2().node());
                }

                if (pending.top() == n) {
                    n->setFlag<1>();
                    processed.push_back(n);

                    EstimateRemainingNodes_Node d;
                    d._var       = -1;
                    d._nand      = n->isNAND();
                    d._p1        = nodeindex[n->parent1().node()];
                    d._p1n       = n->parent1().isInverted();
                    d._p2        = nodeindex[n->parent2().node()];
                    d._p2n       = n->parent2().isInverted();
                    nodeindex[n] = nodes.size();
                    nodes.push_back(d);

                    pending.pop();
                }
            }
        }

        /* check if all variables are in the support (otherwise set their result to
         * |root|) */
        for (auto v = variables.begin(); v != variables.end(); ++v) {
            if (!((*v)->flag<1>())) {
                // std::cout << "WARNING " << __func__ << " variable '" << **v << "' not
                // in support" << std::endl;
                result[v - variables.begin()] = static_cast<long>(processed.size());
            }
        }

        /* clear flags */
        for (const Node* p : processed) {
            p->unsetFlag<1>();
        }

        processed.clear();

        unlockFlags<1>();
    }

    assert(!nodes.empty());
    std::size_t rootindex = nodes.size() - 1;

    std::stack<std::size_t> pending;
    /*
     * sim:
     * bits 0,1: 01X-value of 0-propagation
     * bits 2,3: 01X-value of 1-propagation
     * bit 4: support-flag
     */
    for (auto v = variables.begin(); v != variables.end(); ++v) {
        /* variable is not in the support of root -> size already set */
        if (result[v - variables.begin()] != -1) {
            continue;
        }

        int vi = (*v)->varIndex();

        /* compute propagation results and store them in sim (this also reverts
         * changes to sim done by previous propagations!) */
        for (auto p = nodes.begin(); p != nodes.end(); ++p) {
            EstimateRemainingNodes_Node& n = *p;

            /* variable */
            if (n._var >= 0) {
                if (n._var == vi) {
                    // 1,11,00 = 1C
                    n._sim = 0x1C;
                } else {
                    // 0,01,01 = 0x5
                    n._sim = 0x5;
                }
            }
            /* and */
            else {
                unsigned int s1 = nodes[n._p1]._sim;
                unsigned int s2 = nodes[n._p2]._sim;
                /* compute support flag */
                unsigned int supp = (s1 | s2) & 0x10;

                // negation is swap and inversion
                if (n._p1n) {
                    s1 = ~(((s1 & 0xAAAAAAAA) >> 1) | ((s1 & 0x55555555) << 1));
                }
                if (n._p2n) {
                    s2 = ~(((s2 & 0xAAAAAAAA) >> 1) | ((s2 & 0x55555555) << 1));
                }

                unsigned int s = s1 & s2;
                if (n._nand) {
                    s = ~(((s & 0xAAAAAAAA) >> 1) | ((s & 0x55555555) << 1));
                }

                /* merge with support flag */
                n._sim = (s & 0xF) | supp;

/* if var not in support, the prop values must be X */
#if 0
                if( ( n._sim & 0x10 ) == 0 )
                {
                    assert( ( n._sim & 0xF ) == 0x5 );
                }
#endif
            }
        }

        /* now the sims of all nodes contain
         * 1. the 0-propagation result at bits 0,1
         * 2. the 1-propagation result at bits 2,3
         * 3. the support-flag at bit 4
         */
        std::size_t nodesA      = 0;
        std::size_t nodesB      = 0;
        std::size_t nodesCommon = 0;

        /* count needed nodes for 0-prop.
         * use bit 5 to mark common nodes of 0-prop and 1-prop (0x20)
         * use bit 6 to mark processed nodes by 0-prop (0x40)
         * use bit 7 to mark processed nodes by 1-prop (0x80)
         */

        pending.push(rootindex);
        while (!(pending.empty())) {
            EstimateRemainingNodes_Node& n = nodes[pending.top()];
            pending.pop();

            if ((n._sim & 0x40) == 0) {
                n._sim = n._sim | 0x40;

                /* support node */
                if ((n._sim & 0x10) != 0) {
                    /* check if 0-prop result is X */
                    if ((n._sim & 0x3) == 0x1) {
                        ++nodesA;
                    }
                    /* constant node */
                    else {
                        continue;
                    }
                }
                /* non-support node */
                else {
                    /* node must not be marked as common, yet! */
                    assert((n._sim & 0x20) == 0);

                    /* non-support nodes cannot have a constant propagation result, i.e.
                     * their prop results must be X,X = 01,01 */
                    assert((n._sim & 0xF) == 0x5);

                    ++nodesCommon;

                    /* mark node as common */
                    n._sim = n._sim | 0x20;
                }

                /* push parents if node is "and" */
                if (n._var == -1) {
                    pending.push(n._p1);
                    pending.push(n._p2);
                }
            }
        }

        /* count needed nodes for 1-prop */
        pending.push(rootindex);
        while (!(pending.empty())) {
            EstimateRemainingNodes_Node& n = nodes[pending.top()];
            pending.pop();

            if ((n._sim & 0x80) == 0) {
                n._sim = n._sim | 0x80;

                /* support node */
                if ((n._sim & 0x10) != 0) {
                    /* check if 1-prop result is X */
                    if ((n._sim & 0xC) == 0x4) {
                        ++nodesB;
                    } else {
                        continue;
                    }
                }
                /* non-support node */
                else {
                    /* if node has not been marked as common by 0-prop... */
                    if ((n._sim & 0x20) == 0) {
                        /* non-support nodes cannot have a constant propagation result, i.e.
                         * their prop results must be X,X = 01,01 */
                        assert((n._sim & 0xF) == 0x5);

                        ++nodesCommon;

                        /* mark node as common */
                        n._sim = n._sim | 0x20;
                    } else {
                        /* already counted by 0-prop */
                        continue;
                    }
                }

                /* push parents if node is "and" */
                if (n._var == -1) {
                    pending.push(n._p1);
                    pending.push(n._p2);
                }
            }
        }

        std::size_t est = 1 + nodesA + nodesB + nodesCommon;

        result[v - variables.begin()] = est;
    }

    // double deltaTime = lrabs::cpuTime() - startTime;
    // std::cout << __func__ << " " << deltaTime << std::endl;

    return result;
}

#define CACHE_FLAG_SUPPORT (1 << 10)
#define CACHE_FLAG_COF0 (1 << 11)
#define CACHE_FLAG_COF1 (1 << 12)
#define CACHE_FLAG_POS_EX (1 << 13)
#define CACHE_FLAG_NEG_EX (1 << 14)

aigpp::EdgeRef
aigpp::Manager::quantifyDeep(const aigpp::EdgeRef& root, int varIndex, bool existential)
{
    if (root.isConstant()) {
        return root;
    }

    lockFlags<2 | CACHE_FLAG_SUPPORT | CACHE_FLAG_COF0 | CACHE_FLAG_COF1 | CACHE_FLAG_POS_EX | CACHE_FLAG_NEG_EX>();

#if 0
    {
        /* all flags of all nodes must be cleared */
        for( Node* n = _nodes; n != 0; n = n->next() )
        {
            assert( !( n->flag<2|CACHE_FLAG_SUPPORT|CACHE_FLAG_COF0|CACHE_FLAG_COF1|CACHE_FLAG_POS_EX|CACHE_FLAG_NEG_EX>() ) );
        }
    }
#endif

    assert(!automaticFRAIGing());

    // populate cache
    QuantifyDeepDataMap cache;

    std::stack<Node*> pending;
    pending.push(root.getInternal().node());
    while (!pending.empty()) {
        if (pending.top()->flag<2>()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            QuantifyDeepData d;

            if (pending.top()->varIndex() == varIndex) {
                pending.top()
                    ->setFlag<2 | CACHE_FLAG_SUPPORT | CACHE_FLAG_COF0 | CACHE_FLAG_COF1 | CACHE_FLAG_POS_EX
                              | CACHE_FLAG_NEG_EX>();

                d.cof0   = const0;
                d.cof1   = const1;
                d.pos_ex = const1;
                d.neg_ex = const1;
            } else {
                pending.top()->setFlag<2 | CACHE_FLAG_COF0 | CACHE_FLAG_COF1 | CACHE_FLAG_POS_EX | CACHE_FLAG_NEG_EX>();

                d.cof0   = InternalEdgeRef(pending.top());
                d.cof1   = d.cof0;
                d.pos_ex = d.cof0;
                d.neg_ex = !(d.pos_ex);
            }

            cache[pending.top()] = d;
            pending.pop();
        } else if (!(pending.top()->parent1().node()->flag<2>())) {
            pending.push(pending.top()->parent1().node());
        } else if (!(pending.top()->parent2().node()->flag<2>())) {
            pending.push(pending.top()->parent2().node());
        } else {
            QuantifyDeepData d;

            if (pending.top()->parent1().node()->flag<CACHE_FLAG_SUPPORT>()
                || pending.top()->parent2().node()->flag<CACHE_FLAG_SUPPORT>()) {
                pending.top()->setFlag<CACHE_FLAG_SUPPORT>();
            }

            cache[pending.top()] = d;

            pending.top()->setFlag<2>();
            pending.pop();
        }
    }  // end while

    InternalEdgeRef result;
    if (existential) {
        result = quantifyDeepRecursive(root.getInternal(), cache);
    } else {
        result = !quantifyDeepRecursive(!(root.getInternal()), cache);
    }

    for (QuantifyDeepDataMap::const_iterator p = cache.begin(); p != cache.end(); ++p) {
        p->first->unsetFlag<(2 | CACHE_FLAG_SUPPORT | CACHE_FLAG_COF0 | CACHE_FLAG_COF1 | CACHE_FLAG_POS_EX
                             | CACHE_FLAG_NEG_EX)>();
    }
    cache.clear();

    unlockFlags<2 | CACHE_FLAG_SUPPORT | CACHE_FLAG_COF0 | CACHE_FLAG_COF1 | CACHE_FLAG_POS_EX | CACHE_FLAG_NEG_EX>();

    return EdgeRef(_extRefTable, result);
}

aigpp::InternalEdgeRef
aigpp::Manager::quantifyDeepRecursive(const aigpp::Edge& root, aigpp::Manager::QuantifyDeepDataMap& cache)
{
    Node* rootN = root.node();
    if (!(rootN->flag<CACHE_FLAG_SUPPORT>())) {
        return InternalEdgeRef(root);
    }

    if (rootN->isVar()) {
        QuantifyDeepData& d = cache[rootN];
        if (root.isInverted()) {
            return d.neg_ex;
        } else {
            return d.pos_ex;
        }
    }

    /* ex x f */
    if (!(root.isInverted())) {
        if (rootN->flag<CACHE_FLAG_POS_EX>()) {
            return cache[rootN].pos_ex;
        }
        /* ex x f=a&b */
        else if (!(rootN->isNAND())) {
            /* x occurs in a */
            if (rootN->parent1().node()->flag<CACHE_FLAG_SUPPORT>()) {
                /* x occurs in b */
                if (rootN->parent2().node()->flag<CACHE_FLAG_SUPPORT>()) {
                    return quantifyDeepCofactor(root, cache);
                }
                /* x does not occur in b -> ex x a&b = ( ex x a ) & b */
                else {
                    InternalEdgeRef res
                        = quantifyDeepRecursive(rootN->parent1(), cache) & InternalEdgeRef(rootN->parent2());

                    rootN->setFlag<CACHE_FLAG_POS_EX>();
                    QuantifyDeepData& d2 = cache[rootN];
                    d2.pos_ex            = res;

                    return res;
                }
            }
            /* x does not occur in a -> ex x a&b = a & ( ex x b ) */
            else {
                assert(rootN->parent2().node()->flag<CACHE_FLAG_SUPPORT>());
                InternalEdgeRef res
                    = InternalEdgeRef(rootN->parent1()) & quantifyDeepRecursive(rootN->parent2(), cache);

                rootN->setFlag<CACHE_FLAG_POS_EX>();

                QuantifyDeepData& d2 = cache[rootN];
                d2.pos_ex            = res;

                return res;
            }
        }
        /* ex x f=!(a&b) = ( ex x !a ) | ( ex x !b ) */
        else {
            InternalEdgeRef res
                = quantifyDeepRecursive(!(rootN->parent1()), cache) | quantifyDeepRecursive(!(rootN->parent2()), cache);

            rootN->setFlag<CACHE_FLAG_POS_EX>();

            QuantifyDeepData& d2 = cache[rootN];
            d2.pos_ex            = res;

            return res;
        }
    }
    /* ex x !f */
    else {
        if (rootN->flag<CACHE_FLAG_NEG_EX>()) {
            return cache[rootN].neg_ex;
        }
        /* ex x !f=!(a&b) = ( ex x !a ) | ( ex x !b ) */
        else if (!(rootN->isNAND())) {
            InternalEdgeRef res
                = quantifyDeepRecursive(!(rootN->parent1()), cache) | quantifyDeepRecursive(!(rootN->parent2()), cache);

            rootN->setFlag<CACHE_FLAG_NEG_EX>();

            QuantifyDeepData& d2 = cache[rootN];
            d2.neg_ex            = res;

            return res;
        }
        /* ex x !f=!!(a&b)=a&b */
        else {
            /* x occurs in a */
            if (rootN->parent1().node()->flag<CACHE_FLAG_SUPPORT>()) {
                /* x occurs in b */
                if (rootN->parent2().node()->flag<CACHE_FLAG_SUPPORT>()) {
                    return quantifyDeepCofactor(root, cache);
                }
                /* x does not occur in b -> ex x !!a&b = ( ex x a ) & b */
                else {
                    InternalEdgeRef res
                        = quantifyDeepRecursive(rootN->parent1(), cache) & InternalEdgeRef(rootN->parent2());

                    rootN->setFlag<CACHE_FLAG_NEG_EX>();

                    QuantifyDeepData& d2 = cache[rootN];
                    d2.neg_ex            = res;

                    return res;
                }
            }
            /* x does not occur in a -> ex x !!a&b = a & ( ex x b ) */
            else {
                assert(rootN->parent2().node()->flag<CACHE_FLAG_SUPPORT>());
                InternalEdgeRef res
                    = InternalEdgeRef(rootN->parent1()) & quantifyDeepRecursive(rootN->parent2(), cache);

                rootN->setFlag<CACHE_FLAG_NEG_EX>();

                QuantifyDeepData& d2 = cache[rootN];
                d2.neg_ex            = res;

                return res;
            }
        }
    }
}

aigpp::InternalEdgeRef
aigpp::Manager::quantifyDeepCofactor(const aigpp::Edge& root, aigpp::Manager::QuantifyDeepDataMap& cache)
{
    Node* rootN = root.node();

    // ok
    if (!(rootN->flag<CACHE_FLAG_SUPPORT>())) {
        return InternalEdgeRef(root);
    }

    // ok
    if (root.isInverted()) {
        if (rootN->flag<CACHE_FLAG_NEG_EX>()) {
            return cache[rootN].neg_ex;
        }
    } else {
        if (rootN->flag<CACHE_FLAG_POS_EX>()) {
            return cache[rootN].pos_ex;
        }
    }

    std::stack<Node*> pending;

    /* compute negative cofactor */
    pending.push(rootN);
    while (!pending.empty()) {
        // ok
        if (pending.top()->flag<CACHE_FLAG_COF0>()) {
            pending.pop();
        } else if (!(pending.top()->flag<CACHE_FLAG_SUPPORT>())) {
            pending.top()->setFlag<CACHE_FLAG_COF0>();
            cache[pending.top()].cof0 = InternalEdgeRef(pending.top());
            pending.pop();
        }
        // ok
        else if (!(pending.top()->parent1().node()->flag<CACHE_FLAG_COF0>())) {
            pending.push(pending.top()->parent1().node());
        }
        // ok
        else if (!(pending.top()->parent2().node()->flag<CACHE_FLAG_COF0>())) {
            pending.push(pending.top()->parent2().node());
        }
        // ok
        else {
            aigpp::Manager::QuantifyDeepDataMap::const_iterator p1 = cache.find(pending.top()->parent1().node());
            assert(p1 != cache.end());
            aigpp::Manager::QuantifyDeepDataMap::const_iterator p2 = cache.find(pending.top()->parent2().node());
            assert(p2 != cache.end());

            InternalEdgeRef res = p1->second.cof0.notIf(pending.top()->parent1().isInverted())
                                  & p2->second.cof0.notIf(pending.top()->parent2().isInverted());

            pending.top()->setFlag<CACHE_FLAG_COF0>();
            cache[pending.top()].cof0 = res.notIf(pending.top()->isNAND());

            pending.pop();
        }
    }

    /* compute positive cofactor */
    pending.push(rootN);
    while (!pending.empty()) {
        // ok
        if (pending.top()->flag<CACHE_FLAG_COF1>()) {
            pending.pop();
        } else if (!(pending.top()->flag<CACHE_FLAG_SUPPORT>())) {
            pending.top()->setFlag<CACHE_FLAG_COF1>();
            cache[pending.top()].cof1 = InternalEdgeRef(pending.top());

            pending.pop();
        }
        // ok
        else if (!(pending.top()->parent1().node()->flag<CACHE_FLAG_COF1>())) {
            pending.push(pending.top()->parent1().node());
        }
        // ok
        else if (!(pending.top()->parent2().node()->flag<CACHE_FLAG_COF1>())) {
            pending.push(pending.top()->parent2().node());
        }
        // ok
        else {
            aigpp::Manager::QuantifyDeepDataMap::const_iterator p1 = cache.find(pending.top()->parent1().node());
            assert(p1 != cache.end());
            aigpp::Manager::QuantifyDeepDataMap::const_iterator p2 = cache.find(pending.top()->parent2().node());
            assert(p2 != cache.end());

            InternalEdgeRef res = p1->second.cof1.notIf(pending.top()->parent1().isInverted())
                                  & p2->second.cof1.notIf(pending.top()->parent2().isInverted());

            pending.top()->setFlag<CACHE_FLAG_COF1>();
            cache[pending.top()].cof1 = res.notIf(pending.top()->isNAND());

            pending.pop();
        }
    }

    assert(rootN->allFlags<(CACHE_FLAG_COF0 | CACHE_FLAG_COF1)>());
    auto c = cache.find(rootN);
    assert(c != cache.end());

    if (root.isInverted()) {
        rootN->setFlag<CACHE_FLAG_NEG_EX>();
        c->second.neg_ex = (!(c->second.cof0)) | (!(c->second.cof1));

        return c->second.neg_ex;
    } else {
        rootN->setFlag<CACHE_FLAG_POS_EX>();
        c->second.pos_ex = c->second.cof0 | c->second.cof1;

        return c->second.pos_ex;
    }
}
