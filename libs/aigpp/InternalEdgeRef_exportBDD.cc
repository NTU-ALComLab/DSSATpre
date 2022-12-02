
#include "InternalEdgeRef.hh"

DdNode*
aigpp::InternalEdgeRef::exportBDD(DdManager* bddManager) const
{
    if (isConstant()) {
        if (*this == const0) {
            return Cudd_ReadLogicZero(bddManager);
        } else {
            return Cudd_Not(Cudd_ReadLogicZero(bddManager));
        }
    }

    if (node()->isVar()) {
        DdNode* v = Cudd_NotCond(Cudd_bddIthVar(bddManager, node()->varIndex()), isInverted());
        Cudd_Ref(v);
        return v;
    }

    node()->manager()->lockFlags<1>();

    std::stack<Node*>  pending;
    std::vector<Node*> processed;

    /* count occurences of each node in the cone */
    typedef std::map<Node*, int> CountMap;
    CountMap                     count;

    pending.push(node());
    while (!pending.empty()) {
        Node* n = pending.top();

        assert(!(n->isVar()));

        if (n->flag<1>()) {
            pending.pop();
            continue;
        }

        if (!(n->parent1().node()->isVar()) && !(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
            continue;
        }
        if (!(n->parent2().node()->isVar()) && !(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
            continue;
        }

        if (!(n->parent1().node()->isVar())) {
            ++count[n->parent1().node()];
        }

        if (!(n->parent2().node()->isVar())) {
            ++count[n->parent2().node()];
        }

        count[n] = 0;

        processed.push_back(n);
        n->setFlag<1>();
        pending.pop();
    }

    for (std::vector<Node*>::iterator p = processed.begin(); p != processed.end(); ++p) {
        (*p)->unsetFlag<1>();
    }
    processed.clear();

    /* create BDD */
    typedef std::map<Node*, DdNode*> ABMap;
    ABMap                            aig2bdd;

    pending.push(node());

    bool reachedLimit = false;

    while (!pending.empty()) {
        Node* n = pending.top();

        if (n->flag<1>()) {
            pending.pop();
            continue;
        }

        if (n->isVar()) {
            processed.push_back(n);
            n->setFlag<1>();

            aig2bdd[n] = Cudd_bddIthVar(bddManager, n->varIndex());
            Cudd_Ref(aig2bdd[n]);

            pending.pop();
            continue;
        }

        if (!(n->parent1().node()->flag<1>())) {
            pending.push(n->parent1().node());
            continue;
        }
        if (!(n->parent2().node()->flag<1>())) {
            pending.push(n->parent2().node());
            continue;
        }

        DdNode* p1;
        {
            ABMap::iterator    p = aig2bdd.find(n->parent1().node());
            CountMap::iterator c = count.find(n->parent1().node());
            assert(p != aig2bdd.end());
            assert(c != count.end());

            p1 = p->second;

            --(c->second);
            if (c->second == 0) {
                count.erase(c);
                aig2bdd.erase(p);
            } else {
                Cudd_Ref(p1);
            }
        }
        p1 = Cudd_NotCond(p1, n->parent1().isInverted());

        DdNode* p2;
        {
            ABMap::iterator    p = aig2bdd.find(n->parent2().node());
            CountMap::iterator c = count.find(n->parent2().node());

            p2 = p->second;

            --(c->second);
            if (c->second == 0) {
                count.erase(c);
                aig2bdd.erase(p);
            } else {
                Cudd_Ref(p2);
            }
        }
        p2 = Cudd_NotCond(p2, n->parent2().isInverted());

        DdNode* b = Cudd_bddAnd(_bddManager, p1, p2);

        if (b == 0) {
            Cudd_RecursiveDeref(_bddManager, p1);
            Cudd_RecursiveDeref(_bddManager, p2);
            reachedLimit = true;
            break;
        } else {
            Cudd_Ref(b);
            Cudd_RecursiveDeref(_bddManager, p1);
            Cudd_RecursiveDeref(_bddManager, p2);
        }

        b = Cudd_NotCond(b, n->isNAND());

        aig2bdd[n] = b;

        processed.push_back(n);
        n->setFlag<1>();
    }

    /* clear stack (stack may be non-empty if the bdd limit was reached) */
    while (!pending.empty()) pending.pop();

    /* get result */
    DdNode* bdd = 0;
    if (!reachedLimit) {
        bdd = Cudd_NotCond(aig2bdd[node()], isInverted());
        assert(!Cudd_IsConstant(bdd));
        Cudd_Ref(bdd);
    }

    /* clear cache */
    for (ABMap::iterator p = aig2bdd.begin(); p != aig2bdd.end(); ++p) {
        Cudd_RecursiveDeref(_bddManager, p->second);
    }
    aig2bdd.clear();

    for (std::vector<Node*>::iterator p = processed.begin(); p != processed.end(); ++p) {
        (*p)->unsetFlag<1>();
    }
    processed.clear();

    unlockFlags<1>();

    return bdd;
}
