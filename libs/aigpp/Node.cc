/**************************************************************
 *
 *       AIGPP Package // Node.cc
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
 *
 ***************************************************************/

#include "Node.hh"

/* std */
#include <set>

/* lrabs utilities */
#include <lrabsutil/Resources.hh>

/* local */
#include "Manager.hh"
#include "StatManager.hh"

aigpp::Node::Node(aigpp::Manager* manager, int var) :
    _manager(manager),
    _flags(0),
    _refCount(0),
    _index(-var),
    _parent1(nullptr),
    _parent2(nullptr),
    _sim(nullptr),
#ifdef USE_TEMPSIM
    _tempSim(0ul),
#endif
    _nEqualSimHash(nullptr),
    _pEqualSimHash(nullptr),
    _nEqualSim(nullptr),
    _pEqualSim(nullptr),
    _invertedSimClass(nullptr),
    _prevNode(nullptr),
    _nextNode(nullptr),
    _prevUniqueTableEntry(nullptr),
    _nextUniqueTableEntry(nullptr),
    _satVar(var_Undef)
{
    if (manager->simCreation()) {
        _sim = manager->createSim();
    }
}

aigpp::Node::Node(aigpp::Manager* manager, const aigpp::Edge& p1, const aigpp::Edge& p2) :
    _manager(manager),
    _flags(0),
    _refCount(0),
    _index(0),
    _parent1(p1),
    _parent2(p2),
    _sim(nullptr),
#ifdef USE_TEMPSIM
    _tempSim(0ul),
#endif
    _nEqualSimHash(nullptr),
    _pEqualSimHash(nullptr),
    _nEqualSim(nullptr),
    _pEqualSim(nullptr),
    _invertedSimClass(nullptr),

    _prevNode(nullptr),
    _nextNode(nullptr),
    _prevUniqueTableEntry(nullptr),
    _nextUniqueTableEntry(nullptr),
    _satVar(var_Undef)
{
    if (manager->simCreation()) {
        _sim = manager->createSim();
        updateSim();
    }

#ifdef USE_TEMPSIM
    updateTempSim();
#endif
}

aigpp::Node::~Node()
{
    /* delete _sim if non-zero */
    if (_sim != nullptr) {
        delete _sim;
        _sim = nullptr;
    }
}

bool
aigpp::Node::hasParent(aigpp::Node* n) const
{
    assert(this != n);

    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& flagged = manager()->_universalConstNodeVector;
    assert(flagged.empty());

    bool nIsParent = false;

    pending.push(this);
    while (!pending.empty()) {
        /* n found */
        if (pending.top() == n) {
            nIsParent = true;
            break;
        }

        /* already processed */
        else if (pending.top()->flag<1>()) {
            pending.pop();
        }

        /* parent 1 not processed yet */
        else if (!(pending.top()->isVar()) && !(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        }

        /* parent 2 not processed yet */
        else if (!(pending.top()->isVar()) && !(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        }

        /* var node, or parents processed */
        else {
            flagged.push_back(pending.top());
            pending.top()->setFlag<1>();
            pending.pop();
        }
    }

    /* clear marks of all processed nodes */
    for (const Node* n : flagged) {
        n->unsetFlag<1>();
    }
    flagged.clear();

    while (!pending.empty()) {
        pending.pop();
    }

    manager()->unlockFlags<1>();

    return (nIsParent);
}

aigpp::VarSupport
aigpp::Node::support() const
{
    VarSupport s;

    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    pending.push(this);
    while (!(pending.empty())) {
        const Node* n = pending.top();
        pending.pop();

        if (!(n->flag<1>())) {
            n->setFlag<1>();
            processed.push_back(n);

            if (n->isVar()) {
                s.addVar(n->varIndex());
            } else {
                pending.push(n->parent1().node());
                pending.push(n->parent2().node());
            }
        }
    }

    /* clear flags */
    for (const Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    manager()->unlockFlags<1>();

    return s;
}

std::vector<aigpp::Node*>
aigpp::Node::supportNodes() const
{
    VarSupport supp = support();
    assert(supp.vars() > 0);

    std::vector<Node*> nodes;
    for (std::size_t v = 0; v != manager()->variableCount(); ++v) {
        if (supp.hasVar(v)) {
            nodes.push_back(manager()->_variables[v]);
        }
    }

    return nodes;
}

bool
aigpp::Node::varInSupport(int v) const
{
    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    bool found = false;

    pending.push(this);
    while (!(pending.empty())) {
        const Node* n = pending.top();
        pending.pop();

        if (!(n->flag<1>())) {
            n->setFlag<1>();
            processed.push_back(n);

            if (n->isVar()) {
                if (n->varIndex() == v) {
                    found = true;
                    break;
                }
            } else {
                pending.push(n->parent1().node());
                pending.push(n->parent2().node());
            }
        }
    }

    /* clear stack (may be non-empty if variable was found ) */
    while (!pending.empty()) {
        pending.pop();
    }

    /* clear flags */
    for (const Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    manager()->unlockFlags<1>();

    return found;
}

std::size_t
aigpp::Node::remainingNodesAfterQuantification(int var) const
{
    if (!varInSupport(var)) {
        return 0;
    }

    manager()->lockFlags<1 | 2>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    // propagate const0
    pending.push(this);
    while (!pending.empty()) {
        if (pending.top()->flag<1 | 2>()) {
            pending.pop();
            continue;
        }

        if (pending.top()->isVar()) {
            if (pending.top()->varIndex() != var) {
                pending.top()->_prop0.type = 2;
                pending.top()->_prop0.p1   = nullptr;

                pending.top()->_prop1.type = 2;
                pending.top()->_prop1.p1   = nullptr;
            } else {
                pending.top()->_prop0.type = 0;

                pending.top()->_prop1.type = 1;
            }

            pending.top()->setFlag<1 | 2>();
            processed.push_back(pending.top());

            pending.pop();
            continue;
        }

        if (!(pending.top()->parent1().node()->flag<1 | 2>())) {
            pending.push(pending.top()->parent1().node());
            continue;
        }
        if (!(pending.top()->parent2().node()->flag<1 | 2>())) {
            pending.push(pending.top()->parent2().node());
            continue;
        }

        if ((pending.top()->parent1().node()->_prop0.type == 0 && !(pending.top()->parent1().isInverted()))
            || (pending.top()->parent1().node()->_prop0.type == 1 && (pending.top()->parent1().isInverted()))
            || (pending.top()->parent2().node()->_prop0.type == 0 && !(pending.top()->parent2().isInverted()))
            || (pending.top()->parent2().node()->_prop0.type == 1
                && (pending.top()->parent2().isInverted())))  // 0 NAND ?, 0 AND ?
        {
            pending.top()->_prop0.type = (pending.top()->isNAND() ? 1 : 0);
        } else if ((pending.top()->parent1().node()->_prop0.type == 1 && !(pending.top()->parent1().isInverted()))
                   || (pending.top()->parent1().node()->_prop0.type == 0
                       && (pending.top()->parent1().isInverted())))  // 1 NAND ?, 1 AND ?
        {
            if ((pending.top()->parent2().node()->_prop0.type == 1 && !(pending.top()->parent2().isInverted()))
                || (pending.top()->parent2().node()->_prop0.type == 0
                    && (pending.top()->parent2().isInverted())))  // 1 NAND 1, 1 AND 1
            {
                pending.top()->_prop0.type = (pending.top()->isNAND() ? 0 : 1);
            } else {
                pending.top()->_prop0.type = 3;
                pending.top()->_prop0.p1   = pending.top()->parent2().node();
            }
        } else if ((pending.top()->parent2().node()->_prop0.type == 1 && !(pending.top()->parent2().isInverted()))
                   || (pending.top()->parent2().node()->_prop0.type == 0
                       && (pending.top()->parent2().isInverted())))  // ? NAND 1, ? AND 1
        {
            pending.top()->_prop0.type = 3;
            pending.top()->_prop0.p1   = pending.top()->parent1().node();
        } else  // ? NAND ?, ? AND ?
        {
            pending.top()->_prop0.type = 2;
            pending.top()->_prop0.p1   = pending.top()->parent1().node();
            pending.top()->_prop0.p2   = pending.top()->parent2().node();
        }

        if ((pending.top()->parent1().node()->_prop1.type == 0 && !(pending.top()->parent1().isInverted()))
            || (pending.top()->parent1().node()->_prop1.type == 1 && (pending.top()->parent1().isInverted()))
            || (pending.top()->parent2().node()->_prop1.type == 0 && !(pending.top()->parent2().isInverted()))
            || (pending.top()->parent2().node()->_prop1.type == 1
                && (pending.top()->parent2().isInverted())))  // 0 NAND ?, 0 AND ?
        {
            pending.top()->_prop1.type = (pending.top()->isNAND() ? 1 : 0);
        } else if ((pending.top()->parent1().node()->_prop1.type == 1 && !(pending.top()->parent1().isInverted()))
                   || (pending.top()->parent1().node()->_prop1.type == 0
                       && (pending.top()->parent1().isInverted())))  // 1 NAND ?, 1 AND ?
        {
            if ((pending.top()->parent2().node()->_prop1.type == 1 && !(pending.top()->parent2().isInverted()))
                || (pending.top()->parent2().node()->_prop1.type == 0
                    && (pending.top()->parent2().isInverted())))  // 1 NAND 1, 1 AND 1
            {
                pending.top()->_prop1.type = (pending.top()->isNAND() ? 0 : 1);
            } else {
                pending.top()->_prop1.type = 3;
                pending.top()->_prop1.p1   = pending.top()->parent2().node();
            }
        } else if ((pending.top()->parent2().node()->_prop1.type == 1 && !(pending.top()->parent2().isInverted()))
                   || (pending.top()->parent2().node()->_prop1.type == 0
                       && (pending.top()->parent2().isInverted())))  // ? NAND 1, ? AND 1
        {
            pending.top()->_prop1.type = 3;
            pending.top()->_prop1.p1   = pending.top()->parent1().node();
        } else  // ? NAND ?, ? AND ?
        {
            pending.top()->_prop1.type = 2;
            pending.top()->_prop1.p1   = pending.top()->parent1().node();
            pending.top()->_prop1.p2   = pending.top()->parent2().node();
        }

        pending.top()->setFlag<1 | 2>();
        processed.push_back(pending.top());

        pending.pop();
    }

    const Node* pt;

    // calculate cone0
    std::set<const Node*> cone0;

    pending.push(this);
    while (!pending.empty()) {
        pt = pending.top();
        pending.pop();

        if (!(pt->flag<1>())) {
            continue;
        }

        if (pt->_prop0.type <= 1)  // constant
        {
        } else if (pt->_prop0.type == 2 && pt->_prop0.p1 == nullptr)  // variable
        {
            cone0.insert(pt);
        } else if (pt->_prop0.type == 3)  // forward
        {
            pending.push(pt->_prop0.p1);
        } else  // and
        {
            cone0.insert(pt);

            pending.push(pt->_prop0.p1);
            pending.push(pt->_prop0.p2);
        }

        pt->unsetFlag<1>();
    }

    // calculate cone1
    std::set<const Node*> cone1;

    pending.push(this);
    while (!pending.empty()) {
        pt = pending.top();
        pending.pop();

        if (!(pt->flag<2>())) {
            continue;
        }

        if (pt->_prop1.type <= 1)  // constant
        {
        } else if (pt->_prop1.type == 2 && pt->_prop1.p1 == nullptr)  // variable
        {
            cone1.insert(pt);
        } else if (pt->_prop1.type == 3)  // forward
        {
            pending.push(pt->_prop1.p1);
        } else  // and
        {
            cone1.insert(pt);

            pending.push(pt->_prop1.p1);
            pending.push(pt->_prop1.p2);
        }

        pt->unsetFlag<2>();
    }

    for (const Node* p : processed) {
        p->unsetFlag<1 | 2>();
    }
    processed.clear();

    // calculate common nodes depending on var
    std::size_t common = 0;

    pending.push(this);
    while (!pending.empty()) {
        if (pending.top()->flag<1 | 2>()) {
            pending.pop();
            continue;
        } else if (pending.top()->isVar()) {
            processed.push_back(pending.top());

            if (pending.top()->varIndex() == var) {
                pending.top()->setFlag<1>();
            } else {
                pending.top()->setFlag<2>();
            }
            pending.pop();
            continue;

        } else if (!(pending.top()->parent1().node()->flag<1 | 2>())) {
            pending.push(pending.top()->parent1().node());
            continue;
        } else if (!(pending.top()->parent2().node()->flag<1 | 2>())) {
            pending.push(pending.top()->parent2().node());
            continue;
        } else if (pending.top()->parent1().node()->flag<1>() || pending.top()->parent2().node()->flag<1>()) {
            processed.push_back(pending.top());
            pending.top()->setFlag<1>();

            if (cone0.find(pending.top()) != cone0.end() && cone1.find(pending.top()) != cone1.end()) {
                ++common;
            }
            pending.pop();
            continue;
        } else {
            processed.push_back(pending.top());
            pending.top()->setFlag<2>();
            pending.pop();
            continue;
        }
    }
    for (const Node* p : processed) {
        p->unsetFlag<1 | 2>();
    }
    processed.clear();

    manager()->unlockFlags<1 | 2>();

    //  std::cout << "cone0=" << cone0.size() << " cone1=" << cone1.size() << "
    //  common=" << common << std::endl;

    return cone0.size() + cone1.size() - common;
}

std::size_t
aigpp::Node::nodeCount() const
{
    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& parents = manager()->_universalConstNodeVector;
    assert(parents.empty());

    pending.push(this);

    while (!pending.empty()) {
        /* node was counted/processed already */
        if (pending.top()->flag<1>()) {
            pending.pop();
        }
        /* node is an "and" node and its first parent wasn't processed yet */
        else if (!(pending.top()->isVar()) && !(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        }
        /* node is an "and" node and its second parent wasn't processed yet */
        else if (!(pending.top()->isVar()) && !(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        }
        /* node is an "and" node with both parents processed, or a variable node */
        else {
            parents.push_back(pending.top());
            pending.top()->setFlag<1>();

            pending.pop();
        }
    }

    /* clear marks of all processed nodes */
    for (const Node* n : parents) {
        n->unsetFlag<1>();
    }

    const std::size_t count = parents.size();
    parents.clear();

    manager()->unlockFlags<1>();

    return count;
}

std::set<aigpp::Node*>
aigpp::Node::cone()
{
    manager()->lockFlags<1>();

    std::stack<Node*>& pending = manager()->_universalNodeStack;
    assert(pending.empty());
    std::vector<Node*>& processed = manager()->_universalNodeVector;
    assert(processed.empty());

    std::set<aigpp::Node*> parents;

    pending.push(this);

    while (!pending.empty()) {
        if (pending.top()->flag<1>()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            parents.insert(pending.top());
            processed.push_back(pending.top());
            pending.top()->setFlag<1>();

            pending.pop();
        } else if (!(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        } else if (!(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        } else {
            parents.insert(pending.top());
            processed.push_back(pending.top());
            pending.top()->setFlag<1>();

            pending.pop();
        }
    }

    for (Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    manager()->unlockFlags<1>();

    return parents;
}

std::vector<aigpp::Node*>
aigpp::Node::coneVector()
{
    manager()->lockFlags<1>();

    std::stack<Node*>& pending = manager()->_universalNodeStack;
    assert(pending.empty());

    std::vector<Node*> processed;

    pending.push(this);

    while (!pending.empty()) {
        if (pending.top()->flag<1>()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            processed.push_back(pending.top());
            pending.top()->setFlag<1>();

            pending.pop();
        } else if (!(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        } else if (!(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        } else {
            processed.push_back(pending.top());
            pending.top()->setFlag<1>();

            pending.pop();
        }
    }

    for (Node* p : processed) {
        p->unsetFlag<1>();
    }

    manager()->unlockFlags<1>();

    return processed;
}

std::map<int, std::size_t>
aigpp::Node::maxVarDepth() const
{
    // TODO: get rid of "-1" in tempDepth
    VarSupport sup = support();

    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    std::map<int, std::size_t> depth;
    std::map<const Node*, int> tempDepth;

    for (std::size_t v = 0; v < manager()->variableCount(); ++v) {
        if (!(sup.hasVar(v))) {
            continue;
        }

        tempDepth.clear();

        pending.push(this);
        while (!pending.empty()) {
            if (pending.top()->flag<1>()) {
                pending.pop();
            } else if (pending.top()->isVar()) {
                if ((std::size_t)pending.top()->varIndex() == v) {
                    tempDepth[pending.top()] = 0;
                } else {
                    tempDepth[pending.top()] = -1;
                }

                pending.top()->setFlag<1>();
                processed.push_back(pending.top());
                pending.pop();
            } else if (!(pending.top()->parent1().node()->flag<1>())) {
                pending.push(pending.top()->parent1().node());
            } else if (!(pending.top()->parent2().node()->flag<1>())) {
                pending.push(pending.top()->parent2().node());
            } else {
                int d1 = tempDepth[pending.top()->parent1().node()];
                int d2 = tempDepth[pending.top()->parent2().node()];

                if (d2 > d1) {
                    d1 = d2;
                }

                if (d1 != -1) {
                    ++d1;
                }

                tempDepth[pending.top()] = d1;
                pending.top()->setFlag<1>();
                processed.push_back(pending.top());
                pending.pop();
            }
        }

        int d = tempDepth[this];
        assert(d >= 0);
        depth[static_cast<int>(v)] = static_cast<std::size_t>(d);

        for (const Node* p : processed) {
            p->unsetFlag<1>();
        }
        processed.clear();
    }

    manager()->unlockFlags<1>();
    return depth;
}

std::vector<int>
aigpp::Node::minVarDepth2() const
{
    manager()->lockFlags<1>();

    std::map<const Node*, int>      depth;
    std::vector<int>                result;
    static std::deque<const Node*>  pending;
    static std::vector<const Node*> marked;
    marked.clear();

    depth[this] = 0;
    this->setFlag<1>();
    marked.push_back(this);
    if (this->isVar()) {
        result.push_back(this->varIndex());
    }
    pending.push_back(this);

    while (!pending.empty()) {
        const Node* n = pending.front();
        pending.pop_front();

        if (!(n->isVar())) {
            int d = depth[n];

            const Node* n1 = n->parent1().node();
            const Node* n2 = n->parent2().node();

            if (!(n1->flag<1>())) {
                n1->setFlag<1>();
                depth[n1] = d + 1;
                marked.push_back(n1);
                if (n1->isVar()) {
                    result.push_back(n1->varIndex());
                }

                pending.push_back(n1);
            }

            if (!(n2->flag<1>())) {
                n2->setFlag<1>();
                depth[n2] = d + 1;
                marked.push_back(n2);
                if (n2->isVar()) {
                    result.push_back(n2->varIndex());
                }

                pending.push_back(n2);
            }
        }
    }

    for (const Node* n : marked) {
        n->unsetFlag<1>();
    }

    manager()->unlockFlags<1>();

    return result;
}

std::map<int, std::size_t>
aigpp::Node::minVarDepth() const
{
    // dijkstra's shortest path algorithm

    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<std::pair<const Node*, int>> Q;

    pending.push(this);
    while (!pending.empty()) {
        if (pending.top()->flag<1>()) {
            pending.pop();
        } else {
            const Node* n = pending.top();
            pending.pop();

            n->setFlag<1>();
            Q.emplace_back(n, -1);

            if (!(n->isVar())) {
                if (!(n->parent1().node()->flag<1>())) {
                    pending.push(n->parent1().node());
                }
                if (!(n->parent2().node()->flag<1>())) {
                    pending.push(n->parent2().node());
                }
            }
        }
    }

    for (auto& p : Q) {
        p.first->unsetFlag<1>();
    }

    for (auto p = Q.begin(); p != Q.end(); ++p) {
        if (p->first == this) {
            p->second = 0;

            while ((p != Q.begin()) && (((p - 1)->second > p->second) || (p - 1)->second == -1)) {
                std::swap(*p, *(p - 1));
                --p;
            }

            break;
        }
    }

    std::map<int, std::size_t> depth;

    for (auto current = Q.begin(); current != Q.end(); ++current) {
        if (!(current->first->isVar())) {
            const Node* child = current->first->parent1().node();

            for (auto p = current + 1; p != Q.end(); ++p) {
                if (p->first == child) {
                    if ((p->second == -1) || (p->second > (current->second + 1))) {
                        p->second = current->second + 1;
                        while ((p != current) && (((p - 1)->second > p->second) || (p - 1)->second == -1)) {
                            std::swap(*p, *(p - 1));
                            --p;
                        }
                    }
                    break;
                }
            }

            child = current->first->parent2().node();
            for (auto p = current + 1; p != Q.end(); ++p) {
                if (p->first == child) {
                    if ((p->second == -1) || (p->second > (current->second + 1))) {
                        p->second = current->second + 1;
                        while ((p != current) && (((p - 1)->second > p->second) || (p - 1)->second == -1)) {
                            std::swap(*p, *(p - 1));
                            --p;
                        }
                    }
                    break;
                }
            }
        } else {
            depth[current->first->varIndex()] = (std::size_t)current->second;
        }
    }

    manager()->unlockFlags<1>();

    return depth;
}

int
aigpp::Node::maxDepth(int var) const
{
    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    std::map<const Node*, int> tempDepth;

    pending.push(this);
    while (!pending.empty()) {
        if (pending.top()->flag<1>()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            if (pending.top()->varIndex() == var) {
                tempDepth[pending.top()] = 0;
            } else {
                tempDepth[pending.top()] = -1;
            }

            pending.top()->setFlag<1>();
            processed.push_back(pending.top());
            pending.pop();
        } else if (!(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        } else if (!(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        } else {
            int d1 = tempDepth[pending.top()->parent1().node()];
            int d2 = tempDepth[pending.top()->parent2().node()];

            if (d2 > d1) {
                d1 = d2;
            }

            if (d1 != -1) {
                ++d1;
            }

            tempDepth[pending.top()] = d1;
            pending.top()->setFlag<1>();
            processed.push_back(pending.top());
            pending.pop();
        }
    }

    for (const Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    manager()->unlockFlags<1>();

    return tempDepth[this];
}

int
aigpp::Node::minDepth(int var) const
{
    manager()->lockFlags<1>();

    std::stack<const Node*>& pending = manager()->_universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& processed = manager()->_universalConstNodeVector;
    assert(processed.empty());

    std::map<const Node*, int> tempDepth;

    pending.push(this);
    while (!pending.empty()) {
        if (pending.top()->flag<1>()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            if (pending.top()->varIndex() == var) {
                tempDepth[pending.top()] = 0;
            } else {
                tempDepth[pending.top()] = -1;
            }

            pending.top()->setFlag<1>();
            processed.push_back(pending.top());
            pending.pop();
        } else if (!(pending.top()->parent1().node()->flag<1>())) {
            pending.push(pending.top()->parent1().node());
        } else if (!(pending.top()->parent2().node()->flag<1>())) {
            pending.push(pending.top()->parent2().node());
        } else {
            int d1 = tempDepth[pending.top()->parent1().node()];
            int d2 = tempDepth[pending.top()->parent2().node()];

            if ((d2 != -1) && (d2 < d1 || d1 == -1)) {
                d1 = d2;
            }

            if (d1 != -1) {
                ++d1;
            }

            tempDepth[pending.top()] = d1;
            pending.top()->setFlag<1>();
            processed.push_back(pending.top());
            pending.pop();
        }
    }

    for (const Node* p : processed) {
        p->unsetFlag<1>();
    }
    processed.clear();

    manager()->unlockFlags<1>();

    return tempDepth[this];
}

bool
aigpp::Node::isMuxType() const
{
    if (isVar()) {
        return false;
    }

    Node* p1 = parent1().node();
    Node* p2 = parent2().node();

    if (p1->isVar() || p2->isVar()) {
        return false;
    }

    if (parent1().isInverted() == p1->isNAND()) {
        return false;
    }
    if (parent2().isInverted() == p2->isNAND()) {
        return false;
    }

    return (structurallyEquivalent(p1->parent1(), !(p2->parent1())) && p1->parent2().node() != p2->parent2().node())
           || (structurallyEquivalent(p1->parent1(), !(p2->parent2())) && p1->parent2().node() != p2->parent1().node())
           || (structurallyEquivalent(p1->parent2(), !(p2->parent1())) && p1->parent1().node() != p2->parent2().node())
           || (structurallyEquivalent(p1->parent2(), !(p2->parent2())) && p1->parent1().node() != p2->parent1().node());
}

bool
aigpp::Node::isXorType() const
{
    if (isVar()) {
        return false;
    }

    Node* p1 = parent1().node();
    Node* p2 = parent2().node();

    if (p1->isVar() || p2->isVar()) {
        return false;
    }

    if (parent1().isInverted() == p1->isNAND()) {
        return false;
    }
    if (parent2().isInverted() == p2->isNAND()) {
        return false;
    }

    return ((structurallyEquivalent(p1->parent1(), !(p2->parent1()))
             && structurallyEquivalent(p1->parent2(), !(p2->parent2())))
            || (structurallyEquivalent(p1->parent1(), !(p2->parent2()))
                && structurallyEquivalent(p1->parent2(), !(p2->parent1()))));
}

bool
aigpp::Node::isMuxOrXorType() const
{
    if (isVar()) {
        return false;
    }

    Node* p1 = parent1().node();
    Node* p2 = parent2().node();

    if (p1->isVar() || p2->isVar()) {
        return false;
    }
    if (parent1().isInverted() == p1->isNAND()) {
        return false;
    }
    if (parent2().isInverted() == p2->isNAND()) {
        return false;
    }

    return (structurallyEquivalent(p1->parent1(), !(p2->parent1()))
            || structurallyEquivalent(p1->parent1(), !(p2->parent2()))
            || structurallyEquivalent(p1->parent2(), !(p2->parent1()))
            || structurallyEquivalent(p1->parent2(), !(p2->parent2())));
}

std::vector<aigpp::Node*>
aigpp::Node::findRedundantVars(bool remove)
{
    VarSupport  s     = support();
    std::size_t count = s.vars();

    if (manager()->verbosity() >= 3) {
        std::cout << "checking redundance of " << count << " variables" << std::endl;
    }

    StatManager::instance().incRedundancyChecksCombined();
    StatManager::instance().incRedundancyChecksSingle(count);

    std::vector<Node*> candidates;

    for (std::vector<Node*>::const_iterator v = _manager->_variables.begin(); v != _manager->_variables.end(); ++v) {
        if (s.hasVar((*v)->varIndex())) {
            candidates.push_back(*v);
        }
    }

    double startTime = lrabs::cpuTime();

    std::vector<Node*> redundant = _manager->boolRedundant(this, candidates);

    StatManager::instance().incRedundancyCheckTime(lrabs::cpuTime() - startTime);

    if (manager()->verbosity() >= 3) {
        std::cout << "found " << redundant.size() << " redundant variables" << std::endl;
    }
    StatManager::instance().incRedundanciesDetected(redundant.size());

    if (remove && !(redundant.empty())) {
        bool saveAutoFRAIGing        = _manager->_automaticFRAIGing;
        _manager->_automaticFRAIGing = true;

        int saveVerbosity = _manager->verbosity();
        _manager->setVerbosity(0);

        std::vector<InternalEdgeRef> redundantVars;
        redundantVars.reserve(redundant.size());
        for (Node* n : redundant) {
            redundantVars.emplace_back(n);
        }

        if (saveVerbosity >= 3) {
            std::cout << "trying to eliminate" << std::endl;
            std::cout << "nodecount before: " << nodeCount() << std::endl;
        }
        InternalEdgeRef n(this);
        n.cofactor(redundantVars);
        if (saveVerbosity >= 3) {
            std::cout << "nodecount after: " << nodeCount() << std::endl;
        }

        _manager->setVerbosity(saveVerbosity);

        _manager->_automaticFRAIGing = saveAutoFRAIGing;
    }

    return redundant;
}

void
aigpp::Node::set(aigpp::Manager* manager, const aigpp::Edge& p1, const aigpp::Edge& p2)
{
    _manager  = manager;
    _flags    = 0;
    _refCount = 0;
    _index    = 0;
    _parent1  = p1;
    _parent2  = p2;

    if (manager->simCreation()) {
        updateSim();
    }

#ifdef USE_TEMPSIM
    updateTempSim();
#endif
    _refCount      = 0;
    _nEqualSimHash = nullptr;
    _nEqualSim     = nullptr;
    _prevNode      = nullptr;
    _nextNode      = nullptr;
    _satVar        = var_Undef;
}

bool
aigpp::Node::isMultiAnd(std::vector<aigpp::Edge>& inputs) const
{
    if (isVar()) {
        return false;
    } else if ((parent1().isInverted() != parent1().node()->isNAND() || parent1().node()->isVar())
               && (parent2().isInverted() != parent2().node()->isNAND() || parent2().node()->isVar())) {
        return false;
    }

    std::set<Edge, structurallyLess> inputstemp;

    std::stack<const Node*> pending;
    pending.push(this);

    while (!pending.empty()) {
        const Node* n = pending.top();
        pending.pop();

        assert(!(n->isVar()));

        if (n->parent1().node()->isVar() || (n->parent1().node()->isNAND() != n->parent1().isInverted())) {
            inputstemp.insert(n->parent1());
        } else {
            pending.push(n->parent1().node());
        }

        if (n->parent2().node()->isVar() || (n->parent2().node()->isNAND() != n->parent2().isInverted())) {
            inputstemp.insert(n->parent2());
        } else {
            pending.push(n->parent2().node());
        }
    }

    assert(inputstemp.size() >= 2);
    inputs.assign(inputstemp.begin(), inputstemp.end());

    return true;
}

std::ostream&
aigpp::operator<<(std::ostream& os, const aigpp::Node& n)
{
    if (n.isVar()) {
        std::string s = n._manager->variableName(n.varIndex());
        if (!s.empty()) {
            os << s;
        } else {
            os << "v" << n.varIndex();
        }
    } else {
        if (n.isNAND()) {
            os << "~";
        }
        os << "(" << n.parent1() << " & " << n.parent2() << ")";
    }

    return os;
}
