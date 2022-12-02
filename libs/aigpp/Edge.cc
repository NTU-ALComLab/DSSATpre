/**************************************************************
 *
 *       AIGPP Package // Edge.cc
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

#include "Edge.hh"

#include <stack>

#include "Manager.hh"
#include "Node.hh"

bool
aigpp::Edge::functionallyEquivalent(const aigpp::Edge& n) const
{
    if (node() == n.node()) {
        return (isInverted() == n.isInverted());
    } else if (isConstant() && n.isConstant()) {
        return (_node == n._node);
    } else if (!isConstant() && n.isConstant()) {
        if (node()->isReduced()) {
            return (_node == n._node);
        }

        Manager* m = node()->manager();

        if (isInverted() == n.isInverted()) {
            if (m->simCreation()) {
                if (node()->sim().isZero()) {
                    return (m->satEqualZero(node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqualZero(node()));
            }
        } else {
            if (m->simCreation()) {
                if (node()->sim().isOne()) {
                    return (m->satEqualOne(node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqualOne(node()));
            }
        }
    } else if (isConstant() && !(n.isConstant())) {
        if (n.node()->isReduced()) {
            return (_node == n._node);
        }

        Manager* m = n.node()->manager();

        if (isInverted() == n.isInverted()) {
            if (m->simCreation()) {
                if (n.node()->sim().isZero()) {
                    return (m->satEqualZero(n.node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqualZero(n.node()));
            }
        } else {
            if (m->simCreation()) {
                if (n.node()->sim().isOne()) {
                    return (m->satEqualOne(n.node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqualOne(n.node()));
            }
        }
    } else {
        if (node()->isReduced() && n.node()->isReduced()) {
            return (_node == n._node);
        }

        Manager* m = node()->manager();

        if (isInverted() == n.isInverted()) {
            if (m->simCreation()) {
                if (node()->sim() == n.node()->sim()) {
                    return (m->satEqual(node(), n.node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqual(node(), n.node()));
            }
        } else {
            if (m->simCreation()) {
                if (node()->sim().equalInverted(n.node()->sim())) {
                    return (m->satEqualNegated(node(), n.node()));
                } else {
                    return false;
                }
            } else {
                return (m->satEqualNegated(node(), n.node()));
            }
        }
    }
}

std::size_t
aigpp::Edge::nodeCount() const
{
    if (isConstant()) {
        return 0;
    } else {
        return node()->nodeCount();
    }
}

int
aigpp::Edge::hasSuperNode(const aigpp::Edge& n) const
{
    std::stack<Edge> pending;
    pending.push(Edge(nullptr));
    Edge m = *this;

    do {
        if (m.node() == n.node()) {
            if (m.isInverted() == n.isInverted()) {
                return +1;
            } else {
                return -1;
            }
        } else if (m.node()->isVar() || (m.isInverted() != m.node()->isNAND())) {
            m = pending.top();
            pending.pop();
        } else {
            pending.push(m.node()->parent2());
            m = m.node()->parent1();
        }
    } while (!pending.empty());

    return 0;
}

unsigned long
aigpp::Edge::hash() const
{
    assert(!isConstant());
    return (node()->index() * 2 + (isInverted() ? 1 : 0));
}

std::ostream&
aigpp::operator<<(std::ostream& os, const Edge& n)
{
    if (!n.isConstant()) {
        if (n.isInverted()) {
            os << "~";
        }
        os << *(n.node());
    } else {
        if (n.isInverted()) {
            os << "T";
        } else {
            os << "F";
        }
    }

    return os;
}
