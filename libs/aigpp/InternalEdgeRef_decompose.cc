/**************************************************************
 *
 *       AIGPP // InternalEdgeRef_decompose.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 225 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "InternalEdgeRef.hh"

// std
#include <set>
#include <stack>

// aigpp
#include "Node.hh"

std::vector<aigpp::InternalEdgeRef>
aigpp::InternalEdgeRef::decomposeConjunction() const
{
    if (isConstant()) {
        std::vector<InternalEdgeRef> result;
        result.push_back(*this);
        return result;
    }

    std::set<Edge, structurallyLess> factors;

    std::stack<Edge> pending;
    pending.push(*this);

    while (!(pending.empty())) {
        if (pending.top().node()->isVar()) {
            factors.insert(pending.top());
            pending.pop();
        } else if (pending.top().isInverted() != pending.top().node()->isNAND()) {
            factors.insert(pending.top());
            pending.pop();
        } else {
            Node* n = pending.top().node();
            pending.pop();

            pending.push(n->parent1());
            pending.push(n->parent2());
        }
    }

    std::vector<InternalEdgeRef> result;
    result.reserve(factors.size());
    for (const auto& p : factors) {
        result.emplace_back(p);
    }

    return result;
}

std::vector<aigpp::InternalEdgeRef>
aigpp::InternalEdgeRef::decomposeDisjunction() const
{
    std::vector<InternalEdgeRef> result = (!*this).decomposeConjunction();
    for (auto& p : result) {
        p = !p;
    }
    return result;
}
