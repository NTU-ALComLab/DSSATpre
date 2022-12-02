/**************************************************************
 *
 *       AIGPP Package // InternalEdgeRef_polarity.cc
 *
 *	Author:
 *         Florian Pigorsch
 *	  University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *	Last revision:
 *         $Revision: 205 $
 *	  $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "InternalEdgeRef.hh"

#include "Manager.hh"

void
aigpp::InternalEdgeRef::recursivelySelectBothPolarities(aigpp::Node* root, std::vector<aigpp::Node*>& processed,
                                                        std::stack<aigpp::Node*>& pending)
{
    std::stack<Node*>::size_type initialStackSize = pending.size();

    pending.push(root);
    while (pending.size() > initialStackSize) {
        Node* n = pending.top();
        pending.pop();

        if (!(n->allFlags<(1 | 2)>())) {
            /* n has no polarities yet */
            if (!(n->flag<(1 | 2)>())) {
                processed.push_back(n);
            }
            n->setFlag<(1 | 2)>();

            if (!(n->isVar())) {
                pending.push(n->parent1().node());
                pending.push(n->parent2().node());
            }
        }
    }
}

std::vector<std::pair<aigpp::Node*, bool>>
aigpp::InternalEdgeRef::getVariablePolarities()
{
    if (isConstant()) {
        return std::vector<std::pair<Node*, bool>>();
    }

    /* calculate polarities of all variables occuring in the cone */

    Manager* m = node()->manager();
    m->lockFlags<(1 | 2)>();

    std::stack<Node*>& pending = m->_universalNodeStack;
    assert(pending.empty());

    std::vector<Node*>& processed = m->_universalNodeVector;
    assert(processed.empty());

    if (isInverted()) {
        node()->setFlag<1>();
    } else {
        node()->setFlag<2>();
    }
    processed.push_back(node());

    pending.push(node());

    while (!pending.empty()) {
        Node* n = pending.top();
        pending.pop();

        /* n already has both polarities */
        if (n->allFlags<(1 | 2)>()) {
            /* nothing to do */
        }
        /* n is a terminal node */
        else if (n->isVar()) {
            /* nothing to do */
        }
        /* n is an and node */
        else {
            Node* p = n->parent1().node();

            /* first parent does not have both polarities */
            if (!(p->allFlags<(1 | 2)>())) {
                int ptarget = 0;
                if (n->flag<1>()) {
                    ptarget = -1;
                } else if (n->flag<2>()) {
                    ptarget = +1;
                }
                assert(ptarget != 0);

                if (n->isNAND() != n->parent1().isInverted()) {
                    ptarget *= -1;
                }

                /* p has no polarities set -> first visit */
                if (!(p->flag<(1 | 2)>())) {
                    processed.push_back(p);
                    if (ptarget == -1) {
                        p->setFlag<1>();
                    } else if (ptarget == +1) {
                        p->setFlag<2>();
                    }
                    pending.push(p);
                }
                /* p already has target polarity */
                else if ((ptarget == -1 && p->flag<1>()) || (ptarget == +1 && p->flag<2>())) {
                    /* nothing to do */
                }
                /* p has different polarity than the target polarity */
                else {
                    /* set both polarites in the whole cone */
                    recursivelySelectBothPolarities(p, processed, pending);
                }
            }

            p = n->parent2().node();

            /* second parent does not have both polarities */
            if (!(p->allFlags<(1 | 2)>())) {
                int ptarget = 0;
                if (n->flag<1>()) {
                    ptarget = -1;
                } else if (n->flag<2>()) {
                    ptarget = +1;
                }
                assert(ptarget != 0);

                if (n->isNAND() != n->parent2().isInverted()) {
                    ptarget *= -1;
                }

                /* p has no polarities set -> first visit */
                if (!(p->flag<(1 | 2)>())) {
                    processed.push_back(p);
                    if (ptarget == -1) {
                        p->setFlag<1>();
                    } else if (ptarget == +1) {
                        p->setFlag<2>();
                    }
                    pending.push(p);
                }
                /* p already has target polarity */
                else if ((ptarget == -1 && p->flag<1>()) || (ptarget == +1 && p->flag<2>())) {
                    /* nothing to do */
                }
                /* p has different polarity than the target polarity */
                else {
                    /* set both polarites in the whole cone */
                    recursivelySelectBothPolarities(p, processed, pending);
                }
            }
        }
    }

    std::vector<std::pair<Node*, bool>> polarities;

    for (std::vector<Node*>::const_iterator p = processed.begin(); p != processed.end(); ++p) {
        if ((*p)->isVar()) {
            if ((*p)->flag<1>()) {
                polarities.emplace_back(*p, false);
            }

            if ((*p)->flag<2>()) {
                polarities.emplace_back(*p, true);
            }
        }

        (*p)->unsetFlag<(1 | 2)>();
    }

    processed.clear();

    m->unlockFlags<(1 | 2)>();

    return polarities;
}
