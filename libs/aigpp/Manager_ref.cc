/**************************************************************
 *
 *       AIGPP Package // Manager_ref.cc
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 452 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "Manager.hh"

#include "StatManager.hh"

void
aigpp::Manager::ref(aigpp::Node* n)
{
    StatManager::instance().incRefs();

    if (n != nullptr) {
        n->ref();

        if (n->refCount() == 1 /* && !( n->isVar() ) */) {
            ++_nodeCount;

            StatManager::instance().incCreated();
            StatManager::instance().updateMaxNodes(_nodeCount);
        }
    }
}

void
aigpp::Manager::deref(aigpp::Node* n)
{
    StatManager::instance().incDerefs();

    if (n != nullptr /*&& n->refCount() > 0*/) {
        recursiveDeref(n);
    }
}
