/**************************************************************
 *
 *       AIGPP Package // Manager_fraig.cc
 *
 *       Copyright (C) 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: $
 *         $Author: $
 *
 ***************************************************************/

#include "Manager.hh"

#include <lrabsutil/Hashes.hh>
#include <lrabsutil/Resources.hh>

#include "EquivalenceChecker.hh"

void
aigpp::Manager::toggleAutomaticFRAIGing(bool enable)
{
    _automaticFRAIGing = enable;

    if (_automaticFRAIGing) {
        makeFRAIG();
    }
}

bool
aigpp::Manager::automaticFRAIGing() const
{
    return _automaticFRAIGing;
}

void
aigpp::Manager::makeFRAIG(double timeout)
{
#if 0
    /* setup equivalence checker */
    EquivalenceChecker eqchk;
    {
        std::vector<InternalEdgeRef> roots;
        ExtRefTable* extreftab = getExtRefTable();

        for( unsigned int i = 0; i != extreftab->_targets.size(); ++i )
        {
            if( extreftab->_targets[i]._refCount == 0 ) continue;
            roots.push_back( extreftab->_targets[i]._target );
        }
        eqchk.setOriginalRoots( roots );
    }
#endif

#ifdef FRAIGMANAGER_TIMEOUT
    _fraigManager.makeFRAIG(timeout);
#else
    if (timeout > 0) {
        std::cout << "c ERROR: makeFRAIG called with timeout " << timeout
                  << " but aigpp is compiled without timeout support" << std::endl;
        assert(false);
    } else {
        _fraigManager.makeFRAIG();
    }
#endif

#if 0
    {
        std::vector<InternalEdgeRef> roots;
        ExtRefTable* extreftab = getExtRefTable();
        for( unsigned int i = 0; i != extreftab->_targets.size(); ++i )
        {
            if( extreftab->_targets[i]._refCount == 0 ) continue;
            roots.push_back( extreftab->_targets[i]._target );
        }
        std::cout << "\nfraiging: checking equivalence for " << roots.size() << " node pairs" << std::endl;

        bool ok = eqchk.checkNewRoots( roots );
        std::cout << "eqchk-result: " << ok << std::endl;
    }
#endif
}
