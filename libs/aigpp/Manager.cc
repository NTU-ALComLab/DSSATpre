/**************************************************************
 *
 *       AIGPP Package // Manager.cc
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
 *         $Date$
 *
 ***************************************************************/

#include "Manager.hh"

/* std */
#include <fstream>
#include <sstream>

/* LRABS utilities */
#include <lrabsutil/Assert.hh>
#include <lrabsutil/Resources.hh>

/* cudd */
#include <cuddInt.h>

/* local */
#include "Node.hh"
#include "StatManager.hh"

aigpp::Manager::Manager() :
    _extRefTable(new ExtRefTable),
    _unreducedNodes(0),
    _automaticFRAIGing(true),

    _verbosity(0),
    _nodeCount(0),
    _nextIndex(1),
    _dead(0),
    _reclaimed(0),
    _freeNodes(),

    _freeSim(),
    _simCreation(true),

    _nextUpdateSimVectorBin(0),
#ifdef USE_TEMPSIM
    _tempSimUpToDate(true),
#endif
#ifdef USE_COMPUTEDTABLE
    _computed(),
#endif
    _nodes(nullptr),
    _lastNode(nullptr),

    _toBeRemoved(-1),

    _solver(new Minisat::Solver),
    _noReset(false),
    _satChecksPresent(0),
    _nodesInCNF(0),
    _deletedNodesInCNF(0),
    _equivalentNodesInCNF(0),
    _mitersInCNF(0),
    _openMitersInCNF(0),
    _closedMitersInCNF(0),

    _nextAndSkipSat(false),
    _nextAndSatResult(nullptr),
    _nextAndSatResultIsInverted(false),
    _bddSweepingDuringQuantification(true),

#ifdef ENABLE_BDD
    _bddManager(nullptr),
    _nextReductionSize(100),
    _BDDSimDecay(0.01),
    _bddreducedratio(3),
    _bddReachedNodeLimit(false),
    _bddReachedTimeLimit(false),
    _bddReachedGrowthLimit(false),
#endif

    _skipTestWith(nullptr),
    _inCofactor(false),

    _flagsLocked(0),

    _const0(EdgeRef(getExtRefTable(), aigpp::const0)),
    _const1(EdgeRef(getExtRefTable(), aigpp::const1)),

    _rwManager(this),
    _fraigManager(this),
    _timeoutDuringLastFRAIG(false)
{
#ifdef ENABLE_BDD
    /* enable sifting */
    // Cudd_AutodynEnable( _bddManager, CUDD_REORDER_LAZY_SIFT );
#endif

#ifdef USE_COMPUTEDTABLE
    _computed.setMaxSize(10000);
#endif
}

aigpp::Manager::~Manager()
{
    delete _solver;

#ifdef ENABLE_BDD
    if (_bddManager != nullptr) {
        Cudd_Quit(_bddManager);
    }
#endif

    _const0.invalidate();
    _const1.invalidate();

    _extRefTable->clear();

    for (Node* n : _variables) {
        n->deref();

        if (n->refCount() == 0) {
            --_nodeCount;
        }
        assert(n->refCount() == 0);
    }

    garbageCollect();

    _variables.clear();

    assert(_nodeCount == 0);
    assert(_dead == 0);

    while (!_freeNodes.empty()) {
        delete _freeNodes.top();
        _freeNodes.pop();
    }

    while (!_freeSim.empty()) {
        delete _freeSim.top();
        _freeSim.pop();
    }

    delete _extRefTable;
}

aigpp::Settings&
aigpp::Manager::settings()
{
    return _settings;
}

const aigpp::Settings&
aigpp::Manager::settings() const
{
    return _settings;
}

aigpp::EdgeRef
aigpp::Manager::addVariable(const std::string& name)
{
    return EdgeRef(getExtRefTable(), addVariableInternal(name));
}

aigpp::InternalEdgeRef
aigpp::Manager::addVariableInternal(const std::string& name)
{
    assert(name == "" || _varMap.find(name) == _varMap.end());

    /*
     * create a new variable node.
     */
    Node* n = createNode(variableCount());
    n->setFlag<(Node::FLAG_ISREDUCED | Node::FLAG_UNMODIFIED)>();

    ref(n);

    _variables.push_back(n);
    _variableNames.push_back(name);
    if (!name.empty()) {
        _varMap[name] = n->varIndex();
    }

    if (simCreation()) {
        n->sim().generateRandom(random());
        _simTable.insert(n);
    }

    addToNodesList(n);

#ifdef USE_TEMPSIM
    n->setTempSim(_randomTempSim.getUInt());
#endif

    return InternalEdgeRef(n);
}

aigpp::EdgeRef
aigpp::Manager::lookupVariable(const std::string& name)
{
    auto v = _varMap.find(name);
    assert(v != _varMap.end());

    return variable(v->second);
}

void
aigpp::Manager::rewrite()
{
    const bool autoFraig = automaticFRAIGing();
    const bool simCreate = simCreation();

    toggleAutomaticFRAIGing(false);
    toggleSimCreation(false);

    _rwManager.rewrite();

    toggleSimCreation(simCreate);
    toggleAutomaticFRAIGing(autoFraig);
}

void
aigpp::Manager::rewriteCone(aigpp::Node* root)
{
    const bool autoFraig = automaticFRAIGing();
    const bool simCreate = simCreation();

    toggleAutomaticFRAIGing(false);
    toggleSimCreation(false);

    _rwManager.rewriteCone(root);

    toggleSimCreation(simCreate);
    toggleAutomaticFRAIGing(autoFraig);
}

std::size_t
aigpp::Manager::sharedSize(const std::vector<aigpp::EdgeRef>& roots) const
{
    std::vector<Node*> roots2;
    roots2.reserve(roots.size());

    for (const EdgeRef& r : roots) {
        Node* n = r.getInternal().node();
        if (n != nullptr) {
            roots2.push_back(n);
        }
    }

    return sharedSize(roots2);
}

std::size_t
aigpp::Manager::sharedSize(const std::vector<aigpp::Node*>& roots) const
{
    lockFlags<1>();

    std::stack<const Node*>& pending = _universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& parents = _universalConstNodeVector;
    assert(parents.empty());

    for (Node* r : roots) {
        if (r == nullptr) {
            continue;
        }

        pending.push(r);

        while (!pending.empty()) {
            /* node was counted/processed already */
            if (pending.top()->flag<1>()) {
                pending.pop();
            }
            /* node is an "and" node and its first parent wasn't processed yet */
            else if (!(pending.top()->isVar()) && !(pending.top()->parent1().node()->flag<1>())) {
                pending.push(pending.top()->parent1().node());
            }
            /* node is an "and" node and its second parent wasn't prccocessed yet */
            else if (!(pending.top()->isVar()) && !(pending.top()->parent2().node()->flag<1>())) {
                pending.push(pending.top()->parent2().node());
            }
            /* node is an "and" node with both parents processed, or a variable node
             */
            else {
                parents.push_back(pending.top());
                pending.top()->setFlag<1>();

                pending.pop();
            }
        }
    }

    /* clear marks of all processed nodes */
    for (const Node* n : parents) {
        n->unsetFlag<1>();
    }

    const std::size_t count = parents.size();
    parents.clear();

    unlockFlags<1>();

    return count;
}

std::vector<int>
aigpp::Manager::sharedSupport(const std::vector<aigpp::EdgeRef>& roots) const
{
    if (roots.empty()) {
        return std::vector<int>();
    }

    std::vector<Node*> roots2;
    roots2.reserve(roots.size());

    for (const EdgeRef& r : roots) {
        Node* n = r.getInternal().node();
        if (n != nullptr) {
            roots2.push_back(n);
        }
    }

    return sharedSupport(roots2);
}

std::vector<int>
aigpp::Manager::sharedSupport(const std::vector<aigpp::Edge>& roots) const
{
    if (roots.empty()) {
        return std::vector<int>();
    }

    std::vector<Node*> roots2;
    roots2.reserve(roots.size());

    for (const Edge& r : roots) {
        Node* n = r.node();
        if (n != nullptr) {
            roots2.push_back(n);
        }
    }

    return sharedSupport(roots2);
}

std::vector<int>
aigpp::Manager::sharedSupport(const std::vector<aigpp::Node*>& roots) const
{
    if (roots.empty()) {
        return std::vector<int>();
    }

    lockFlags<1>();

    std::stack<const Node*>& pending = _universalConstNodeStack;
    assert(pending.empty());

    std::vector<const Node*>& parents = _universalConstNodeVector;
    assert(parents.empty());

    std::vector<int> result;

    for (Node* r : roots) {
        if (r == nullptr) {
            continue;
        }

        pending.push(r);

        while (!pending.empty()) {
            /* node was counted/processed already */
            if (pending.top()->flag<1>()) {
                pending.pop();
            }
            /* node is an unprocessed variable */
            else if (pending.top()->isVar()) {
                parents.push_back(pending.top());
                pending.top()->setFlag<1>();

                result.push_back(pending.top()->varIndex());

                pending.pop();
            }
            /* node is an "and" node and its first parent wasn't processed yet */
            else if (!(pending.top()->parent1().node()->flag<1>())) {
                pending.push(pending.top()->parent1().node());
            }
            /* node is an "and" node and its second parent wasn't processed yet */
            else if (!(pending.top()->parent2().node()->flag<1>())) {
                pending.push(pending.top()->parent2().node());
            }
            /* node is an "and" node with both parents processed */
            else {
                parents.push_back(pending.top());
                pending.top()->setFlag<1>();

                pending.pop();
            }
        }
    }

    /* clear marks of all processed nodes */
    for (const Node* n : parents) {
        n->unsetFlag<1>();
    }
    parents.clear();
    unlockFlags<1>();

    return result;
}

long
aigpp::Manager::compareQuality(aigpp::Node* n1, aigpp::Node* n2) const
{
    if (n1->isVar() && !(n2->isVar())) {
        return -1;
    }

    if (n2->isVar() && !(n1->isVar())) {
        return +1;
    }

    if (settings().getNodeReplacement() == Settings::NODEREPLACEMENT_NONE) {
        return 0;
    } else if (settings().getNodeReplacement() == Settings::NODEREPLACEMENT_CONESIZE) {
        if (_toBeRemoved != -1) {
            bool v1 = n1->varInSupport(_toBeRemoved);
            bool v2 = n2->varInSupport(_toBeRemoved);

            if (v1 && !v2) {
                return +1;
            } else if (!v1 && v2) {
                return -1;
            }
        }

        return (long)n1->nodeCount() - (long)n2->nodeCount();
    } else {
        NEVER_GET_HERE;
        return 0;
    }
}

void
aigpp::Manager::printStats(std::ostream& os) const
{
    os << "#############################################" << std::endl;

    os << "nodes            " << _nodeCount << std::endl
       << "dead             " << _dead << std::endl
       << "variables        " << variableCount() << std::endl
       << StatManager::instance() << std::endl
       << "fraig_manager:" << std::endl;

    os << "\nfraig_manager:\n";
    _fraigManager.printStats(os);
    os << "\nrewriting_manager:\n";
    _rwManager.printStats(os);
}

void
aigpp::Manager::checkIntegrity() const
{
    std::cout << "integrity check" << std::endl;

    _unique.checkIntegrity(_nodes);
    _simTable.checkIntegrity(_nodes);

    /* check that nodecounts... are correct */
    std::size_t mynodecount      = 0;
    std::size_t mydeadcount      = 0;
    std::size_t myunreducedcount = 0;
    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->refCount() == 0) {
            ++mydeadcount;
        } else {
            ++mynodecount;
        }

        if (!(n->isReduced())) {
            ++myunreducedcount;
        }
    }
    if (mynodecount != nodeCount() || mydeadcount != deadNodeCount() || myunreducedcount != unreducedNodes()) {
        std::cout << "bad counts:" << std::endl
                  << "nodes:     " << mynodecount << " " << nodeCount() << std::endl
                  << "dead:      " << mydeadcount << " " << deadNodeCount() << std::endl
                  << "unreduced: " << myunreducedcount << " " << unreducedNodes() << std::endl;
        assert(false);
    }

    /* check that reference counts are correct */
    std::map<Node*, int> countrefs;
    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        countrefs[n] = 0;
    }

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            ++countrefs[n];
            continue;
        }

        if (n->refCount() == 0) {
            continue;
        }

        ++(countrefs[n->parent1().node()]);
        ++(countrefs[n->parent2().node()]);
    }

    for (const ExtRefTable::Item& p : _extRefTable->_targets) {
        if (!(p._target.isConstant())) {
            ++(countrefs[p._target.node()]);
        }
    }

    //    bool valid = true;

    //    for( const std::pair<Node*, int>& p: countrefs )
    //    {
    //        assert( p.first->refCount() >= p.second );
    //        if( p.first->refCount() < p.second )
    //        {
    //            std::cout << "node " << p.first << " has a refcount of " <<
    //            p.first->refCount() << " but is referenced at least " <<
    //            p.second << " times" << std::endl; valid = false;
    //        }
    //        else if( p.first->refCount() < p.second )
    //        {
    //            std::cout << "node " << p.first << " has a refcount of " <<
    //            p.first->refCount() << " but is referenced " << p.second << "
    //            times" << std::endl;
    //        }
    //    }

    //    assert( valid );

    std::cout << "everything is ok!" << std::endl;
}

int
aigpp::Manager::modifiedNodes() const
{
    int modified = 0;

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (!(n->flag<Node::FLAG_UNMODIFIED>())) {
            assert(!(n->isVar()));

            ++modified;
        }
    }

    return modified;
}

void
aigpp::Manager::resetModified()
{
    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        n->setFlag<Node::FLAG_UNMODIFIED>();
    }
}

void
aigpp::Manager::toggleSimCreation(bool on)
{
    /* fraiging -> simcreation */
    assert(!automaticFRAIGing() || simCreation());

    if (automaticFRAIGing() && !on) {
        std::cout << "WARNING: You are trying to disable sim vector creation while "
                     "automatic\n"
                     "fraiging is enabled. This is not possible!\n"
                     "Leaving sim vector creation enabled!"
                  << std::endl;
        return;
    }

    if (on == simCreation()) {
        return;
    }

    _simCreation = on;

    /* recreate sim of all nodes that don't have a sim */
    if (simCreation()) {
        for (Node* n = _nodes; n != nullptr; n = n->next()) {
            if (n->_sim == nullptr) {
                n->_sim = createSim();
                if (n->isVar()) {
                    n->sim().generateRandom(random());
                } else {
                    n->updateSim();
                }
                _simTable.insert(n);
            }
        }
    }
}

bool
aigpp::Manager::simCreation() const
{
    return _simCreation;
}

aigpp::Node*
aigpp::Manager::createNode(const aigpp::Edge& p1, const aigpp::Edge& p2)
{
    Node* newnode;

    /* no free nodes on the "freenodes" stack -> create new object */
    if (_freeNodes.empty()) {
        newnode = new Node(this, p1, p2);
    }
    /* there is a free node object -> take and initialize it */
    else {
        newnode = _freeNodes.top();
        _freeNodes.pop();

        if (simCreation()) {
            newnode->_sim = createSim();
        }

        newnode->set(this, p1, p2);
    }

    return newnode;
}

aigpp::Node*
aigpp::Manager::createNode(int varindex)
{
    Node* newnode;

    newnode = new Node(this, varindex);

    return newnode;
}

void
aigpp::Manager::releaseNode(aigpp::Node* n)
{
    if (n->_sim != nullptr) {
        releaseSim(n->_sim);
        n->_sim = nullptr;
    }

    _freeNodes.push(n);
}

aigpp::SimVector*
aigpp::Manager::createSim()
{
    SimVector* newsim;

    if (_freeSim.empty()) {
        newsim = new SimVector;
    } else {
        newsim = _freeSim.top();
        _freeSim.pop();
    }

    return newsim;
}

void
aigpp::Manager::releaseSim(aigpp::SimVector* s)
{
    _freeSim.push(s);
}

void
aigpp::Manager::clearNodeCaches()
{
    _cache.clear();
}

void
aigpp::Manager::clearNodeCachesZ()
{
    _cacheZ.clear();
}

void
aigpp::Manager::replace(aigpp::Node* oldNode, aigpp::Node* newNode, bool inverted)
{
    assert(simCreation() == true);

    assert(!(oldNode == 0) && !(oldNode->isVar()));
    assert(!(newNode == 0) && !(newNode->isVar()));
    assert(oldNode->manager() == newNode->manager());
    assert(newNode->refCount() == 0);

    // assert( !( newNode->isNAND() ) );

    double startTime = lrabs::cpuTime();

    if (oldNode->manager()->verbosity() >= 3) {
        std::cout << "replace node " << oldNode->index() << (inverted ? " by complement of node " : " by node ")
                  << newNode->index() << std::endl
                  << (oldNode->parent1().isInverted() ? "!" : "") << oldNode->parent1().node()->index() << " "
                  << (oldNode->isNAND() ? "NAND" : "AND") << " " << (oldNode->parent2().isInverted() ? "!" : "")
                  << oldNode->parent2().node()->index() << " <- " << (newNode->parent1().isInverted() ? "!" : "")
                  << newNode->parent1().node()->index() << " " << (newNode->isNAND() ? "NAND" : "AND") << " "
                  << (newNode->parent2().isInverted() ? "!" : "") << newNode->parent2().node()->index() << std::endl;
    }

#if 0
  if( !inverted )
  {
    if( !( oldNode->sim() == newNode->sim() ) )
    {
      std::cout << "!inverted and sim != sim" << std::endl;

      std::cout << "sim equivalent: " << ( oldNode->sim() == newNode->sim() ) << std::endl
                << "sim antivalent: " << ( oldNode->sim().equalInverted( newNode->sim() ) ) << std::endl;
      /*
      std::cout << "sat equivalent: " << satEqual( oldNode, newNode ) << std::endl
                << "sat antivalent: " << satEqualNegated( oldNode, newNode ) << std::endl;
      */
      assert( false );
    }
  }
  else
  {
    if( !( oldNode->sim().equalInverted( newNode->sim() ) ) )
    {
      std::cout << "inverted and sim != !sim" << std::endl;

      std::cout << "sim equivalent: " << ( oldNode->sim() == newNode->sim() ) << std::endl
                << "sim antivalent: " << ( oldNode->sim().equalInverted( newNode->sim() ) ) << std::endl;
      /*
      std::cout << "sat equivalent: " << satEqual( oldNode, newNode ) << std::endl
                << "sat antivalent: " << satEqualNegated( oldNode, newNode ) << std::endl;
      */
      assert( false );
    }
  }
#endif

    if (oldNode->isSatVarValid()) {
        oldNode->invalidateSatVar();

        if (newNode->isSatVarValid()) {
            newNode->invalidateSatVar();
            /*
              oldNode->setSatVar( newNode->_satVar );
              _nodesWithSATVars.erase( newNode );
            */
            _nodesWithSATVars.erase(newNode);
            _nodesWithSATVars.erase(oldNode);
        } else {
            _nodesWithSATVars.erase(oldNode);
            assert(!inverted);
        }
    } else {
        if (newNode->isSatVarValid()) {
            newNode->invalidateSatVar();
            /*
              oldNode->setSatVar( newNode->_satVar );
              _nodesWithSATVars.erase( newNode );
              _nodesWithSATVars.insert( oldNode );
            */
            _nodesWithSATVars.erase(newNode);
        } else {
            // nothing
        }
    }

    _unique.remove(oldNode);

    if (inverted != newNode->isNAND()) {
        StatManager::instance().incReplacementsInverted();
        oldNode->setFlag<Node::FLAG_ISNAND>();
    } else {
        StatManager::instance().incReplacements();
        oldNode->unsetFlag<Node::FLAG_ISNAND>();
    }

    ref(newNode->parent1().node());
    ref(newNode->parent2().node());

    deref(oldNode->parent1().node());
    deref(oldNode->parent2().node());

    oldNode->_parent1 = newNode->parent1();
    oldNode->_parent2 = newNode->parent2();

    fixNodesList();

    _unique.insert(oldNode);

    StatManager::instance().incReplacementTime(lrabs::cpuTime() - startTime);
}

void
aigpp::Manager::fixNodesList()
{
    std::vector<Node*> nodesList;
    nodesList.reserve(_nodeCount);

    lockFlags<1>();

    std::stack<Node*>& pending = _universalNodeStack;
    assert(pending.empty());

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        nodesList.push_back(n);
        n->setFlag<1>();
    }

    _nodes    = nullptr;
    _lastNode = nullptr;

    for (Node* n : nodesList) {
        if (!(n->flag<1>())) {
            continue;
        }

        pending.push(n);
        while (!pending.empty()) {
            // node was already processed -> skip it
            if (!(pending.top()->flag<1>())) {
                pending.pop();
            }
            // variable node
            else if (pending.top()->isVar()) {
                pending.top()->unsetFlag<1>();
                addToNodesList(pending.top());
                pending.pop();
            }
            // "normal" and node
            // is first parent already processed
            else if (pending.top()->parent1().node()->flag<1>()) {
                pending.push(pending.top()->parent1().node());
            }
            // is second parent already processed
            else if (pending.top()->parent2().node()->flag<1>()) {
                pending.push(pending.top()->parent2().node());
            }
            // all parent nodes processed
            else {
                pending.top()->unsetFlag<1>();
                addToNodesList(pending.top());
                pending.pop();
            }
        }
    }

    unlockFlags<1>();
}

void
aigpp::Manager::checkNodesList() const
{
    lockFlags<1>();

    int i = 0;

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            n->setFlag<1>();
        } else {
            if (!(n->parent1().node()->flag<1>()) || !(n->parent2().node()->flag<1>())) {
                std::cout << "check failed on node #" << i << " (id=" << n->index() << ")" << std::endl
                          << *n << std::endl;

                if (!(n->parent1().node()->flag<1>())) {
                    std::cout << "first parent: " << n->parent1() << std::endl;
                }
                if (!(n->parent2().node()->flag<1>())) {
                    std::cout << "second parent: " << n->parent2() << std::endl;
                }

                abort();
            }
            n->setFlag<1>();
        }
        ++i;
    }

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        n->unsetFlag<1>();
    }

    unlockFlags<1>();
}

void
aigpp::Manager::existentialQInplace(aigpp::InternalEdgeRef& root, const aigpp::InternalEdgeRef& var, bool internalCall)
{
    assert(!var.isConstant());
    assert(var.node()->isVar());

    if (root.isConstant()) {
        return;
    }

    StatManager::instance().incQuantification();

    if (verbosity() >= 1) {
        std::cout << "quantifying variable " << var << std::endl;
    }

    if (!internalCall) {
        if (!(root.node()->varInSupport(var.node()->varIndex()))) {
            return;
        }
    }

    if (root.node()->refCount() == 1) {
#if 0
        std::cout << "quantification root refcount=1" << std::endl;
#endif
        _skipTestWith = root.node();
    }

    InternalEdgeRef cofPos = root.cofactor(var, false);
    InternalEdgeRef cofNeg = root.cofactor(!var, false);

    root.clear();
    root = cofPos + cofNeg;

    _skipTestWith = nullptr;

    if (!internalCall && !(root.isConstant())) {
        /* QUIET version */
#if 0
        std::cout << "size after: " << root.nodeCount() << std::endl;
#endif

        BDDSweepInternal(root);

        /* QUIET version */
#if 0
        std::cout << "size after (BDD): " << root.nodeCount() << std::endl;
#endif

        if (settings().getRedundancyRemoval() == Settings::REDUNDANCYREMOVAL_QUANTIFICATION) {
            root.node()->findRedundantVars(true);
        }
    }
}

void
aigpp::Manager::existentialQInplace(aigpp::InternalEdgeRef& root, const std::vector<aigpp::InternalEdgeRef>& vars)
{
    if (root.isConstant()) {
        return;
    }

    for (const aigpp::InternalEdgeRef& var : vars) {
        if (root.isConstant()) {
            break;
        }
        if (!(root.node()->varInSupport(var.node()->varIndex()))) {
            continue;
        }

        existentialQInplace(root, var, /*internalCall=*/true);

        BDDSweepInternal(root);

        if (!(root.isConstant()) && settings().getRedundancyRemoval() == Settings::REDUNDANCYREMOVAL_QUANTIFICATION) {
            root.node()->findRedundantVars(true);
        }
    }
}

void
aigpp::Manager::universalQInplace(aigpp::InternalEdgeRef& root, const std::vector<aigpp::InternalEdgeRef>& vars)
{
    if (root.isConstant()) {
        return;
    }

    root.toggleInverted();
    existentialQInplace(root, vars);
    root.toggleInverted();
}
