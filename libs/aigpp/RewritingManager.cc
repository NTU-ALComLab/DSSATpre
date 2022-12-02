/**************************************************************
 *
 *       AIGPP // RewritingManager.cc
 *
 *       Copyright (C) 2008 Florian Pigorsch
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

#include "RewritingManager.hh"

/* stdlib */
#include <fstream>

/* lrabs */
#include <lrabsutil/Math.hh>
#include <lrabsutil/Resources.hh>

/* local */
#include "EquivalenceChecker.hh"
#include "Manager.hh"

#define _unused(x) ((void)(x))

aigpp::RewritingManager::RewritingManager(aigpp::Manager* manager) :
    _manager(manager),
    _nodes(nullptr),
    _nodesCapacity(0),
    _nodesUsed(0),
    _rewritten(nullptr),
    _modified(nullptr),
    _sim(nullptr),
    _hasSim(nullptr),
    _cuts(nullptr),
    _mark(0),
    _deadNodes(0),
#if 0
     _modifiedRWTable(false),
#endif
    _totalApplications(0),
    _totalNodesRewritten(0),
    _totalCutsComputed(0),
    _totalNodesChecked(0),
    _totalImplementationsChecked(0),
    _totalDelta(0),
    _totalEstimatedDelta(0),
    _totalRewritingTime(0),

    _totalApplicationsLocal(0),
    _totalNodesRewrittenLocal(0),
    _totalDeltaLocal(0),
    _totalEstimatedDeltaLocal(0),
    _totalRewritingTimeLocal(0)
{
    _npnClassImpl.assign(NumberOfNPNClasses, NPNClassImplementation());

    loadStaticRWDB();
}

aigpp::RewritingManager::~RewritingManager()
{
#if 0
    if( _modifiedRWTable )
    {
        std::cout << "rw table modified -> need to save" << std::endl;
    }
#endif

    /* cleanup rw db */
    for (const auto& i : _npnClassImpl) {
        for (const auto& p : i.implementations) {
            _rwTable.deref(p);
        }
    }

    if (_nodes != nullptr) {
        delete[] _nodes;
        delete[] _rewritten;
        delete[] _modified;
        delete[] _sim;
        delete[] _hasSim;
        for (unsigned int i = 0; i != _nodesCapacity; ++i) {
            if (_cuts[i] != nullptr) {
                releaseCuts(_cuts[i]);
            }
        }
        delete[] _cuts;
    }

    for (unsigned int i = 0; i != _freeCuts.size(); ++i) {
        delete _freeCuts[i];
    }
}

void
aigpp::RewritingManager::rewrite()
{
    if (_manager->nodeCount() == 0) {
        return;
    }

    ++_totalApplications;

    double startTime = lrabs::cpuTime();

    _manager->garbageCollect();
    _manager->resetSolver(true);

    assert(_manager->deadNodeCount() == 0);

#if 0
    /* setup equivalence checker */
    EquivalenceChecker eqchk;
    {
        std::vector<InternalEdgeRef> roots;
        ExtRefTable* extreftab = _manager->getExtRefTable();

        for( unsigned int i = 0; i != extreftab->_targets.size(); ++i )
        {
            if( extreftab->_targets[i]._refCount == 0 ) continue;
            roots.push_back( extreftab->_targets[i]._target );
        }
        eqchk.setOriginalRoots( roots );
    }
#endif

    auto         totalNodes  = static_cast<unsigned int>(_manager->nodeCount());
    unsigned int newCapacity = 2 * totalNodes;

    /* clear old structures */
    for (unsigned int i = 0; i != _nodesCapacity; ++i) {
        releaseCuts(_cuts[i]);
    }
    _id.clear();
    _andMap.clear();
    _backMap.clear();
    _impl2node.assign(_impl2node.size(), -1000000000);
    _impl2nodeUsed.clear();
    _markedNodes.clear();
    _mark      = 0;
    _deadNodes = 0;
    _nodesUsed = totalNodes;

    /* reuse or reallocate arrays */
    if (_nodes != nullptr) {
        /* arrays have to be reallocated */
        if (_nodesCapacity < newCapacity) {
            delete[] _nodes;
            delete[] _rewritten;
            delete[] _modified;
            delete[] _sim;
            delete[] _hasSim;
            delete[] _cuts;

            _nodesCapacity = newCapacity;

            _nodes     = new NodeStruct[_nodesCapacity];
            _rewritten = new unsigned int[_nodesCapacity];
            _modified  = new unsigned int[_nodesCapacity];
            _sim       = new unsigned int[_nodesCapacity];
            _hasSim    = new unsigned int[_nodesCapacity];
            _cuts      = new CutStruct*[_nodesCapacity];
        }
    }
    /* allocate new arrays */
    else {
        _nodesCapacity = newCapacity;

        _nodes     = new NodeStruct[_nodesCapacity];
        _rewritten = new unsigned int[_nodesCapacity];
        _modified  = new unsigned int[_nodesCapacity];
        _sim       = new unsigned int[_nodesCapacity];
        _hasSim    = new unsigned int[_nodesCapacity];
        _cuts      = new CutStruct*[_nodesCapacity];
    }

    /* clear arrays */
    std::fill(_rewritten, _rewritten + _nodesCapacity, 0);
    std::fill(_modified, _modified + _nodesCapacity, 0);
    std::fill(_hasSim, _hasSim + _nodesCapacity, 0);
    std::fill(_cuts, _cuts + _nodesCapacity, (CutStruct*)nullptr);

    /* copy AIG */
    unsigned int index = 0;
    for (Node* n = _manager->_nodes; n != nullptr; n = n->next(), ++index) {
        _id[n] = index;

        NodeStruct& ns = _nodes[index];

        if (n->isVar()) {
            ns.var  = true;
            ns.nand = false;
            assert(n->refCount() < (1 << 27));
            ns.ref = n->refCount();
        } else {
            ns.var  = false;
            ns.nand = n->isNAND();
            ns.neg1 = n->parent1().isInverted();
            ns.neg2 = n->parent2().isInverted();
            assert(n->refCount() < (1 << 27));
            ns.ref = n->refCount();

            ns.p1 = _id[n->parent1().node()];
            ns.p2 = _id[n->parent2().node()];

            std::pair<int, int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0), 2 * ns.p2 + (ns.neg2 ? 1 : 0));
            if (ip.first > ip.second) {
                std::swap(ip.first, ip.second);
            }

            int andindex = 2 * index + (ns.nand ? 1 : 0);

            _andMap[ip] = andindex;
        }
    }

    assert(index == totalNodes);

    int estimatedDelta = 0;

    for (unsigned int nodeindex = 0; nodeindex != totalNodes; ++nodeindex) {
        if (_rewritten[nodeindex] != 0) {
            continue;
        }

        NodeStruct& node = _nodes[nodeindex];
        if (node.ref == 0) {
            continue;
        }

        CutStruct* cuts = computeCuts(nodeindex);

        if (node.var) {
            continue;
        }
        if (_nodes[node.p1].ref >= 2 && _nodes[node.p2].ref >= 2) {
            continue;
        }

        _totalNodesChecked += 1;

        NPNImplInfo bestimpl = getBestImplementation(nodeindex, cuts);

        if (bestimpl.cutIndex != -1) {
            _totalNodesRewritten += 1;

            estimatedDelta += bestimpl.delta;

            /* find selected cut in cutlist */
            CutStruct* cutit = cuts;
            for (int i = 0; i != bestimpl.cutIndex; ++i) {
                cutit = cutit->next;
            }

            if (bestimpl.npnclass.npnClass == 0) {
                FourInputFunction f0 = getFunction(nodeindex, *cutit);
                assert(f0 == 0 || f0 == 0xFFFF);
                bool isOne = (f0 == 0xFFFF);

                buildNewConeConstant(nodeindex, isOne);
            } else {
                buildNewCone(nodeindex, *cutit, bestimpl);
            }
        }
    }

    auto realNodesBefore = static_cast<int>(_manager->nodeCount());

    {
        unsigned int index = 0;

        _backMap.resize(totalNodes);

        for (Node* n = _manager->_nodes; n != nullptr; n = n->next(), ++index) {
            _backMap[index] = InternalEdgeRef(n);
        }
        assert(index == totalNodes);
    }

    std::vector<bool>            processed(_nodesUsed, false);
    std::vector<InternalEdgeRef> cache(_nodesUsed);

    for (const auto& p : _manager->getExtRefTable()->_indexMap) {
        Node* node = Edge(p.first).node();

        int rwid = _id[node];

        if (processed[rwid]) {
            continue;
        }

        recursiveRebuild(rwid, processed, cache);
    }

    cache.clear();
    _backMap.clear();

    _manager->fixNodesList();
    _manager->garbageCollect();

    auto realNodesAfter = static_cast<int>(_manager->nodeCount());

    _totalEstimatedDelta += estimatedDelta;
    _totalDelta += (realNodesAfter - realNodesBefore);
    _totalRewritingTime += (lrabs::cpuTime() - startTime);

#if 0
    {
        std::vector<InternalEdgeRef> roots;
        ExtRefTable* extreftab = _manager->getExtRefTable();
        for( unsigned int i = 0; i != extreftab->_targets.size(); ++i )
        {
            if( extreftab->_targets[i]._refCount == 0 ) continue;
            roots.push_back( extreftab->_targets[i]._target );
        }
        std::cout << "\nrewriting: checking equivalence for " << roots.size() << " node pairs" << std::endl;

        bool ok = eqchk.checkNewRoots( roots );
        std::cout << "eqchk-result: " << ok << std::endl;
    }
#endif

#if 0
    if( _modifiedRWTable )
    {
        saveRWDB( "rwdb.txt", _npnClassImpl, _rwTable );
    }
#endif
}

void
aigpp::RewritingManager::rewriteCone(aigpp::Node* root)
{
    /* skip constants, variables, and dead nodes */
    if (root == nullptr) {
        return;
    }
    if (root->isVar()) {
        return;
    }
    if (root->refCount() == 0) {
        return;
    }

    ++_totalApplications;

    double startTime = lrabs::cpuTime();

    _manager->garbageCollect();
    _manager->resetSolver(true);

    assert(_manager->deadNodeCount() == 0);
    auto         totalNodes  = static_cast<unsigned int>(_manager->nodeCount());
    unsigned int newCapacity = 2 * totalNodes;

    /* clear old structures */
    for (unsigned int i = 0; i != _nodesCapacity; ++i) {
        releaseCuts(_cuts[i]);
    }
    _id.clear();
    _andMap.clear();
    _backMap.clear();
    _impl2node.assign(_impl2node.size(), -1000000000);
    _impl2nodeUsed.clear();
    _markedNodes.clear();
    _mark      = 0;
    _deadNodes = 0;
    _nodesUsed = totalNodes;

    /* reuse or reallocate arrays */
    if (_nodes != nullptr) {
        /* arrays have to be reallocated */
        if (_nodesCapacity < newCapacity) {
            delete[] _nodes;
            delete[] _rewritten;
            delete[] _modified;
            delete[] _sim;
            delete[] _hasSim;
            delete[] _cuts;

            _nodesCapacity = newCapacity;

            _nodes     = new NodeStruct[_nodesCapacity];
            _rewritten = new unsigned int[_nodesCapacity];
            _modified  = new unsigned int[_nodesCapacity];
            _sim       = new unsigned int[_nodesCapacity];
            _hasSim    = new unsigned int[_nodesCapacity];
            _cuts      = new CutStruct*[_nodesCapacity];
        }
    }
    /* allocate new arrays */
    else {
        _nodesCapacity = newCapacity;

        _nodes     = new NodeStruct[_nodesCapacity];
        _rewritten = new unsigned int[_nodesCapacity];
        _modified  = new unsigned int[_nodesCapacity];
        _sim       = new unsigned int[_nodesCapacity];
        _hasSim    = new unsigned int[_nodesCapacity];
        _cuts      = new CutStruct*[_nodesCapacity];
    }

    /* clear arrays */
    std::fill(_rewritten, _rewritten + _nodesCapacity, 0);
    std::fill(_modified, _modified + _nodesCapacity, 0);
    std::fill(_hasSim, _hasSim + _nodesCapacity, 0);
    std::fill(_cuts, _cuts + _nodesCapacity, (CutStruct*)nullptr);

    /* copy AIG */
    unsigned int index = 0;
    for (Node* n = _manager->_nodes; n != nullptr; n = n->next(), ++index) {
        _id[n] = index;

        NodeStruct& ns = _nodes[index];

        if (n->isVar()) {
            ns.var  = true;
            ns.nand = false;
            assert(n->refCount() < (1 << 27));
            ns.ref = n->refCount();
        } else {
            ns.var  = false;
            ns.nand = n->isNAND();
            ns.neg1 = n->parent1().isInverted();
            ns.neg2 = n->parent2().isInverted();
            assert(n->refCount() < (1 << 27));
            ns.ref = n->refCount();

            ns.p1 = _id[n->parent1().node()];
            ns.p2 = _id[n->parent2().node()];

            std::pair<int, int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0), 2 * ns.p2 + (ns.neg2 ? 1 : 0));
            if (ip.first > ip.second) {
                std::swap(ip.first, ip.second);
            }

            int andindex = 2 * index + (ns.nand ? 1 : 0);

            _andMap[ip] = andindex;
        }
    }
    assert(index == totalNodes);

    /* mark "non-cone" nodes as rewritten, such that they are skipped by the
     * rewriting procedure */
    std::set<Node*> cone = root->cone();
    index                = 0;
    for (Node* n = _manager->_nodes; n != nullptr; n = n->next(), ++index) {
        if (cone.find(n) != cone.end()) {
            _rewritten[index] = 1;
        }
    }

    int estimatedDelta = 0;

    for (unsigned int nodeindex = 0; nodeindex != totalNodes; ++nodeindex) {
        if (_rewritten[nodeindex] != 0) {
            continue;
        }

        NodeStruct& node = _nodes[nodeindex];
        if (node.ref == 0) {
            continue;
        }

        CutStruct* cuts = computeCuts(nodeindex);

        if (node.var) {
            continue;
        }
        if (_nodes[node.p1].ref >= 2 && _nodes[node.p2].ref >= 2) {
            continue;
        }

        _totalNodesChecked += 1;

        NPNImplInfo bestimpl = getBestImplementation(nodeindex, cuts);

        if (bestimpl.cutIndex != -1) {
            _totalNodesRewritten += 1;

            estimatedDelta += bestimpl.delta;

            /* find selected cut in cutlist */
            CutStruct* cutit = cuts;
            for (int i = 0; i != bestimpl.cutIndex; ++i) {
                cutit = cutit->next;
            }

            if (bestimpl.npnclass.npnClass == 0) {
                FourInputFunction f0 = getFunction(nodeindex, *cutit);
                assert(f0 == 0 || f0 == 0xFFFF);
                bool isOne = (f0 == 0xFFFF);

                buildNewConeConstant(nodeindex, isOne);
            } else {
                buildNewCone(nodeindex, *cutit, bestimpl);
            }
        }
    }

    auto realNodesBefore = static_cast<int>(_manager->nodeCount());

    {
        unsigned int index = 0;

        _backMap.resize(totalNodes);

        for (Node* n = _manager->_nodes; n != nullptr; n = n->next(), ++index) {
            _backMap[index] = InternalEdgeRef(n);
        }
        assert(index == totalNodes);
    }

    std::vector<bool>            processed(_nodesUsed, false);
    std::vector<InternalEdgeRef> cache(_nodesUsed);

    for (const auto& p : _manager->getExtRefTable()->_indexMap) {
        Node* node = Edge(p.first).node();
        int   rwid = _id[node];

        if (processed[rwid]) {
            continue;
        }

        recursiveRebuild(rwid, processed, cache);
    }

    cache.clear();
    _backMap.clear();

    _manager->fixNodesList();
    _manager->garbageCollect();

    auto realNodesAfter = static_cast<int>(_manager->nodeCount());

    _totalEstimatedDelta += estimatedDelta;
    _totalDelta += (realNodesAfter - realNodesBefore);
    _totalRewritingTime += (lrabs::cpuTime() - startTime);
}

aigpp::EdgeRef
aigpp::RewritingManager::localRewriteCone(const aigpp::EdgeRef& root)
{
    /* skip constants, variables, and dead nodes */
    if (root.isConstant()) {
        return root;
    }
    if (root.isVariable()) {
        return root;
    }

    const InternalEdgeRef& iroot = root.getInternal();
    assert(iroot.node()->refCount() > 0);

    ++_totalApplicationsLocal;

    double startTime = lrabs::cpuTime();

    std::vector<Node*> cone = iroot.node()->coneVector();

    auto         totalNodes  = static_cast<unsigned int>(cone.size());
    unsigned int newCapacity = 2 * totalNodes;

    /* clear old structures */
    for (unsigned int i = 0; i != _nodesCapacity; ++i) {
        releaseCuts(_cuts[i]);
    }
    _id.clear();
    _andMap.clear();
    _backMap.clear();
    _impl2node.assign(_impl2node.size(), -1000000000);
    _impl2nodeUsed.clear();
    _markedNodes.clear();
    _mark      = 0;
    _deadNodes = 0;
    _nodesUsed = totalNodes;

    /* reuse or reallocate arrays */
    if (_nodes != nullptr) {
        /* arrays have to be reallocated */
        if (_nodesCapacity < newCapacity) {
            delete[] _nodes;
            delete[] _rewritten;
            delete[] _modified;
            delete[] _sim;
            delete[] _hasSim;
            delete[] _cuts;

            _nodesCapacity = newCapacity;

            _nodes     = new NodeStruct[_nodesCapacity];
            _rewritten = new unsigned int[_nodesCapacity];
            _modified  = new unsigned int[_nodesCapacity];
            _sim       = new unsigned int[_nodesCapacity];
            _hasSim    = new unsigned int[_nodesCapacity];
            _cuts      = new CutStruct*[_nodesCapacity];
        }
    }
    /* allocate new arrays */
    else {
        _nodesCapacity = newCapacity;

        _nodes     = new NodeStruct[_nodesCapacity];
        _rewritten = new unsigned int[_nodesCapacity];
        _modified  = new unsigned int[_nodesCapacity];
        _sim       = new unsigned int[_nodesCapacity];
        _hasSim    = new unsigned int[_nodesCapacity];
        _cuts      = new CutStruct*[_nodesCapacity];
    }

    /* clear arrays */
    std::fill(_rewritten, _rewritten + _nodesCapacity, 0);
    std::fill(_modified, _modified + _nodesCapacity, 0);
    std::fill(_hasSim, _hasSim + _nodesCapacity, 0);
    std::fill(_cuts, _cuts + _nodesCapacity, (CutStruct*)nullptr);

    /* copy aig */
    unsigned int index = 0;
    for (std::vector<Node*>::const_iterator it = cone.begin(); it != cone.end(); ++it, ++index) {
        Node* n = *it;
        _id[n]  = index;

        NodeStruct& ns = _nodes[index];
        assert(n->refCount() < (1 << 27));
        ns.ref = n->refCount();

        if (n->isVar()) {
            ns.var  = true;
            ns.nand = false;
        } else {
            ns.var  = false;
            ns.nand = n->isNAND();
            ns.neg1 = n->parent1().isInverted();
            ns.neg2 = n->parent2().isInverted();

            assert(_id.find(n->parent1().node()) != _id.end());
            assert(_id.find(n->parent2().node()) != _id.end());

            ns.p1 = _id[n->parent1().node()];
            ns.p2 = _id[n->parent2().node()];

            std::pair<int, int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0), 2 * ns.p2 + (ns.neg2 ? 1 : 0));
            if (ip.first > ip.second) {
                std::swap(ip.first, ip.second);
            }

            int andindex = 2 * index + (ns.nand ? 1 : 0);

            _andMap[ip] = andindex;
        }
    }
    assert(index == totalNodes);

    int estimatedDelta = 0;

    for (unsigned int nodeindex = 0; nodeindex != totalNodes; ++nodeindex) {
        if (_rewritten[nodeindex] != 0) {
            continue;
        }

        NodeStruct& node = _nodes[nodeindex];
        if (node.ref == 0) {
            continue;
        }

        CutStruct* cuts = computeCuts(nodeindex);

        if (node.var) {
            continue;
        }
        if (_nodes[node.p1].ref >= 2 && _nodes[node.p2].ref >= 2) {
            continue;
        }

        _totalNodesChecked += 1;

        NPNImplInfo bestimpl = getBestImplementation(nodeindex, cuts);

        if (bestimpl.cutIndex != -1) {
            _totalNodesRewrittenLocal += 1;

            estimatedDelta += bestimpl.delta;

            if (bestimpl.npnclass.npnClass == 0) {
                buildNewConeConstant(nodeindex, isOutputNegated(bestimpl.npnclass.t2)
                                                    != DynamicRewritingTable::isInverted(bestimpl.implIndex));
            } else {
                /* find selected cut int cutlist */
                CutStruct* cutit = cuts;
                for (int i = 0; i != bestimpl.cutIndex; ++i) {
                    cutit = cutit->next;
                }

                buildNewCone(nodeindex, *cutit, bestimpl);
            }
        }
    }

    auto realNodesBefore = static_cast<int>(_manager->nodeCount());

    /*** build new cone ***/
    std::vector<bool>            processed(_nodesUsed, false);
    std::vector<InternalEdgeRef> cache(_nodesUsed);

    /* initialize cache with variables */
    for (std::vector<Node*>::const_iterator it = cone.begin(); it != cone.end(); ++it) {
        Node* n = *it;

        if (n->isVar()) {
            int index        = _id[n];
            processed[index] = true;
            cache[index]     = InternalEdgeRef(n);
        }
    }

    std::stack<int> pending;
    pending.push(_id[iroot.node()]);

    while (!pending.empty()) {
        int n = pending.top();

        if (processed[n]) {
            pending.pop();
            continue;
        }

        const NodeStruct& ns = _nodes[n];
        assert(!ns.var);

        if (!processed[ns.p1]) {
            pending.push(ns.p1);
        } else if (!processed[ns.p2]) {
            pending.push(ns.p2);
        } else {
            InternalEdgeRef p1 = cache[ns.p1].notIf(ns.neg1);
            InternalEdgeRef p2 = cache[ns.p2].notIf(ns.neg2);

            InternalEdgeRef result = (p1 & p2).notIf(ns.nand);

            cache[n]     = result;
            processed[n] = true;
            pending.pop();
        }
    }

    InternalEdgeRef iresult     = cache[_id[iroot.node()]].notIf(iroot.isInverted());
    auto            newConeSize = static_cast<int>(iresult.nodeCount());
    int             deltaN      = newConeSize - static_cast<int>(cone.size());

    cache.clear();

    auto realNodesAfter = static_cast<int>(_manager->nodeCount());

    _totalEstimatedDeltaLocal += estimatedDelta;
    _totalDeltaLocal += deltaN;
    double deltaT = (lrabs::cpuTime() - startTime);
    _totalRewritingTimeLocal += deltaT;

    std::cout << "local_rw time=" << deltaT << " delta=" << deltaN
              << " globaldelta=" << (realNodesAfter - realNodesBefore) << std::endl;

    return _manager->DebugGetExternal(iresult);
}

void
aigpp::RewritingManager::printStats(std::ostream& os) const
{
    os << "rw_applications " << _totalApplications << "\n"
       << "rw_nodeschecked " << _totalNodesChecked << "\n"
       << "rw_cutscomputed " << _totalCutsComputed << "\n"
       << "rw_implschecked " << _totalImplementationsChecked << "\n"
       << "rw_rewritten    " << _totalNodesRewritten << "\n"
       << "rw_estdelta     " << _totalEstimatedDelta << "\n"
       << "rw_delta        " << _totalDelta << "\n"
       << "rw_time         " << _totalRewritingTime << "\n"
       << "rw_local_applications " << _totalApplicationsLocal << "\n"
       << "rw_local_rewritten    " << _totalNodesRewrittenLocal << "\n"
       << "rw_local_estdelta     " << _totalEstimatedDeltaLocal << "\n"
       << "rw_local_delta        " << _totalDeltaLocal << "\n"
       << "rw_local_time         " << _totalRewritingTimeLocal << "\n"
       << std::flush;
}

int
aigpp::RewritingManager::lookupAnd(unsigned int p1, unsigned int p2) const
{
    if (p1 > p2) {
        std::swap(p1, p2);
    }

    std::pair<unsigned int, unsigned int> ip(p1, p2);

    auto p = _andMap.find(ip);
    if (p != _andMap.end()) {
        return p->second;
    } else {
        return -1;
    }
}

void
aigpp::RewritingManager::growArrays()
{
    unsigned int newcapacity = 2 * _nodesCapacity;

    auto* new_nodes = new NodeStruct[newcapacity];
    for (unsigned int i = 0; i != _nodesCapacity; ++i) {
        std::swap(_nodes[i], new_nodes[i]);
    }
    delete[] _nodes;
    _nodes = new_nodes;

    auto* new_rewritten = new unsigned int[newcapacity];
    std::fill(new_rewritten, new_rewritten + newcapacity, 0);
    std::copy(_rewritten, _rewritten + _nodesCapacity, new_rewritten);
    delete[] _rewritten;
    _rewritten = new_rewritten;

    auto* new_modified = new unsigned int[newcapacity];
    std::fill(new_modified, new_modified + newcapacity, 0);
    std::copy(_modified, _modified + _nodesCapacity, new_modified);
    delete[] _modified;
    _modified = new_modified;

    auto* new_sim = new unsigned int[newcapacity];
    std::copy(_sim, _sim + _nodesCapacity, new_sim);
    delete[] _sim;
    _sim = new_sim;

    auto* new_hasSim = new unsigned int[newcapacity];
    std::fill(new_hasSim, new_hasSim + newcapacity, 0);
    std::copy(_hasSim, _hasSim + _nodesCapacity, new_hasSim);
    delete[] _hasSim;
    _hasSim = new_hasSim;

    auto** new_cuts = new CutStruct*[newcapacity];
    std::fill(new_cuts, new_cuts + newcapacity, (CutStruct*)nullptr);
    std::copy(_cuts, _cuts + _nodesCapacity, new_cuts);
    delete[] _cuts;
    _cuts = new_cuts;

    _markedNodes.reserve(newcapacity);
    for (unsigned int i = _nodesCapacity; i != newcapacity; ++i) {
        _markedNodes.push_back(0);
    }

    _nodesCapacity = newcapacity;
}

void
aigpp::RewritingManager::refNode(unsigned int nodeindex)
{
    NodeStruct& ns = _nodes[nodeindex];
    ++ns.ref;
}

void
aigpp::RewritingManager::derefNode(unsigned int nodeindex)
{
    --_nodes[nodeindex].ref;
    if (_nodes[nodeindex].ref > 0) {
        return;
    }

    static lrabs::SimpleStack<unsigned int> pending;
    pending.push(nodeindex);

    while (!pending.empty()) {
        const NodeStruct& ns = _nodes[pending.top()];
        pending.pop();

        assert(ns.ref == 0);

        if (!ns.var) {
            ++_deadNodes;

            --_nodes[ns.p1].ref;
            --_nodes[ns.p2].ref;

            if (_nodes[ns.p1].ref == 0) {
                pending.push(ns.p1);
            }
            if (_nodes[ns.p2].ref == 0) {
                pending.push(ns.p2);
            }
        }
    }
}

void
aigpp::RewritingManager::CutStruct::print(std::ostream& os) const
{
    for (unsigned int i = 0; i != size; ++i) {
        os << ((i > 0) ? " " : "") << nodes[i];
    }
}

bool
aigpp::RewritingManager::CutStruct::isSubsetOf(const aigpp::RewritingManager::CutStruct& bigger) const
{
    /* dont check size, signature, etc.! */
    unsigned int i1 = 0;
    unsigned int i2 = 0;

    while (i1 != size && i2 != bigger.size) {
        if (nodes[i1] == bigger.nodes[i2]) {
            ++i1;
            ++i2;
        } else if (nodes[i1] < bigger.nodes[i2]) {
            return false;
        } else {
            ++i2;
        }
    }

    return (i1 == size);
}

aigpp::RewritingManager::CutStruct*
aigpp::RewritingManager::computeCuts(unsigned int nodeindex)
{
    const NodeStruct& node = _nodes[nodeindex];

    assert(_cuts[nodeindex] == 0);
    CutStruct* cuts = nullptr;

    CutStruct* ucut = getNewCut();

    if (!node.var) {
        if (_cuts[node.p1] == nullptr) {
            computeCuts(node.p1);
        }

        if (_cuts[node.p2] == nullptr) {
            computeCuts(node.p2);
        }

        CutStruct* cuts1 = _cuts[node.p1];
        CutStruct* cuts2 = _cuts[node.p2];

        for (CutStruct* c1 = cuts1; c1 != nullptr; c1 = c1->next) {
            for (CutStruct* c2 = cuts2; c2 != nullptr; c2 = c2->next) {
                assert(c1->size <= 4 && c2->size <= 4);

                /* signature check */
                ucut->signature = c1->signature | c2->signature;
                if (!lessEqual4(ucut->signature)) {
                    continue;
                }

                /* compute union */
                ucut->size = 0;
                {
                    unsigned int i1 = 0;
                    unsigned int i2 = 0;
                    while (i1 < c1->size && i2 < c2->size) {
                        if (c1->nodes[i1] == c2->nodes[i2]) {
                            ucut->nodes[ucut->size] = c1->nodes[i1];
                            ++i1;
                            ++i2;
                            ++ucut->size;
                        } else if (c1->nodes[i1] < c2->nodes[i2]) {
                            ucut->nodes[ucut->size] = c1->nodes[i1];
                            ++i1;
                            ++ucut->size;
                        } else {
                            ucut->nodes[ucut->size] = c2->nodes[i2];
                            ++i2;
                            ++ucut->size;
                        }
                    }

                    if (ucut->size > 4) {
                        continue;
                    }

                    while (i1 < c1->size) {
                        ucut->nodes[ucut->size] = c1->nodes[i1];
                        ++i1;
                        ++ucut->size;
                    }

                    while (i2 < c2->size) {
                        ucut->nodes[ucut->size] = c2->nodes[i2];
                        ++i2;
                        ++ucut->size;
                    }
                }

                /* create cut object if number of nodes <=4 */
                if (ucut->size <= 4) {
                    if (insertCut(&cuts, ucut)) {
                        ucut = getNewCut();
                    }
                }
            }
        }
    }

    /* add singleton cut */
    ucut->signature = (unsigned int)1 << (nodeindex % 32);
    ucut->size      = 1;
    ucut->nodes[0]  = nodeindex;

    ucut->next = cuts;
    cuts       = ucut;

    _cuts[nodeindex] = cuts;

    return cuts;
}

bool
aigpp::RewritingManager::insertCut(aigpp::RewritingManager::CutStruct** cuts,
                                   aigpp::RewritingManager::CutStruct*  newcut)
{
    /* check if newcut is superset of some existing cut */
    for (CutStruct* c = *cuts; c != nullptr; c = c->next) {
        if (c->size <= newcut->size && (c->signature & newcut->signature) == c->signature) {
            if (c->isSubsetOf(*newcut)) {
                return false;
            }
        }
    }

    /* check if newcut is a subset of some existing cuts */
    CutStruct* last = nullptr;
    for (CutStruct* c = *cuts; c != nullptr; /**/) {
        if (c->size >= newcut->size && (c->signature & newcut->signature) == newcut->signature) {
            if (newcut->isSubsetOf(*c)) {
                if (last == nullptr) {
                    *cuts = c->next;
                    delete c;
                    c = *cuts;
                } else {
                    last->next = c->next;
                    delete c;
                    c = last->next;
                }
            } else {
                last = c;
                c    = c->next;
            }
        } else {
            last = c;
            c    = c->next;
        }
    }

    /* add new cut */
    newcut->next = *cuts;
    *cuts        = newcut;

    return true;
}

aigpp::FourInputFunction
aigpp::RewritingManager::getFunction(unsigned int nodeindex, const aigpp::RewritingManager::CutStruct& cut)
{
    static lrabs::SimpleStack<unsigned int> processed(_nodesCapacity);
    static unsigned int                     tempSims[] = {0xAAAA, 0xCCCC, 0xF0F0, 0xFF00};

    for (unsigned int i = 0; i != cut.size; ++i) {
        _hasSim[cut.nodes[i]] = 1;
        processed.push(cut.nodes[i]);
        _sim[cut.nodes[i]] = (_nodes[cut.nodes[i]].nand ? ~tempSims[i] : tempSims[i]);
    }

    assert(_pending.empty());
    _pending.push(nodeindex);

    while (!_pending.empty()) {
        if (_hasSim[_pending.top()] != 0) {
            _pending.pop();
            continue;
        }

        const NodeStruct& node = _nodes[_pending.top()];

        /* push unprocessed parents */
        if (_hasSim[node.p1] == 0) {
            _pending.push(node.p1);
        } else if (_hasSim[node.p2] == 0) {
            _pending.push(node.p2);
        }
        /* process "and" node */
        else {
            /* get (inverted) parent sims */
            unsigned int s1 = (node.neg1 ? ~_sim[node.p1] : _sim[node.p1]);
            unsigned int s2 = (node.neg2 ? ~_sim[node.p2] : _sim[node.p2]);

            /* build (inverted) and of sims */
            s1 &= s2;
            _sim[_pending.top()] = (node.nand ? ~s1 : s1);
            processed.push(_pending.top());
            _hasSim[_pending.top()] = 1;

            _pending.pop();
        }
    }

    assert(_hasSim[nodeindex] != 0);

    for (unsigned int i = 0; i != processed.size(); ++i) {
        _hasSim[processed[i]] = 0;
    }
    processed.clear();

    return _sim[nodeindex] & 0xFFFF;
}

const aigpp::NPNClassAndTransformations&
aigpp::RewritingManager::getNPNClass(const aigpp::FourInputFunction& f) const
{
    return _staticF2N[f];
}

aigpp::RewritingManager::NPNImplInfo
aigpp::RewritingManager::getBestImplementation(unsigned int nodeindex, CutStruct* cuts)
{
    NPNImplInfo bestImpl;
    bestImpl.cutIndex = -1;

    if (_markedNodes.size() < _nodesCapacity) {
        _mark = 1;
        _markedNodes.assign(_nodesCapacity, 0);
    }

    int cutindex = 0;
    for (CutStruct* cut = cuts; cut != nullptr; cut = cut->next, ++cutindex) {
        /* skip single input cuts */
        if (cut->size <= 1) {
            continue;
        }

        /* get function table */
        FourInputFunction f = getFunction(nodeindex, *cut);

        /* get npn class */
        const NPNClassAndTransformations& npnclass = _staticF2N[f];

        /* special handling for constant function cuts */
        if (npnclass.npnClass == 0) {
            bestImpl.cutIndex  = cutindex;
            bestImpl.npnclass  = npnclass;
            bestImpl.implIndex = 0;
            bestImpl.delta     = -countDeletedNodesConstant(nodeindex) + 1;

            continue;
        }

#if 0
        /* add new implemenation if it is better than any known implementations
           continue, if new impl was added */
        if( addNewImplementation( nodeindex, *cut, npnclass ) ) continue;
#endif

        /* get implementations */
        const NPNClassImplementation& impls = _npnClassImpl[npnclass.npnClass];

        /* initialize implementation mapping */
        /* if >= 0: index of matching node,
           if == -1000000000: uninitialized,
           if < 0: number of needed nodes */

        for (unsigned int i = 0; i != _impl2nodeUsed.size(); ++i) {
            _impl2node[_impl2nodeUsed[i]] = -1000000000;
        }
        _impl2nodeUsed.clear();

        if (impls.usedInputs < cut->size) {
            selectBestImplSmallerInput(nodeindex, *cut, cutindex, npnclass, bestImpl);
        } else if (impls.usedInputs == cut->size) {
            selectBestImpl(nodeindex, *cut, cutindex, npnclass, bestImpl);
        } else {
            /* QUIET version */
#if 0
            std::cout << "minimum implementation with more inputs" << std::endl;
#endif
        }
    }

    return bestImpl;
}

int
aigpp::RewritingManager::countNewNodes(int implindex)
{
    /* return if already computed */
    if (_impl2node[implindex] != -1000000000) {
        return _impl2node[implindex];
    }

    assert(_pending.empty());
    _pending.push(implindex);

    while (!_pending.empty()) {
        if (_impl2node[_pending.top()] == -1000000000) {
            /* no input! */
            assert(_pending.top() >= 4);

            const DynamicRewritingTable::IndexPair& ip = _rwTable.getAnd(_pending.top() * 2);

            int needed1 = _impl2node[ip.first / 2];
            if (needed1 == -1000000000) {
                _pending.push(ip.first / 2);
                continue;
            }

            int needed2 = _impl2node[ip.second / 2];
            if (needed2 == -1000000000) {
                _pending.push(ip.second / 2);
                continue;
            }

            /* if one of the parents doesn't exist -> return and report number of
             * needed new nodes */
            if (needed1 < 0 || needed2 < 0) {
                /* need at least on node (for the and) */
                int counter = -1;
                if (needed1 < 0) {
                    counter += needed1;
                }
                if (needed2 < 0) {
                    counter += needed2;
                }

                _impl2node[_pending.top()] = counter;
                _impl2nodeUsed.push(_pending.top());
            } else {
                int andindex = lookupAnd(needed1 ^ (ip.first % 2), needed2 ^ (ip.second % 2));
                if (andindex <= 0) {
                    _impl2node[_pending.top()] = -1;
                } else if (_nodes[andindex / 2].ref > 0) {
                    _impl2node[_pending.top()] = andindex;
                } else {
                    _impl2node[_pending.top()] = -1000;
                }

                _impl2nodeUsed.push(_pending.top());
            }
        }

        _pending.pop();
    }

    return _impl2node[implindex];
}

int
aigpp::RewritingManager::countDeletedNodes(unsigned int nodeindex, int implindex)
{
    assert(_impl2node[implindex] > -1000000000);
    assert(!_nodes[nodeindex].var);

    assert(_pending.empty());

    /* mark nodes needed by implemenation graph */
    _pending.push(implindex);
    while (!_pending.empty()) {
        if (_impl2node[_pending.top()] >= 0) {
            markNode(_impl2node[_pending.top()] / 2);
            _pending.pop();
        } else {
            assert(_pending.top() >= 4);

            const DynamicRewritingTable::IndexPair& ip = _rwTable.getAnd(_pending.top() * 2);
            _pending.pop();

            _pending.push(ip.first / 2);
            _pending.push(ip.second / 2);
        }
    }

    if (isNodeMarked(nodeindex)) {
        clearMarks();

        /* QUIET version */
#if 0
        std::cout << __func__ << ": node rewritten by itself" << std::endl;
#endif
        return 0;
    }

    int deletedNodes = 0;

    static std::vector<int>        refcount(_nodesCapacity, -1);
    static lrabs::SimpleStack<int> usedrefs(_nodesCapacity);

    if (refcount.size() < _nodesCapacity) {
        refcount.assign(_nodesCapacity, -1);
    }

    usedrefs.push(nodeindex);
    refcount[nodeindex] = 1;

    _pending.push(nodeindex);
    while (!_pending.empty()) {
        int n = _pending.top();
        _pending.pop();

        int r;

        /* not initialized, yet */
        if (refcount[n] < 0) {
            r = _nodes[n].ref;

            /* node is not needed by implemenation graph -> decrement its refcount */
            if (!isNodeMarked(n)) {
                --r;
            }

            refcount[n] = r;
            usedrefs.push(n);
        } else {
            r = --(refcount[n]);
        }

        if (r > 0) {
            continue;
        }

        assert(r == 0);
        assert(!_nodes[n].var);

        ++deletedNodes;

        _pending.push(_nodes[n].p1);
        _pending.push(_nodes[n].p2);
    }

    for (unsigned int i = 0; i != usedrefs.size(); ++i) {
        refcount[usedrefs[i]] = -1;
    }
    usedrefs.clear();

    clearMarks();

    return deletedNodes;
}

int
aigpp::RewritingManager::countDeletedNodesSmallerInput(unsigned int nodeindex, int implindex)
{
    assert(_impl2node[implindex] > -1000000000);
    assert(!_nodes[nodeindex].var);

    assert(_pending.empty());
    static lrabs::SimpleStack<unsigned int> marked;
    marked.clear();

    /* mark nodes needed by implemenation graph */
    _pending.push(implindex);
    while (!_pending.empty()) {
        if (_impl2node[_pending.top()] >= 0) {
            markNode(_impl2node[_pending.top()] / 2);
            marked.push(_impl2node[_pending.top()] / 2);
            _pending.pop();
        } else {
            assert(_pending.top() >= 4);

            const DynamicRewritingTable::IndexPair& ip = _rwTable.getAnd(_pending.top() * 2);
            _pending.pop();

            _pending.push(ip.first / 2);
            _pending.push(ip.second / 2);
        }
    }

    if (isNodeMarked(nodeindex)) {
        clearMarks();

        /* QUIET version */
#if 0
        std::cout << __func__ << ": node rewritten by itself" << std::endl;
#endif
        return 0;
    }

    /* push marks towards inputs */
    for (unsigned int i = 0; i != marked.size(); ++i) {
        if (!_nodes[marked[i]].var) {
            _pending.push(_nodes[marked[i]].p1);
            _pending.push(_nodes[marked[i]].p2);
        }

        while (!_pending.empty()) {
            unsigned int n = _pending.top();
            _pending.pop();

            if (!isNodeMarked(n)) {
                markNode(n);

                if (!_nodes[n].var) {
                    _pending.push(_nodes[n].p1);
                    _pending.push(_nodes[n].p2);
                }
            }
        }
    }

    int deletedNodes = 0;

    static std::vector<int>        refcount(_nodesCapacity, -1);
    static lrabs::SimpleStack<int> usedrefs(_nodesCapacity);

    if (refcount.size() < _nodesCapacity) {
        refcount.assign(_nodesCapacity, -1);
    }

    usedrefs.push(nodeindex);
    refcount[nodeindex] = 1;

    _pending.push(nodeindex);
    while (!_pending.empty()) {
        unsigned int n = _pending.top();
        int          r;

        if (refcount[n] < 0) {
            r = _nodes[n].ref;

            /* node is not needed by implemenation graph -> decrement its refcount */
            if (!isNodeMarked(n)) {
                --r;
            }

            refcount[n] = r;
            usedrefs.push(n);
        } else {
            r = --(refcount[n]);
        }

        assert(r >= 0);

        _pending.pop();

        if (r > 0) {
            continue;
        }

        ++deletedNodes;

        assert(!_nodes[n].var);

        _pending.push(_nodes[n].p1);
        _pending.push(_nodes[n].p2);
    }

    for (unsigned int i = 0; i != usedrefs.size(); ++i) {
        refcount[usedrefs[i]] = -1;
    }
    usedrefs.clear();

    clearMarks();

    return deletedNodes;
}

int
aigpp::RewritingManager::countDeletedNodesConstant(unsigned int nodeindex)
{
    assert(!_nodes[nodeindex].var);

    assert(_pending.empty());

    int deletedNodes = 0;

    static std::vector<int>        refcount(_nodesCapacity, -1);
    static lrabs::SimpleStack<int> usedrefs(_nodesCapacity);

    if (refcount.size() < _nodesCapacity) {
        refcount.assign(_nodesCapacity, -1);
    }

    usedrefs.push(nodeindex);
    refcount[nodeindex] = 1;

    _pending.push(nodeindex);
    while (!_pending.empty()) {
        unsigned int n = _pending.top();
        int          r;

        if (refcount[n] < 0) {
            r = _nodes[n].ref - 1;

            refcount[n] = r;
            usedrefs.push(n);
        } else {
            r = --(refcount[n]);
        }

        assert(r >= 0);

        _pending.pop();

        if (r > 0) {
            continue;
        }

        ++deletedNodes;

        assert(!_nodes[n].var);

        if (!_nodes[_nodes[n].p1].var) {
            _pending.push(_nodes[n].p1);
        }
        if (!_nodes[_nodes[n].p2].var) {
            _pending.push(_nodes[n].p2);
        }
    }

    for (unsigned int i = 0; i != usedrefs.size(); ++i) {
        refcount[usedrefs[i]] = -1;
    }
    usedrefs.clear();

    return deletedNodes;
}

void
aigpp::RewritingManager::markRewritten(unsigned int nodeindex)
{
    assert(_pending.empty());

    _pending.push(nodeindex);
    while (!_pending.empty()) {
        if (_rewritten[_pending.top()]) {
            _pending.pop();
            continue;
        }

        const NodeStruct& ns = _nodes[_pending.top()];

        _rewritten[_pending.top()] = 1;

        _pending.pop();

        if (!ns.var) {
            _pending.push(ns.p1);
            _pending.push(ns.p2);
        }
    }
}

void
aigpp::RewritingManager::buildNewCone(unsigned int nodeindex, const aigpp::RewritingManager::CutStruct& cut,
                                      const aigpp::RewritingManager::NPNImplInfo& bestImpl)
{
    assert(bestImpl.npnclass.npnClass != 0);

    for (unsigned int i = 0; i != _impl2nodeUsed.size(); ++i) {
        _impl2node[_impl2nodeUsed[i]] = -1000000000;
    }
    _impl2nodeUsed.clear();

    for (unsigned int implinput = 0; implinput != 4; ++implinput) {
        unsigned int cutindex = getInputPermutation(bestImpl.npnclass.t2, implinput);
        if (cutindex >= cut.size) {
            continue;
        }

        /* get cut node with index "cutindex" */
        unsigned int cutnode = cut.nodes[cutindex];

        if (_nodes[cutnode].nand != isInputNegated(bestImpl.npnclass.t2, implinput)) {
            cutnode = 1 + cutnode * 2;
        } else {
            cutnode = cutnode * 2;
        }

        _impl2node[implinput] = cutnode;
        _impl2nodeUsed.push(implinput);
    }

    assert(_pending.empty());
    _pending.push(bestImpl.implIndex / 2);

    static lrabs::SimpleStack<unsigned int> toMarkRewritten;

    while (!_pending.empty()) {
        if (_impl2node[_pending.top()] < 0) {
            const DynamicRewritingTable::IndexPair& ip = _rwTable.getAnd(_pending.top() * 2);

            int needed1 = _impl2node[ip.first / 2];
            if (needed1 < 0) {
                _pending.push(ip.first / 2);
                continue;
            }
            needed1 ^= (ip.first % 2);

            int needed2 = _impl2node[ip.second / 2];
            if (needed2 < 0) {
                _pending.push(ip.second / 2);
                continue;
            }
            needed2 ^= (ip.second % 2);

            /* not top node */
            if (_pending.size() != 1) {
                int andnode = lookupAnd(needed1, needed2);

                /* not found */
                if (andnode < 0) {
                    if (_nodesUsed == _nodesCapacity) {
                        growArrays();
                    }

                    NodeStruct& ns        = _nodes[_nodesUsed];
                    _modified[_nodesUsed] = 1;
                    andnode               = 2 * _nodesUsed;
                    ++_nodesUsed;

                    ns.var  = false;
                    ns.nand = false;
                    ns.neg1 = ((unsigned int)needed1 % 2) != 0;
                    ns.neg2 = ((unsigned int)needed2 % 2) != 0;
                    ns.ref  = 0;
                    ns.p1   = needed1 / 2;
                    ns.p2   = needed2 / 2;

                    refNode(ns.p1);
                    refNode(ns.p2);

                    _impl2node[_pending.top()] = andnode;
                    _impl2nodeUsed.push(_pending.top());

                    /* insert new node into mapping */
                    std::pair<unsigned int, unsigned int> ip((unsigned int)needed1, (unsigned int)needed2);
                    if (ip.first > ip.second) {
                        std::swap(ip.first, ip.second);
                    }
                    _andMap[ip] = andnode;

                    computeCuts(andnode / 2);

                    toMarkRewritten.push(andnode / 2);
                }
                /* matching node found */
                else {
#ifndef NDEBUG
                    const NodeStruct& ns = _nodes[andnode / 2];
                    assert(andnode % 2 == ns.nand);
#endif
                    _impl2node[_pending.top()] = andnode;
                    _impl2nodeUsed.push(_pending.top());

                    toMarkRewritten.push(andnode / 2);
                }
            }
            /* top node */
            else {
                NodeStruct& ns = _nodes[nodeindex];
                assert(!ns.var);

                /*
                int andnode = lookupAnd( needed1, needed2 );
                if( andnode >= 0 )
                {
                    std::cout << "INFO: non-top node already exists" << std::endl;
                }
                */

#if 0
                bool targetinv = ( isOutputNegated( bestImpl.npnclass.t2 ) != DynamicRewritingTable::isInverted( bestImpl.implIndex ) );
                checkEquivalence( 2 * nodeindex + ( ( targetinv != ns.nand ) ? 1 : 0 ), needed1, needed2 );
#endif

                refNode(needed1 / 2);
                refNode(needed2 / 2);

                _modified[nodeindex] = 1;

                /* erase old representation from _andMap */
                std::pair<unsigned int, unsigned int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0), 2 * ns.p2 + (ns.neg2 ? 1 : 0));
                if (ip.first > ip.second) {
                    std::swap(ip.first, ip.second);
                }
                _andMap.erase(ip);

                derefNode(ns.p1);
                derefNode(ns.p2);

                ns.neg1 = ((unsigned int)needed1 & 1) != 0;
                ns.neg2 = ((unsigned int)needed2 & 1) != 0;

                ns.p1 = needed1 / 2;
                ns.p2 = needed2 / 2;

                ns.nand
                    = isOutputNegated(bestImpl.npnclass.t2) != DynamicRewritingTable::isInverted(bestImpl.implIndex);

                /* add new representation */
                std::pair<unsigned int, unsigned int> newip((unsigned int)needed1, (unsigned int)needed2);
                if (newip.first > newip.second) {
                    std::swap(newip.first, newip.second);
                }
                _andMap[newip] = 2 * nodeindex + (ns.nand ? 1 : 0);

                releaseCuts(_cuts[nodeindex]);
                _cuts[nodeindex] = nullptr;
                computeCuts(nodeindex);

                toMarkRewritten.push(nodeindex);
            }
        } else {
            /* top node implementation already exists! */

            /* three possibilities:
             * 1. we built the input cone itself
             *    -> just skip replacement
             * 2. the input cone is equivalent to a cut node c
             *    -> replace cone by the cut node (or by "c & c")
             * 3. the input cone is equivalent to some existing cone, and we
                  found it
             *    -> replace old cone by new cone since we have already
                     proven, that this replacement will lead to a node decrease
             */
            if (_pending.size() == 1) {
                int i = _impl2node[_pending.top()];
                toMarkRewritten.push(i / 2);
                // std::cout << "top node known " << nodeindex << " " << i << std::endl;

                /* case 1: built cone itself ( i/2 == nodeindex ) */
                if (i / 2 == (int)nodeindex) {
                    /* QUIET version */
#if 0
                    std::cout << "INFO: built cone itself -> no rewriting"
                              << std::endl;
#endif
                }
                /* case 2: build cone input ( i < 8 ) */
                else if (i < 8) {
                    /* QUIET version */
#if 0
                    std::cout << "INFO: built existing node (input)" << std::endl;
#endif
                    refNode(i / 2);
                    refNode(i / 2);

                    NodeStruct& ns = _nodes[nodeindex];
                    assert(!ns.var);
                    _modified[nodeindex] = 1;

                    /* erase old representation from _andMap */
                    std::pair<unsigned int, unsigned int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0),
                                                             2 * ns.p2 + (ns.neg2 ? 1 : 0));
                    if (ip.first > ip.second) {
                        std::swap(ip.first, ip.second);
                    }
                    _andMap.erase(ip);

                    derefNode(ns.p1);
                    derefNode(ns.p2);

                    ns.neg1 = ((unsigned int)i & 1) != 0;
                    ns.neg2 = ns.neg1;

                    ns.p1 = i / 2;
                    ns.p2 = ns.p1;

                    ns.nand = isOutputNegated(bestImpl.npnclass.t2)
                              != DynamicRewritingTable::isInverted(bestImpl.implIndex);

                    /* add new representation */
                    std::pair<unsigned int, unsigned int> newip((unsigned int)i, (unsigned int)i);
                    _andMap[newip] = 2 * nodeindex + (ns.nand ? 1 : 0);

                    releaseCuts(_cuts[nodeindex]);
                    _cuts[nodeindex] = nullptr;
                    computeCuts(nodeindex);
                }
                /* case 3: build equivalent, existing cone which is smaller */
#if 0
                else
                {
    /* QUIET version */
#    if 0
                    std::cout << "INFO: built existing node" << std::endl;
#    endif
                    /* get existing nodestruct */
                    NodeStruct& nsold = _nodes[i/2];
                    if( nsold.var )
                    {
    /* QUIET version */
#    if 0
                        std::cout << "var: " << nsold.p1 << std::endl;
#    endif
                    }

                    assert( !nsold.var );

                    /* reference old node's parents */
                    refNode( nsold.p1 );
                    refNode( nsold.p2 );

                    /* get current nodestruct */
                    NodeStruct& ns = _nodes[nodeindex];
                    assert( !ns.var );
                    _modified[nodeindex] = 1;

                    /* erase old representation from _andMap */
                    std::pair<unsigned int, unsigned int> ip( 2 * ns.p1 + ( ns.neg1 ? 1 : 0 ) , 2 * ns.p2 + ( ns.neg2 ? 1 : 0 ) );
                    if( ip.first > ip.second ) std::swap( ip.first, ip.second );
                    _andMap.erase( ip );

                    derefNode( ns.p1 );
                    derefNode( ns.p2 );

                    ns.neg1 = nsold.neg1;
                    ns.neg2 = nsold.neg2;

                    ns.p1 = nsold.p1;
                    ns.p2 = nsold.p1;

                    ns.nand = isOutputNegated( bestImpl.npnclass.t2 ) != DynamicRewritingTable::isInverted( bestImpl.implIndex );

                    releaseCuts( _cuts[nodeindex] );
                    _cuts[nodeindex] = 0;
                    computeCuts( nodeindex );
                }
#endif
            }
        }

        _pending.pop();
    }

    markRewritten(nodeindex);

    for (unsigned int i = 0; i != toMarkRewritten.size(); ++i) {
        markRewritten(toMarkRewritten[i]);
    }
    toMarkRewritten.clear();
}

void
aigpp::RewritingManager::buildNewConeConstant(unsigned int nodeindex, bool invert)
{
    // std::cout << "[INFO] rewriting with constant" << std::endl;

    NodeStruct& ns = _nodes[nodeindex];
    assert(!ns.var);

    unsigned int someVariable = 0;
    refNode(someVariable);
    refNode(someVariable);

#if 0
    checkEquivalence( 2 * nodeindex + ( invert ? 1 : 0 ), 2 * someVariable, 2 * someVariable + 1 );
#endif

    _modified[nodeindex] = 1;

    /* erase old representation from _andMap */
    std::pair<unsigned int, unsigned int> ip(2 * ns.p1 + (ns.neg1 ? 1 : 0), 2 * ns.p2 + (ns.neg2 ? 1 : 0));
    if (ip.first > ip.second) {
        std::swap(ip.first, ip.second);
    }
    _andMap.erase(ip);

    derefNode(ns.p1);
    derefNode(ns.p2);

    ns.neg1 = false;
    ns.neg2 = true;

    ns.p1 = someVariable;
    ns.p2 = someVariable;

    ns.nand = invert;

    /* add new representation */
    std::pair<unsigned int, unsigned int> newip(2 * someVariable, 2 * someVariable + 1);
    if (newip.first > newip.second) {
        std::swap(newip.first, newip.second);
    }
    _andMap[newip] = 2 * nodeindex + (ns.nand ? 1 : 0);

    CutStruct* cut = getNewCut();
    cut->signature = (unsigned int)1 << (nodeindex % 32);
    cut->size      = 1;
    cut->nodes[0]  = nodeindex;
    cut->next      = nullptr;

    releaseCuts(_cuts[nodeindex]);
    _cuts[nodeindex] = cut;

    markRewritten(nodeindex);
}

void
aigpp::RewritingManager::loadStaticRWDB()
{
#if 0
    _modifiedRWTable = false;
#endif
    assert(_rwTable.size() == 0);

    for (int i = 0; i != _initialImplNodesCount; ++i) {
#ifndef NDEBUG
        const int nodeid = 2 * (i + 4);
#endif

        const int input1 = _initialImplNodes[2 * i];
        const int input2 = _initialImplNodes[2 * i + 1];

        const int realid = _rwTable.insertAnd(input1, input2);
        _unused(realid);

#ifndef NDEBUG
        assert(nodeid == realid);
#endif
    }

    /* load npn classes */
    const int* it = _initialNPNImplementations;
    for (unsigned int npnclass = 0; npnclass != NumberOfNPNClasses; ++npnclass) {
        int size = *it;
        assert(size == -1 || size >= 0);
        ++it;

        unsigned int usedMask = *it;
        assert(usedMask < (1 << 4));
        ++it;

        unsigned int count = *it;
        ++it;

        std::list<DynamicRewritingTable::Index> implementations;

        for (unsigned int i = 0; i != count; ++i) {
            int index = *it;
            ++it;

            _rwTable.ref(index);
            implementations.push_back(index);
        }

        _npnClassImpl[npnclass].size            = size;
        _npnClassImpl[npnclass].usedInputsMask  = usedMask;
        _npnClassImpl[npnclass].usedInputs      = lrabs::countOnes(usedMask);
        _npnClassImpl[npnclass].implementations = implementations;
    }

    /* setup impl2node lookup map */
    _impl2node.assign(_rwTable.size() + 4, -1000000000);
    _impl2nodeUsed.clear();
}

#if 0
void
aigpp::RewritingManager::
loadRWDB( std::string filename,
          std::vector<aigpp::NPNClassImplementation>& npn2impl,
          aigpp::DynamicRewritingTable& rwDB )
{
    assert( _rwTable.size() == 0 );

    std::ifstream file( filename.c_str() );
    if( !file )
    {
        std::cout << "WARNING: cannot open rewriting db '" << filename << "'" << std::endl;
    }
    else
    {

        std::string identifier;
        file >> identifier;

        if( identifier != "rwdb_2" )
        {
            std::cout << "WARNING: cannot load rewriting db '" << filename << "': wrong identifier '" << identifier << "', expected 'rwdb_2'" << std::endl;
            return;
        }

        _modifiedRWTable = false;

        int numberOfNodes = 0;
        file >> numberOfNodes;

        /* load implementation nodes */
        std::map<int, DynamicRewritingTable::Index> old2new;
        old2new[0] = 0;
        old2new[2] = 2;
        old2new[4] = 4;
        old2new[6] = 6;

        for( int i = 0; i != numberOfNodes; ++i )
        {
            int x, y, z;
            file >> x >> y >> z;

            assert( x % 2 == 0 );
            assert( old2new.find( x ) == old2new.end() );
            assert( old2new.find( y & ~1 ) != old2new.end() );
            assert( old2new.find( z & ~1 ) != old2new.end() );

            DynamicRewritingTable::Index yi = old2new[ y & ~1 ];
            if( ( y & 1 ) != 0 ) yi = DynamicRewritingTable::invert( yi );

            DynamicRewritingTable::Index zi = old2new[ z & ~1 ];
            if( ( z & 1 ) != 0 ) zi = DynamicRewritingTable::invert( zi );

            DynamicRewritingTable::Index xi = rwDB.lookupAnd( yi, zi );
            assert( xi == DynamicRewritingTable::InvalidIndex );

            xi = rwDB.insertAnd( yi, zi );
            old2new[ x ] = xi;
        }

        /* load npn classes */
        for( unsigned int npnclass = 0; npnclass != NumberOfNPNClasses; ++npnclass )
        {
            int size = 0;
            file >> size;
            assert( size == -1 || size >= 0 );

            unsigned int inputsMask = 0;
            file >> inputsMask;
            assert( inputsMask < ( 1 << 4 ) );

            unsigned int count = 0;
            file >> count;

            std::list<DynamicRewritingTable::Index> implementations;

            for( unsigned int i = 0; i != count; ++i )
            {
                int oldindex = 0;
                file >> oldindex;

                assert( old2new.find( oldindex & ~1 ) != old2new.end() );

                DynamicRewritingTable::Index newindex = old2new[ oldindex & ~1 ];
                rwDB.ref( newindex );
                if( ( oldindex & 1 ) != 0 ) newindex = DynamicRewritingTable::invert( newindex );

                implementations.push_back( newindex );
            }

            npn2impl[npnclass].size = size;
            npn2impl[npnclass].usedInputsMask = inputsMask;
            npn2impl[npnclass].usedInputs = lrabs::countOnes( inputsMask );
            npn2impl[npnclass].implementations = implementations;
        }
    }


    /* check equivalences */
    for( unsigned int n = 0; n != NumberOfNPNClasses; ++n )
    {
        std::list<DynamicRewritingTable::Index>& implementations = npn2impl[n].implementations;
        if( implementations.empty() ) continue;

        FourInputFunction ffirst = rwDB.getFunction( *( implementations.begin() ) );

        for( std::list<DynamicRewritingTable::Index>::iterator impl = implementations.begin();
             impl != implementations.end(); ++impl )
        {
            FourInputFunction f2 = rwDB.getFunction( *impl );
            assert( ffirst == f2 );
            assert( _staticF2N[f2].npnClass == n );
        }
    }

    /* setup impl2node lookup map */
    _impl2node.assign( _rwTable.size() + 4, -1000000000 );
    _impl2nodeUsed.clear();
}

void
aigpp::RewritingManager::
saveRWDB( std::string filename,
          const std::vector<aigpp::NPNClassImplementation>& npn2impl,
          const aigpp::DynamicRewritingTable& rwDB ) const
{
    std::ofstream file( filename.c_str() );
    if( !file )
    {
        std::cout << "WARNING: cannot open rewriting db '" << filename << "'" << std::endl;
    }
    else
    {
        _modifiedRWTable = false;

        file << "rwdb_2\n";

        /* count used nodes */
        {
            std::stack<DynamicRewritingTable::Index> pending;
            std::set<DynamicRewritingTable::Index> counted;
            counted.insert( 0 );
            counted.insert( 2 );
            counted.insert( 4 );
            counted.insert( 6 );

            for( unsigned int npnclass = 0; npnclass != NumberOfNPNClasses; ++npnclass )
            {
                const NPNClassImplementation& impl = npn2impl[npnclass];
                for( std::list<DynamicRewritingTable::Index>::const_iterator i = impl.implementations.begin(); i != impl.implementations.end(); ++i )
                {
                    pending.push( DynamicRewritingTable::getNonInvertedIndex( *i ) );
                }
            }

            while( !pending.empty() )
            {
                if( counted.find( pending.top() ) != counted.end() )
                {
                    pending.pop();
                    continue;
                }

                assert( !DynamicRewritingTable::isInput( pending.top() ) );


                DynamicRewritingTable::IndexPair ip = rwDB.getAnd( pending.top() );

                if( counted.find( DynamicRewritingTable::getNonInvertedIndex( ip.first ) ) == counted.end() )
                {
                    pending.push( DynamicRewritingTable::getNonInvertedIndex( ip.first ) );
                    continue;
                }

                if( counted.find( DynamicRewritingTable::getNonInvertedIndex( ip.second ) ) == counted.end() )
                {
                    pending.push( DynamicRewritingTable::getNonInvertedIndex( ip.second ) );
                    continue;
                }

                counted.insert( DynamicRewritingTable::getNonInvertedIndex( pending.top() ) );
                pending.pop();
            }

            file << counted.size() - 4 << std::endl;
        }

        /* save implementation nodes */
        std::map<DynamicRewritingTable::Index, int> old2new;
        old2new[0] = 0;
        old2new[2] = 2;
        old2new[4] = 4;
        old2new[6] = 6;
        int nextIndex = 8;

        std::stack<DynamicRewritingTable::Index> pending;

        for( unsigned int npnclass = 0; npnclass != NumberOfNPNClasses; ++npnclass )
        {
            const NPNClassImplementation& impl = npn2impl[npnclass];
            for( std::list<DynamicRewritingTable::Index>::const_iterator i = impl.implementations.begin(); i != impl.implementations.end(); ++i )
            {
                pending.push( DynamicRewritingTable::getNonInvertedIndex( *i ) );
            }
        }

        while( !pending.empty() )
        {
            if( old2new.find( pending.top() ) != old2new.end() )
            {
                pending.pop();
                continue;
            }

            assert( !DynamicRewritingTable::isInput( pending.top() ) );

            DynamicRewritingTable::IndexPair ip = rwDB.getAnd( pending.top() );

            if( old2new.find( DynamicRewritingTable::getNonInvertedIndex( ip.first ) ) == old2new.end() )
            {
                pending.push( DynamicRewritingTable::getNonInvertedIndex( ip.first ) );
                continue;
            }

            if( old2new.find( DynamicRewritingTable::getNonInvertedIndex( ip.second ) ) == old2new.end() )
            {
                pending.push( DynamicRewritingTable::getNonInvertedIndex( ip.second ) );
                continue;
            }


            old2new[ pending.top() ] = nextIndex;
            pending.pop();

            file << nextIndex;
            nextIndex += 2;

            int i1 = old2new[ DynamicRewritingTable::getNonInvertedIndex( ip.first ) ];
            if( DynamicRewritingTable::isInverted( ip.first ) ) i1 = i1 | 1;
            file << " " << i1;

            int i2 = old2new[ DynamicRewritingTable::getNonInvertedIndex( ip.second ) ];
            if( DynamicRewritingTable::isInverted( ip.second ) ) i2 = i2 | 1;
            file << " " << i2;

            file << std::endl;
        }

        /* save npn classes */
        for( unsigned int npnclass = 0; npnclass != NumberOfNPNClasses; ++npnclass )
        {
            const NPNClassImplementation& impl = npn2impl[npnclass];
            file << impl.size;
            file << " " << impl.usedInputsMask;
            file << " " << impl.implementations.size();

            for( std::list<DynamicRewritingTable::Index>::const_iterator i = impl.implementations.begin(); i != impl.implementations.end(); ++i )
            {
                DynamicRewritingTable::Index oldindex = *i;
                assert( old2new.find( oldindex & ~1 ) != old2new.end() );

                int newindex = old2new[ oldindex & ~1 ];
                if( ( oldindex & 1 ) != 0 ) newindex |= 1;

                file << " " << newindex;
            }

            file << std::endl;
        }
    }
}
#endif

void
aigpp::RewritingManager::checkEquivalence(unsigned int nodeindex, unsigned int parent1, unsigned int parent2) const
{
    /*
    std::cout << "structure before equivalence check" << std::endl;
    print( std::cout );

    std::cout << "equivalent " << ( ( nodeindex % 2 == 0 ) ? "" : "!" ) <<
    nodeindex/2 << " "
              << ( ( ( parent1 % 2 ) == 0 ) ? "" : "!" ) << parent1/2
              << " and "
              << ( ( ( parent2 % 2 ) == 0 ) ? "" : "!" ) << parent2/2
              << std::endl;
    */

    Minisat::Solver           solver;
    std::vector<Minisat::Var> satVars(_nodesUsed, var_Undef);

    static lrabs::SimpleStack<unsigned int> pending;

    pending.push(nodeindex / 2);
    pending.push(parent1 / 2);
    pending.push(parent2 / 2);

    while (!pending.empty()) {
        if (satVars[pending.top()] != var_Undef) {
            pending.pop();
            continue;
        }

        const NodeStruct& ns = _nodes[pending.top()];
        if (ns.var) {
            satVars[pending.top()] = solver.newVar();
            pending.pop();
            continue;
        } else if (satVars[ns.p1] == var_Undef) {
            pending.push(ns.p1);
            continue;
        } else if (satVars[ns.p2] == var_Undef) {
            pending.push(ns.p2);
            continue;
        } else {
            satVars[pending.top()] = solver.newVar();

            Minisat::Lit lit1 = Minisat::mkLit(satVars[pending.top()], ns.nand);
            Minisat::Lit lit2 = Minisat::mkLit(satVars[ns.p1], ns.neg1);
            Minisat::Lit lit3 = Minisat::mkLit(satVars[ns.p2], ns.neg2);

            Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);
            clause1[0] = ~lit1;
            clause1[1] = lit2;

            clause2[0] = ~lit1;
            clause2[1] = lit3;

            clause3[0] = lit1;
            clause3[1] = ~lit2;
            clause3[2] = ~lit3;

            solver.addClause(clause1);
            solver.addClause(clause2);
            solver.addClause(clause3);

            pending.pop();
        }
    }

    Minisat::Lit lit1 = Minisat::mkLit(satVars[nodeindex / 2], (nodeindex & 1) == 0);
    Minisat::Lit lit2 = Minisat::mkLit(satVars[parent1 / 2], (parent1 & 1) != 0);
    Minisat::Lit lit3 = Minisat::mkLit(satVars[parent2 / 2], (parent2 & 1) != 0);

    Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

    clause1[0] = ~lit1;
    clause1[1] = lit2;

    clause2[0] = ~lit1;
    clause2[1] = lit3;

    clause3[0] = lit1;
    clause3[1] = ~lit2;
    clause3[2] = ~lit3;

    solver.addClause(clause1);
    solver.addClause(clause2);
    solver.addClause(clause3);

    if (solver.solve()) {
        std::cout << "not equivalent!" << std::endl;
        abort();
    }
}

unsigned int
aigpp::RewritingManager::usedNodes() const
{
    unsigned int used = 0;

    for (unsigned int i = 0; i != _nodesUsed; ++i) {
        if (_nodes[i].ref > 0) {
            ++used;
        }
    }

    return used;
}

void
aigpp::RewritingManager::print(std::ostream& os) const
{
    for (unsigned int i = 0; i != _nodesUsed; ++i) {
        const NodeStruct& ns = _nodes[i];

        if (ns.var) {
            os << i << ": variable "
               << " ref=" << ns.ref << std::endl;
        } else {
            os << i << ": " << (ns.nand ? "nand" : "and") << " " << (ns.neg1 ? "!" : "") << ns.p1 << " "
               << (ns.neg2 ? "!" : "") << ns.p2 << " ref=" << ns.ref << std::endl;
        }
    }
}

void
aigpp::RewritingManager::recursiveRebuild(unsigned int nodeindex, std::vector<bool>& processed,
                                          std::vector<aigpp::InternalEdgeRef>& cache)
{
    if (processed[nodeindex]) {
        return;
    }

    const NodeStruct& ns = _nodes[nodeindex];
    // assert( ns.ref > 0 );

    if (ns.var) {
        assert(nodeindex < _backMap.size());
        assert(_backMap[nodeindex].node()->isVar());
        cache[nodeindex]     = _backMap[nodeindex];
        processed[nodeindex] = true;
        return;
    }

    if (!processed[ns.p1]) {
        recursiveRebuild(ns.p1, processed, cache);
    }

    if (!processed[ns.p2]) {
        recursiveRebuild(ns.p2, processed, cache);
    }

    /* old node */
    if (nodeindex < _backMap.size()) {
        Node* node = _backMap[nodeindex].node();

        InternalEdgeRef p1 = cache[ns.p1].notIf(ns.neg1);
        InternalEdgeRef p2 = cache[ns.p2].notIf(ns.neg2);

        if (node->_parent1.structurallyEquivalent(p1) && node->_parent2.structurallyEquivalent(p2)) {
            assert(node->isNAND() == ns.nand);

            cache[nodeindex]     = InternalEdgeRef(node);
            processed[nodeindex] = true;
        } else {
            _manager->_unique.remove(node);

            if (ns.nand) {
                node->setFlag<Node::FLAG_ISNAND>();
            } else {
                node->unsetFlag<Node::FLAG_ISNAND>();
            }

            _manager->ref(p1.node());
            _manager->ref(p2.node());

            _manager->deref(node->_parent1.node());
            _manager->deref(node->_parent2.node());

            node->_parent1 = p1;
            node->_parent2 = p2;

            _manager->_unique.insert(node);

            cache[nodeindex]     = InternalEdgeRef(node);
            processed[nodeindex] = true;
        }
    }
    /* new node */
    else {
        InternalEdgeRef p1 = cache[ns.p1].notIf(ns.neg1);
        InternalEdgeRef p2 = cache[ns.p2].notIf(ns.neg2);

        InternalEdgeRef result = _manager->SimpleAnd(p1, p2).notIf(ns.nand);

        cache[nodeindex]     = result;
        processed[nodeindex] = true;
    }
}

unsigned int
aigpp::RewritingManager::getConeSize(unsigned int nodeindex, const aigpp::RewritingManager::CutStruct& cut) const
{
    assert(_pending.empty());

    static int processedMark = 0;
    ++processedMark;

    static std::vector<int> processed;
    if (processed.size() != _nodesCapacity) {
        processed.assign(_nodesCapacity, 0);
    }

    for (unsigned int i = 0; i != cut.size; ++i) {
        processed[cut.nodes[i]] = processedMark;
    }

    unsigned int count = 0;
    _pending.push(nodeindex);

    while (!_pending.empty()) {
        if (processed[_pending.top()] == processedMark) {
            _pending.pop();
        } else {
            const NodeStruct& ns = _nodes[_pending.top()];
            assert(!ns.var);

            ++count;

            processed[_pending.top()] = processedMark;
            _pending.pop();

            _pending.push(ns.p1);
            _pending.push(ns.p2);
        }
    }

    return count;
}

unsigned int
aigpp::RewritingManager::getImplInputs(int implIndex) const
{
    static lrabs::SimpleStack<int> pending;
    bool                           used[4] = {false, false, false, false};

    pending.push(implIndex);
    while (!pending.empty()) {
        if (pending.top() < 4) {
            used[pending.top()] = true;
            pending.pop();
        } else {
            const DynamicRewritingTable::IndexPair& ip = _rwTable.getAnd(pending.top() * 2);
            pending.pop();

            pending.push(ip.first / 2);
            pending.push(ip.second / 2);
        }
    }

    unsigned int usedCount = 0;
    if (used[0]) {
        ++usedCount;
    }
    if (used[1]) {
        ++usedCount;
    }
    if (used[2]) {
        ++usedCount;
    }
    if (used[3]) {
        ++usedCount;
    }

    return usedCount;
}

void
aigpp::RewritingManager::selectBestImpl(unsigned int nodeindex, const aigpp::RewritingManager::CutStruct& cut,
                                        int cutindex, const aigpp::NPNClassAndTransformations& npntrans,
                                        aigpp::RewritingManager::NPNImplInfo& bestImpl)
{
    /* map implementation variables to cut nodes */
    for (int implinput = 0; implinput != 4; ++implinput) {
        unsigned int cutinput = getInputPermutation(npntrans.t2, implinput);
        if (cutinput >= cut.size) {
            continue;
        }

        /* get cut node with index "cutindex" */
        int cutnode = cut.nodes[cutinput];

        if (_nodes[cutnode].nand != isInputNegated(npntrans.t2, implinput)) {
            cutnode = 1 + cutnode * 2;
        } else {
            cutnode = cutnode * 2;
        }

        _impl2node[implinput] = cutnode;
        _impl2nodeUsed.push(implinput);
    }

    const NPNClassImplementation& impls = _npnClassImpl[npntrans.npnClass];

    int implindex = 0;
    for (auto impl = impls.implementations.begin(); impl != impls.implementations.end(); ++impl, ++implindex) {
        _totalImplementationsChecked += 1;

        int newNodes = countNewNodes(*impl / 2);
        assert(newNodes > -1000000000);

        /* skip impl if newNodes >= 0, i.e. a matching toplevel node was found */
        if (newNodes >= 0) {
            /* node itself */
            if ((unsigned int)newNodes / 2 == nodeindex) {
                continue;
            }
        }

        int deletedNodes = countDeletedNodes(nodeindex, *impl / 2);

        /* skip impl if no nodes would be deleted */
        if (deletedNodes == 0) {
            continue;
        }

        int delta = -deletedNodes;

        if (newNodes < 0) {
            delta -= newNodes;
        } else {
            delta += 1;
        }

        if (delta < 0) {
            if (bestImpl.cutIndex == -1 || bestImpl.delta > delta) {
                bestImpl.cutIndex  = cutindex;
                bestImpl.npnclass  = npntrans;
                bestImpl.implIndex = *impl;
                bestImpl.delta     = delta;
            }
        }
    }
}

void
aigpp::RewritingManager::selectBestImplSmallerInput(unsigned int                              nodeindex,
                                                    const aigpp::RewritingManager::CutStruct& cut, int cutindex,
                                                    const aigpp::NPNClassAndTransformations& npntrans,
                                                    aigpp::RewritingManager::NPNImplInfo&    bestImpl)
{
    const NPNClassImplementation& impls = _npnClassImpl[npntrans.npnClass];

    /* map implementation variables to cut nodes */
    for (unsigned int implinput = 0; implinput != 4; ++implinput) {
        unsigned int cutinput = getInputPermutation(npntrans.t2, implinput);
        if (cutinput >= cut.size) {
            continue;
        }

        if ((impls.usedInputsMask & (1 << implinput)) != 0) {
            /* get cut node with index "cutindex" */
            int cutnode = cut.nodes[cutinput];

            if (_nodes[cutnode].nand != isInputNegated(npntrans.t2, implinput)) {
                cutnode = 1 + cutnode * 2;
            } else {
                cutnode = cutnode * 2;
            }

            _impl2node[implinput] = cutnode;
            _impl2nodeUsed.push(implinput);
        }
    }

    int implindex = 0;
    for (auto impl = impls.implementations.begin(); impl != impls.implementations.end(); ++impl, ++implindex) {
        _totalImplementationsChecked += 1;

        int newNodes = countNewNodes(*impl / 2);
        assert(newNodes > -1000000000);

        /* skip impl if newNodes >= 0, i.e. a matching toplevel node was found */
        if (newNodes >= 0) {
            /* node itself */
            if ((unsigned int)newNodes / 2 == nodeindex) {
                continue;
            }
        }

        int deletedNodes = countDeletedNodesSmallerInput(nodeindex, *impl / 2);

        /* skip impl if no nodes would be deleted */
        if (deletedNodes == 0) {
            continue;
        }

        int delta = -deletedNodes;

        if (newNodes < 0) {
            delta -= newNodes;
        } else {
            delta += 1;
        }

        if (delta <= 0) {
            if (bestImpl.cutIndex == -1 || bestImpl.delta > delta) {
                bestImpl.cutIndex  = cutindex;
                bestImpl.npnclass  = npntrans;
                bestImpl.implIndex = *impl;
                bestImpl.delta     = delta;
            }
        }
    }
}

#if 0
bool
aigpp::RewritingManager::
addNewImplementation(
    unsigned int nodeindex,
    const aigpp::RewritingManager::CutStruct& cut,
    const aigpp::NPNClassAndTransformations& npntrans )
{
    unsigned int size = getConeSize( nodeindex, cut ) + cut.size;

    NPNClassImplementation& impls = _npnClassImpl[npntrans.npnClass];

    /* return if there already exists an implementation with has fewer inputs or
       equal inputs and smaller nodes */

    if( impls.size != -1 &&
        ( ( cut.size > impls.usedInputs ) ||
          ( cut.size == impls.usedInputs && size > (unsigned int)impls.size ) ) )
    {
        return false;
    }

    bool betterSize = false;

    if( impls.implementations.empty() || impls.size == -1 )
    {
        std::cout << "no known implementations so far for npn class " << npntrans.npnClass << std::endl;
    }
    else if( cut.size < impls.usedInputs )
    {
        betterSize = true;

        std::cout << "best known implementation of npn class " << npntrans.npnClass << " has " << impls.usedInputs << " inputs. new implementation has " << cut.size << " inputs" << std::endl;
    }
    else if( size < (unsigned int)impls.size )
    {
        betterSize = true;

        std::cout << "best known implementation of npn class " << npntrans.npnClass << " has size " << impls.size << ". new implementation has size " << size << std::endl;
    }

    /* delete old impls if we found a smaller implementation */
    if( betterSize )
    {
        if( !impls.implementations.empty() )
        {
            for( std::list<DynamicRewritingTable::Index>::const_iterator i = impls.implementations.begin();
                 i != impls.implementations.end(); ++i )
            {
                _rwTable.deref( *i );
            }

            impls.implementations.clear();
        }

        impls.size = size;
        impls.usedInputs = cut.size;
    }

    /* add new impl */
    std::map<int, DynamicRewritingTable::Index> mapping;

    int usedInputsMask = 0;

    for( unsigned int i = 0; i != cut.size; ++i )
    {
        DynamicRewritingTable::Index idx = _rwTable.lookupInput( getInputPermutation( npntrans.t1, i ) );

        if( _nodes[cut.nodes[i]].nand != isInputNegated( npntrans.t1, i ) )
        {
            idx = DynamicRewritingTable::invert( idx );
        }

        mapping[cut.nodes[i]] = idx;

        usedInputsMask |= ( 1 << ( idx/2 ) );
    }

    assert( cut.size == lrabs::countOnes( usedInputsMask ) );
    assert( impls.usedInputsMask == 0 ||
            ( ( usedInputsMask & ~( impls.usedInputsMask ) ) == 0 ) );
    impls.usedInputsMask = usedInputsMask;

    static lrabs::SimpleStack<int> pending;
    pending.push( nodeindex );

    int newNodesCreated = 0;

    while( !pending.empty() )
    {
        if( mapping.find( pending.top() ) != mapping.end() )
        {
            pending.pop();
        }
        else
        {
            const NodeStruct& ns = _nodes[pending.top()];

            assert( !ns.var );

            if( mapping.find( ns.p1 ) == mapping.end() )
            {
                pending.push( ns.p1 );
                continue;
            }
            if( mapping.find( ns.p2 ) == mapping.end() )
            {
                pending.push( ns.p2 );
                continue;
            }

            DynamicRewritingTable::Index i1 = mapping[ns.p1];
            if( ns.neg1 )
            {
                i1 = DynamicRewritingTable::invert( i1 );
            }

            DynamicRewritingTable::Index i2 = mapping[ns.p2];
            if( ns.neg2 )
            {
                i2 = DynamicRewritingTable::invert( i2 );
            }

            DynamicRewritingTable::Index iand = _rwTable.lookupAnd( i1, i2 );
            if( iand == DynamicRewritingTable::InvalidIndex )
            {
                ++newNodesCreated;

                iand = _rwTable.insertAnd( i1, i2 );
                _impl2node.push_back( -1000000000 );
            }

            if( ns.nand )
            {
                iand = DynamicRewritingTable::invert( iand );
            }

            mapping[pending.top()] = iand;
            pending.pop();
        }
    }

    assert( mapping.find( nodeindex ) != mapping.end() );

    DynamicRewritingTable::Index nindex = mapping[nodeindex];
    if( isOutputNegated( npntrans.t1 ) )
    {
        nindex = DynamicRewritingTable::invert( nindex );
    }

    bool implFound = false;
    for( std::list<DynamicRewritingTable::Index>::const_iterator i = impls.implementations.begin();
         i != impls.implementations.end(); ++i )
    {
        if( *i == nindex )
        {
            implFound = true;
            break;
        }
    }

    if( !implFound )
    {
        std::cout << "adding new impl" << std::endl;

        impls.implementations.push_back( nindex );
        _rwTable.ref( nindex );

        _modifiedRWTable = true;
        return true;
    }
    else
    {
        return false;
    }
}
#endif

aigpp::RewritingManager::CutStruct*
aigpp::RewritingManager::getNewCut()
{
    _totalCutsComputed += 1;

    if (_freeCuts.empty()) {
        return new CutStruct;
    } else {
        CutStruct* cut = _freeCuts.top();
        _freeCuts.pop();
        return cut;
    }
}

void
aigpp::RewritingManager::releaseCuts(aigpp::RewritingManager::CutStruct* c)
{
    while (c != nullptr) {
        _freeCuts.push(c);
        c = c->next;
    }
}
