/**************************************************************
 *
 *       AIGPP // FRAIGManager.cc
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

#include "FRAIGManager.hh"

#include <lrabsutil/Resources.hh>

#include "Manager.hh"
#include "Node.hh"

#define _unused(x) ((void)(x))

aigpp::FRAIGManager::FRAIGManager(aigpp::Manager* manager) :
    _manager(manager),
    _varCount(0),
    _size(0),
    _capacity(1024),
    _solver(nullptr),

    _initTime(0),
    _cnfCreationTime(0),
    _satTime(0),
    _satMaxTime(0),
    _propTime(0),
    _splitTime(0),
    _rebuildTime(0),
    _fraigTime(0),
    _satChecks(0),
    _satEquiv(0),
    _satChecksConst(0),
    _satEquivConst(0),

    _totalInitTime(0),
    _totalCnfCreationTime(0),
    _totalSatTime(0),
    _totalSatMaxTime(0),
    _totalPropTime(0),
    _totalSplitTime(0),
    _totalRebuildTime(0),
    _totalFraigTime(0),
    _totalSatChecks(0),
    _totalSatEquiv(0),
    _totalSatChecksConst(0),
    _totalSatEquivConst(0),

    _totalApplications(0)
{
    _nodes     = new NodeStruct[_capacity];
    _unreduced = new bool[_capacity];
#ifdef FRAIGMANAGER_TIMEOUT
    _reduced = new bool[_capacity];
#endif
    _sim          = new unsigned long[_capacity];
    _satVars      = new Minisat::Var[_capacity];
    _depths       = new int[_capacity];
    _equivalences = new std::size_t[_capacity];
    _aigNodes     = new aigpp::Node*[_capacity];

#ifdef FRAIGMANAGER_RESET_ENABLE
    _newEquivalencesSinceLastReset = 0;
    _checksSinceLastReset          = 0;
#endif
}

aigpp::FRAIGManager::~FRAIGManager()
{
    assert(_solver == 0);

    for (std::size_t i = 0; i != _freeClasses.size(); ++i) {
        delete _freeClasses[i];
    }

    delete[] _nodes;
    delete[] _unreduced;
#ifdef FRAIGMANAGER_TIMEOUT
    delete[] _reduced;
#endif
    delete[] _sim;
    delete[] _satVars;
    delete[] _depths;
    delete[] _equivalences;
    delete[] _aigNodes;
}

void
aigpp::FRAIGManager::makeFRAIG(
#ifdef FRAIGMANAGER_TIMEOUT
    double timeout
#endif
)
{
    assert(_manager->simCreation() == true);

    _manager->_timeoutDuringLastFRAIG = false;

    _initTime        = 0;
    _cnfCreationTime = 0;
    _satTime         = 0;
    _satMaxTime      = 0;
    _propTime        = 0;
    _splitTime       = 0;
    _fraigTime       = 0;
    _rebuildTime     = 0;
    _satChecks       = 0;
    _satEquiv        = 0;
    _satChecksConst  = 0;
    _satEquivConst   = 0;

#ifdef FRAIGMANAGER_RESET_ENABLE
    _newEquivalencesSinceLastReset = 0;
    _checksSinceLastReset          = 0;
#endif

    double startTime = lrabs::cpuTime();

#ifdef FRAIGMANAGER_TIMEOUT
    bool   hasTimeoutTime = (timeout > 0);
    double timeoutTime    = startTime + (hasTimeoutTime ? timeout : 0);
//    bool timeoutReached = false;
//    _unused(timeoutReached);
#endif

    if (_manager->_unreducedNodes == 0) {
        return;
    }
    _manager->garbageCollect();
    if (_manager->_unreducedNodes == 0) {
        return;
    }

    /* create data structures */
    initStructures();
    initClasses();

    Minisat::vec<Minisat::Lit> assumptions1(1), assumptions2(2);

    /* perform sat checks until all equivalence classes are singletons */
    /* todo: stop if all classes are composed completely of previously reduced
     * nodes */
    while (!_nonEmptyClasses.empty()) {
#ifdef FRAIGMANAGER_TIMEOUT
        if (hasTimeoutTime) {
            if (lrabs::cpuTime() > timeoutTime) {
                //                timeoutReached = true;
                _manager->_timeoutDuringLastFRAIG = true;

                // std::cout << "c fraiging: timeout reached!" << std::endl;

                break;
            }
        }
#endif

        auto currentClass = getEquivalenceClassToSplit();

        /* get non-singleton equivalence class */
        EquivalenceClass* e = *currentClass;

        /* representant */
        std::size_t n2 = e->_repr >> 1;

        /* find pair with at least one unreduced node */
        auto candidate = e->_nodes.begin();

        if (n2 == 0ul || !_unreduced[n2]) {
            while (candidate != e->_nodes.end()) {
                if (_unreduced[*candidate >> 1u]) {
                    break;
                }
                ++candidate;
            }
        }

        /* no candidate found -> remove class and continue */
        if (candidate == e->_nodes.end()) {
#ifdef FRAIGMANAGER_TIMEOUT
            _reduced[n2] = true;
            for (std::size_t p : e->_nodes) {
                _reduced[p >> 1u] = true;
            }
#endif

#ifdef FRAIGMANAGER_SIMP_ENABLE
            markNodeAsDone(e->_repr >> 1);
            for (std::size_t p : e->_nodes) {
                markNodeAsDone(p >> 1u);
            }
#endif

            _nonEmptyClasses.erase(currentClass);
            e->_nodes.clear();
            releaseEquivalenceClassObject(e);

            continue;
        }

        /* equivalence check of candidate with the representant */
        std::size_t n1        = *candidate >> 1u;
        bool        isEqCheck = ((e->_repr & 1ul) == (*candidate & 1ul));

        /* create necessary clause */
        addClauses(n1);
        /* n2 may be 0, in case e is the "zero" class */
        if (n2 != 0ul) {
            addClauses(n2);
        }

        /* zero or one check */
        if (n2 == 0ul) {
            ++_satChecksConst;

#ifdef FRAIGMANAGER_RESET_ENABLE
            ++_checksSinceLastReset;
#endif

            assumptions1[0] = Minisat::mkLit(_satVars[n1], !isEqCheck);

            double startCheckTime = lrabs::cpuTime();
            bool   satisfiable    = _solver->solve(assumptions1);
            double deltaT         = lrabs::cpuTime() - startCheckTime;
            _satTime += deltaT;
            if (deltaT > _satMaxTime) {
                _satMaxTime = deltaT;
            }

            /* equivalence found */
            if (!satisfiable) {
                ++_satEquivConst;

#ifdef FRAIGMANAGER_RESET_ENABLE
                ++_newEquivalencesSinceLastReset;
#endif

                _equivalences[n1] = (0 << 1) + (!isEqCheck ? 1 : 0);

#ifdef FRAIGMANAGER_SIMP_ENABLE
                incRefCount(0);
                markNodeAsDone(n1);
#endif

                /* remove node from its equivalence class */
                e->_nodes.erase(candidate);
                if (e->_nodes.empty()) {
                    _nonEmptyClasses.erase(currentClass);
                    releaseEquivalenceClassObject(e);
                }
#ifdef FRAIGMANAGER_DEPTHSORT
                else {
                    /* update depth */

                    std::size_t firstCandidate = ~0ul;

                    for (std::size_t n : e->_nodes) {
                        if (_unreduced[e->_repr >> 1u] || _unreduced[n >> 1u]) {
                            firstCandidate = n >> 1u;
                            break;
                        }
                    }

                    if (firstCandidate == ~0ul) {
                        e->_depth = 0;
                    } else {
                        assert(e->_depth <= _depths[firstCandidate]);
                        e->_depth = _depths[firstCandidate];
                    }
                }
#endif

#ifdef FRAIGMANAGER_RESET_ENABLE
                /*
                std::cout << "eq found: " << _size << " " << _checksSinceLastReset << "
                " << _newEquivalencesSinceLastReset << std::endl; if( ( 100 *
                _checksSinceLastReset >= _size ) && ( 10 *
                _newEquivalencesSinceLastReset >= _checksSinceLastReset ) )
                {
                    resetSolver();
                }
                */
                resetSolver();
#endif
            } else /* no equivalence -> split classes */
            {
                propagate();
                split();
            }
        }

        /* equivalence or antivalence check */
        else {
            ++_satChecks;

#ifdef FRAIGMANAGER_RESET_ENABLE
            ++_checksSinceLastReset;
#endif

            /* check if n1==0, n2==1 is satisfiable (equivalence) */
            /* check if n1==1, n2==1 is satisfiable (antivalence) */
            assumptions2[0] = Minisat::mkLit(_satVars[n1], isEqCheck);
            assumptions2[1] = Minisat::mkLit(_satVars[n2]);

            double startCheckTime = lrabs::cpuTime();
            bool   satisfiable    = _solver->solve(assumptions2);
            double deltaT         = lrabs::cpuTime() - startCheckTime;
            _satTime += deltaT;
            if (deltaT > _satMaxTime) {
                _satMaxTime = deltaT;
            }

            if (!satisfiable) {
                /* add impl clause */
                addBinary(Minisat::mkLit(_satVars[n1], !isEqCheck), Minisat::mkLit(_satVars[n2], true));

                /* check if n1==1, n2==0 is satisfiable (equivalence) */
                /* check if n1==0, n2==0 is satisfiable (antivalence) */
                assumptions2[0] = Minisat::mkLit(_satVars[n1], !isEqCheck);
                assumptions2[1] = ~Minisat::mkLit(_satVars[n2]);

                startCheckTime = lrabs::cpuTime();
                satisfiable    = _solver->solve(assumptions2);
                deltaT         = lrabs::cpuTime() - startCheckTime;
                _satTime += deltaT;
                if (deltaT > _satMaxTime) {
                    _satMaxTime = deltaT;
                }
            }

            /* equivalence found */
            if (!satisfiable) {
                ++_satEquiv;

#ifdef FRAIGMANAGER_RESET_ENABLE
                ++_newEquivalencesSinceLastReset;
#endif
                _equivalences[n1] = (n2 << 1) + (!isEqCheck ? 1 : 0);

#ifdef FRAIGMANAGER_SIMP_ENABLE
                incRefCount(n2);
                markNodeAsDone(n1);
#endif

                if (!_unreduced[n1] && _unreduced[n2]) {
                    _unreduced[n2] = false;
                }

                e->_nodes.erase(candidate);
                if (e->_nodes.empty()) {
#ifdef FRAIGMANAGER_SIMP_ENABLE
                    markNodeAsDone(e->_repr >> 1);
#endif

                    _nonEmptyClasses.erase(currentClass);
                    releaseEquivalenceClassObject(e);
                }
#ifdef FRAIGMANAGER_DEPTHSORT
                else {
                    /* update depth */

                    std::size_t firstCandidate1 = e->_repr >> 1;
                    std::size_t firstCandidate2 = ~0ul;

                    for (std::list<std::size_t>::const_iterator n = e->_nodes.begin(); n != e->_nodes.end(); ++n) {
                        if (_unreduced[firstCandidate1] || _unreduced[*n >> 1]) {
                            firstCandidate2 = *n >> 1;
                            break;
                        }
                    }

                    if (firstCandidate2 == ~0ul) {
                        e->_depth = _depths[firstCandidate1];
                    } else {
                        e->_depth = _depths[firstCandidate2];
                    }
                }
#endif

#ifdef FRAIGMANAGER_RESET_ENABLE
                /*
                std::cout << "eq found: " << _size << " " << _checksSinceLastReset << "
                " << _newEquivalencesSinceLastReset << std::endl; if( ( 100 *
                _checksSinceLastReset >= _size ) && ( 10 *
                _newEquivalencesSinceLastReset >= _checksSinceLastReset ) )
                {
                    resetSolver();
                }
                */
                resetSolver();
#endif
            } else /* no equivalence -> split classes */
            {
                propagate();
                split();
            }
        }
    }

#ifdef FRAIGMANAGER_TIMEOUT
    // std::cout << "FRAIGing: checks=" << ( _satChecks + _satChecksConst ) << "
    // equiv=" << ( _satEquiv + _satEquivConst ) << std::endl;
#endif

#if 1

    double startRebuildTime = lrabs::cpuTime();

    /* save old nodes */
    std::vector<Node*> nodes;
    nodes.reserve(_manager->nodeCount());
    for (Node* n = _manager->_nodes; n != nullptr; n = n->next()) {
        nodes.push_back(n);
        n->_refCount = 0;
    }

    /* clear counters and data structures */
    _manager->_nodes     = nullptr;
    _manager->_lastNode  = nullptr;
    _manager->_nodeCount = 0;
    _manager->_unique.clear();
#    ifdef USE_COMPUTEDTABLE
    _manager->_computed.clear();
#    endif
    _manager->_simTable.clear();
    _manager->_unreducedNodes = 0;

    std::map<Node*, Edge> cache;

    /* add variables again */
    for (Node* v : _manager->_variables) {
        _manager->addToNodesList(v);
        _manager->_simTable.insert(v);

        cache[v] = Edge(v);
#    ifdef FRAIGMANAGER_TIMEOUT
        _manager->ref(v);
#    endif
        _manager->ref(v);
    }

    /* rebuild AIG */
    static lrabs::SimpleStack<Node*> pending;
    for (const auto& e : _manager->_extRefTable->_targets) {
        if (e._refCount == 0) {
            continue;
        }
        if (e._target.isConstant()) {
            continue;
        }
        pending.push(e._target.node());
    }

    while (!pending.empty()) {
        Node* pt = pending.top();

        /* a variable node, already processed in "variable loop" -> nothing to do */
        if (pt->isVar()) {
            pending.pop();
            continue;
        }
        /* already present in cache -> nothing to do */
        else if (cache.find(pt) != cache.end()) {
            pending.pop();
            continue;
        }
        /* node is in equivalence mapping */
        else if (_equivalences[_id[pt]] != ~0ul) {
            Edge eq(_aigNodes[_equivalences[_id[pt]] >> 1], (_equivalences[_id[pt]] & 1) != 0);
            if (eq.isConstant()) {
                pending.pop();
                cache[pt] = eq;
#    ifdef FRAIGMANAGER_TIMEOUT
                _manager->ref(eq.node());
#    endif
            } else if (cache.find(eq.node()) != cache.end()) {
                pending.pop();
                cache[pt] = cache[eq.node()].notIf(eq.isInverted());
#    ifdef FRAIGMANAGER_TIMEOUT
                _manager->ref(cache[pt].node());
#    endif
            } else {
                pending.push(eq.node());
            }
            continue;
        }

        /* make sure first parent has been processed */
        auto p1 = cache.find(pt->parent1().node());
        if (p1 == cache.end()) {
            pending.push(pt->parent1().node());
            continue;
        }

        /* make sure second parent has been processed */
        auto p2 = cache.find(pt->parent2().node());
        if (p2 == cache.cend()) {
            pending.push(pt->parent2().node());
            continue;
        }

        Edge e1 = p1->second.notIf(pt->parent1().isInverted());
        Edge e2 = p2->second.notIf(pt->parent2().isInverted());

#    ifdef FRAIGMANAGER_TIMEOUT
        if (e1.isConstant() || e2.isConstant()) {
            //            assert( timeoutReached );

            if ((e1.structurallyFalse()) || (e2.structurallyFalse())) {
                cache[pt] = const0.notIf(pt->isNAND());
            } else if (e1.structurallyTrue()) {
                cache[pt] = e2.notIf(pt->isNAND());
                _manager->ref(e2.node());
            } else if (e2.structurallyTrue()) {
                cache[pt] = e1.notIf(pt->isNAND());
                _manager->ref(e1.node());
            }

            pending.pop();
            continue;
        }

        if (e1.structurallyEquivalent(e2)) {
            //            assert( timeoutReached );

            cache[pt] = e1.notIf(pt->isNAND());
            _manager->ref(e1.node());

            pending.pop();
            continue;
        } else if (e1.structurallyAntivalent(e2)) {
            //            assert( timeoutReached );

            cache[pt] = const0.notIf(pt->isNAND());

            pending.pop();
            continue;
        }
#    else
        if (e1.isConstant() || e2.isConstant()) {
            /* QUIET version */
#        if 0
            std::cout << "node " << pt->isNAND() << " " << pt << " " << _id[pt] << std::endl;
            std::cout << "p1 " << e1.isInverted() << " " << e1.node() << " " << _id[e1.node()] << std::endl;
            std::cout << "p2 " << e2.isInverted() << " " << e2.node() << " " << _id[e2.node()] << std::endl;
            std::cout << "at least one parent is constant, but node is not equivalent to anything" << std::endl;
            std::cout << e1.isConstant() << " " << e2.isConstant() << std::endl;

            std::cout << "reduced " << pt->flag<Node::FLAG_ISREDUCED>() << " " << e1.node()->flag<Node::FLAG_ISREDUCED>() << std::endl;
#        endif
            abort();
        }
        assert(e1.node() != e2.node());
#    endif

        if (e1.node()->index() > e2.node()->index()) {
            e1.swap(e2);
        }

        {
            Edge res;
#    ifdef FRAIGMANAGER_TIMEOUT
            if (_manager->_unique.lookup(e1, e2, res)) {
                cache[pt] = res.notIf(pt->isNAND());
                _manager->ref(res.node());

                pending.pop();
                continue;
            }
#    else
            assert(!_manager->_unique.lookup(e1, e2, res));
#    endif
        }

        pt->_parent1 = e1;
        pt->_parent2 = e2;

        _manager->ref(e1.node());
        _manager->ref(e2.node());

        _manager->_unique.insert(pt);
        _manager->_simTable.insert(pt);
        _manager->addToNodesList(pt);
        cache[pt] = Edge(pt);
#    ifdef FRAIGMANAGER_TIMEOUT
        _manager->ref(pt);
#    endif

        pending.pop();
    }

    /* modify external references */
    for (auto& e : _manager->_extRefTable->_targets) {
        if (e._refCount == 0) {
            continue;
        }
        if (e._target.isConstant()) {
            continue;
        }

        assert(cache.find(e._target.node()) != cache.end());

        e._target._node = cache[e._target.node()].notIf(e._target.isInverted())._node;
        _manager->ref(e._target.node());
    }
    _manager->_extRefTable->rebuildMapping();

#    ifdef FRAIGMANAGER_TIMEOUT
    for (const auto& p : cache) {
        if (p.second.node() != nullptr) {
            _manager->deref(p.second.node());
        }
    }

    std::set<Node*> deadNodes;
    for (Node* n = _manager->nodesList(); n != nullptr; n = n->next()) {
        if (n->refCount() == 0) {
            deadNodes.insert(n);
        }
    }

    _manager->garbageCollect();
#    endif

    /* delete remaining nodes */
    for (Node* n : nodes) {
#    ifdef FRAIGMANAGER_TIMEOUT
        if (deadNodes.find(n) != deadNodes.end()) {
            // std::cout << "trying double-release dead node " << *n << std::endl;
            continue;
        }
#    endif
        if (n->refCount() == 0) {
            _manager->releaseNode(n);
        } else {
#    ifdef FRAIGMANAGER_TIMEOUT
            if (n->flag<Node::FLAG_ISREDUCED>()) {
                /**/
            } else if (_reduced[_id[n]]) {
                n->setFlag<Node::FLAG_ISREDUCED>();
            } else {
                //                assert( timeoutReached );
                ++(_manager->_unreducedNodes);
            }
#    else
            n->setFlag<Node::FLAG_ISREDUCED>();
#    endif
        }
    }

    _manager->updateSimTable();

    _rebuildTime += lrabs::cpuTime() - startRebuildTime;

#endif

#ifdef FRAIGMANAGER_TIMEOUT
    /*
    if( _manager->_unreducedNodes != 0 )
    {
        std::cout << "c unreduced nodes " << _manager->_unreducedNodes <<
    std::endl;
    }
    */
#endif

    /* delete data structures */
    clearStructures();

    _fraigTime = (lrabs::cpuTime() - startTime);
    // std::cout << "fraig_stats " << _fraigTime << " " << ( _satEquiv +
    // _satEquivConst )  << "/" << ( _satChecks + _satChecksConst ) << std::endl;

    _totalInitTime += _initTime;
    _totalCnfCreationTime += _cnfCreationTime;
    _totalSatTime += _satTime;
    if (_satMaxTime > _totalSatMaxTime) {
        _totalSatMaxTime = _satMaxTime;
    }
    _totalPropTime += _propTime;
    _totalSplitTime += _splitTime;
    _totalRebuildTime += _rebuildTime;
    _totalFraigTime += _fraigTime;
    _totalSatChecks += _satChecks;
    _totalSatEquiv += _satEquiv;
    _totalSatChecksConst += _satChecksConst;
    _totalSatEquivConst += _satEquivConst;

    ++_totalApplications;
}

void
aigpp::FRAIGManager::printStats(std::ostream& os) const
{
    os << "fr_applications " << _totalApplications << "\n"
       << "fr_satequiv     " << _totalSatChecks << " " << _totalSatEquiv << "\n"
       << "fr_satconst     " << _totalSatChecksConst << " " << _totalSatEquivConst << "\n"
       << "fr_time         " << _totalFraigTime << "\n"
       << "fr_init_time    " << _totalInitTime << "\n"
       << "fr_cnf_time     " << _totalCnfCreationTime << "\n"
       << "fr_sat_time     " << _totalSatTime << "\n"
       << "fr_prop_time    " << _totalPropTime << "\n"
       << "fr_split_time   " << _totalSplitTime << "\n"
       << "fr_rebuild_time " << _totalRebuildTime << "\n"
       << std::flush;
}

void
aigpp::FRAIGManager::initStructures()
{
    double startTime = lrabs::cpuTime();

    assert(_solver == 0);

#ifdef FRAIGMANAGER_SIMP_ENABLE
    _solver             = new Minisat::SimpSolver();
    _solver->clause_lim = 8;
#else
    _solver = new Minisat::Solver();
#endif

    _solver->ccmin_mode = 1;

    _size     = _manager->nodeCount() + 1;
    _varCount = _manager->variableCount();

    /* resize arrays if needed */
    if (_capacity < _size) {
        while (_capacity < _size) {
            _capacity <<= 1;
        }

        assert(_capacity >= _size);

        delete[] _nodes;
        delete[] _unreduced;
#ifdef FRAIGMANAGER_TIMEOUT
        delete[] _reduced;
#endif
        delete[] _sim;
        delete[] _satVars;
        delete[] _depths;
        delete[] _equivalences;
        delete[] _aigNodes;

        _nodes     = new NodeStruct[_capacity];
        _unreduced = new bool[_capacity];
#ifdef FRAIGMANAGER_TIMEOUT
        _reduced = new bool[_capacity];
#endif
        _sim          = new unsigned long[_capacity];
        _satVars      = new Minisat::Var[_capacity];
        _depths       = new int[_capacity];
        _equivalences = new std::size_t[_capacity];
        _aigNodes     = new aigpp::Node*[_capacity];
    }

    for (std::size_t i = 0; i != _size; ++i) {
        _satVars[i]      = var_Undef;
        _equivalences[i] = ~0ul;
    }

    /* init infos for zero node */
    _unreduced[0] = false;
#ifdef FRAIGMANAGER_TIMEOUT
    _reduced[0] = true;
#endif
    _sim[0] = 0ul;

    _nodes[0]._var = true;
#ifdef FRAIGMANAGER_SIMP_ENABLE
    _nodes[0]._originalRef = 1;
    _nodes[0]._ref         = 1;

    _nodes[0]._done = false;
#endif

    _aigNodes[0] = nullptr;
    _id[nullptr] = 0;
    _depths[0]   = 0;

    /* init variables */
    std::size_t nextid = 1;
    for (std::size_t i = 0; i != _varCount; ++i) {
        _unreduced[nextid] = false;
#ifdef FRAIGMANAGER_TIMEOUT
        _reduced[nextid] = true;
#endif
        _id[_manager->_variables[i]] = nextid;

        NodeStruct& ns = _nodes[nextid];
        ns._var        = true;
#ifdef FRAIGMANAGER_SIMP_ENABLE
        ns._originalRef = 1;
        ns._ref         = 1;
        ns._done        = false;
#endif

        _depths[nextid]   = 1;
        _aigNodes[nextid] = _manager->_variables[i];

        ++nextid;
    }

    /* init remaining nodes */
    for (aigpp::Node* n = _manager->_nodes; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            continue;
        }

        _unreduced[nextid] = !(n->flag<aigpp::Node::FLAG_ISREDUCED>());
#ifdef FRAIGMANAGER_TIMEOUT
        _reduced[nextid] = n->flag<aigpp::Node::FLAG_ISREDUCED>();
#endif

        _id[n]            = nextid;
        _aigNodes[nextid] = n;

        NodeStruct& ns = _nodes[nextid];
        ns._var        = false;
        ns._nand       = n->isNAND();
        ns._neg1       = n->parent1().isInverted();
        ns._neg2       = n->parent2().isInverted();
        ns._p1         = _id[n->parent1().node()];
        ns._p2         = _id[n->parent2().node()];
#ifdef FRAIGMANAGER_SIMP_ENABLE
        ns._ref  = 0;
        ns._done = false;

        incRefCount(ns._p1);
        incRefCount(ns._p2);
#endif
        _depths[nextid] = 1 + std::max(_depths[ns._p1], _depths[ns._p2]);

        ++nextid;
    }

    assert(nextid == _size);

    _initTime += lrabs::cpuTime() - startTime;
}

struct sortByDepth : public std::binary_function<int, int, bool>
{
    const int* _depths;

    explicit sortByDepth(const int* depths) : _depths(depths) {}

    bool operator()(const int& x, const int& y) { return (_depths[x >> 1] <= _depths[y >> 1]); };
};

void
aigpp::FRAIGManager::initClasses()
{
    double startTime = lrabs::cpuTime();

    std::set<Node*> processed;

    /* create zero equivalence class */
    {
        EquivalenceClass* e = allocateEquivalenceClassObject(0);

        static SimVector zerosim;
        assert(zerosim.isZero());

        /* lookup (possibly) equivalent nodes */
        Node* nodeclass = _manager->_simTable.lookup(zerosim);
        if (nodeclass != nullptr) {
            Node* simnode = nodeclass;
            do {
                processed.insert(simnode);
                e->_nodes.push_back(_id[simnode] << 1);

                simnode = simnode->_nEqualSim;
            } while (simnode != nodeclass);
        }

        /* lookup (possibly) antivalent nodes */
        nodeclass = _manager->_simTable.lookupInverted(zerosim);
        if (nodeclass != nullptr) {
            Node* simnode = nodeclass;
            do {
                processed.insert(simnode);
                e->_nodes.push_back((_id[simnode] << 1) | 1);

                simnode = simnode->_nEqualSim;
            } while (simnode != nodeclass);
        }

        /* >= 1 class members (in addition to the "0" node) -> sort and add class */
        if (!e->_nodes.empty()) {
            /* sort members by depth */
            e->_nodes.sort(sortByDepth(_depths));

#ifdef FRAIGMANAGER_DEPTHSORT
            /* find first unreduced node */
            std::size_t firstUnreduced = ~0ul;

            for (std::list<std::size_t>::const_iterator n = e->_nodes.begin(); n != e->_nodes.end(); ++n) {
                if (_unreduced[*n >> 1]) {
                    firstUnreduced = *n >> 1;
                    break;
                }
            }

            if (firstUnreduced != ~0ul) {
                e->_depth = 0;
                _nonEmptyClasses.push_back(e);
            } else {
                e->_depth = _depths[firstUnreduced];
                _nonEmptyClasses.push_back(e);
            }
#else  /* FRAIGMANAGER_DEPTHSORT */
            _nonEmptyClasses.push_back(e);
#endif /* FRAIGMANAGER_DEPTHSORT */
        }
        /* no class members -> delete class */
        else {
            releaseEquivalenceClassObject(e);
        }
    }

    /* create remaining equivalence classes based on simulation vectors */
    for (Node* n = _manager->_nodes; n != nullptr; n = n->next()) {
        if (processed.find(n) != processed.end()) {
            continue;
        }

        EquivalenceClass* e = allocateEquivalenceClassObject();

        /* lookup (possibly) equivalent nodes */
        // Node* nodeclass = _manager->_simTable.lookup( n->sim() );
        Node* nodeclass = _manager->_simTable.getClass(n);
        if (nodeclass != nullptr) {
            Node* simnode = nodeclass;
            do {
                processed.insert(simnode);
                e->_nodes.push_back(_id[simnode] << 1);

                simnode = simnode->_nEqualSim;
            } while (simnode != nodeclass);
        }

        /* lookup (possibly) antivalent nodes */
        // nodeclass = _manager->_simTable.lookupInverted( n->sim() );
        if (nodeclass != nullptr) {
            nodeclass = _manager->_simTable.getInvertedClass(nodeclass);
        } else {
            nodeclass = _manager->_simTable.lookupInverted(n->sim());
        }

        if (nodeclass != nullptr) {
            Node* simnode = nodeclass;
            do {
                processed.insert(simnode);
                e->_nodes.push_back((_id[simnode] << 1) | 1ul);

                simnode = simnode->_nEqualSim;
            } while (simnode != nodeclass);
        }

        /* >= 2 class members -> sort and add class */
        if (e->_nodes.size() >= 2) {
            /* sort members by depth */
            e->_nodes.sort(sortByDepth(_depths));

#ifdef FRAIGMANAGER_DEPTHSORT
            /* find first candidate pair for eq checks */
            std::size_t firstCandidate1 = ~0ul;
            std::size_t firstCandidate2 = ~0ul;

            for (std::list<std::size_t>::const_iterator n = e->_nodes.begin(); n != e->_nodes.end(); ++n) {
                if (firstCandidate1 == ~0ul) {
                    firstCandidate1 = *n >> 1u;
                } else if (_unreduced[firstCandidate1]) {
                    firstCandidate2 = *n >> 1u;
                    break;
                } else {
                    if (_unreduced[*n >> 1u]) {
                        firstCandidate2 = *n >> 1u;
                        break;
                    }
                }
            }

            if (firstCandidate2 == ~0ul) {
                e->_depth = _depths[firstCandidate1];

                e->_repr = *e->_nodes.begin();
                e->_nodes.erase(e->_nodes.begin());

                _nonEmptyClasses.push_back(e);

#    ifdef FRAIGMANAGER_TIMEOUT
                _reduced[e->_repr >> 1u] = true;
#    endif
            } else {
                assert(_depths[firstCandidate2] >= _depths[firstCandidate1]);

                e->_depth = _depths[firstCandidate2];

                /* select shallowest member as representant */
                e->_repr = *e->_nodes.begin();
                e->_nodes.erase(e->_nodes.begin());

                _nonEmptyClasses.push_back(e);

#    ifdef FRAIGMANAGER_TIMEOUT
                _reduced[e->_repr >> 1] = true;
#    endif
            }
#else /* FRAIGMANAGER_DEPTHSORT */
            /* select shallowest member as representant */
            e->_repr = *e->_nodes.begin();
            e->_nodes.erase(e->_nodes.begin());

            _nonEmptyClasses.push_back(e);

#    ifdef FRAIGMANAGER_TIMEOUT
            _reduced[e->_repr >> 1u] = true;
#    endif
#endif /* FRAIGMANAGER_DEPTHSORT */
        }
        /* <= 1 class member -> delete class */
        else {
            std::size_t id_n = _id[n];
#ifdef FRAIGMANAGER_TIMEOUT
            _reduced[id_n] = true;
#endif

#ifdef FRAIGMANAGER_SIMP_ENABLE
            markNodeAsDone(id_n);
#endif

            e->_nodes.clear();
            releaseEquivalenceClassObject(e);
        }
    }

    _initTime += lrabs::cpuTime() - startTime;
}

void
aigpp::FRAIGManager::clearStructures()
{
    assert(_solver != nullptr);

    delete _solver;
    _solver = nullptr;

    for (EquivalenceClass* eq : _nonEmptyClasses) {
        eq->_nodes.clear();
        releaseEquivalenceClassObject(eq);
    }
    _nonEmptyClasses.clear();
}

#ifdef FRAIGMANAGER_SIMP_ENABLE
void
aigpp::FRAIGManager::markNodeAsDone(std::size_t index)
{
    NodeStruct& ns = _nodes[index];

    // std::cout << __func__ << " " << index << " " << ns._ref << " " << ns._done
    // << std::endl;

    assert(!ns._done);
    ns._done = true;

    if (!ns._var && ns._ref == 0) {
        Minisat::Var& satvar_n = _satVars[index];
        if (satvar_n != var_Undef) {
            // std::cout << "unfreeze1 " << index << std::endl;
            _solver->setFrozen(satvar_n, false);
        }
    }
}

void
aigpp::FRAIGManager::decRefCount(std::size_t index)
{
    NodeStruct& ns = _nodes[index];
    // std::cout << __func__ << " " << index << " " << ns._ref << " " << ns._done
    // << std::endl;

    assert(ns._ref > 0);

    --ns._ref;

    if (ns._ref == 0 && !ns._var && ns._done) {
        Minisat::Var& satvar_n = _satVars[index];
        if (satvar_n != var_Undef) {
            // std::cout << "unfreeze2 " << index << std::endl;
            _solver->setFrozen(satvar_n, false);
        }
    }
}

void
aigpp::FRAIGManager::incRefCount(std::size_t index)
{
    NodeStruct& ns = _nodes[index];

    ++ns._originalRef;
    ++ns._ref;
}
#endif

aigpp::FRAIGManager::EquivalenceClass*
aigpp::FRAIGManager::allocateEquivalenceClassObject(std::size_t repr)
{
    if (!_freeClasses.empty()) {
        EquivalenceClass* e = _freeClasses.top();
        _freeClasses.pop();
        e->_repr = repr;
#ifdef FRAIGMANAGER_DEPTHSORT
        e->_depth = 0;
#endif
        return e;
    } else {
        return new EquivalenceClass(repr);
    }
}

void
aigpp::FRAIGManager::releaseEquivalenceClassObject(aigpp::FRAIGManager::EquivalenceClass* e)
{
    assert(e->_nodes.empty());
    _freeClasses.push(e);
}

std::list<aigpp::FRAIGManager::EquivalenceClass*>::iterator
aigpp::FRAIGManager::getEquivalenceClassToSplit()
{
    assert(!_nonEmptyClasses.empty());

#ifdef FRAIGMANAGER_DEPTHSORT
    int  minDepth = -1;
    auto minClass = _nonEmptyClasses.end();

    for (auto c = _nonEmptyClasses.begin(); c != _nonEmptyClasses.end(); ++c) {
        if (minClass == _nonEmptyClasses.end() || (*c)->_depth < minDepth) {
            minClass = c;
            minDepth = (*c)->_depth;
        }
    }

    return minClass;
#else
    return _nonEmptyClasses.begin();
#endif
}

void
aigpp::FRAIGManager::resetSolver()
{
    /* clear variable mapping */
    for (std::size_t i = 0; i != _size; ++i) {
        _satVars[i] = var_Undef;
    }

    /* create new solver object */
    delete _solver;

#ifdef FRAIGMANAGER_SIMP_ENABLE
    _solver             = new Minisat::SimpSolver();
    _solver->clause_lim = 8;
#else
    _solver = new Minisat::Solver();
#endif

    _solver->ccmin_mode = 1;

#if defined(FRAIGMANAGER_SIMP_ENABLE) and defined(FRAIGMANAGER_RESET_ENABLE)
    for (std::size_t i = 0; i != _size; ++i) {
        _nodes[i]._ref = _nodes[i]._originalRef;
    }
    for (std::size_t i = 0; i != _size; ++i) {
        if (_equivalences[i] != ~0ul) {
            ++_nodes[_equivalences[i] >> 1]._ref;
        }
    }
#endif

#if defined(FRAIGMANAGER_RESET_ENABLE)
    _newEquivalencesSinceLastReset = 0;
    _checksSinceLastReset          = 0;
#endif
}

void
aigpp::FRAIGManager::addUnit(const Minisat::Lit lit1)
{
    Minisat::vec<Minisat::Lit> c(1);
    c[0] = lit1;
    _solver->addClause(c);
}

void
aigpp::FRAIGManager::addBinary(const Minisat::Lit lit1, const Minisat::Lit lit2)
{
    Minisat::vec<Minisat::Lit> c(2);
    c[0] = lit1;
    c[1] = lit2;
    _solver->addClause(c);
}

void
aigpp::FRAIGManager::addTernary(const Minisat::Lit lit1, const Minisat::Lit lit2, const Minisat::Lit lit3)
{
    Minisat::vec<Minisat::Lit> c(3);
    c[0] = lit1;
    c[1] = lit2;
    c[2] = lit3;
    _solver->addClause(c);
}

void
aigpp::FRAIGManager::addClauses(std::size_t node)
{
    double startTime = lrabs::cpuTime();

    static lrabs::SimpleStack<std::size_t> pending;

    pending.push(0);
    std::size_t n = node;
    while (!pending.empty()) {
        Minisat::Var& satvar_n = _satVars[n];

        if (satvar_n != var_Undef) {
            n = pending.top();
            pending.pop();
        } else if (_nodes[n]._var) {
            satvar_n = _solver->newVar();

#ifdef FRAIGMANAGER_SIMP_ENABLE
            _solver->setFrozen(satvar_n, true);
#endif

            n = pending.top();
            pending.pop();
        } else if (_equivalences[n] != ~0ul) {
            std::size_t e  = _equivalences[n] >> 1;
            bool        en = (_equivalences[n] & 1) != 0;

            /* constant */
            if (e == 0) {
                satvar_n = _solver->newVar();

#ifdef FRAIGMANAGER_SIMP_ENABLE
                _solver->setFrozen(satvar_n, true);
#endif

                addUnit(Minisat::mkLit(satvar_n, !en));

                n = pending.top();
                pending.pop();
            }
            /* non constant, eq node not processed yet */
            else if (_satVars[e] == var_Undef) {
                pending.push(n);
                n = e;
            }
            /* non constant, eq node processed */
            else {
                satvar_n = _solver->newVar();

#ifdef FRAIGMANAGER_SIMP_ENABLE
                _solver->setFrozen(satvar_n, true);
#endif

                addBinary(Minisat::mkLit(satvar_n, true), Minisat::mkLit(_satVars[e], en));
                addBinary(Minisat::mkLit(satvar_n), Minisat::mkLit(_satVars[e], !en));

#ifdef FRAIGMANAGER_SIMP_ENABLE
                // decRefCount( e );
#endif

                n = pending.top();
                pending.pop();
            }
        } else if (_satVars[_nodes[n]._p1] == var_Undef) {
            pending.push(n);
            n = _nodes[n]._p1;
        } else if (_satVars[_nodes[n]._p2] == var_Undef) {
            pending.push(n);
            n = _nodes[n]._p2;
        } else {
            satvar_n = _solver->newVar();

#ifdef FRAIGMANAGER_SIMP_ENABLE
            _solver->setFrozen(satvar_n, true);
#endif

            Minisat::Lit vn  = Minisat::mkLit(satvar_n, _nodes[n]._nand);
            Minisat::Lit vn1 = Minisat::mkLit(_satVars[_nodes[n]._p1], _nodes[n]._neg1);
            Minisat::Lit vn2 = Minisat::mkLit(_satVars[_nodes[n]._p2], _nodes[n]._neg2);

            addBinary(~vn, vn1);
            addBinary(~vn, vn2);
            addTernary(vn, ~vn1, ~vn2);

#ifdef FRAIGMANAGER_SIMP_ENABLE
            decRefCount(_nodes[n]._p1);
            decRefCount(_nodes[n]._p2);
#endif

            n = pending.top();
            pending.pop();
        }
    }

    _cnfCreationTime += lrabs::cpuTime() - startTime;
}

void
aigpp::FRAIGManager::propagate()
{
    double startTime = lrabs::cpuTime();

    aigpp::VariableAssignment assignment;

    static lrabs::Random rng;

    /* extract counter example from solver */
    for (std::size_t i = 0; i != _varCount; ++i) {
        unsigned long& s = _sim[i + 1];

        const Minisat::Var& v = _satVars[i + 1];

        Minisat::lbool b = Minisat::l_Undef;
        if (v != var_Undef) {
            b = _solver->model[v];
        }

        // is this variable present in the counter example?
        if (b == Minisat::l_True) {
            s = ~0ul;
            assignment.setPositive(_aigNodes[i + 1]->varIndex());
        } else if (b == Minisat::l_False) {
            s = 0ul;
            assignment.setNegative(_aigNodes[i + 1]->varIndex());
        } else {
            s = rng.getULong();
        }

        /* flip one random bit (except most significant bit (bit #31 for 32 bit
         * numbers) -> used to store the original counter example) */
        s ^= (1ul << rng.getUInt(sizeof(unsigned int) * CHAR_BIT - 1));
    }

    /* learn counter example */
    _manager->addCounterExample(assignment);

    /* propagate counter example */
    for (std::size_t i = _varCount + 1; i != _size; ++i) {
        const NodeStruct& n = _nodes[i];

        unsigned long s1 = _sim[n._p1] ^ (n._neg1 ? ~0ul : 0ul);
        unsigned long s2 = _sim[n._p2] ^ (n._neg2 ? ~0ul : 0ul);
        s1 &= s2;

        _sim[i] = s1 ^ (n._nand ? ~0ul : 0ul);
    }

    _propTime += lrabs::cpuTime() - startTime;
}

void
aigpp::FRAIGManager::split()
{
    double startTime = lrabs::cpuTime();

    std::list<EquivalenceClass*>::const_iterator oldend = _nonEmptyClasses.end();
    static std::list<EquivalenceClass*>          addedClass;

    for (auto pe = _nonEmptyClasses.begin(); pe != oldend;
         /**/) {
        EquivalenceClass* e = *pe;

        /* init splitmap with representant */
        unsigned long repr_s = _sim[e->_repr >> 1];
        if ((e->_repr & 1) != 0) {
            repr_s = ~repr_s;
        }

        /* total members of eq class >= 3 -> general handling*/
        if (e->_nodes.size() >= 2) {
            _splitMap.clear();
            _splitMap.insert(SplitMap::value_type(repr_s, e));

            for (auto member = e->_nodes.begin(); member != e->_nodes.end();
                 /**/) {
                unsigned long s = _sim[(*member) >> 1u];
                if (((*member) & 1) != 0) {
                    s = ~s;
                }

                if (s == repr_s) {
                    ++member;
                } else {
                    auto eqclassit = _splitMap.find(s);
                    if (eqclassit != _splitMap.end()) {
                        if (eqclassit->second->_nodes.empty()) {
                            _nonEmptyClasses.push_back(eqclassit->second);
                        }

                        eqclassit->second->_nodes.push_back(*member);
                    } else {
                        EquivalenceClass* eqclass = allocateEquivalenceClassObject(*member);
                        addedClass.push_back(eqclass);

                        _splitMap.insert(SplitMap::value_type(s, eqclass));
#ifdef FRAIGMANAGER_TIMEOUT
                        _reduced[*member >> 1] = true;
#endif
                    }

                    e->_nodes.erase(member++);
                }
            }

            if (!e->_nodes.empty()) {
#ifdef FRAIGMANAGER_DEPTHSORT
                /* update depth */

                if ((e->_repr >> 1u) == 0u) {
                    /* zero class */
                    std::size_t firstCandidate = ~0ul;

                    for (std::size_t n : e->_nodes) {
                        if (_unreduced[n >> 1u]) {
                            firstCandidate = n >> 1u;
                            break;
                        }
                    }

                    if (firstCandidate == ~0ul) {
                        e->_depth = 0;
                        ++pe;
                    } else {
                        e->_depth = _depths[firstCandidate];
                        ++pe;
                    }
                } else {
                    /* eq (non-zero)  class */
                    std::size_t firstCandidate1 = e->_repr >> 1;
                    std::size_t firstCandidate2 = ~0ul;

                    for (std::size_t n : e->_nodes) {
                        if (_unreduced[firstCandidate1]) {
                            firstCandidate2 = n >> 1u;
                            break;
                        } else {
                            if (_unreduced[n >> 1u]) {
                                firstCandidate2 = n >> 1u;
                                break;
                            }
                        }
                    }

                    if (firstCandidate2 == ~0ul) {
                        e->_depth = _depths[firstCandidate1];
                        ++pe;
                    } else {
                        e->_depth = _depths[firstCandidate2];
                        ++pe;
                    }
                }
#else  /* FRAIGMANAGER_DEPTHSORT */
                ++pe;
#endif /* FRAIGMANAGER_DEPTHSORT */
            } else {
#ifdef FRAIGMANAGER_SIMP_ENABLE
                markNodeAsDone(e->_repr >> 1u);
#endif

                releaseEquivalenceClassObject(e);
                _nonEmptyClasses.erase(pe++);
            }

            for (EquivalenceClass* eq : addedClass) {
                assert((eq->_repr >> 1u) != 0u);

#ifdef FRAIGMANAGER_DEPTHSORT
                if (!eq->_nodes.empty()) {
                    /* update depth */
                    /* find first candidate pair for eq checks */

                    std::size_t firstCandidate1 = eq->_repr >> 1u;
                    std::size_t firstCandidate2 = ~0ul;

                    for (std::size_t n : e->_nodes) {
                        if (_unreduced[firstCandidate1]) {
                            firstCandidate2 = n >> 1u;
                            break;
                        } else {
                            if (_unreduced[n >> 1u]) {
                                firstCandidate2 = n >> 1u;
                                break;
                            }
                        }
                    }

                    if (firstCandidate2 == ~0ul) {
                        eq->_depth = _depths[firstCandidate1];
                    } else {
                        eq->_depth = _depths[firstCandidate2];
                    }
                } else {
#    ifdef FRAIGMANAGER_SIMP_ENABLE
                    markNodeAsDone(eq->_repr >> 1u);
#    endif

                    releaseEquivalenceClassObject(eq);
                }
#else  /* FRAIGMANAGER_DEPTHSORT */
                if (eq->_nodes.empty()) {
                    releaseEquivalenceClassObject(eq);
                }
#endif /* FRAIGMANAGER_DEPTHSORT */
            }
            addedClass.clear();
        }
        /* total members of eq class = 2 -> special handling */
        else {
            unsigned long s = _sim[*(e->_nodes.begin()) >> 1u];
            if ((*(e->_nodes.begin()) & 1u) != 0) {
                s = ~s;
            }

            /* split! -> two one-element classes */
            if (s != repr_s) {
#ifdef FRAIGMANAGER_TIMEOUT
                _reduced[*(e->_nodes.begin()) >> 1u] = true;
#endif

#ifdef FRAIGMANAGER_SIMP_ENABLE
                markNodeAsDone(e->_repr >> 1u);
                markNodeAsDone(*(e->_nodes.begin()) >> 1u);
#endif

                e->_nodes.clear();
                releaseEquivalenceClassObject(e);
                _nonEmptyClasses.erase(pe++);
            }
            /* no split */
            else {
#ifdef FRAIGMANAGER_DEPTHSORT
                e->_depth = _depths[*e->_nodes.begin() >> 1u];
#endif
                ++pe;
            }
        }
    }

    _splitTime += lrabs::cpuTime() - startTime;
}
