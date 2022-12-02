/**************************************************************
 *
 *       AIGPP // FRAIGManager.hh
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

#ifndef AIGPP_FRAIGMANAGER_HH
#define AIGPP_FRAIGMANAGER_HH

#include <iostream>
#include <list>
#include <stack>
#include <unordered_map>

#include <lrabsutil/Hashes.hh>
#include <lrabsutil/SimpleStack.hh>

/* minisat2 */
#include "core/Solver.h"
#include "core/SolverTypes.h"
#include "simp/SimpSolver.h"

#define FRAIGMANAGER_TIMEOUT
#define FRAIGMANAGER_DEPTHSORT

//#define FRAIGMANAGER_SIMP_ENABLE
#define FRAIGMANAGER_SIMP_DISABLE

#if defined(FRAIGMANAGER_SIMP_ENABLE) and defined(FRAIGMANAGER_SIMP_DISABLE)
#    error "FRAIGManager: simplification both enabled and disabled"
#elif not defined(FRAIGMANAGER_SIMP_ENABLE) and not defined(FRAIGMANAGER_SIMP_DISABLE)
#    error "FRAIGManager: simplification neither enabled nor disabled"
#endif

#define FRAIGMANAGER_RESET_ENABLE
//#define FRAIGMANAGER_RESET_DISABLE

#if defined(FRAIGMANAGER_RESET_ENABLE) and defined(FRAIGMANAGER_RESET_DISABLE)
#    error "FRAIGManager: resets both enabled and disabled"
#elif not defined(FRAIGMANAGER_RESET_ENABLE) and not defined(FRAIGMANAGER_RESET_DISABLE)
#    error "FRAIGManager: resets neither enabled nor disabled"
#endif

namespace aigpp {
class Manager;
class Node;
}  // namespace aigpp

namespace aigpp {
class FRAIGManager
{
   public:
    explicit FRAIGManager(Manager* manager);
    FRAIGManager(const FRAIGManager&) = delete;
    FRAIGManager(FRAIGManager&&)      = delete;
    FRAIGManager& operator=(const FRAIGManager&) = delete;
    FRAIGManager& operator=(FRAIGManager&&) = delete;
    ~FRAIGManager();

    void makeFRAIG(
#ifdef FRAIGMANAGER_TIMEOUT
        double timeout = -1
#endif
    );

    void printStats(std::ostream& os = std::cout) const;

   private:
    void initStructures();
    void initClasses();
    void clearStructures();

    struct NodeStruct
    {
        bool _var : 1;
        bool _nand : 1;

        bool _neg1 : 1;
        bool _neg2 : 1;

        std::size_t _p1;
        std::size_t _p2;
#ifdef FRAIGMANAGER_SIMP_ENABLE
        int  _ref;
        int  _originalRef;
        bool _done;
#endif
    };

#ifdef FRAIGMANAGER_SIMP_ENABLE
    void markNodeAsDone(std::size_t index);
    void decRefCount(std::size_t index);
    void incRefCount(std::size_t index);
#endif

    struct EquivalenceClass
    {
        EquivalenceClass() :
            _nodes(),
            _repr(-1)
#ifdef FRAIGMANAGER_DEPTHSORT
            ,
            _depth(0)
#endif
        {}

        explicit EquivalenceClass(std::size_t n) :
            _nodes(),
            _repr(n)
#ifdef FRAIGMANAGER_DEPTHSORT
            ,
            _depth(0)
#endif
        {}

        std::list<std::size_t> _nodes;
        std::size_t            _repr;
#ifdef FRAIGMANAGER_DEPTHSORT
        int _depth;
#endif
    };

    EquivalenceClass* allocateEquivalenceClassObject(std::size_t repr = ~0ul);
    void              releaseEquivalenceClassObject(EquivalenceClass* e);

    std::list<EquivalenceClass*>::iterator getEquivalenceClassToSplit();

    /* sat solver stuff */
    void resetSolver();
    void addUnit(const Minisat::Lit lit1);
    void addBinary(const Minisat::Lit lit1, const Minisat::Lit lit2);
    void addTernary(const Minisat::Lit lit1, const Minisat::Lit lit2, const Minisat::Lit lit3);
    void addClauses(std::size_t node);

    void propagate();
    void split();

   private:
    aigpp::Manager* _manager;

    /* mapping from aig nodes to node indices */
    typedef std::unordered_map<Node*, std::size_t, lrabs::hashPointer<Node>> IDMap;
    IDMap                                                                    _id;
    std::size_t                                                              _varCount;

    /* all of the following arrays have size "_capacity" */
    std::size_t _size;
    std::size_t _capacity;
    NodeStruct* _nodes;
    bool*       _unreduced;
#ifdef FRAIGMANAGER_TIMEOUT
    bool* _reduced;
#endif

    unsigned long* _sim;
    Minisat::Var*  _satVars;
    int*           _depths;
    std::size_t*   _equivalences;
    Node**         _aigNodes;

    std::list<EquivalenceClass*>                                 _nonEmptyClasses;
    lrabs::SimpleStack<EquivalenceClass*>                        _freeClasses;
    typedef std::unordered_map<unsigned long, EquivalenceClass*> SplitMap;
    SplitMap                                                     _splitMap;

#ifdef FRAIGMANAGER_SIMP_ENABLE
    Minisat::SimpSolver* _solver;
#else
    Minisat::Solver* _solver;
#endif

#ifdef FRAIGMANAGER_RESET_ENABLE
    std::size_t _checksSinceLastReset;
    std::size_t _newEquivalencesSinceLastReset;
#endif

    double      _initTime, _cnfCreationTime, _satTime, _satMaxTime, _propTime, _splitTime, _rebuildTime, _fraigTime;
    std::size_t _satChecks, _satEquiv, _satChecksConst, _satEquivConst;

    double _totalInitTime, _totalCnfCreationTime, _totalSatTime, _totalSatMaxTime, _totalPropTime, _totalSplitTime,
        _totalRebuildTime, _totalFraigTime;
    std::size_t _totalSatChecks, _totalSatEquiv, _totalSatChecksConst, _totalSatEquivConst;

    std::size_t _totalApplications;
};
}  // namespace aigpp

#endif /* AIGPP_FRAIGMANAGER_HH */
