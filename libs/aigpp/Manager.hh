/**************************************************************
 *
 *       AIGPP Package // Manager.hh
 *
 *       Copyright (C) 2006, 2007, 2008 Florian Pigorsch
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

#ifndef AIGPP_MANAGER_HH
#define AIGPP_MANAGER_HH

/* std */
#include <cassert>
#include <climits>
#include <deque>
#include <fstream>
#include <map>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/* CUDD */
#include <cudd.h>

/* Minisat */
#include <core/Solver.h>

/* LRABS utilities */
#include <lrabsutil/Random.hh>

/* local */
#include "ComputedTable.hh"
#include "Edge.hh"
#include "EdgeRef.hh"
#include "ExtRefTable.hh"
#include "FRAIGManager.hh"
#include "InternalEdgeRef.hh"
#include "Node.hh"
#include "NodeHashMap.hh"
#include "RewritingManager.hh"
#include "Settings.hh"
#include "SimVector.hh"
#include "SimulationTable.hh"
#include "UniqueTable.hh"
#include "VariableAssignment.hh"

#define ENABLE_BDD

/*
#if defined(USE_MINICRAIG)
    #if defined(USE_MINICRAIG2)
        #error defined both USE_MINICRAIG and USE_MINICRAIG2
    #endif
#else
    #if not defined(USE_MINICRAIG2)
        #error defined neither USE_MINICRAIG nor USE_MINICRAIG2
    #endif
#endif
*/

namespace aigpp {
/*!
 * \class Manager
 * \brief The Manager class is used for managing the creation and deletion of
 * AIG nodes, and holds some statisticis about the represented AIG.
 */
class Manager
{
   public:
    /*!
     * \ingroup ext
     * Constructor. Creates an empty Manager.
     */
    Manager();

    Manager(const Manager&) = delete;
    Manager(Manager&&)      = delete;
    Manager& operator=(const Manager&) = delete;
    Manager& operator=(Manager&&) = delete;

    /*!
     * \ingroup ext
     * Destructor.
     */
    ~Manager();

    /*!
     * \ingroup ext
     * Returns the Settings object of this Manager.
     */
    Settings& settings();
    /*!
     * \ingroup ext
     * Returns the Settings object of this Manager.
     */
    const Settings& settings() const;

    /*!
     * \ingroup ext
     * Sets the verbosity to level \a v. A level of 0 represents the "quiet" mode.
     * \sa verbosity() const
     */
    void setVerbosity(int v);
    /*!
     * \ingroup ext
     * Returns the verbosity level.
     * \sa setVerbosity( int )
     */
    int verbosity() const;

    /*** constant access ****************/
    /************************************/

    /*!
     * \ingroup ext
     * Get an EdgeRef representing the constant 0 function, i.e. the contradiction
     * function. \sa getConst1() const
     */
    const EdgeRef& getConst0() const;

    /*!
     * \ingroup ext
     * Get an EdgeRef representing the constant 1 function, i.e. the tautology
     * function. \sa getConst0() const
     */
    const EdgeRef& getConst1() const;

    /**Rewriting Procedures**/
    void rewrite();
    void rewriteCone(Node* root);

    /*** variable access ****************/
    /************************************/

    /*!
     * \ingroup ext
     * Create a new variable node and return an uncomplemented EdgeRef pointing to
     * the new variable node. Assign the \a name to the new variable. \a name must
     * either be empty or unique. If \a name is empty, the Manager automatically
     * creates a generic name for the variable.
     * \sa lookupVariable( const std::string& )
     */
    EdgeRef addVariable(const std::string& name = "");

    /*!
     * \ingroup ext
     * Return an EdgeRef pointing to the variable node with the given \a name.
     * A variable node with the \a name must exist!
     * \sa addVariable( const std::string& )
     */
    EdgeRef lookupVariable(const std::string& name);

    /*!
     * \ingroup ext
     * Checks if a variable with the given name exists.
     * \sa addVariable( const std::string& )
     * \sa lookupVariable(const std::string& )
     */
    bool hasVariableName(const std::string& name) { return _varMap.find(name) != _varMap.end(); }

    /*!
     * \ingroup ext
     * Returns the number of variables present in the AIG Manager.
     */
    std::size_t variableCount() const;

    /*!
     * \ingroup ext
     * Checks if the given integer is a valid variable index, i.e.
     * positive and less than variableCount().
     *
     * \param index
     *
     * \return true if valid
     */
    bool isVarIndexValid(std::size_t index) const;

    /*!
     * \ingroup ext
     * Returns a non-inverted EdgeRef pointing to the variable with the given \a
     * index. \a index must be a valid variable index! \sa isVarIndexValid( int )
     * const, EdgeRef::variableIndex() const
     */
    EdgeRef variable(std::size_t index) const;

    /*!
     * \ingroup ext
     * Returns the name of the variable with the given \a index.
     * \a index must be a valid variable index!
     *
     * \sa EdgeRef::variableName() const, EdgeRef::variableIndex() const,
     * isVarIndexValid( int ) const,
     */
    const std::string& variableName(std::size_t index) const;

    /*** statistics on the geometry *****/
    /************************************/
    /*!
     * \ingroup ext
     * Return the current number of existing nodes.
     */
    std::size_t nodeCount() const;
    std::size_t deadNodeCount() const;

    /*!
     * \ingroup ext
     * Returns the number of nodes in the shared cones of the given \a roots.
     *
     * \sa EdgeRef::nodeCount() const
     */
    std::size_t sharedSize(const std::vector<EdgeRef>& roots) const;

    /*!
     * \ingroup ext
     * Returns the shared support of the given \a roots in form of a vector of
     * variable indices.
     *
     * \sa EdgeRef::support() const
     */
    std::vector<int> sharedSupport(const std::vector<EdgeRef>& roots) const;
    std::vector<int> sharedSupport(const std::vector<Edge>& roots) const;
    std::vector<int> sharedSupport(const std::vector<Node*>& roots) const;

    /*!
     * Compares the structural quality of the two given nodes. Returns
     * 0 if both nodes have equal quality,
     * a value < 0 if \a n1 is better than \a n2,
     * a value > 0 if \a n2 is better than \a n1.
     */
    long compareQuality(Node* n1, Node* n2) const;

    /*!
     * \ingroup ext
     * Print some statistics to the given output stream \a os.
     */
    void printStats(std::ostream& os = std::cout) const;

    /*!
     * Access the list of nodes
     */
    Node* nodesList();

    /*** operations *********************/
    /************************************/
    /*!
     * \ingroup ext
     * Substitute the \a variables by the corresponding \a replacements in
     * the given \a functions, and return the resulting functions. The
     * \a variables vector must only contain references to distinct,
     * unnegated variable nodes.  The \a variables and \a replacements
     * vectors must have equal size.
     * \sa EdgeRef::compose( const std::vector<EdgeRef>&, const
     * std::vector<EdgeRef>& ) const
     */
    std::vector<EdgeRef> vectorCompose(const std::vector<EdgeRef>& functions, const std::vector<EdgeRef>& variables,
                                       const std::vector<EdgeRef>& replacements);

    /*** fraig **************************/
    /************************************/
    /*!
     * @ingroup ext
     * Enable or disable automatic FRAIGing. If it is \a enabled, equivalent AIG
     * nodes are merged on-thy-fly during construction.
     */
    void toggleAutomaticFRAIGing(bool enable);

    /*!
     * @ingroup ext
     * Return true if automatic FRAIGing is enabled.
     */
    bool automaticFRAIGing() const;

    /*!
     * @ingroup ext
     * Merge equivalent nodes in the AIG.
     */
    void makeFRAIG(double timeout = -1);

    /*!
     * Enable/disable the creation of SimVectors for new nodes
     *
     * If disabled, simulation vectors are only created, when actually used.
     * If enabled, new nodes are created with SimVectors.
     *
     * If newly enabled, all nodes that do not have a SimVector are updated.
     */
    void toggleSimCreation(bool on);
    bool simCreation() const;

    /*!
     * Recompute simulation vectors for all nodes (except input nodes)
     */
    void recomputeSimulation(bool warnings = false);

    std::vector<DdNode*> computeBDD(const std::vector<EdgeRef>& roots, DdManager* bddMan,
                                    const std::map<int, DdNode*>& varMapping) const;

    enum BDDSweepingResult
    {
        BDDSweeping_Skipped,
        BDDSweeping_Success,
        BDDSweeping_BDDCouldNotBeConstructed,
        BDDSweeping_BDDTooBig
    };

    /*!
     * Perform BDD sweeping based AIG restructuring.
     *
     * \param bigNode root of the AIG cone to restructure
     *
     * \return status of the operation
     */
    BDDSweepingResult BDDSweep(const EdgeRef& bigNode, bool force = false);
    BDDSweepingResult BDDSweepInternal(const InternalEdgeRef& bigNode, bool force = false);

    void toggleBDDSweepingDuringQuantification(bool enable);

    void existentialQInplace(InternalEdgeRef& root, const InternalEdgeRef& var, bool internalCall = false);
    void existentialQInplace(InternalEdgeRef& root, const std::vector<InternalEdgeRef>& vars);
    void universalQInplace(InternalEdgeRef& root, const std::vector<InternalEdgeRef>& vars);

    EdgeRef bddQuantify(const EdgeRef& root, const std::vector<EdgeRef>& variables, bool existential, bool preferAIG,
                        double timeout = -1, bool* timeoutReached = nullptr);
    void    bddQuantifyInplace(EdgeRef& root, const std::vector<EdgeRef>& variables, bool existential, bool preferAIG,
                               bool optimizeLast, double timeout = -1, bool* timeoutReached = nullptr);
    std::vector<long> estimateRemainingNodesAfterQuantification(Node* root, const std::vector<Node*>& variables) const;

    void derefBDD(DdNode* bdd)
    {
        assert(_bddManager != nullptr);
        assert(bdd != nullptr);
        Cudd_RecursiveDeref(_bddManager, bdd);
    }

    struct HybridQuantifyResult
    {
        HybridQuantifyResult() : aig(), bdd(nullptr) {}

        ~HybridQuantifyResult()
        {
            if (bdd != nullptr) {
                std::cerr << "WARNING: HybridQuantifyResult has bdd != 0 on deletion\n";
                assert(false);
            }
        }

        EdgeRef aig;
        DdNode* bdd;
    };

    void hybridQuantify(HybridQuantifyResult& baseFunction, EdgeRef& conjunctiveFunction,
                        const std::vector<EdgeRef>& variables, bool existential, bool optimizeLast);
    void hybridQuantify(HybridQuantifyResult& baseFunction, const std::vector<EdgeRef>& variables, bool existential,
                        bool optimizeLast);

    int              selectBestVariableToQuantify(const std::set<int>& variables, const EdgeRef& root);
    std::vector<int> computeVariableOrdering(const std::set<int>& variables, const EdgeRef& root);
    EdgeRef          hybridQuantify2(EdgeRef root, const std::vector<EdgeRef>& variables);

    EdgeRef quantifyDeep(const EdgeRef& root, int varIndex, bool existential);
    struct QuantifyDeepData
    {
        InternalEdgeRef cof0, cof1, pos_ex, neg_ex;
    };
    typedef std::unordered_map<Node*, QuantifyDeepData, lrabs::hashPointer<Node>> QuantifyDeepDataMap;
    InternalEdgeRef quantifyDeepRecursive(const Edge& root, QuantifyDeepDataMap& cache);
    InternalEdgeRef quantifyDeepCofactor(const Edge& root, QuantifyDeepDataMap& cache);

    void collectMultiAndInputs(Node* root, std::vector<Edge>& inputs) const;
    bool collectXORInputs(Node* root, std::vector<Edge>& inputs) const;

#if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
    std::pair<InternalEdgeRef, InternalEdgeRef> quantifierEliminationBySubstitution(const InternalEdgeRef& root,
                                                                                    const InternalEdgeRef& var,
                                                                                    std::size_t wantedNodes = 0,
                                                                                    double      timeout     = -1,
                                                                                    bool*       failed      = nullptr);
    InternalEdgeRef dcmin(const InternalEdgeRef& function, const InternalEdgeRef& care_set, std::size_t wantedNodes = 0,
                          double timeout = -1, bool* failed = nullptr);
#endif

    void cleanBDD();

    /*!
     *   Returns a satisfying assignment to the variables of \e e,
     *   i.e. an assignment such that \e e's function evaluates to
     *   true. If \e minimize is true, the assignment is minimized.
     *
     *   \param e
     *   \param minimize
     *
     *   \return satisfying assignment
     */
    std::pair<bool, VariableAssignment> satisfyingAssignment(const Edge& e, bool minimize);

    /* export methods */
    enum AIGERMode
    {
        AIGER_BinaryMode = 0,
        AIGER_ASCIIMode  = 1
    };

    /*!
     * Exports the AIG to AIGER format (http://fmv.jku.at/aiger/).
     *
     * \param os stream to export to
     * \param mode AIGER export mode
     */
    void exportAIGER(std::ostream& os, AIGERMode mode = AIGER_ASCIIMode) const;
    void exportRootsAIGER(std::ostream& os) const;
    void exportAIGER(std::ostream& os, const std::vector<EdgeRef>& outputs, const std::vector<std::string>& outputNames,
                     AIGERMode mode = AIGER_ASCIIMode) const;
    void exportAIGER(std::ostream& os, const std::vector<Edge>& outputs,
                     const std::vector<std::string>& outputNames) const;
    void exportAIGERBinary(std::ostream& os, const std::vector<Edge>& outputs,
                           const std::vector<std::string>& outputNames) const;

    void writeBinary(std::ostream& os, size_t delta) const;

    void exportSimulation(std::ostream& os) const;

    void exportSeparationMatrix(std::ostream& os, int simvalues = -1);

    /*!
     * Exports the whole AIG as a BLIF model. The model's inputs are the AIGs
     * variables, the model's outputs are the external references (EdgeRef) to the
     * AIG. The model doesn't contain any latches.
     */
    void exportBLIF(std::ostream& os) const;

    /*!
     * \ingroup ext
     * Dumps the complete graph structure onto the output stream \a os.
     * The output format is the <a href="http://www.graphviz.org/">GraphViz</a>
     * format. Variable nodes are displayed as boxes, AND nodes are displayed as
     * circles (dotted outline in case of a NAND node), external references are
     * displayed as "houses". Edges are represented by arrows between nodes,
     * inverted edges are drawn in dotted style.
     */
    void dumpGraphViz(std::ostream& os, bool displayNodeIndexes = false) const;
    /*!
     * \ingroup ext
     * Dumps the graph structure rooted in the given \link EdgeRef EdgeRefs
     * \endlink (\a roots) onto the output stream \a os. The output format is the
     * <a href="http://www.graphviz.org/>GraphViz</a> format. Variable nodes are
     * displayed as boxes, AND nodes are displayed as circles (dotted outline in
     * case of a NAND node), the \a roots are displayed as "houses". Edges are
     * represented by arrows between nodes, inverted edges are drawn in dotted
     * style.
     */
    void dumpGraphViz(std::ostream& os, const std::vector<EdgeRef>& roots, bool displayNodeIndexes = false) const;
    void dumpGraphViz(std::ostream& os, const std::vector<InternalEdgeRef>& roots,
                      bool displayNodeIndexes = false) const;
    /*!
     * \ingroup ext
     * Dumps the graph structure rooted in the given EdgeRef (\a root) onto the
     * output stream \a os. The output format is the <a
     * href="http://www.graphviz.org/>GraphViz</a> format. Variable nodes are
     * displayed as boxes, AND nodes are displayed as circles (dotted outline in
     * case of a NAND node), the \a root is displayed as a "house". Edges are
     * represented by arrows between nodes, inverted edges are drawn in dotted
     * style.
     */
    void dumpGraphViz(std::ostream& os, const EdgeRef& root, bool displayNodeIndexes = false) const;
    void dumpGraphViz(std::ostream& os, const InternalEdgeRef& root, bool displayNodeIndexes = false) const;

    void exportGraphML(std::ostream& os, const std::vector<EdgeRef>& roots) const;
    void exportGML(std::ostream& os, const std::vector<EdgeRef>& roots) const;

    /*!
     * Creates an CNF represenation of the function given by the AIG edge "root"
     * by applying the tseitin transformation to the cone of "root". "mapping" is
     * an unsigned int -> EdgeRef map, which contains an entry for at least each
     * variable occuring in the cone of "root", "rootLiteral" is set to the
     * literal which is equivalent to the given "root", "nextFreeIndex" specifies
     * the next index to be used for enconding (tseitin) variables. setting
     * "updateMapping" to true, adds encoding variables to the mapping. if "root"
     * is (structurally) equivalent to "true", an empty clauseset is returned, if
     * "root" is (structurally) equivalent to "false", a clauseset containing the
     * empty clause only is returned.
     */
    template <typename T>
    std::vector<std::vector<T>> createCNF(const EdgeRef& root, std::map<T, EdgeRef>& mapping, T& rootLiteral,
                                          T& nextFreeIndex, bool updateMapping = false) const;

    std::vector<std::vector<int>> createCNF(const EdgeRef& root, int& rootLiteral) const;

    void checkIntegrity() const;

    void garbageCollect();

    int  modifiedNodes() const;
    void resetModified();

    ExtRefTable* getExtRefTable() const { return _extRefTable; }

    void addCounterExample(const VariableAssignment& counterex);

    InternalEdgeRef DebugAnd(const InternalEdgeRef& n1, const InternalEdgeRef& n2, bool createNAND, bool doTrivial,
                             bool doUnique, bool doSAT);
    EdgeRef         DebugGetExternal(const InternalEdgeRef& e) { return EdgeRef(_extRefTable, e); }

    std::size_t unreducedNodes() const { return _unreducedNodes; }

   private:
    void ref(Node* n);
    void deref(Node* n);

    bool AndTrivialCase(const InternalEdgeRef& n1, const InternalEdgeRef& n2, InternalEdgeRef& result);

    /*!
     * Creates an edge to an and-node representing the conjunction of
     * \e n1 and \e n2.
     * The method first checks the \e unique table for an existing
     * isomorphic node and immediately returns if this check is successful.
     * Then some trivial rewriting rules are applied (see AndTrivialCase()).
     * If automatic FRAIGing is enabled (see toggleAutomaticFRAIGing()),
     * a set of equivalent candidate nodes is build using simulation vectors.
     * For each candidate a SAT checkis run to prove equivalence with the new
     * node. If an equivalent node is found, the method returns an edge to this
     * node. Finally, after all checks failed, a new node is created and an edge
     * to this node is returned.
     */
    InternalEdgeRef And(const InternalEdgeRef& n1, const InternalEdgeRef& n2);

    /*!
     * Creates an edge to an and-node representing the conjunction of
     * \e n1 and \e n2.
     * Unlike the And() method, this method does neither perform rewriting nor
     * SAT checking, but only checks for isomorphic nodes.
     */
    InternalEdgeRef SimpleAnd(const InternalEdgeRef& n1, const InternalEdgeRef& n2);

    /*!
     * Create an And-node object (or take a free node from the _freeNodes stack)
     */
    Node* createNode(const Edge& p1, const Edge& p2);

    /*!
     * Create a Variable-node object (or take a free node from the _freeNodes
     * stack)
     */
    Node* createNode(int varindex);

    /*!
     * Push the Node n onto the _freeNodes stack
     */
    void releaseNode(Node* n);

    /*!
     * Create a SimVector (or take a free SimVector from the _freeSim stack)
     */
    SimVector* createSim();

    /*!
     * Push the SimVector s onto the _freeSim stack
     */
    void releaseSim(SimVector* s);

    bool               satEqual(Node* n1, Node* n2);
    bool               satEqualNegated(Node* n1, Node* n2);
    bool               satEqualZero(Node* n1);
    bool               satEqualOne(Node* n1);
    bool               isBoolRedundant(Node* node, Node* var) const;
    std::vector<Node*> boolRedundant(Node* node, const std::vector<Node*>& candidates) const;

    void resetSolver(bool force = false);

    void createClauses(Node* n1);
    void createClauses(Node* n1, Node* n2);

    lrabs::Random& random() { return _random; };

    const SimulationTable& simtable() const;
    SimulationTable&       simtable();

    void updateSimTable(bool force = false);
    void updateSimTable(const VariableAssignment& assignment);

#ifdef USE_TEMPSIM
    /*
     * set the temp sim of the variables as induced by the assignment and
     * propagate through the aig. also update the additionalNode.
     */
    void updateTempSim(const VariableAssignment& assignment, Node* additionalNode = nullptr);

    /*
     * propagate temp sim through the aig if the aig's temp sim is not uptodate.
     * also update the additionalNode.
     */
    void propagateTempSim(Node* additionalNode = nullptr);
    void removeUnequalTempSim(std::set<Node*>& list, unsigned long tempSim) const;
#endif

    void initializeCandidates(std::set<Node*>& list, aigpp::Node* simclass
#ifdef USE_TEMPSIM
                              ,
                              unsigned long tempSim
#endif
                              ) const;
    void updateCandidates(std::set<Node*>& list, Node* newsimclass) const;

    VariableAssignment getCounterExample() const;
    void               addRandomCounterExample();
    void               mergeCounterExamples();

    void clearNodeCaches();
    void clearNodeCachesZ();

    void buildRecursiveUsingCache(Node* n);
    void buildCofactorRecursiveUsingCache(Node* n, const Edge& var);
    void buildRecursiveZ(Node* n);

    DdNode*         buildBDDfromAIG(Node* aig);
    DdNode*         buildBDDfromAIGlimited(Node* aig, double growthLimit);
    void            buildAIGfromBDD(DdNode* bdd, Node* originalAIG);
    InternalEdgeRef createAIG(DdNode* bdd);

   private:
    std::vector<InternalEdgeRef> vectorCompose(const std::vector<InternalEdgeRef>& functions,
                                               const std::vector<InternalEdgeRef>& variables,
                                               const std::vector<InternalEdgeRef>& replacements);

    std::size_t     sharedSize(const std::vector<Node*>& roots) const;
    InternalEdgeRef variableInternal(std::size_t index) const;
    InternalEdgeRef addVariableInternal(const std::string& name = "");

    /* !
     *  Reclaim node, i.e. bring a dead node with refcount=0 back to life.
     *  Must not be called with a NULL pointer.
     *  The node's refCount() must be 0.
     */
    void reclaim(Node* n);

    /* !
     *  Reclaim node if it is dead. May be called with a NULL pointer.
     */
    void reclaimIfDead(Node* n);

    void recursiveDeref(Node* n);

    void replace(Node* oldNode, Node* newNode, bool inverted);

    void addToNodesList(Node* n);
    void removeFromNodesList(Node* n);
    void fixNodesList();
    void checkNodesList() const;

    template <unsigned int FLAGS>
    bool flagsLocked() const;

    template <unsigned int FLAGS>
    void lockFlags() const;

    template <unsigned int FLAGS>
    void unlockFlags() const;

    int nextIndex();

   private:
    ExtRefTable* _extRefTable;

    std::size_t _unreducedNodes;
    bool        _automaticFRAIGing;

    int _verbosity;

    std::size_t _nodeCount;
    int         _nextIndex;

    std::size_t _dead, _reclaimed;

    std::stack<Node*>      _freeNodes;
    std::stack<SimVector*> _freeSim;
    bool                   _simCreation;

   public:
    std::vector<Node*>         _variables;
    std::map<std::string, int> _varMap;
    std::vector<std::string>   _variableNames;

   private:
    Settings _settings;

    SimulationTable _simTable;
    int             _nextUpdateSimVectorBin;

    std::deque<VariableAssignment> _newCounterExamples, _counterExamples;

    lrabs::Random _random;
#ifdef USE_TEMPSIM
    lrabs::Random _randomTempSim;   /* rng for temp sim */
    bool          _tempSimUpToDate; /* temp sim contains valid simulation data */
#endif

    UniqueTable _unique;
#ifdef USE_COMPUTEDTABLE
    ComputedTable _computed;
#endif

    Node* _nodes;
    Node* _lastNode;

    typedef NodeHashMap<InternalEdgeRef, Node::FLAG_CACHEVALID, Node::FLAG_CACHECONST0 | Node::FLAG_CACHECONST1>
             CacheMap;
    CacheMap _cache;

    typedef NodeHashMap<std::pair<InternalEdgeRef, InternalEdgeRef>, Node::FLAG_ZCACHEVALID> PairCacheMap;
    PairCacheMap                                                                             _cacheZ;

    int _toBeRemoved;

    Minisat::Solver* _solver;
    // minisat::Solver* _solver;

    bool            _noReset;
    int             _satChecksPresent;
    int             _nodesInCNF;
    int             _deletedNodesInCNF;
    int             _equivalentNodesInCNF;
    int             _mitersInCNF;
    int             _openMitersInCNF;
    int             _closedMitersInCNF;
    std::set<Node*> _nodesWithSATVars;

    mutable std::stack<const Node*>  _universalConstNodeStack;
    mutable std::vector<const Node*> _universalConstNodeVector;
    mutable std::stack<Node*>        _universalNodeStack;
    mutable std::vector<Node*>       _universalNodeVector;

    /* bdd sweeping stuff */
    bool  _nextAndSkipSat;
    Node* _nextAndSatResult;
    bool  _nextAndSatResultIsInverted;
    bool  _bddSweepingDuringQuantification;

#ifdef ENABLE_BDD
    DdManager*  _bddManager;
    std::size_t _nextReductionSize;
    double      _BDDSimDecay;
    double      _bddreducedratio;
    bool        _bddReachedNodeLimit;
    bool        _bddReachedTimeLimit;
    bool        _bddReachedGrowthLimit;
#endif

    std::stack<Node*> _recursiveNodeStack;

    Node* _skipTestWith;

    bool _inCofactor;

    mutable unsigned int _flagsLocked;

    EdgeRef _const0, _const1;

    RewritingManager _rwManager;
    FRAIGManager     _fraigManager;
    bool             _timeoutDuringLastFRAIG;

    friend class Node;
    friend struct Edge;
    friend class InternalEdgeRef;
    friend class EdgeRef;
    friend class RewritingManager;
    friend class FRAIGManager;
};

}  // namespace aigpp

#include "Manager.icc"
#include "Manager_export.icc"

#endif /* AIGPP_MANAGER_HH */
