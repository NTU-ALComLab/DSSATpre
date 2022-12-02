/**************************************************************
 *
 *       AIGPP // RewritingManager.hh
 *
 *       Copyright (C) 2008 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 649 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#ifndef AIGPP_REWRITINGMANAGER_HH
#define AIGPP_REWRITINGMANAGER_HH

/* stdlib/stl/.... */
#include <map>
#include <unordered_map>

/* lrabs */
#include <lrabsutil/Hashes.hh>
#include <lrabsutil/SimpleStack.hh>

/* local */
namespace aigpp {
class Manager;
}

#include "DynamicRewritingTable.hh"
#include "FourInputFunction.hh"
#include "NPNClass.hh"
#include "Node.hh"

namespace aigpp {
class RewritingManager
{
   public:
    explicit RewritingManager(Manager* m);
    ~RewritingManager();

    void    rewrite();
    void    rewriteCone(aigpp::Node* root);
    EdgeRef localRewriteCone(const EdgeRef& root);

    void printStats(std::ostream& os = std::cout) const;

   private:
    struct NodeStruct
    {
        NodeStruct() : var(false), nand(false), neg1(false), neg2(false), ref(0), p1(0), p2(0) {}

        bool var : 1;
        bool nand : 1;
        bool neg1 : 1;
        bool neg2 : 1;
        int  ref : (32 - 4);

        unsigned int p1;
        unsigned int p2;
    };

    int  lookupAnd(unsigned int p1, unsigned int p2) const;
    void growArrays();
    void refNode(unsigned int nodeindex);
    void derefNode(unsigned int nodeindex);

    struct CutStruct
    {
        unsigned int signature;
        unsigned int size;
        unsigned int nodes[8];
        CutStruct*   next;

        void print(std::ostream& os = std::cout) const;
        bool isSubsetOf(const CutStruct& bigger) const;
    };

    CutStruct*  computeCuts(unsigned int nodeindex);
    static bool insertCut(CutStruct** cuts, CutStruct* newcut);

    FourInputFunction                 getFunction(unsigned int nodeindex, const CutStruct& cut);
    const NPNClassAndTransformations& getNPNClass(const FourInputFunction& f) const;

    struct NPNImplInfo
    {
        int                        cutIndex;
        NPNClassAndTransformations npnclass;
        int                        implIndex;
        int                        delta;
    };

    NPNImplInfo getBestImplementation(unsigned int nodeindex, CutStruct* cuts);
    int         countNewNodes(int implindex);
    int         countDeletedNodes(unsigned int nodeindex, int implindex);
    int         countDeletedNodesSmallerInput(unsigned int nodeindex, int implindex);
    int         countDeletedNodesConstant(unsigned int nodeindex);

    void markRewritten(unsigned int nodeindex);

    void buildNewCone(unsigned int nodeindex, const CutStruct& cut, const NPNImplInfo& bestImpl);

    void buildNewConeConstant(unsigned int nodeindex, bool invert);

    void loadStaticRWDB();

#if 0
        void loadRWDB( std::string filename,
                       std::vector<NPNClassImplementation>& npn2impl,
                       DynamicRewritingTable& rwdb );

        void saveRWDB( std::string filename,
                       const std::vector<NPNClassImplementation>& npn2impl,
                       const DynamicRewritingTable& rwdb ) const;
#endif

    inline bool lessEqual4(unsigned int v)
    {
        unsigned int c = 0;

        for (; v; c++) {
            v &= v - 1;
        }

        return (c <= 4);
    }

    void         checkEquivalence(unsigned int nodeindex, unsigned int parent1, unsigned int parent2) const;
    unsigned int usedNodes() const;

    void print(std::ostream& os = std::cout) const;

    void recursiveRebuild(unsigned int nodeindex, std::vector<bool>& processed, std::vector<InternalEdgeRef>& cache);
    unsigned int getConeSize(unsigned int nodeindex, const CutStruct& cut) const;
    unsigned int getImplInputs(int implindex) const;

    void selectBestImpl(unsigned int nodeindex, const CutStruct& cut, int cutindex,
                        const NPNClassAndTransformations& npntrans, NPNImplInfo& bestImpl);

    void selectBestImplSmallerInput(unsigned int nodeindex, const CutStruct& cut, int cutindex,
                                    const NPNClassAndTransformations& npntrans, NPNImplInfo& bestImpl);

#if 0
        bool addNewImplementation(
            unsigned int nodeindex, const CutStruct& cut,
            const NPNClassAndTransformations& npntrans );
#endif

    void markNode(unsigned int nodeindex);
    bool isNodeMarked(unsigned int nodeindex) const;
    void clearMarks();

    CutStruct* getNewCut();
    void       releaseCuts(CutStruct* c);

    lrabs::SimpleStack<CutStruct*> _freeCuts;

    Manager* _manager;

    typedef std::unordered_map<Node*, unsigned int, lrabs::hashPointer<Node>> IDMap;
    IDMap                                                                     _id;

    NodeStruct*   _nodes;
    unsigned      _nodesCapacity, _nodesUsed;
    unsigned int* _rewritten;
    unsigned int* _modified;
    unsigned int* _sim;
    unsigned int* _hasSim;
    CutStruct**   _cuts;

    std::map<std::pair<unsigned int, unsigned int>, unsigned int> _andMap;

    std::vector<InternalEdgeRef> _backMap;

    std::vector<int>        _impl2node;
    lrabs::SimpleStack<int> _impl2nodeUsed;

    std::vector<int> _markedNodes;
    int              _mark;

    mutable lrabs::SimpleStack<unsigned int> _pending;

    DynamicRewritingTable                   _rwTable;
    std::vector<NPNClassImplementation>     _npnClassImpl;
    static const NPNClassAndTransformations _staticF2N[65536];
    static const int                        _initialImplNodesCount;
    static const int                        _initialImplNodes[];
    static const int                        _initialNPNImplementations[];

    int _deadNodes;
#if 0
        mutable bool _modifiedRWTable;
#endif

    int _totalApplications;
    int _totalNodesRewritten, _totalCutsComputed, _totalNodesChecked, _totalImplementationsChecked, _totalDelta,
        _totalEstimatedDelta;
    double _totalRewritingTime;

    int    _totalApplicationsLocal;
    int    _totalNodesRewrittenLocal, _totalDeltaLocal, _totalEstimatedDeltaLocal;
    double _totalRewritingTimeLocal;
};
}  // namespace aigpp

#include "RewritingManager.icc"

#endif /* AIGPP_REWRITINGMANAGER_HH */
