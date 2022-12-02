/**************************************************************
 *
 *       AIGPP Package // Node.hh
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

#ifndef AIGPP_NODE_HH
#define AIGPP_NODE_HH

#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>

#include <core/Solver.h>

namespace aigpp {
class Manager;
}

#define USE_TEMPSIM

#include "Edge.hh"
#include "InternalEdgeRef.hh"
#include "SimVector.hh"
#include "VarSupport.hh"

namespace aigpp {

/*!
 * \class Node
 * \brief The Node class represents an AIG node.
 *
 * The Node class represents either an And/Nand Node or a variable Node (see
 * isVar()). Nodes are created internally by a Manager object.
 */
class Node
{
   public:
    ~Node();

    /*!
     * Returns the Node's responsible Manager.
     *
     * \return manager
     */
    Manager* manager() const;

    /*!
     * Returns the Node's reference count.
     *
     * \return reference count
     */
    int refCount() const;

    int index() const { return _index; };

    /*!
     * Returns true if this Node represents a Boolean variable.
     *
     * \return true if variable
     */
    bool isVar() const;

    /*!
     * Returns the index of the associated variable. This method assumes that
     * the Node actually represents a Boolean variable (see isVar()).
     *
     * \return variable index
     */
    int varIndex() const;

    /*!
     * Returns true if the node represents a NAND node.
     *
     * \return true if NAND
     */
    bool isNAND() const;

    /*!
     * Returns true if the node is functionally reduced.
     */
    bool isReduced() const;

    /*!
     * Returns an edge to the first parent of the node. This method assumes that
     * the Node actually is an And/Nand node (see isVar()).
     *
     * \return first parent
     */
    const Edge& parent1() const;

    /*!
     * Returns an edge to the second parent of the node. This method assumes that
     * the Node actually is an And/Nand node (see isVar()).
     *
     * \return second parent
     */
    const Edge& parent2() const;

    /*!
     * Returns true if the given Node \e n occurs in the cone of this Node.
     *
     * \param n
     *
     * \return true if \e n is in the cone
     */
    bool hasParent(Node* n) const;

    /*!
     * Returns the simulation vector of this Node.
     *
     * \return simulation vector
     */
    const SimVector& sim() const;
    void             updateSim();
    void             updateSim(unsigned int bin);

#ifdef USE_TEMPSIM
    unsigned long tempSim() const;
    void          updateTempSim();
    void          setTempSim(unsigned long s);
#endif

    /*!
     * Returns the support set of this Node.
     *
     * \return support set
     */
    VarSupport support() const;

    /*!
     * Returns the support set of this Node, i.e. the vector of variable Nodes
     * in the support of this Node
     *
     * \return support vector
     */
    std::vector<Node*> supportNodes() const;

    /*!
     * Checks if the support of *this contains the variable with index \e v.
     * \param v
     * \return true if v is in the support
     */
    bool varInSupport(int v) const;

    /*!
     * Returns an approximation of the number of nodes remaining after a
     * quantification of variable \e var.
     *
     * \param var variable index
     *
     * \return number of nodes remaining
     */
    std::size_t remainingNodesAfterQuantification(int var) const;

    /*!
     * Returns the size of this Node, e.g. the size of the Node's cone.
     *
     * \return node size
     */
    std::size_t nodeCount() const;

    /*!
     * Returns the set of nodes in the cone of this Node.
     *
     * \return cone
     */
    std::set<Node*>    cone();
    std::vector<Node*> coneVector();

    /*!
     * Returns a map containing the index and the maximum depth of each variable
     * in this node. Variables not occuring in the cone have no entry in the map.
     *
     * \return depth map
     */
    std::map<int, std::size_t> maxVarDepth() const;

    /*!
     * Returns a map containing the index and the minimum depth of each variable
     * in this node. Variables not occuring in the cone have no entry in the map.
     *
     * \return depth map
     */
    std::map<int, std::size_t> minVarDepth() const;
    std::vector<int>           minVarDepth2() const;

    int maxDepth(int var) const;
    int minDepth(int var) const;

    bool isMuxType() const;
    bool isXorType() const;
    bool isMuxOrXorType() const;
    bool isMultiAnd(std::vector<Edge>& inputs) const;

    std::vector<Node*> findRedundantVars(bool remove = false);

    inline Node* next() const;

    inline bool isCacheValid() const;
    inline bool isZCacheValid() const;

    enum Flags
    {
        FLAG_ISNAND      = 4,
        FLAG_CACHEVALID  = 8,
        FLAG_ZCACHEVALID = 16,
        FLAG_ISREDUCED   = 32,
        FLAG_CACHECONST0 = 64,
        FLAG_CACHECONST1 = 128,
        FLAG_REWRITTEN   = 256,
        FLAG_UNMODIFIED  = 512
    };

    template <unsigned int FLAGS>
    inline bool flag() const;
    template <unsigned int FLAGS>
    inline bool allFlags() const;
    template <unsigned int FLAGS>
    inline void setFlag() const;
    template <unsigned int FLAGS>
    inline void unsetFlag() const;
    void        clearFlags() const;

    template <unsigned int FLAGS>
    inline void setFlagRecursive() const;

   protected:
    /*
     * create an aig-leaf representing the variable with index var
     */
    Node(Manager* manager, int var);

    /*
     * create an inner aig-node with parents parent1, parent2
     */
    Node(Manager* manager, const Edge& parent1, const Edge& parent2);

    Node(const Node& n);

    Node& operator=(const Node& n) = delete;

    void set(Manager* m, const Edge& p1, const Edge& p2);

    void ref();
    void deref();

    SimVector& sim();

    Minisat::Var satVar() const;
    bool         isSatVarValid() const;
    void         setSatVar(Minisat::Var v);
    void         invalidateSatVar();

    struct PropagationData
    {
        int          type;
        aigpp::Node* p1;
        aigpp::Node* p2;
    };

    mutable PropagationData _prop0, _prop1;

   protected:
    Manager* _manager;

   private:
    mutable unsigned int _flags;

    int _refCount;
    int _index;

    Edge _parent1, _parent2;

    SimVector* _sim;

#ifdef USE_TEMPSIM
    unsigned long _tempSim;
#endif

    Node* _nEqualSimHash;
    Node* _pEqualSimHash;
    Node* _nEqualSim;
    Node* _pEqualSim;
    Node* _invertedSimClass;

    Node* _prevNode;
    Node* _nextNode;

    Node* _prevUniqueTableEntry;
    Node* _nextUniqueTableEntry;

    Minisat::Var _satVar;

    friend class Manager;
    friend struct Edge;
    friend class InternalEdgeRef;
    friend class SimulationTable;
    friend class UniqueTable;
    friend class RewritingManager;
    friend class FRAIGManager;

    friend std::ostream& operator<<(std::ostream& os, const Node& n);
};

std::ostream& operator<<(std::ostream& os, const Node& n);

}  // namespace aigpp

#include "Node.icc"

#endif /* AIGPP_NODE_HH */
