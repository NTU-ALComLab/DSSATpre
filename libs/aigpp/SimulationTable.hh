/**************************************************************
 *
 *       AIGPP Package // SimulationTable.hh
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

#ifndef SIMULATIONTABLE_HH
#define SIMULATIONTABLE_HH

#include "SimVector.hh"

namespace aigpp {
class Node;

/*!
 *  \brief The SimulationTable class is a hash table that maps SimVectors
 *         to lists containing pointers to all AIG nodes with the
 *         corresponding SimVector, e.g. the equivalence class wrt. to the
 *         SimVector.
 *
 *
 *   Collisions are resolved by chaining.
 *   The hash table grows each time it is filled at least 75% but it
 *   will never shrink.
 */
class SimulationTable
{
   public:
    SimulationTable();
    SimulationTable(const SimulationTable&) = delete;
    SimulationTable(SimulationTable&&)      = delete;
    SimulationTable& operator=(const SimulationTable&) = delete;
    SimulationTable& operator=(SimulationTable&&) = delete;
    ~SimulationTable();

    Node* getClass(Node* n) const;
    Node* getInvertedClass(Node* nonInvertedClass) const;

    /*!
     *   Returns the linked list containing all Nodes whose SimVector is equal
     *   to \e sim, or 0 if no such Nodes exist.
     *
     *   \param sim SimVector to find
     *
     *   \return list of Nodes
     */
    Node* lookup(const SimVector& sim) const;

    /*!
     *   Returns the linked list containing all Nodes whose SimVector is equal
     *   to the negation of \e sim, or 0 if no such Nodes exist.
     *
     *   \param sim SimVector to find
     *
     *   \return list of Nodes
     */
    Node* lookupInverted(const SimVector& sim, Node* nonInvertedClass = nullptr) const;

    /*!
     *  Insert Node \n into the appropriate list.
     *
     *   \param n Node to insert
     */
    void insert(Node* n);

    /*!
     *  Clears the SimulationTable.
     */
    void clear();
    void clearRememberClasses();
    void rebuild(Node* nodes, bool useRememberedClasses = false);

    /*!
     *  Grows the hash table to the size given by \e Primes[newPrimesIndex].
     *  If the new size is less or equal to the current size, nothing is done.
     *
     *   \param newPrimesIndex index for the desired size
     */
    void resize(int newPrimesIndex);

    /*!
     * Remove dead nodes from the table.
     */
    void garbageCollect();

    /*!
     *   Returns the size of the linked list of simulation-equivalent Nodes
     *   starting in \e simclass
     *
     *   \param simclass head of the list
     *
     *   \return size of the list
     */
    static std::size_t simClassSize(Node* simclass);

    /*!
     * returns the maximum number of nodes in a simclass of the current simulation
     * table.
     */
    std::size_t maxSimClassSize() const;

    /*!
     *   Checks the integrity of the SimulationTable and aborts if integrity
     *   is violated.
     */
    void checkIntegrity(Node* nodeslist) const;

    std::size_t entries() const { return _entries; };

   public:
    Node**       _table;
    std::size_t  _capacity;
    std::size_t  _entries;
    unsigned int _primesIndex;
};

}  // namespace aigpp

#endif /* SIMULATIONTABLE_HH */
