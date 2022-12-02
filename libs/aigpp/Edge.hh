/**************************************************************
 *
 *       AIGPP Package // Edge.hh
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
 *
 ***************************************************************/

#ifndef AIGPP_EDGE_HH
#define AIGPP_EDGE_HH

#include <cstdint>
#include <iostream>

namespace aigpp {
class Node;
}

namespace aigpp {
/*!
 *  \brief The Edge class represents an AIG edge possibly
 *         containing a complementation marker.
 */
struct Edge
{
    /*!
     *   Constructor. Creates a Edge representing the constant 0
     *   function.
     */
    Edge();

    /*!
     *   Constructor. Creates a Edge based on the Node* \e n and the
     *   complementation flag \e i.
     *
     *   \param n source Node*
     *   \param i complementation flag
     */
    explicit Edge(Node* n, bool i = false);

    /*!
     *   Returns a pointer to the referenced Node.
     *
     *   \return pointer to the referenced Node
     */
    inline Node* node() const;

    /*!
     *   Returns the value of complementation marker.
     *
     *   \return complementation marker
     */
    inline bool isInverted() const;

    /*!
     *   Returns true, if the Edge represents a constant function.
     *
     *   \return true, if constant
     */
    inline bool isConstant() const;

    /*!
     *   Sets the referenced Node pointer.
     *
     *   \param n source Node pointer
     */
    void setNode(Node* n);

    /*!
     *   Sets the complementation marker.
     *
     *   \param i source complementation marker
     */
    void setInverted(bool i);

    /*!
     *   Toggles the complementation marker.
     */
    void toggleInverted();

    /*!
     *   Toggles the complementation marker if \e invert is true.
     */
    void toggleInverted(bool invert);

    /*!
     *   Swaps *this Edge and \e e.
     */
    void swap(Edge& e);

    /*!
     *   Clears the Edge, e.g. set it to constant 0.
     */
    void clear();

    /*!
     *   Returns the Edge representing the complement of \e *this, e.g.
     *   <em>!*this</em>.
     *
     *   \return complement of \e *this
     */
    Edge operator!() const;

    /*!
     *   Returns the Edge representing the complement of \e *this if
     *   invert is true, otherwise return \e *this.
     *
     *   \param invert complementation flag
     *
     *   \return complement of \e *this if \e invert is true, \e *this otherwise
     */
    Edge notIf(bool invert) const;

    /*!
     * Returns true, if \e *this and \e n are structurally equal.
     *
     * \param n
     *
     * \return true if equal
     */
    bool structurallyEquivalent(const Edge& n) const { return (_node == n._node); };
    bool structurallyAntivalent(const Edge& n) const { return (_node == (!n)._node); };
    bool structurallyTrue() const { return (_node == reinterpret_cast<Node*>(0x1)); };
    bool structurallyFalse() const { return (_node == reinterpret_cast<Node*>(0x0)); };

    bool functionallyEquivalent(const Edge& n) const;
    bool functionallyAntivalent(const Edge& n) const { return functionallyEquivalent(!n); };
    bool functionallyTrue() const { return functionallyEquivalent(Edge(nullptr, true)); };
    bool functionallyFalse() const { return functionallyEquivalent(Edge(nullptr, false)); };

#if 0
        /*!
         * Returns true, if \e *this and \e n are structurally unequal.
         *
         * \param n
         *
         * \return true if unequal
         */
        bool structurallyUnequal( const Edge& n ) const
        {
            return( _node != n._node );
        };

        /*!
         *   Returns true if \e *this and \e n are functionally equivalent.
         *
         *   \param n second operand
         *
         *   \return true, if equivalent
         */
        bool operator==( const Edge& n ) const;

        /*!
         *   Returns true if \e *this and \e n are NOT functionally equivalent.
         *
         *   \param n second operand
         *
         *   \return true, if NOT equivalent
         */
        bool operator!=( const Edge& n ) const
        {
            return !operator==( n );
        };
#endif

    /*!
     *  Returns the size of the cone of the referenced Node.
     *
     *   \return cone size
     */
    std::size_t nodeCount() const;

    /*!
     *  Returns +1, if \e n is directly implied by \e *this, -1, if \e !n is
     * directly implied by \e *this, 0 otherwise, e.g. \n +1, if \e n occurs on a
     * path rooted in \e *this containing no complementation markers. \n -1, if \e
     * n occurs in complemented form on a path rooted in \e *this containing no
     * complementation markers. \n 0, if \e n does not occur on any path rooted in
     * \e *this containing no complementation markers.
     *
     *   \return
     */
    int hasSuperNode(const Edge& n) const;

    unsigned long hash() const;

   protected:
    Node* _node;

    friend class ExtRefTable;
    friend class Manager;
    friend class FRAIGManager;
};

/*!
 *   Pretty printer for Edge.
 *
 *   \param os output stream
 *   \param n Edge to be printed
 *
 *   \return modified stream
 */
std::ostream& operator<<(std::ostream& os, const Edge& n);

inline bool
structurallyEquivalent(const Edge& e1, const Edge& e2)
{
    return e1.structurallyEquivalent(e2);
}

inline bool
structurallyAntivalent(const Edge& e1, const Edge& e2)
{
    return e1.structurallyAntivalent(e2);
}

inline bool
structurallyTrue(const Edge& e)
{
    return e.structurallyTrue();
}

inline bool
structurallyFalse(const Edge& e)
{
    return e.structurallyFalse();
}

inline bool
functionallyEquivalent(const Edge& e1, const Edge& e2)
{
    return e1.functionallyEquivalent(e2);
}

inline bool
functionallyAntivalent(const Edge& e1, const Edge& e2)
{
    return e1.functionallyAntivalent(e2);
}

inline bool
functionallyTrue(const Edge& e)
{
    return e.functionallyTrue();
}

inline bool
functionallyFalse(const Edge& e)
{
    return e.functionallyFalse();
}

struct structurallyEqual
{
    bool operator()(const Edge& e1, const Edge& e2) const { return e1.structurallyEquivalent(e2); };
};

struct structurallyLess
{
    bool operator()(const Edge& e1, const Edge& e2) const
    {
        if (e1.isInverted() && !(e2.isInverted())) {
            return true;
        } else if (!(e1.isInverted()) && e2.isInverted()) {
            return false;
        } else {
            return (e1.node() < e2.node());
        }
    };
};

}  // namespace aigpp

#include "Edge.icc"

#endif /* AIGPP_EDGE_HH */
