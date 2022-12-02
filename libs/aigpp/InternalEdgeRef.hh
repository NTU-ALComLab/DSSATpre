/**************************************************************
 *
 *       AIGPP Package // InternalEdgeRef.hh
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 373 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#ifndef AIGPP_INTERNALNODEREF_HH
#define AIGPP_INTERNALNODEREF_HH

// stdlib/stl includes
#include <cassert>
#include <iostream>
#include <map>
#include <stack>
#include <utility>
#include <vector>

#include <cudd.h>

// aigpp forward declarations
namespace aigpp {
class Node;
class VarSupport;
class EdgeRef;
}  // namespace aigpp

// aigpp includes
#include "Edge.hh"
#include "VariableAssignment.hh"

namespace aigpp {
/*!
 *  \brief The InternalEdgeRef class represents an AIG edge, is responsible for
 *         the garbage collection of the referenced Node, and provides
 *         methods for the usual boolean function operations.
 */
class InternalEdgeRef : public Edge
{
   public:
    /*!
     *   Constructor. Creates a InternalEdgeRef representing the 0 constant.
     */
    InternalEdgeRef();

    /*!
     *   Copy constructor. Creates a copy of n.
     *
     *   \param n source InternalEdgeRef
     */
    InternalEdgeRef(const InternalEdgeRef& n);

    /*!
     *   Constructor. Creates a InternalEdgeRef based on the Edge n.
     *
     *   \param n source Edge
     */
    explicit InternalEdgeRef(const Edge& n);

    /*!
     *   Constructor. Creates a InternalEdgeRef based on the Node* n and the
     *   complementation flag i.
     *
     *   \param n source Node*
     *   \param i complementation flag
     */
    explicit InternalEdgeRef(Node* n, bool i = false);

    /*!
     *   Destructor. Destroys the InternalEdgeRef object and decrements the
     *   referenced Node's reference counter.
     */
    ~InternalEdgeRef();

    /*!
     *   Clears the InternalEdgeRef, e.g. set it to constant 0.
     */
    void clear();

    /*!
     *   Assignment operator. Created a copy of n and returns a reference
     *   to the created copy.
     *
     *   \param n source InternalEdgeRef
     *
     *   \return copy
     */
    InternalEdgeRef& operator=(const InternalEdgeRef& n);

    static void recursivelySelectBothPolarities(Node* root, std::vector<Node*>& processed, std::stack<Node*>& pending);
    /*!
     *  Returns a vector of variables along with the polarity in this AIG.
     */
    std::vector<std::pair<Node*, bool>> getVariablePolarities();

    /*!
     *   Creates a InternalEdgeRef representing the mux/ITE operation
     *   <em>(*this * nodeThen) + (!*this * nodeElse)</em>
     *
     *   \param nodeThen then branch
     *   \param nodeElse else branch
     *
     *   \return InternalEdgeRef representing the mux operation
     */
    InternalEdgeRef Mux(const InternalEdgeRef& nodeThen, const InternalEdgeRef& nodeElse) const;

    /*!
     *   Creates a InternalEdgeRef representing the xnor operation, e.g.
     *   <em>!(*this xor n )</em>
     *
     *   \param n the second operator of the operation
     *
     *   \return InternalEdgeRef representing the xnor operation
     */
    InternalEdgeRef Xnor(const InternalEdgeRef& n) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of an existential
     *   quantification operation, e.g. <em>exists var. *this</em>
     *   Quantification is expressed using the positive and negative
     *   cofactors of *this wrt. var.
     *   \n
     *   The parameter \e var must a InternalEdgeRef pointing to a variable node!
     *
     *   \param var the quantification variable
     *   \param internalCall indicate that this method call came from within
     *          the library (-> skip unneccessary checks...)
     *
     *   \return InternalEdgeRef representing existential quantification wrt. \e
     * var
     */
    InternalEdgeRef existentialQ(const InternalEdgeRef& var, bool internalCall = false) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of a universal
     *   quantification operation, e.g. <em>forall var. *this</em>
     *   Quantification is expressed using the positive and negative
     *   cofactors of *this wrt. var.
     *   \n
     *   The parameter \e var must a InternalEdgeRef pointing to a variable node!
     *
     *   \param var the quantification variable
     *
     *   \return InternalEdgeRef representing universal quantification wrt. \e var
     */
    InternalEdgeRef universalQ(const InternalEdgeRef& var) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of multiple existential
     *   quantification operations, e.g.
     *   <em>exists var[0]...var[n]. *this</em>
     *   \n
     *   The vector \e vars must contain InternalEdgeRefs pointing to variable
     * nodes only!
     *
     *   \param vars the quantification variables
     *
     *   \return InternalEdgeRef representing existential quantification wrt. \e
     * vars
     */
    InternalEdgeRef existentialQ(const std::vector<InternalEdgeRef>& vars) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of multiple universal
     *   quantification operations, e.g.
     *   <em>exists var[0]...var[n]. *this</em>
     *   \n
     *   The vector \e vars must contain InternalEdgeRefs pointing to variable
     * nodes only!
     *
     *   \param vars the quantification variables
     *
     *   \return InternalEdgeRef representing universal quantification wrt. \e
     * vars
     */
    InternalEdgeRef universalQ(const std::vector<InternalEdgeRef>& vars) const;

    /*!
     *   Creates a InternalEdgeRef representing the original AIG with all
     * variables \e vars1[i] replaced by \e vars2[i] and vice versa. \n The
     * vectors \e vars1 and \e vars2 must contain InternalEdgeRefs pointing to
     *   variable nodes only,
     *   no variable must occur twice in \e vars1 and \e vars2,
     *   and the sets of variables contained in \e vars1 and \e vars2 must be
     *   disjunct.
     *
     *   \param vars1 the first variable vector
     *   \param vars2 the second variable vector
     *
     *   \return InternalEdgeRef representing the result of the swap operation
     */
    InternalEdgeRef swapVariables(const std::vector<InternalEdgeRef>& vars1,
                                  const std::vector<InternalEdgeRef>& vars2) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of a cofactor operation
     *   wrt. var, e.g. <em>*this[var<-1]</em>.
     *   If \e var is complemented the negative cofactor is created, otherwise
     *   the positive cofactor is built.
     *   \n
     *   The parameter \e var must a InternalEdgeRef pointing to a variable node!
     *
     *   \param var the cofactor variable
     *   \param checkSupport check if var occurs in the support before cofactoring
     *
     *   \return InternalEdgeRef representing the cofactor of \e *this wrt. \e var
     */
    InternalEdgeRef cofactor(const InternalEdgeRef& var, bool checkSupport = true) const;

    /*!
     *   Creates a InternalEdgeRef representing the result of multiple cofactor
     *   operations, e.g. <em>*this[vars[0]<-1, ..., vars[n]<-1]</em>.
     *   \n
     *   The vector \e vars must contain InternalEdgeRefs pointing to variable
     * nodes only!
     *
     *   \param vars the cofactor variables
     *
     *   \return InternalEdgeRef representing the cofactor of \e *this wrt. \e
     * vars
     */
    InternalEdgeRef cofactor(const std::vector<InternalEdgeRef>& vars) const;

    static std::vector<InternalEdgeRef> cofactors(const std::vector<InternalEdgeRef>& roots,
                                                  const InternalEdgeRef&              var);

    /*!
     *   Creates a InternalEdgeRef representing the result of a substitution
     *   operation, e.g.
     *   <em>*this[vars[0]<-replacements[0], ..., vars[n]<-replacements[n]]"</em>.
     *   \n
     *   The set \e vars must contain non-complemented InternalEdgeRefs pointing
     * to variable nodes only, no variable must occur twice in \e vars, and \e
     * vars and \e replacements must be of the same size.
     *
     *   \param vars the variables to be replaced
     *   \param replacements the replacement nodes
     *
     *   \return InternalEdgeRef representing the substitution operation
     */
    InternalEdgeRef compose(const std::vector<InternalEdgeRef>& vars,
                            const std::vector<InternalEdgeRef>& replacements) const;

    InternalEdgeRef composeZ(const std::vector<InternalEdgeRef>& vars, const std::vector<InternalEdgeRef>& replacements,
                             const InternalEdgeRef& Z) const;

    /*!
     *   Returns the InternalEdgeRef representing the complement of \e *this, e.g.
     *   <em>!*this</em>.
     *
     *   \return complement of \e *this
     */
    inline InternalEdgeRef operator!() const;

    /*!
     *   Returns the InternalEdgeRef representing the complement of \e *this if
     *   invert is true, otherwise return \e *this.
     *
     *   \param invert complementation flag
     *
     *   \return complement of \e *this if \e invert is true, \e *this otherwise
     */
    InternalEdgeRef notIf(bool invert) const;

    /*!
     *   Returns the InternalEdgeRef representing the conjunction of \e *this and
     * \e n2, e.g. <em>*this * n2</em>.
     *
     *   \param n2 second operand
     *
     *   \return conjunction
     */
    InternalEdgeRef operator&(const InternalEdgeRef& n2) const;

    /*!
     *   Returns the InternalEdgeRef representing the disjunction of \e *this and
     * \e n2, .e.g. <em>*this + n2</em>.
     *
     *   \param n2 second operand
     *
     *   \return disjunction
     */
    inline InternalEdgeRef operator|(const InternalEdgeRef& n2) const;

    /*!
     *   Returns the InternalEdgeRef representing the exclusive or of \e *this and
     * \e n2, .e.g. <em>*this ^ n2</em>.
     *
     *   \param n2 second operand
     *
     *   \return exclusive or
     */
    inline InternalEdgeRef operator^(const InternalEdgeRef& n2) const;

    /*!
     *   Returns the InternalEdgeRef representing the conjunction of \e *this and
     * \e n2, e.g. <em>*this * n2</em>.
     *
     *   \param n2 second operand
     *
     *   \return conjunction
     */
    inline InternalEdgeRef operator*(const InternalEdgeRef& n2) const;

    /*!
     *   Returns the InternalEdgeRef representing the exclusive or of \e *this and
     * \e n2, .e.g. <em>*this ^ n2</em>.
     *
     *   \param n2 second operand
     *
     *   \return exclusive or
     */
    inline InternalEdgeRef operator+(const InternalEdgeRef& n2) const;

    /*!
     *   Returns the InternalEdgeRef representing the conjunction of \e *this and
     * \e n, e.g. <em>*this * n</em>, and assigns the result to \e *this
     *
     *   \param n second operand
     *
     *   \return conjunction
     */
    InternalEdgeRef operator&=(const InternalEdgeRef& n);

    /*!
     *   Returns the InternalEdgeRef representing the disjunction of \e *this and
     * \e n, e.g. <em>*this + n</em>, and assigns the result to \e *this
     *
     *   \param n second operand
     *
     *   \return disjunction
     */
    InternalEdgeRef operator|=(const InternalEdgeRef& n);

    /*!
     *   Returns the InternalEdgeRef representing the exclusive or of \e *this and
     * \e n, e.g. <em>*this ^ n</em>, and assigns the result to \e *this
     *
     *   \param n second operand
     *
     *   \return exclusive or
     */
    InternalEdgeRef operator^=(const InternalEdgeRef& n);

    /*!
     *   Returns the InternalEdgeRef representing the conjunction of \e *this and
     * \e n, e.g. <em>*this * n</em>, and assigns the result to \e *this
     *
     *   \param n second operand
     *
     *   \return conjunction
     */
    InternalEdgeRef operator*=(const InternalEdgeRef& n);

    /*!
     *   Returns the InternalEdgeRef representing the disjunction of \e *this and
     * \e n, e.g. <em>*this + n</em>, and assigns the result to \e *this
     *
     *   \param n second operand
     *
     *   \return disjunction
     */
    InternalEdgeRef operator+=(const InternalEdgeRef& n);

    /*!
     *   Returns the InternalEdgeRef representing the conjunction of all \e nodes,
     *   e.g. <em>nodes[0] * nodes[1] * ... * nodes[n]</em>.
     *
     *   \param nodes operands
     *
     *   \return conjunction
     */
    static InternalEdgeRef multiAnd(std::vector<InternalEdgeRef> nodes);

    /*!
     *   Returns the InternalEdgeRef representing the disjunction of all \e nodes,
     *   e.g. <em>nodes[0] + nodes[1] + ... + nodes[n]</em>.
     *
     *   \param nodes operands
     *
     *   \return conjunction
     */
    static InternalEdgeRef multiOr(std::vector<InternalEdgeRef> nodes);

    /*!
     * Create a BDD for this node using the given BDD manager and the given
     * initial variable mapping. The BDD is created recursively following
     * the AIG structure. If an AIG variable is encountered,
     * the BDD corresponding to the variable's index from the varMap is used.
     * If the map doesn't contain a BDD for this variable index, a new
     * BDD is created and inserted into the map.
     * If the BDD creation failes (e.g. because of a memory out) the null
     * pointer is returned.
     */
    DdNode* exportBDD(DdManager* bddManager, std::map<int, DdNode*>& varMap) const;

    /*!
     * Decompose a "big" conjunction to the single factors, e.g. create a list of
     * edges such that the conjunction of these edges ist equivalent to *this.
     */
    std::vector<InternalEdgeRef> decomposeConjunction() const;
    std::vector<InternalEdgeRef> decomposeDisjunction() const;

   protected:
    void setNode(Node* n);

   private:
    friend class Manager;
};

/*!
 * InternalEdgeRef representing the constant 0 function.
 */
static const InternalEdgeRef const0 = InternalEdgeRef(nullptr, false);

/*!
 * InternalEdgeRef representing the constant 1 function.
 */
static const InternalEdgeRef const1 = InternalEdgeRef(nullptr, true);
}  // namespace aigpp

#include "InternalEdgeRef.icc"

#endif /* AIGPP_INTERNALNODEREF_HH */
