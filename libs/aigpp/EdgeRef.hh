/**************************************************************
 *
 *       AIGPP Package // EdgeRef.hh
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

#ifndef AIGPP_EDGEREF_HH
#define AIGPP_EDGEREF_HH

/* stdlib + stl */
#include <iostream>
#include <vector>

/* aigpp forwards */
namespace aigpp {
class ExtRefTable;
class InternalEdgeRef;
class Manager;
}  // namespace aigpp

#include "VariableAssignment.hh"

namespace aigpp {
/*!
 * \ingroup ext
 * The EdgeRef class represents an AIG cone. Together with the Manager class,
 * EdgeRef defines the main (user visible) interface for AIG manipulation.
 * The EdgeRef class provides methods to retrieve structural information about
 * the represented AIG, as well as methods to create new \link EdgeRef EdgeRefs
 * \endlink based on simple and complex Boolean operations.
 */
class EdgeRef
{
   public:
    /*!
     *
     * Creates an uninitialized EdgeRef object. You have to initialize the object
     * using the operator=() before any other method calls!
     * \sa isValid, operator=
     */
    EdgeRef();

    /*!
     *
     * Copy constructor. Creates a copy of \a base.
     */
    EdgeRef(const EdgeRef& base);
    EdgeRef(EdgeRef&& base) noexcept;

    /*!
     *
     * Destructor.
     */
    ~EdgeRef();

    /*!
     *
     * Assigns \a base to this EdgeRef and returns a reference to this EdgeRef.
     */
    EdgeRef& operator=(const EdgeRef& base);
    EdgeRef& operator=(EdgeRef&& base);

    /*!
     *
     * Invalidates the EdgeRef object. You have to initialize the object
     * using the operator=() before any other method calls!
     * \sa isValid(), operator=
     */
    void invalidate();

    /*!
     *
     * Returns true, if the EdgeRef is initialized and points to an AIG.
     * The default constructor (EdgeRef()) creates an uninitialized and thus
     * invalid EdgeRef. \sa invalidate
     */
    bool isValid() const;

    /*!
     *
     * Returns true if the referenced AIG is inverted; otherwise returns false.
     */
    bool isInverted() const;

    /*!
     *
     * Returns true if the referenced AIG is a constant AIG; otherwise returns
     * false. The test for constantness is a purely structural one, i.e. this
     * function checks whether the EdgeRef points to the constant node. No SAT
     * calls are involved here. \sa aigpp::structurallyTrue( const EdgeRef& ),
     * aigpp::structurallyFalse( const EdgeRef& )
     */
    bool isConstant() const;

    /*!
     *
     * Returns true if the referenced AIG is a variable node; otherwise returns
     * false;
     */
    bool isVariable() const;

    /*!
     *
     * Returns the index of the referenced AIG variable. If the AIG is not a
     * variable the method will fail. \sa isVariable
     */
    int variableIndex() const;

    /*!
     *
     * Returns the name of the referenced AIG variable. If the AIG is not a
     * variable the method will fail. \sa isVariable, Manager::variableName
     */
    const std::string& variableName() const;

    /*!
     *
     * Returns the number of nodes in the referenced AIG cone.
     * Returns 0 if the referenced AIG is constant.
     * \sa Manager::sharedSize( const std::vector<EdgeRef>& ) const
     */
    std::size_t nodeCount() const;

    /*!
     *
     * Returns a vector of indices of the variables that occur in the cone of the
     * referenced AIG.
     * \sa supportVariables() const, Manager::variable( int ) const
     */
    std::vector<int> support() const;

    /*!
     *
     * Returns a vector of references to the variables that occur in the cone of
     * the referenced AIG. \sa support() const
     */
    std::vector<EdgeRef> supportVariables() const;

    /*!
     *
     * Returns a pair of \a bool and VariableAssignment.
     * If the boolean function represented by the
     * referenced AIG is satisfiable, the pair's bool is true and the pair's
     * VariableAssignment holds a satisfying assignment to the variables in the
     * support. If \a minimize is true a minimal satisfying assignment is
     * returned. If the boolean function represented by the referenced AIG is
     * unsatisfiable, the pair's bool is false and the VariableAssignment is
     * empty.
     *
     * \sa support() const
     */
    std::pair<bool, VariableAssignment> satisfyingAssignment(bool minimize = false) const;

    /*!
     *
     * Creates an EdgeRef pointing to the logically complemented AIG of this
     * EdgeRef if \i ist true; if \a i is false it returns a copy of this EdgeRef.
     * \sa Not
     */
    EdgeRef notIf(bool i) const;

    /*!
     *
     * Creates an EdgeRef pointing to the cofactor of this \link EdgeRef EdgeRef
     * \endlink's AIG with respect to the \a variable, i.e. the function creates a
     * copy of this \link EdgeRef EdgeRef \endlink's AIG with all non-inverted
     * occurences of \a variable replaced by \a true and all inverted occurences
     * of \a variable replaced by \a false, and returns an EdgeRef pointing to the
     * resulting AIG.
     *
     * This method creates an AIG representing the function \f$f_{|x=true|}\f$ if
     * the variable pointed to by \a variable is non-inverted and
     * \f$f_{|x=false|}\f$ if the variable is inverted (\f$f\f$ is the function
     * represented by this \link EdgeRef EdgeRef \endlink's AIG and \f$x\f$
     * corresponds to the \a variable).
     *
     * The AIG pointed to by \a variable must be a variable node (isVariable).
     */
    EdgeRef cofactor(const EdgeRef& variable) const;
    /*!
     *
     * Creates an EdgeRef pointing to the cofactor of this \link EdgeRef EdgeRef
     * \endlink's AIG with respect to the \a variables, i.e. the function creates
     * a copy of this \link EdgeRef EdgeRef \endlink's AIG with all non-inverted
     * occurences of any variable in the vector \a variables replaced by \a true
     * and all inverted occurences replaced by \a false, and returns an EdgeRef
     * pointing to the resulting AIG.
     *
     * The AIGs pointed to by the \link EdgeRef EdgeRefs \endlink in the \a
     * variables vector must be a variable nodes (EdgeRef::isVariable). There must
     * not be any duplicates in the \a variables vector.
     *
     * \sa cofactor( const EdgeRef& ) const
     */
    EdgeRef cofactor(const std::vector<EdgeRef>& variables) const;

    /*!
     *
     * Creates an EdgeRef pointing to a copy of this \link EdgeRef EdgeRef
     * \endlink's AIG in which all occurences of the \a variables are replaced by
     * the corresponing \a replacements.
     *
     * The AIGs pointed to by the \link EdgeRef EdgeRefs \endlink in the \a
     * variables vector must be a variable nodes (EdgeRef::isVariable), there must
     * not be any duplicates in the \a variables vector, and the \a variables and
     * \a replacements vectors must have the same size.
     *
     * \sa Manager::vectorCompose( const std::vector<EdgeRef>&, const
     * std::vector<EdgeRef>&, const std::vector<EdgeRef>& )
     */
    EdgeRef compose(const std::vector<EdgeRef>& variables, const std::vector<EdgeRef>& replacements) const;

    /*!
     *
     * Creates an EdgeRef pointing to the existential quantification of this \link
     * EdgeRef EdgeRef \endlink's AIG with respect to the \a variable, i.e. the
     * function creates an AIG representing the function \f$\exists x . f\f$ where
     * \f$x\f$ is the variable and \f$f\f$ is the boolean function represented by
     * this \link EdgeRef EdgeRef \endlink's AIG. The formula is equivalent to the
     * expression \f$f_{|x=0} \vee f_{|x=1}\f$.
     *
     * The AIG pointed to by the EdgeRef \a variable must be a variable node
     * (EdgeRef::isVariable).
     *
     * \sa cofactor( const EdgeRef& ) const
     */
    EdgeRef existentialQ(const EdgeRef& variable) const;
    /*!
     *
     * Creates an EdgeRef pointing to the existential quantification of this \link
     * EdgeRef EdgeRef \endlink's AIG with respect to the variables in the \a
     * variables vector, i.e. the function creates an AIG representing the
     * function \f$\exists x_0 \exists x_1 \ldots \exists x_{n-1} . f\f$ where
     * \f$x_i\f$ are the variables in the \a variables vector and \f$f\f$ is the
     * boolean function represented by this \link EdgeRef EdgeRef \endlink's AIG.
     *
     * The AIGs pointed to by the \link EdgeRef EdgeRefs \endlink in the \a
     * variables vector must be a variables node (EdgeRef::isVariable).
     *
     * \sa existentialQ( const EdgeRef& ) const
     */
    EdgeRef existentialQ(const std::vector<EdgeRef>& variables) const;
    /*!
     *
     * Creates an EdgeRef pointing to the universal quantification of this \link
     * EdgeRef EdgeRef \endlink's AIG with respect to the \a variable, i.e. the
     * function creates an AIG representing the function \f$\forall x . f\f$ where
     * \f$x\f$ is the variable and \f$f\f$ is the boolean function represented by
     * this \link EdgeRef EdgeRef \endlink's AIG. The formula is equivalent to the
     * expression \f$f_{|x=0} \wedge f_{|x=1}\f$.
     *
     * The AIG pointed to by the EdgeRef \a variable must be a variable node
     * (EdgeRef::isVariable).
     *
     * \sa cofactor( const EdgeRef& ) const
     */
    EdgeRef universalQ(const EdgeRef& variable) const;
    /*!
     *
     * Creates an EdgeRef pointing to the universal quantification of this \link
     * EdgeRef EdgeRef \endlink's AIG with respect to the variables in the \a
     * variables vector, i.e. the function creates an AIG representing the
     * function \f$\forall x_0 \forall x_1 \ldots \forall x_{n-1} . f\f$ where
     * \f$x_i\f$ are the variables in the \a variables vector and \a f is the
     * boolean function represented by this \link EdgeRef EdgeRef \endlink's AIG.
     *
     * The AIGs pointed to by the \link EdgeRef EdgeRefs \endlink in the \a
     * variables vector must be a variables node (EdgeRef::isVariable).
     *
     * \sa universalQ( const EdgeRef& ) const
     */
    EdgeRef universalQ(const std::vector<EdgeRef>& variables) const;

    /*!
     * Returns the internal representation of the referenced AIG.
     */
    const InternalEdgeRef& getInternal() const;

   private:
    EdgeRef(ExtRefTable* extRefTable, const InternalEdgeRef& internal);

    static std::vector<InternalEdgeRef> toInternal(const std::vector<EdgeRef>& v);
    static std::vector<EdgeRef>         fromInternal(ExtRefTable* extRefTable, const std::vector<InternalEdgeRef>& v);

   private:
    ExtRefTable* _extRefTable;
    int          _index;

    /* friends */
    friend class Manager;

    friend std::ostream& operator<<(std::ostream& os, const EdgeRef& e);

    friend EdgeRef Not(const EdgeRef& e);
    friend EdgeRef And(const EdgeRef& e1, const EdgeRef& e2);
    friend EdgeRef Or(const EdgeRef& e1, const EdgeRef& e2);
    friend EdgeRef Xor(const EdgeRef& e1, const EdgeRef& e2);
    friend EdgeRef Mux(const EdgeRef& select, const EdgeRef& eThen, const EdgeRef& eElse);
    friend EdgeRef Equiv(const EdgeRef& e1, const EdgeRef& e2);

    friend EdgeRef multiAnd(Manager& manager, const std::vector<EdgeRef>& e);
    friend EdgeRef multiOr(Manager& manager, const std::vector<EdgeRef>& e);
};

/*!
 * \ingroup ext
 * Pretty prints the EdgeRef \a e to the output stream \a os.
 */
std::ostream& operator<<(std::ostream& os, const EdgeRef& e);

/*!
 * \ingroup ext
 * Creates the complement of the EdgeRef \a e, i.e. creates an EdgeRef that
 * points to a complemented version of \a e's AIG. \sa notIf
 */
EdgeRef Not(const EdgeRef& e);
/*!
 * \ingroup ext
 * Creates the complement of the EdgeRef \a e, i.e. creates an EdgeRef that
 * points to a complemented version of \a e's AIG. \sa notIf
 */
EdgeRef operator!(const EdgeRef& e);

/*!
 * \ingroup ext
 * Creates the conjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the conjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef And(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Creates the conjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the conjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef operator&(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Creates the conjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the conjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef operator*(const EdgeRef& e1, const EdgeRef& e2);

/*!
 * \ingroup ext
 * Creates the disjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the disjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef Or(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Creates the disjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the disjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef operator|(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Creates the disjunction of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the disjunction of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef operator+(const EdgeRef& e1, const EdgeRef& e2);

/*!
 * \ingroup ext
 * Creates the exclusive or of the two \link EdgeRef EdgeRefs \endlink \a e1 and
 * \a e2, i.e. creates an EdgeRef that points to the exclusive or of \a e1's and
 * \a e2's AIGs.
 */
EdgeRef Xor(const EdgeRef& e1, const EdgeRef& e2);

/*!
 * \ingroup ext
 * Creates the negated exclusive or of the two \link EdgeRef EdgeRefs \endlink
 * \a e1 and \a e2, i.e. creates an EdgeRef that points to the negated exclusive
 * or of \a e1's and \a e2's AIGs.
 */
EdgeRef Equiv(const EdgeRef& e1, const EdgeRef& e2);

/*!
 * \ingroup ext
 * Creates a multiplexer of the two \link EdgeRef EdgeRefs \endlink \a eThen and
 * \a eElse using the select signal \a select, i.e. creates an EdgeRef that
 * points to an AIG representation of \f$( select \rightarrow eThen ) \wedge (
 * \overline{select} \rightarrow eElse)\f$.
 *
 * Example:
 * \code
 * Manager m;
 *
 * EdgeRef a = m.addVariable( "a" );
 * EdgeRef b = m.addVariable( "b" );
 * EdgeRef c = m.addVariable( "c" );
 *
 * EdgeRef mux1 = Mux( a, b, c );
 * EdgeRef mux2 = (!a | b ) & ( a | c );
 *
 * functionallyEquivalent( mux1, mux2 ); // -> true
 * \endcode
 */
EdgeRef Mux(const EdgeRef& select, const EdgeRef& eThen, const EdgeRef& eElse);

/*!
 * \ingroup ext
 * Creates the conjunction of the \link EdgeRef EdgeRefs \endlink is the vector
 * \a e. \sa And
 */
EdgeRef multiAnd(Manager& mananger, const std::vector<EdgeRef>& e);
/*!
 * \ingroup ext
 * Creates the disjunction of the \link EdgeRef EdgeRefs \endlink is the vector
 * \a e. \sa Or
 */
EdgeRef multiOr(Manager& manager, const std::vector<EdgeRef>& e);

/*!
 * \ingroup ext
 * Returns true if the AIG's pointed to by \a e1 and \a e2 are structurally
 * identical. \sa functionallyEquivalent
 */
bool structurallyEquivalent(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Returns true if the AIG's pointed to by \a e1 and \a e2 are structurally
 * antivalent, i.e. one AIG is the structural complement of the other. \sa
 * functionallyAntivalent
 */
bool structurallyAntivalent(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Returns true if \a e points to the constant \a true AIG.
 * \sa functionallyTrue, EdgeRef::isConstant
 */
bool structurallyTrue(const EdgeRef& e);
/*!
 * \ingroup ext
 * Returns true if \a e points to the constant \a false AIG.
 * \sa functionallyFalse, EdgeRef::isConstant
 */
bool structurallyFalse(const EdgeRef& e);

/*!
 * \ingroup ext
 * Returns true if the AIG's pointed to by \a e1 and \a e2 are functionally
 * equivalent, i.e. the represented boolean functions are equal. This function
 * several methods to prove/disprove the equivalence:
 * - purely structural checks
 * - simulation vector based disproval
 * - and finally a SAT based equivalence check (if both AIGs are not
 * functionally reduced)
 *
 * \sa structurallyEquivalent
 */
bool functionallyEquivalent(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Returns true if the AIG's pointed to by \a e1 and \a e2 are functionally
 * antivalent, i.e. the represented boolean functions are complementary. This
 * function several methods to prove/disprove the antivalence:
 * - purely structural checks
 * - simulation vector based disproval
 * - and finally a SAT based antivalence check (if both AIGs are not
 * functionally reduced)
 *
 * \sa structurallyAntivalent
 */
bool functionallyAntivalent(const EdgeRef& e1, const EdgeRef& e2);
/*!
 * \ingroup ext
 * Returns true if the AIG pointed to by \a e is a tautology (always true), i.e.
 * the represented boolean function is a tautology.
 * This function several methods to prove/disprove the tautology:
 * - purely structural checks
 * - simulation vector based disproval
 * - and finally a SAT based tautology check (if the AIG is not functionally
 * reduced)
 *
 * \sa structurallyTrue
 */
bool functionallyTrue(const EdgeRef& e);
/*!
 * \ingroup ext
 * Returns true if the AIG pointed to by \a e is a contradiction (always false),
 * i.e. the represented boolean function is a contradiction. This function
 * several methods to prove/disprove the contradiction:
 * - purely structural checks
 * - simulation vector based disproval
 * - and finally a SAT based contradiction check (if the AIG is not functionally
 * reduced)
 *
 * \sa structurallyTrue
 */
bool functionallyFalse(const EdgeRef& e);
}  // namespace aigpp

#include "EdgeRef_operations.icc"

#endif /* AIGPP_EDGEREF_HH */
