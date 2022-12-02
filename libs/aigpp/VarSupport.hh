/**************************************************************
 *
 *       AIGPP Package // VarSupport.hh
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

#ifndef AIGPP_VARSUPPORT_HH
#define AIGPP_VARSUPPORT_HH

#include <cassert>
#include <climits>
#include <iostream>
#include <vector>

namespace aigpp {

/*!
 * \class VarSupport
 * \brief The VarSupport class represents the support of a Node, e.g.
 *        the set of variables found in the cone of a Node object.
 */
class VarSupport
{
   public:
    typedef unsigned long Mask;

    static constexpr unsigned int MaskSize = sizeof(Mask) * CHAR_BIT;

    /*!
     * Constructor. Creates an empty support set.
     */
    VarSupport();

    /*!
     * Constructor. Creates a support set contatind the variable index \e var.
     *
     * \param var variable index
     */
    explicit VarSupport(std::size_t var);

    /*!
     * Constructor. Creates a support set equal to the union of the two given
     * support sets \e vs1 and \e vs2.
     *
     * \param vs1 first parent VarSupport
     * \param vs2 second parent VarSupport
     */
    VarSupport(const VarSupport& vs1, const VarSupport& vs2);

    /*!
     * Resets \e *this to the support set equal to the union of the two given
     * support sets \e vs1 and \e vs2.
     *
     * \param vs1 first parent VarSupport
     * \param vs2 second parent VarSupport
     */
    void set(const VarSupport& vs1, const VarSupport& vs2);

    /*!
     * Returns the size of the vector which stores the represented support set.
     *
     * \return size of the vector storing the support set
     */
    std::size_t size() const;

    /*!
     * Returns the number of variables in the represented support set.
     *
     * \return snumber of variables
     */
    std::size_t vars() const;

    /*!
     * Returns true if the set contains the variable with the given \e index.
     *
     * \param index variable index
     *
     * \return true if variable is found
     */
    bool hasVar(std::size_t index) const;

    /*!
     * Inserts the variable with the given \e index. If the variable already
     * exists in the set, nothing happens.
     *
     * \param index variable index
     */
    void addVar(std::size_t index);

    /*!
     * Inserts all variables from the given support set \e s.
     *
     * \param s support set
     */
    void addVars(const VarSupport& s);

    /*!
     * Create the intersection of both given support sets \e vs1 and  \e vs2.
     *
     * \param vs1 first support set
     * \param vs2 second support set
     *
     * \return intersection
     */
    static VarSupport intersect(const VarSupport& vs1, const VarSupport& vs2);

    /*!
     * Returns true if \e s is equal to this set; otherwise returns false.
     * Two sets are considered equal if they contain the same elements.
     *
     * \param s other support set
     *
     * \return true if equal
     */
    bool operator==(const VarSupport& s) const;

    /*!
     * Returns true if \e s is NOT equal to this set; otherwise returns false.
     * Two sets are considered equal if they contain the same elements.
     *
     * \param s other support set
     *
     * \return true if unequal
     */
    bool operator!=(const VarSupport& s) const;

    /*!
     * Returns true if \e s and this set intersect, e.g. both contain at least
     * one common item.
     *
     * \param s other support set
     *
     * \return true if intersection is non-empty
     */
    bool intersects(const VarSupport& s) const;

   protected:
    std::vector<Mask> _support;

    friend bool isSubsetOf(const VarSupport& subset, const VarSupport& superset);

    friend std::ostream& operator<<(std::ostream& os, const VarSupport& s);
};

/*!
 * Returns true if \e subset is a subset of \e superset.
 *
 * \param subset
 * \param superset
 *
 * \return true if subset
 */
bool isSubsetOf(const VarSupport& subset, const VarSupport& superset);

/*!
 * Pretty printer for VarSupport.
 *
 * \param os output stream
 * \param s VarSupport to be printed
 *
 * \return modified stream
 */
std::ostream& operator<<(std::ostream& os, const VarSupport& s);
}  // namespace aigpp

#include "VarSupport.icc"

#endif /* AIGPP_VARSUPPORT_HH */
