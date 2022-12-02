/**************************************************************
 *
 *       AIGPP Package // VariableAssignment.hh
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *      Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#ifndef VARIABLEASSIGNMENT_HH
#define VARIABLEASSIGNMENT_HH

/* stdlib + stl */
#include <cassert>
#include <climits>
#include <iostream>
#include <vector>

/* LRABS utilities */
#include <lrabsutil/Assert.hh>

namespace aigpp {
/*!
 * \ingroup ext
 * \class VariableAssignment
 * \brief The VariableAssignment class is used to store an assignment
 *        to Boolean variables.
 */
class VariableAssignment
{
   public:
    typedef unsigned long BinType;
    enum
    {
        BinSize         = sizeof(BinType) * CHAR_BIT,
        VariablesPerBin = BinSize / 2
    };

    /*!
     * Type of assigment.
     */
    enum Assignment
    {
        Unassigned = 0ul, /*!< Not assigned any value */
        Positive   = 1ul, /*!< Positive assignment */
        Negative   = 2ul, /*!< Negative assignment */
        Mask       = 3ul  /*!< Mask value. Only used internally. */
    };

    /*!
     * Constructor. Creates an empty VariableAssignment.
     * Optionally one can specify the maximum variable index to be assigned
     * in the VariableAssignment. Giving \a maxVar reserves enough space to
     * store variable assignments with indices up to \a maxVar (see reserve()
     * for details).
     */
    explicit VariableAssignment(std::size_t maxVar = 0);

    VariableAssignment(const VariableAssignment& other)     = default;
    VariableAssignment(VariableAssignment&& other) noexcept = default;
    VariableAssignment& operator=(const VariableAssignment& other) = default;
    VariableAssignment& operator=(VariableAssignment&& other) noexcept = default;

    /*!
     * Assign \a a to the variable index \a var.
     */
    void set(std::size_t var, Assignment a);

    /*!
     * Assign \a Positive to the variable index \a var.
     */
    inline void setPositive(std::size_t var);

    /*!
     * Assign \a Negative to the variable index \a var.
     */
    inline void setNegative(std::size_t var);

    /*!
     * Remove the assignment from variable index \a var.
     */
    inline void setUnassigned(std::size_t var);

    /*!
     * Get the assignment of the variable index \a var.
     */
    Assignment get(std::size_t var) const;

    /*!
     * Reserve enough space to store variable assignments with indices up
     * to \a maxVar. Using this method is only optional and may speed up
     * the set() method.
     */
    inline void reserve(std::size_t maxVar);

    /*!
     * Compare this VariableAssignment and \a v. Return true, if both
     * VariableAssignments are equal; otherwise return false.
     */
    bool operator==(const VariableAssignment& v) const;

    /*!
     * Compare this VariableAssignment and \a v. Return true, if both
     * VariableAssignments are not equal; otherwise return false.
     */
    bool operator!=(const VariableAssignment& v) const;

    /*!
     * Compare this VariableAssignment and \a v. Return true, if this
     * VariableAssignment is lexicographically smaller than \a v; otherwise return
     * false.
     */
    bool operator<(const VariableAssignment& v) const;

    /*!
     * Return true, if this VariableAssignment implies \a v, i.e. if a variable is
     * assigned a value by this VariableAssignment the same value or nothing is
     * assigned by \a v.
     */
    bool isSuperSetOf(const VariableAssignment& v) const;

    /*!
     * Returns true, if this VariableAssignment and \a v are compatibe, i.e. if a
     * variable is assigned a value by this VariableAssignment the same value or
     * nothing is assigned by \a v, and vice versa.
     */
    bool compatible(const VariableAssignment& v) const;

    void mergeWith(const VariableAssignment& v);

   private:
    std::vector<BinType> _bins;
    std::size_t          _maxUnassignedBin;

    friend std::ostream& operator<<(std::ostream& os, const VariableAssignment& va);
};

/*!
 * \ingroup ext
 * Pretty prints the VariableAssignment \a va to the output stream \a os.
 */
std::ostream& operator<<(std::ostream& os, const VariableAssignment& va);

}  // namespace aigpp

#include "VariableAssignment.icc"

#endif /* VARIABLEASSIGNMENT_HH */
