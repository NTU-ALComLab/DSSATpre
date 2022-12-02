/**************************************************************
 *
 *       AIGPP Package // SimVector.hh
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

#ifndef AIGPP_SIMVECTOR_HH
#define AIGPP_SIMVECTOR_HH

/* std */
#include <cassert>
#include <climits>
#include <iostream>
#include <vector>

/* LRABS utilities */
#include <lrabsutil/Random.hh>

/* local */
#include "Primes.hh"

#define SIMVECSIZE 4096

namespace aigpp {
/*!
 * \class SimVector
 * \brief The SimVector class represents a simulation vector of fixed size.
 */
class SimVector
{
   public:
    /*!
     * \typedef BinType
     * BinType defines the type used for simulation vector bins.
     */
    typedef unsigned long BinType;

    /*!
     * Compile-time constants used in SimVector.
     */
    enum Constants
    {
        BinSize  = CHAR_BIT * sizeof(BinType), /*!< Bit-size of one bin */
        BinCount = SIMVECSIZE / BinSize        /*!< Number of bins */
    };

    /*!
     * Constructor. Creates a SimVector containing 0 only.
     */
    SimVector();

    /*!
     * Copy Constructor
     */
    SimVector(const SimVector& other);

    /*!
     * Constructor. Creates a random SimVector using the supplied random
     * number generator. Bins 0 to maxrandbins-1 are filled with random
     * values, bins maxrandbins to BinCount-1 are filled with 0.
     *
     * \param random
     * \param maxrandbins maximum number of random bins
     */
    explicit SimVector(lrabs::Random& random, std::size_t maxrandbins = BinCount);

    /*!
     * Constructor. Creates a SimVector representing the bitwise conjunction
     * of the two given SimVectors \e s1 and \e s2. If \e c1 (\e c2) is true,
     * \e s1 (\e s2) is complemented before calculating the conjunction.
     *
     * \param s1 first source SimVector
     * \param c1 first complementation flag
     * \param s2 second source SimVector
     * \param c2 second complementation flag
     */
    SimVector(const SimVector& s1, bool c1, const SimVector& s2, bool c2);

    /*!
     * Assignment operator
     */
    SimVector& operator=(const SimVector& other);

    /*!
     * Returns the (fixed) size (number of bins) of the SimVector.
     *
     * \return number of bins
     */
    std::size_t size() const;

    /*!
     * Return the hash-value of the SimVector.
     *
     * \return hash-value
     */
    unsigned long hash() const;

    /*!
     * Return the inverted hash-value of the SimVector.
     *
     * \return inverted hash-value
     */
    unsigned long hashInverted() const;

    void updateHash();

    /*!
     * Creates a random SimVector using the supplied random
     * number generator. Bins 0 to maxrandbins-1 are filled with random
     * values, bins maxrandbins to BinCount-1 are filled with 0.
     *
     * \param random
     * \param maxrandbins maximum number of random bins
     */
    void generateRandom(lrabs::Random& random, std::size_t maxrandbins = BinCount);

    /*!
     * Sets the value of the SimVector bin with index \e bin to \e value.
     *
     * \param bin
     * \param value
     */
    void setBin(std::size_t bin, BinType value);

    /*!
     * Return the bin at position \e index. The access is unchecked,
     * e.g. the caller must ensure that \e index has a valid value.
     *
     * \param index
     *
     * \return bin at position \e index
     */
    const BinType& operator[](std::size_t index) const;

    /*!
     * Returns the bit at \a index. The \a index must be valid!.
     */
    bool getBit(std::size_t index) const;

    /*!
     * Compares \e *this and \e s. Returns true, if \e *this is smaller than
     * \e s in terms of the lexicographical order defined on SimVectors.
     *
     * \param s SimVector to compare with
     *
     * \return true, if \e *this is smaller than \e s
     */
    bool operator<(const SimVector& s) const;

    /*!
     * Compares \e *this and \e s. Returns true, if \e *this is bitwise equal
     * to \e s.
     *
     * \param s SimVector to compare with
     *
     * \return true, if \e *this is equal to \e s
     */
    bool operator==(const SimVector& s) const;

    /*!
     * Compares \e *this and \e s. Returns true, if \e *this is NOT bitwise
     * equal to \e s.
     *
     * \param s SimVector to compare with
     *
     * \return true, if \e *this is NOT equal to \e s
     */
    bool operator!=(const SimVector& s) const;

    /*!
     * Compares \e *this and \e s. Returns true, if \e *this is bitwise equal
     * to the inversion of \e s.
     *
     * \param s SimVector to compare with
     *
     * \return true, if \e *this is equal to the inversion of \e s
     */
    bool equalInverted(const SimVector& s) const;

    /*!
     * Returns a bitwise complemented copy of \e *this.
     *
     * \return complemented copy
     */
    SimVector invert() const;

    /*!
     * Checks if the bits of all bins are 0.
     *
     * \return true, if all bits are 0
     */
    bool isZero() const;

    /*!
     * Checks if the bits of all bins are 1.
     *
     * \return true, if all bits are 1
     */
    bool isOne() const;

    /*!
     * Update all bins using the bitwise conjunction of the
     * (possibly complemented) SimVectors \e s1 and \e s2.
     *
     * \param s1 first source SimVector
     * \param c1 first complementation flag
     * \param s2 second source SimVector
     * \param c2 second complementation flag
     */
    void update(const SimVector& s1, bool c1, const SimVector& s2, bool c2);

    void update(const SimVector& s1, bool c1, const SimVector& s2, bool c2, bool invert);

    /*!
     * Update the <em>bin</em>th bin using the bitwise conjunction of the
     * (possibly complemented) <em>bin</em>th bins of \e s1 and \e s2. If
     * \e invert is true, the result is complemented.
     *
     * \param s1 first source SimVector
     * \param c1 first complementation flag
     * \param s2 second source SimVector
     * \param c2 second complementation flag
     * \param invert result complementation flag
     * \param bin index of the bins to use
     */
    void updateBin(const SimVector& s1, bool c1, const SimVector& s2, bool c2, bool invert, std::size_t bin);

    /*!
     * Computes the bitwise exclusive or of \a s1 and \a s2.
     */
    static SimVector Xor(const SimVector& s1, const SimVector& s2);

   private:
    BinType _bins[BinCount];
    BinType _hash;

    friend std::ostream& operator<<(std::ostream& os, const SimVector& s);
};

/*!
 *   Pretty printer for SimVector.
 *
 *   \param os output stream
 *   \param s SimVector to be printed
 *
 *   \return modified stream
 */
std::ostream& operator<<(std::ostream& os, const SimVector& s);
}  // namespace aigpp

#include "SimVector.icc"

#endif /* AIGPP_SIMVECTOR_HH */
