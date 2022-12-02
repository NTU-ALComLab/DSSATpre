/**************************************************************
 *
 *       AIGPP Package // SimVector.cc
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

#include "SimVector.hh"

/* std */
#include <algorithm>
#include <cassert>

/* local */
#include "Primes.hh"

aigpp::SimVector::SimVector() : _hash(0ul)
{
    std::fill(_bins, _bins + BinCount, 0ul);
}

aigpp::SimVector::SimVector(const aigpp::SimVector& other) : _hash(other._hash)
{
    std::copy(other._bins, other._bins + BinCount, _bins);
}

aigpp::SimVector::SimVector(lrabs::Random& random, std::size_t maxrandbins) : _hash(0ul)
{
    generateRandom(random, maxrandbins);
}

aigpp::SimVector::SimVector(const aigpp::SimVector& s1, bool c1, const aigpp::SimVector& s2, bool c2) : _hash(0ul)
{
    update(s1, c1, s2, c2);
}

aigpp::SimVector&
aigpp::SimVector::operator=(const aigpp::SimVector& other)
{
    if (this != &other) {
        _hash = other._hash;
        std::copy(other._bins, other._bins + BinCount, _bins);
    }

    return *this;
}

void
aigpp::SimVector::updateHash()
{
    _hash = 0ul;

    const BinType*       b = _bins;
    const unsigned long* p = Primes;

    for (std::size_t i = size(); i != 0; --i, ++b, ++p) {
        _hash ^= (*b) * (*p);
    }
}

void
aigpp::SimVector::generateRandom(lrabs::Random& random, std::size_t maxrandbins)
{
    assert(maxrandbins <= BinCount);

    _hash = 0ul;

    for (std::size_t i = 0; i != maxrandbins; ++i) {
        _bins[i] = random.getULong();
        _hash ^= Primes[i] * _bins[i];
    }

    if (maxrandbins < BinCount) {
        std::fill(_bins + maxrandbins, _bins + BinCount, 0ul);
    }
}

aigpp::SimVector
aigpp::SimVector::invert() const
{
    SimVector s = *this;

    s._hash = 0ul;

    for (std::size_t i = 0; i != size(); ++i) {
        s._bins[i] = ~s._bins[i];
        s._hash ^= Primes[i] * s._bins[i];
    }

    return s;
}

void
aigpp::SimVector::update(const aigpp::SimVector& s1, bool c1, const aigpp::SimVector& s2, bool c2)
{
    _hash = 0ul;

    for (std::size_t i = 0; i != size(); ++i) {
        _bins[i] = (c1 ? ~s1[i] : s1[i]) & (c2 ? ~s2[i] : s2[i]);
        _hash ^= Primes[i] * _bins[i];
    }
}

void
aigpp::SimVector::update(const aigpp::SimVector& s1, bool c1, const aigpp::SimVector& s2, bool c2, bool invert)
{
    _hash = 0ul;

    for (std::size_t i = 0; i != size(); ++i) {
        _bins[i] = (c1 ? ~s1[i] : s1[i]) & (c2 ? ~s2[i] : s2[i]);
        if (invert) {
            _bins[i] = ~_bins[i];
        }

        _hash ^= Primes[i] * _bins[i];
    }
}

void
aigpp::SimVector::updateBin(const aigpp::SimVector& s1, bool c1, const aigpp::SimVector& s2, bool c2, bool invert,
                            std::size_t bin)
{
    _hash ^= Primes[bin] * _bins[bin];

    _bins[bin] = (c1 ? ~s1[bin] : s1[bin]) & (c2 ? ~s2[bin] : s2[bin]);

    if (invert) {
        _bins[bin] = ~_bins[bin];
    }

    _hash ^= Primes[bin] * _bins[bin];
}

bool
aigpp::SimVector::equalInverted(const aigpp::SimVector& s) const
{
    const BinType* p1 = _bins;
    const BinType* p2 = s._bins;

    for (std::size_t i = 0; i != size(); ++i, ++p1, ++p2) {
        if (*p1 != ~*p2) {
            return false;
        }
    }

    return true;
}

bool
aigpp::SimVector::isZero() const
{
    const BinType* p = _bins;

    for (std::size_t i = 0; i != size(); ++i, ++p) {
        if (*p != 0ul) {
            return false;
        }
    }

    return true;
}

bool
aigpp::SimVector::isOne() const
{
    const BinType* p = _bins;

    for (std::size_t i = 0; i != size(); ++i, ++p) {
        if (*p != ~0ul) {
            return false;
        }
    }

    return true;
}

aigpp::SimVector
aigpp::SimVector::Xor(const aigpp::SimVector& s1, const aigpp::SimVector& s2)
{
    SimVector result;
    result._hash = 0ul;

    for (std::size_t i = 0; i != BinCount; ++i) {
        result._bins[i] = s1[i] ^ s2[i];
        result._hash ^= Primes[i] * result._bins[i];
    }

    return result;
}

std::ostream&
aigpp::operator<<(std::ostream& os, const aigpp::SimVector& s)
{
    static char charmap[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

    for (std::size_t i = 0; i != s.size(); ++i) {
        SimVector::BinType bin = s[i];

        for (int b = 0; b < SimVector::BinSize / 4; ++b) {
            os << charmap[bin & 0xFul];
            bin >>= 4;
        }
    }
    return os;
}
