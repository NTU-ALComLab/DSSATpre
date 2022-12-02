// 64ok
/**************************************************************
 *
 *       AIGPP Package // VarSupport.cc
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

#include "VarSupport.hh"

#include <algorithm>

aigpp::VarSupport::VarSupport(const aigpp::VarSupport& vs1, const aigpp::VarSupport& vs2) :
    _support(((vs1.size() >= vs2.size()) ? vs1._support : vs2._support))
{
    auto p = _support.begin();

    if (vs1.size() >= vs2.size()) {
        auto p2 = vs2._support.begin();

        while (p2 != vs2._support.end()) {
            *p |= *p2;
            ++p2;
            ++p;
        }
    } else {
        auto p1 = vs1._support.begin();

        while (p1 != vs1._support.end()) {
            *p |= *p1;
            ++p1;
            ++p;
        }
    }
}

void
aigpp::VarSupport::set(const aigpp::VarSupport& vs1, const aigpp::VarSupport& vs2)
{
    if (vs1.size() >= vs2.size()) {
        _support = vs1._support;

        auto p  = _support.begin();
        auto p2 = vs2._support.begin();

        while (p2 != vs2._support.end()) {
            *p |= *p2;
            ++p2;
            ++p;
        }
    } else {
        _support = vs2._support;

        auto p  = _support.begin();
        auto p1 = vs1._support.begin();

        while (p1 != vs1._support.end()) {
            *p |= *p1;
            ++p1;
            ++p;
        }
    }
}

std::size_t
aigpp::VarSupport::vars() const
{
    std::size_t count = 0;

    for (Mask p : _support) {
        if (p == 0ul) {
            continue;
        }

        for (std::size_t i = 0; i < MaskSize; ++i) {
            if (p & (1ul << i)) {
                ++count;
            }
        }
    }

    return count;
}

void
aigpp::VarSupport::addVars(const aigpp::VarSupport& s)
{
    if (size() >= s.size()) {
        auto p  = _support.begin();
        auto p2 = s._support.begin();

        while (p2 != s._support.end()) {
            *p |= *p2;
            ++p;
            ++p2;
        }
    } else {
        auto p  = _support.begin();
        auto p2 = s._support.begin();

        while (p != _support.end()) {
            *p |= *p2;
            ++p;
            ++p2;
        }

        while (p2 != s._support.end()) {
            _support.push_back(*p2);
            ++p2;
        }
    }
}

aigpp::VarSupport
aigpp::VarSupport::intersect(const aigpp::VarSupport& vs1, const aigpp::VarSupport& vs2)
{
    VarSupport v;

    v._support.resize(std::min(vs1.size(), vs2.size()));

    auto p1 = vs1._support.begin();
    auto p2 = vs2._support.begin();
    auto p  = v._support.begin();

    while (p != v._support.end()) {
        *p = *p1 & *p2;
        ++p1;
        ++p2;
        ++p;
    }

    return v;
}

bool
aigpp::VarSupport::intersects(const aigpp::VarSupport& v) const
{
    auto p1 = _support.begin();
    auto p2 = v._support.begin();

    if (size() <= v.size()) {
        while (p1 != _support.end()) {
            if (((*p1) & (*p2)) != 0ul) {
                return true;
            }
            ++p1;
            ++p2;
        }
        return false;
    } else {
        while (p2 != v._support.end()) {
            if (((*p1) & (*p2)) != 0) {
                return true;
            }
            ++p1;
            ++p2;
        }
        return false;
    }
}

bool
aigpp::isSubsetOf(const aigpp::VarSupport& subset, const aigpp::VarSupport& superset)
{
    auto p1 = subset._support.begin();
    auto p2 = superset._support.begin();

    if (subset.size() <= superset.size()) {
        while (p1 != subset._support.end()) {
            if (((*p1) & (*p2)) != *p1) {
                return false;
            }
            ++p1;
            ++p2;
        }
        return true;
    } else {
        while (p2 != superset._support.end()) {
            if (((*p1) & (*p2)) != *p1) {
                return false;
            }
            ++p1;
            ++p2;
        }
        while (p1 != subset._support.end()) {
            if (*p1 != 0ul) {
                return false;
            }
            ++p1;
        }

        return true;
    }
}

std::ostream&
aigpp::operator<<(std::ostream& os, const aigpp::VarSupport& s)
{
    for (std::size_t i = 0; i != s._support.size(); ++i) {
        std::vector<aigpp::VarSupport::Mask>::const_reference m = s._support[i];

        if (m != 0ul) {
            for (std::size_t j = 0; j != aigpp::VarSupport::MaskSize; ++j) {
                if (m & (1ul << j)) {
                    os << (i * aigpp::VarSupport::MaskSize + j) << " ";
                }
            }
        }
    }

    return os;
}
