/**************************************************************
 *
 *       AIGPP Package // VariableAssignment.cc
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

#include "VariableAssignment.hh"

#include <cassert>
#include <iostream>

std::ostream&
aigpp::operator<<(std::ostream& os, const aigpp::VariableAssignment& va)
{
    std::size_t vars = va._bins.size() * VariableAssignment::VariablesPerBin;

    for (std::size_t i = 0; i < vars; ++i) {
        VariableAssignment::Assignment a = va.get(i);

        if (a == VariableAssignment::Positive) {
            os << "+" << i;
        } else if (a == VariableAssignment::Negative) {
            os << "-" << i;
        }
    }

    return os;
}
