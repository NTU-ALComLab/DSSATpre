/*
 * This file is part of HQSpre.
 *
 * Copyright 2016-19 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>  // Modified: std::setprecision
#include <limits>  // Modified: std::numeric_limits

#include <easylogging++.hpp>
#include "aux.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"

/**
 * \file formula_inout.cpp
 * \brief Implementation of methods and functions for reading and writing formulas.
 * \author Ralf Wimmer \date 01/2016
 */

namespace hqspre {

static inline std::uint32_t
str2uint(const std::string& token, const std::int32_t min, const std::int32_t max)
{
    std::size_t   pos   = 0;
    unsigned long value = 0ul;

    try {
        value = std::stoul(token, &pos, 10);
    } catch (std::invalid_argument& e) {
        throw FileFormatException("Invalid integer");
    } catch (std::out_of_range& e) {
        throw FileFormatException("Integer out of range.");
    }

    if (value < static_cast<unsigned long>(min) || value > static_cast<unsigned long>(max)) {
        throw FileFormatException("Integer out of range.");
    }

    if (pos < token.size()) {
        throw FileFormatException("Invalid integer");
    }

    return static_cast<std::uint32_t>(value);
}

static inline std::int32_t
str2int(const std::string& token, const std::int32_t min, const std::int32_t max)
{
    std::size_t pos   = 0;
    long        value = 0l;

    try {
        value = std::stol(token, &pos, 10);
    } catch (std::invalid_argument& e) {
        throw FileFormatException("Invalid integer");
    } catch (std::out_of_range& e) {
        throw FileFormatException("Integer out of range.");
    }

    if (value < min || value > max) {
        throw FileFormatException("Integer out of range.");
    }

    if (pos < token.size()) {
        throw FileFormatException("Invalid integer");
    }

    return static_cast<std::int32_t>(value);
}

// Modified: added function to read probability (double within 0 1)
static inline double
str2prob(const std::string& token)
{
    std::size_t   pos   = 0;
    double        value = 0;

    try {
        value = std::stod(token, &pos);
    } catch (std::invalid_argument& e) {
        throw FileFormatException("Invalid double");
    } catch (std::out_of_range& e) {
        throw FileFormatException("Value out of range.");
    }

    if (value < 0 || value > 1) {
        throw FileFormatException("Probability out of range.");
    }

    if (pos < token.size()) {
        throw FileFormatException("Invalid double");
    }

    return value;
}

/**
 * \brief Reads a formula in DQDIMACS format from an input stream
 *
 * If the current formula is not empty, it is cleared first by calling
 * Formula::reset().
 */
void
Formula::read(std::istream& stream)
{
    val_assert(_prefix);

    reset();
    // Modified: recover stochastic after reset
    _prefix->setStochastic(_enforce_stochastic);

    std::string        token;
    std::uint32_t      num_clauses = 0;
    std::set<Variable> universal_vars;

    // Read the (DQ)DIMACS file header.
    while (!stream.eof()) {
        stream >> token;
        if (token == "c") {
            std::getline(stream, token);  // consume the rest of the line
        } else if (token == "p") {
            // DIMACS file header
            std::string s;
            stream >> s;
            if (s != "cnf") {
                throw FileFormatException("Input is not a (DQ)DIMACS file");
            }

            stream >> token;
            const std::uint32_t num_vars = str2uint(token, 0, _settings.max_num_vars);

            stream >> token;
            num_clauses = str2uint(token, 0, _settings.max_num_clauses);

            _clauses.reserve(num_clauses);
            setMaxVarIndex(static_cast<Variable>(num_vars));
            break;
        } else {
            throw FileFormatException("File header missing.");
        }
    }

    if (stream.eof() && num_clauses != 0) {
        throw FileFormatException("File header missing.");
    }

    // Read the variable declarations */
    while (!stream.eof()) {
        stream >> token;

        if (token == "a") {
            // list of universally quantified variables
            while (stream) {
                stream >> token;
                const Variable var = str2uint(token, 0, _settings.max_num_vars);
                if (var == 0) {
                    break;
                }

                setUVar(var);
                universal_vars.insert(var);
                // Modified: let probability be 0.5 if --stochastic is used
                if (isStochastic()) {
                    _prefix->setProb(var, 0.5);
                }
            }
        } else if (token == "e") {
            // existential variable depending on all universal vars
            // declared so far.
            while (stream) {
                stream >> token;
                const Variable exist_var = str2uint(token, 0, _settings.max_num_vars);
                if (exist_var == 0) {
                    break;
                }
                if (_prefix->type() == PrefixType::DQBF) {
                    setEVar(exist_var, universal_vars);
                } else {
                    setEVar(exist_var);
                }
            }
        } else if (token == "d") {
            // If we have a QBF prefix and encounter a "d"-line,
            // we have to switch to DQBF.
            if (_prefix->type() == PrefixType::QBF) {
                VLOG(2) << "Switching to DQBF prefix representation.";

                _dqbf_prefix = _qbf_prefix->convertToDQBF();
                delete _prefix;
                _prefix     = _dqbf_prefix;
                _qbf_prefix = nullptr;
            }
            // existential variable with dependencies
            stream >> token;
            const Variable     exist_var = str2uint(token, 1, _settings.max_num_vars);
            std::set<Variable> deps;

            while (stream) {
                stream >> token;
                const Variable all_var = str2uint(token, 0, static_cast<std::int32_t>(maxVarIndex()));
                if (all_var == 0) {
                    break;
                }
                if (all_var > maxVarIndex()) {
                    throw FileFormatException("[ERROR] Undeclared variable in dependency set");
                }
                if (!isUniversal(all_var)) {
                    throw FileFormatException("[ERROR] Variable in the dependency set that is not universal");
                }
                deps.insert(all_var);
            }
            setEVar(exist_var, std::move(deps));
        } else if (token == "c") {
            std::getline(stream, token);  // consume the rest of the line
        } else if (token == "r") {  // Modified: handle stochastic instances
            val_assert(isStochastic() || numUVars() == 0);  // Reject case where both a and r are used
            if (!isStochastic()) {
                VLOG(2) << "Switching to (D)SSAT prefix representation.";
                setStochastic();
            }
            // Read probability
            stream >> token;
            const double prob = str2prob(token);
            // list of random quantified variables
            while (stream) {
                stream >> token;
                const Variable var = str2uint(token, 0, _settings.max_num_vars);
                if (var == 0) {
                    break;
                }

                setUVar(var);
                universal_vars.insert(var);
                _prefix->setProb(var, prob);
            }
        } else {
            break;
        }
    }

    // We have read all variable declearations. Re-generate the list of deleted variables.
    _deleted_var_numbers.clear();
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) _deleted_var_numbers.push_back(var);
    }

    // Now read clauses
    std::uint32_t      clauses_read = 0;
    Clause::ClauseData clause;

    while (stream.good()) {
        clause.clear();

        do {
            const std::int32_t lit = str2int(token, -maxVarIndex(), +maxVarIndex());
            if (lit == 0) {
                break;
            }
            const Variable var = static_cast<Variable>(abs(lit));
            clause.push_back(var2lit(var, lit < 0));

            if (clause.back() < minLitIndex() || clause.back() > maxLitIndex()) {
                throw FileFormatException("Invalid literal in clause");
            }
            if (stream.eof()) {
                throw FileFormatException("Unexpected end of file.");
            }
            stream >> token;
        } while (true);

        addClause(std::move(clause));
        ++clauses_read;

        if (clauses_read == num_clauses) {
            break;
        }
        stream >> token;
    }

    // Check if we can switch from DQBF to QBF
    if (!_enforce_dqbf && _prefix->type() == PrefixType::DQBF && isQBF()) {
        val_assert(!_qbf_prefix && _dqbf_prefix);

        VLOG(2) << "Switching to QBF prefix representation.";
        _qbf_prefix = _dqbf_prefix->convertToQBF();
        delete _dqbf_prefix;
        _dqbf_prefix = nullptr;
        _prefix      = _qbf_prefix;
    }

    VLOG(1) << "Read formula with " << numEVars() << " existential vars, " << numUVars() << " univ. vars, and "
            << numClauses() << " clauses.";
}


/**
 * \brief Writes the current formula in DQDIMACS format to the given output
 * stream.
 *
 * If the formula is a QBF, then the output file is actually in QDIMACS format.
 * \param stream output stream to which the formula should be written
 * \param compact rename the variables such that there are no deleted variables
 * in between
 */
void
Formula::write(std::ostream& stream, bool compact) const
{
    val_assert(_prefix);

    if (!stream) {
        LOG(ERROR) << "Could not write to output stream!\n";
        std::exit(-1);
    }

    const_cast<Formula*>(this)->updateVars();
    val_assert(checkConsistency());

    std::vector<Variable> translation_table(maxVarIndex() + 1, 0);
    Variable              current = 0;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) {
            continue;
        }
        if (compact) {
            translation_table[var] = ++current;
        } else {
            translation_table[var] = var;
        }
    }

    // If necessary, print the translation table
    if (compact) {
        for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
            if (!varDeleted(var)) {
                stream << "c ren " << var << " -> " << translation_table[var] << '\n';
            }
        }
    }
    // Modified: Insert one more random if pure random eliminated

    if (_elim_prob >= 0) {
        VLOG(1) << "Unit random variables eliminated. Dummy random variable inserted to keep probability.";
        // print the header
        stream << "p cnf " << (compact ? current + 1 : maxVarIndex() + 1) << ' ' << numClauses() + 1 << '\n';
        _prefix->write(stream, &translation_table);
        stream << "r " << std::setprecision(std::numeric_limits<double>::max_digits10) << _elim_prob << " " << (compact ? current + 1 : maxVarIndex() + 1) << " 0\n";
    } else {
        // print the header
        stream << "p cnf " << (compact ? current : maxVarIndex()) << ' ' << numClauses() << '\n';
        _prefix->write(stream, &translation_table);
    }
    writeClauses(stream, &translation_table);
    if (_elim_prob >= 0)
        stream << (compact ? current + 1 : maxVarIndex() + 1) << " 0\n";
}

/**
 * \brief Print all clauses in DIMACS format to the given stream.
 */
void
Formula::writeClauses(std::ostream& stream, std::vector<Variable>* translation_table) const
{
    auto trans = [translation_table](const Literal lit) -> Literal {
        if (translation_table == nullptr) {
            return lit;
        } else {
            return var2lit((*translation_table)[lit2var(lit)], isNegative(lit));
        }
    };

    // Print the clauses.
    for (ClauseID c_nr = 0; c_nr <= maxClauseIndex(); ++c_nr) {
        if (clauseDeleted(c_nr)) {
            continue;
        }
        for (const Literal lit : _clauses[c_nr]) {
            stream << lit2dimacs(trans(lit)) << ' ';
        }
        stream << "0\n";
    }

    for (const Literal lit : _unit_stack) {
        stream << lit2dimacs(lit) << " 0\n";
    }
}

std::ostream&
operator<<(std::ostream& stream, const Formula& formula)
{
    formula.write(stream);
    return stream;
}

std::istream&
operator>>(std::istream& stream, Formula& formula)
{
    formula.read(stream);
    return stream;
}

}  // end namespace hqspre
