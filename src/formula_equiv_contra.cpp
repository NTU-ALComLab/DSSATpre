/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/17 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
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
#include <stack>
#include <utility>
#include <vector>

#include <easylogging++.hpp>

#include "aux.hpp"
#include "bool_vector.hpp"
#include "clause.hpp"
#include "exceptions.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"

/**
 * \file formula_equiv_contra.cpp
 * \brief Detection of equivalent variables
 * \author Ralf Wimmer
 * \date 02/2016
 */

namespace hqspre {

/**
 * \brief Tarjan's algorithm for computing strongly connected components.
 */
static void
visitSCC(const std::vector<Formula::ImplicationSet>& implications, const Literal n, int& calls, std::vector<int>& scc,
         std::vector<int>& root, int& nextscc, std::stack<Literal>& stk)
{
    if (root[n] != -1) {
        return;
    }

    root[n] = calls;
    scc[n]  = -1;
    stk.push(n);

    const int savecalls = calls;
    ++calls;

    for (const BinaryClause& clause : implications[n]) {
        const Literal impl = clause.getLiteral();
        visitSCC(implications, impl, calls, scc, root, nextscc, stk);

        if (scc[impl] == -1 && root[impl] != -1 && root[impl] < root[n]) {
            root[n] = root[impl];
        }
    }

    if (root[n] == savecalls) {
        Literal m = 0;
        do {
            m = stk.top();
            stk.pop();

            scc[m] = nextscc;
        } while (m != n);

        ++nextscc;
    }
}

/**
 * \brief Tries to find equivalent literals by SCC decomposition of the
 * implication graph.
 *
 * This function analyzes the graph formed by the implications.
 * All literals in the same strongly connected component are equivalent.
 * They are replaced by a single representative.
 * \throws UNSATException if an SCC contains more than one universal literal.
 * \throws UNSATException if an SCC contains a universal literal and an
 * existential one that does not depend upon it. \return true if the formula was
 * modified.
 */
bool
Formula::findEquivalences()
{
    if (!_implications_added) {
        return false;
    }

    VLOG(3) << __FUNCTION__;

    _implications_added = false;

    ScopeTimer  equiv(getTimer(WhichTimer::EQUIV_LITS));
    std::size_t count = 0;

    // First determine the strongly connected components of the implication graph
    std::vector<int>    root(maxLitIndex() + 1, -1);
    std::vector<int>    scc(maxLitIndex() + 1, -1);
    int                 calls   = 0;
    int                 nextscc = 0;
    std::stack<Literal> stk;

    std::vector<Literal> candidates;
    candidates.reserve(maxLitIndex() + 1);
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (!varDeleted(var) && !_implications[var2lit(var, false)].empty()
            && !_implications[var2lit(var, true)].empty()) {
            candidates.push_back(var2lit(var, false));
            candidates.push_back(var2lit(var, true));
        }
    }

    for (const Literal lit : candidates) {
        visitSCC(_implications, lit, calls, scc, root, nextscc, stk);
        if (_interrupt) {
            _implications_added = true;
            return false;
        }
    }

    // Check if there is an SCC which contains a variable both
    // positively and negatively. Then the formula is UNSAT.
    std::vector<std::vector<Literal>> sccs(nextscc);

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        const Literal pos = var2lit(var, false);
        const Literal neg = var2lit(var, true);

        if (scc[pos] != -1 && scc[pos] == scc[neg]) {
            // var and ~var in the same SCC -> UNSAT
            throw UNSATException("A literal and its negation in one SCC");
        }
        if (scc[pos] != -1) {
            sccs[scc[pos]].push_back(pos);
        }
        if (scc[neg] != -1) {
            sccs[scc[neg]].push_back(neg);
        }
    }

    // Now 'sccs' is a vector of all SCCs.
    for (const auto& p : sccs) {
        if (_interrupt) {
            break;
        }
        if (p.size() >= 2) {
            count += replaceSCC(p);
        }
    }

    stat(Statistics::EQUIV_LITS) += count;

    VLOG_IF(count > 0, 2) << __FUNCTION__ << " found " << count << " equivalences";

    return count > 0;
}

/**
 * \brief Replaces an SCC in the implication graph by a single literal.
 *
 * \throws UNSATException if an SCC contains more than one universal literal.
 * \throws UNSATException if an SCC contains a universal literal and an
 * existential one that does not depend upon it. \return the number of
 * eliminated variables.
 */
std::size_t
Formula::replaceSCC(const std::vector<Literal>& scc)
{
    val_assert(_prefix);

    if (scc.size() <= 1) {
        return 0;
    }

    // Moved here so that it's defined
    std::size_t count = 0;

    // Modified: For (D)SSAT, multiple universals doesn't imply UNSAT
    if (isStochastic()) {
        // Store exists and randoms separately
        std::vector<Literal> random_lits;
        std::vector<Literal> exist_lits;
        for (const Literal q : scc) {
            if (isUniversal(lit2var(q))) {
                random_lits.push_back(q);
            }
            else {
                exist_lits.push_back(q);
            }
        }  // end loop over scc

        // Handle exists:
        // -Cannot merge exists with incomparable dependencies
        if (exist_lits.size() > 1) {
            // Handling exists:
            // - For SSAT, replace all literals with the leftmost one
            // - For DSSAT, replace a literal l1 by l2 if dep(l1) is a subset of dep(l2)
            for (auto it1 = exist_lits.begin(); it1 != exist_lits.end(); ++it1) {
                    const Literal l1 = *it1;
                    const Variable v1 = lit2var(l1);
                    if (varDeleted(v1)) {
                        continue;
                    }
                for (auto it2 = it1 + 1; it2 != exist_lits.end(); ++it2) {
                    const Literal l2 = *it2;
                    const Variable v2 = lit2var(l2);
                    if (varDeleted(v2)) {
                        continue;
                    }
                    // If dep(v2) <= dep(v1), we can replace l1 with l2,
                    // and skip to next l1 since v1 will become deleted
                    if (dependenciesSubset(v2, v1)) {
                        VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(l1) << " by exist. literal "
                                << lit2dimacs(l2);
                        replaceLiteral(l1, l2);
                        ++count;
                        break;
                    }
                    // If dep(v1) <= dep(v2), we can replace l2 with l1.
                    // Since v1 still exists, we continue from the next l2
                    if (dependenciesSubset(v1, v2)) {
                        VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(l2) << " by exist. literal "
                                << lit2dimacs(l1);
                        replaceLiteral(l2, l1);
                        ++count;
                    }
                } // end inner loop
            } // end outer loop
        }  // end if (exist_lits.size() > 1)

        // Handling randoms
        Literal rand_replacement;
        Variable rand_replacement_var;

        if (!random_lits.empty()) {
            rand_replacement = random_lits.front();
            rand_replacement_var = lit2var(rand_replacement);
        }

        if (random_lits.size() > 1) {
            // More than one random variable in an SCC -> introduce x_eq and x_same:
            // Pr[Psi] = Pr[x_same] * Pr[Psi[x_eq/x_1,...,x_n]],
            // where Pr[x_same]=Pr[x_1=x_2=...=x_n],
            // and Pr[x_eq]=Pr[x_1=x_2=..=x_n=1 | x_1=x_2=..=x_n]
            // For the implementation,
            // x_same is directly recorded using prob_elim
            if (varDeleted(rand_replacement_var)) {
                return count;
            }
            VLOG(2) << "Multiple random literals in an equivalence class.";
            // Handling dependencies:
            // For SSAT, we chose the leftmost literal as representative
            if (_prefix->type() == PrefixType::QBF) {
                for (const Literal q : random_lits) {
                    if (q == rand_replacement) {
                        continue;
                    }
                    const Variable v = lit2var(q);
                    if (varDeleted(v)) {
                        continue;
                    }
                    if (getLevel(v) < getLevel(rand_replacement_var)) {
                        rand_replacement = q;
                        rand_replacement_var = v;
                    }
                }
            }  // end if SSAT
            // For DSSAT, we have to update the dependency of each exist that
            // depends on any q in random_lits to depend on the representative as well.
            else {
                val_assert(_prefix->type() == PrefixType::DQBF);
                const std::vector<Variable>& existVars = _prefix->getExistVars();
                for (const Variable e : existVars) {
                    if (depends(rand_replacement_var, e)) {
                        continue;
                    }
                    for (const Literal q : random_lits) {
                        if (varDeleted(lit2var(q))) {
                            continue;
                        }
                        if (depends(lit2var(q), e)) {
                            addDependency(rand_replacement_var, e);
                            break;
                        }
                    }
                }
            }  // end if DSSAT
            // End handling dependencies

            double prob_pos = (isPositive(rand_replacement)) ? _prefix->getProb(rand_replacement_var) : (1 - _prefix->getProb(rand_replacement_var));
            double prob_neg = 1 - prob_pos;

            // Replace all randoms with replacement
            for (const Literal q : random_lits) {
                if (q == rand_replacement) {
                    continue;
                }
                const Variable v = lit2var(q);
                if (varDeleted(v)) {
                    continue;
                }

                VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(q) << " by rand. literal "
                        << lit2dimacs(rand_replacement);
                double prob = (isPositive(q)) ? _prefix->getProb(v) : (1 - _prefix->getProb(v));
                prob_pos *= prob;
                prob_neg *= (1-prob);
                replaceLiteral(q, rand_replacement);
                ++count;
            } // End replacing randoms

            // Handling probability
            double prob_same = prob_pos + prob_neg;
            double prob_eq = isPositive(rand_replacement) ? (prob_pos / prob_same) : (prob_neg / prob_same);
            val_assert(prob_same >=0 && prob_same <= 1);
            val_assert(prob_eq >=0 && prob_eq <= 1);
            // -1: no random variable eliminated
            if (_elim_prob < 0) {
                _elim_prob = prob_same;
            }
            else {
                _elim_prob *= prob_same;
            }

            _prefix->setProb(rand_replacement_var, prob_eq);
        }  // end if (random_lits.size() > 1)

        if (!random_lits.empty() && !exist_lits.empty()) {
            val_assert(!varDeleted(rand_replacement_var));
            val_assert(isUniversal(rand_replacement_var));

            // For each survived exist literal, check whether it depends on the random representative.
            // If yes: the exist literal can be replaced by the random representative.
            for (const Literal q : exist_lits) {
                const Variable v = lit2var(q);
                if (varDeleted(v)) {
                    continue;
                }
                if (depends(rand_replacement_var, v)) {
                    VLOG(2) << "Exist literal equivalent to a random one in its dependency.";
                    replaceLiteral(q, rand_replacement);
                    val_assert(_occ_list[q].empty() && _occ_list[negate(q)].empty());

                    VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(q) << " by rand. literal "
                            << lit2dimacs(rand_replacement);
                    ++count;
                }
            }
        }  // end if (hasRand && hasExist)
    } else {  // original implementation for (D)QBF
        // count universals on scc
        unsigned int countUniversal = 0;
        Literal universal_literal = 0;

        for (const Literal q : scc) {
            if (isUniversal(lit2var(q))) {
                ++countUniversal;
                universal_literal = q;
            }
        }  // end loop over scc

        // std::size_t count = 0;
        
        if (countUniversal > 1) {
            // More than one universal variable in an SCC -> UNSAT
            throw UNSATException("More than one universal variable in one SCC");
        } else if (countUniversal == 1) {
            // Check if all existential variables in this SCC depend on the universal
            // one. If no: UNSAT If yes: all existential variables can be replaced by
            // the universal one.

            const Variable univ_var = lit2var(universal_literal);
            val_assert(isUniversal(univ_var));

            for (const Literal q : scc) {
                if (q == universal_literal) {
                    continue;
                }
                if (varDeleted(lit2var(q))) {
                    continue;
                }

                const Variable current_var = lit2var(q);
                val_assert(isExistential(current_var));

                if (!depends(current_var, univ_var)) {
                    // The current SCC contains a universal variable and
                    // an existential one that does not depend on the universal.
                    // -> UNSAT
                    throw UNSATException(
                        "Existential variable equivalent to a universal "
                        "one, but independent of it");
                } else {
                    // Replace q by universal_literal:
                    replaceLiteral(q, universal_literal);
                    val_assert(_occ_list[q].empty() && _occ_list[negate(q)].empty());

                    VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(q) << " by univ. literal "
                            << lit2dimacs(universal_literal);

                    ++count;
                }
            }  // end for
        }      // end if (countUniversal == 1)
        else {
            val_assert(countUniversal == 0);

            // At this position the following holds:
            // - SCC contains more than one variable
            // - SCC contains only existential variables
            // - SCC does not contain x and -x at the same time
            // To do:
            // - replace all variables in the SCC by an arbitrary representative
            // - the dependency set of the representative is the intersection of
            //   all dependency sets

            const Literal replacement = scc.front();
            const Variable replacement_var = lit2var(replacement);
            if (varDeleted(replacement_var)) {
                return count;
            }
            // new dependencies for the equivalent variables. At the end it should
            // contain the intersection of the dependency sets of all equivalent
            // variables.

            for (const Literal q : scc) {
                if (q == replacement) {
                    continue;
                }
                if (varDeleted(lit2var(q))) {
                    continue;
                }

                VLOG(3) << __FUNCTION__ << "(): replacing " << lit2dimacs(q) << " by exist. literal "
                        << lit2dimacs(replacement);

                _prefix->intersectDependencies(replacement_var, lit2var(q));
                replaceLiteral(q, replacement);
                count++;
            }
        }  // end if (countUniversal == 0)
    }

    return count;
}

/**
 * \brief Tries to find binary clauses which complete an equivalence or
 * contradiction definition
 *
 * If the clause (a * b), exists, the method tries to add (~a * ~b), (~a * b)
 * and (a * b) Therefore it is checked whether these clauses are (hidden)
 * blocked within the formula. Also performs the variable replacement and/ unit
 * propagation if the search was successful \return true if the formula was
 * modified
 */
bool
Formula::findHiddenEquivAndContraDefinitions()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer hidden_ec(getTimer(WhichTimer::HIDDEN_EQUIV_CONTRA));
    _process_limit.setLimit(PreproMethod::HEC);

    const std::size_t old_stat_hidden_unit  = stat(Statistics::HIDDEN_UNIT);
    const std::size_t old_stat_hidden_equiv = stat(Statistics::HIDDEN_EQUIV_LITS);

    Clause::ClauseData   binary_clause(2, 0);
    std::vector<Literal> scc;
    std::uint64_t        sign = 0ul;

    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        if (_interrupt) {
            break;
        }
        binary_clause[0] = lit;
        // reset equivalence container
        scc.clear();
        scc.push_back(lit);

        if (_process_limit.reachedLimit()) {
            VLOG(2) << "Terminate " << __FUNCTION__ << " due to process limit.";
            break;
        }

        for (BinaryClause bin : _implications[lit]) {
            // remove hidden literals from previous run
            binary_clause.resize(2);
            // first search for equivalence
            binary_clause[1] = negate(bin.getLiteral());
            val_assert(!_seen[binary_clause[0]]);
            val_assert(!_seen[binary_clause[1]]);
            _seen[binary_clause[0]] = 1u;
            _seen[binary_clause[1]] = 1u;
            // Add hidden literals and check for tautology
            bool found = addHiddenLiterals(-1, binary_clause, sign);

            if (!found) {
                // Clause is no hidden tautology => check for hidden blocked clause
                found = (clauseBlocked(binary_clause) != 0u);
            }
            clearSeen(binary_clause);

            if (found) {
                VLOG(3) << __FUNCTION__ << " found equiv binary " << lit2dimacs(binary_clause[0]) << " "
                        << lit2dimacs(binary_clause[1]) << " => " << lit2dimacs(lit) << " and "
                        << lit2dimacs(bin.getLiteral()) << " are hidden equivalent.";

                // We have to add the clause temporarily, otherwise upcoming checks may
                // become unsound
                Clause::ClauseData added_clause{binary_clause[0], binary_clause[1]};
                addClause(std::move(added_clause), true, false, ClauseStatus::OPTIONAL);
                scc.push_back(bin.getLiteral());
            }

            // remove hidden literals
            binary_clause.resize(2);
            // now search for contradiction
            binary_clause[1] = bin.getLiteral();
            val_assert(!_seen[binary_clause[0]]);
            val_assert(!_seen[binary_clause[1]]);
            _seen[binary_clause[0]] = 1u;
            _seen[binary_clause[1]] = 1u;
            // Add hidden literals and check for tautology
            found = addHiddenLiteralsBinary(-1, binary_clause, sign);

            if (!found) {
                // Clause is no hidden tautology => check for hidden blocked clause
                found = (clauseBlocked(binary_clause) != 0u);
            }
            clearSeen(binary_clause);

            if (found) {
                VLOG(3) << __FUNCTION__ << " found contradiction binary (hidden tautology) "
                        << lit2dimacs(binary_clause[0]) << " " << lit2dimacs(binary_clause[1]) << " => "
                        << lit2dimacs(binary_clause[1]) << " is a hidden unit";

                // We have to add the clause temporarily, otherwise upcoming checks may
                // become unsound
                Clause::ClauseData added_clause{binary_clause[0], binary_clause[1]};
                addClause(std::move(added_clause), true, false, ClauseStatus::OPTIONAL);
                pushUnit(binary_clause[1], PureStatus::UNIT);
                ++stat(Statistics::HIDDEN_UNIT);
            }
        }

        // After all checks for a literal have been done, we can now propagate the
        // found units and equivalences It's possibly unsafe to perform these
        // operations immediately after a unit/equivalence was detected
        unitPropagation();
        if (scc.size() > 1) {
            const std::size_t new_equiv = replaceSCC(scc);
            stat(Statistics::EQUIV_LITS) += new_equiv;
            stat(Statistics::HIDDEN_EQUIV_LITS) += new_equiv;
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << (stat(Statistics::HIDDEN_UNIT) - old_stat_hidden_unit)
            << " new units via contradictions and " << (stat(Statistics::HIDDEN_EQUIV_LITS) - old_stat_hidden_equiv)
            << " new equivalences";

    return (stat(Statistics::HIDDEN_UNIT) > old_stat_hidden_unit)
           || (stat(Statistics::HIDDEN_EQUIV_LITS) > old_stat_hidden_equiv);
}

/**
 * \brief Tries to find unit literals by analyzing contradicting implication
 * chains.
 *
 * If a literal \f$a\f$ transitively implies \f$\neg a\f$, then \f$\neg a\f$ is
 * a unit literal. This is done by investigating the implication graph defined
 * by the data structure 'implications'. \return true if the formula was
 * modified
 */
bool
Formula::findContradictions()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer contra(getTimer(WhichTimer::CONTRADICTIONS));

    _process_limit.setLimit(PreproMethod::CONTRA);

    unsigned int       count    = 0;
    const unsigned int litcount = maxLitIndex() + 1;

    std::vector<int> usedSource;
    usedSource.reserve(litcount);
    std::vector<int> usedTarget;
    usedTarget.reserve(litcount);
    std::vector<int> usedSourceMap(litcount, -1);
    std::vector<int> usedTargetMap(litcount, -1);

    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        if (_interrupt) {
            return false;
        }

        if (!_implications[lit].empty()) {
            usedSourceMap[lit] = static_cast<int>(usedSource.size());
            usedSource.push_back(lit);

            usedTargetMap[negate(lit)] = static_cast<int>(usedTarget.size());
            usedTarget.push_back(negate(lit));
        }
    }

    using ReachMatrix = std::vector<BoolVector>;
    ReachMatrix reach(usedSource.size());

    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        if (_interrupt) {
            return false;
        }
        if (usedSourceMap[lit] == -1) {
            continue;
        }

        BoolVector& reach_n = reach[usedSourceMap[lit]];
        val_assert(reach_n.uninitialized());
        reach_n.initialize(usedTarget.size(), false);

        _process_limit.decreaseLimitBy(_implications[lit].size());

        for (BinaryClause impl : _implications[lit]) {
            reach_n.set(usedTargetMap[impl.getLiteral()], true);
        }
    }

    for (unsigned int via = 0; via != usedSource.size(); ++via) {
        if (_interrupt) {
            return false;
        }

        const int viaTarget = usedTargetMap[usedSource[via]];
        if (viaTarget == -1) {
            continue;
        }

        if (_process_limit.reachedLimit()) {
            VLOG(2) << __FUNCTION__ << " Building data structure terminated due to process limit.";

            // Deactivate this method if we run into our limits.
            // In this case, constructing the matrix is too expensive.
            _settings.contradictions = false;
            break;
        }

        _process_limit.decreaseLimitBy(2, usedSource.size());

        const auto& reach_via = reach[via];
        for (std::size_t from = 0; from != usedSource.size(); ++from) {
            BoolVector& reach_from = reach[from];

            if (!reach_from.get(viaTarget)) {
                continue;
            }

            reach_from |= reach_via;
        }
    }

    // check for i -> !i (i.e., we have the unit literal !i)
    for (std::size_t i = 0; i != usedSource.size(); ++i) {
        const int j = usedTargetMap[negate(usedSource[i])];

        if (reach[i].get(j)) {
            VLOG(3) << __FUNCTION__ << ": " << lit2dimacs(usedSource[i]) << " -> " << lit2dimacs(usedTarget[j]);

            ++count;
            pushUnit(usedTarget[j], PureStatus::UNIT);
        }
    }

    unitPropagation();

    stat(Statistics::CONTRADICTIONS) += count;
    VLOG(2) << __FUNCTION__ << " found " << count << " contradictions.";

    return count > 0;
}

}  // end namespace hqspre
