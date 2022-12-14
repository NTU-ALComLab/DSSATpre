/*
 * This file is part of HQSpre.
 *
 * Copyright 2016-2021 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
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

#ifndef HQSPRE_FORMULA_IPP_
#define HQSPRE_FORMULA_IPP_

/*
 * \file formula.ipp
 * \brief Inline functions for variable, clause, and dependency management
 * \author Ralf Wimmer
 * \date 12/2015-09/2021
 */

namespace hqspre {

/**
 * \brief Constructs an empty formula without variables and clauses.
 */
inline Formula::Formula() :
    _prefix(nullptr),
    _dqbf_prefix(nullptr),
    _qbf_prefix(nullptr),
    _clauses(),
    _gates(),
    _occ_list(),
    _implications(),
    _deleted_clause_numbers(),
    _deleted_var_numbers(),
    _unit_stack(),
    _assignment(),
    _dont_touch(),
    _variable_score(),
    _clause_sizes(),
    _candidates(_variable_score),
    _blocked_candidates(_clause_sizes),
    _removed_lits(),
    _seen(),
    _seen2(),
    _settings(),
    _elim_prob(-1),
    _process_limit(),
    _statistics(),
    _timers()
{
    _qbf_prefix = new QBFPrefix();
    _prefix = _qbf_prefix;
}

/**
 * \brief Moves a formula into another new formula object.
 */
inline Formula::Formula(Formula&& other) noexcept :
    _prefix(other._prefix),
    _dqbf_prefix(other._dqbf_prefix),
    _qbf_prefix(other._qbf_prefix),
    _clauses(std::move(other._clauses)),
    _gates(std::move(other._gates)),
    _occ_list(std::move(other._occ_list)),
    _implications(std::move(other._implications)),
    _implications_added(other._implications_added),
    _deleted_clause_numbers(std::move(other._deleted_clause_numbers)),
    _deleted_var_numbers(std::move(other._deleted_var_numbers)),
    _unit_stack(std::move(other._unit_stack)),
    _assignment(std::move(other._assignment)),
    _dont_touch(std::move(other._dont_touch)),
    _variable_score(std::move(other._variable_score)),
    _clause_sizes(std::move(other._clause_sizes)),
    _candidates(_variable_score, std::move(other._candidates)),
    _blocked_candidates(_clause_sizes, std::move(other._blocked_candidates)),
    _removed_lits(std::move(other._removed_lits)),
    _seen(std::move(other._seen)),
    _seen2(std::move(other._seen2)),
    _settings(std::move(other._settings)),
    _elim_prob(std::move(other._elim_prob)),
    _process_limit(std::move(other._process_limit)),
    _interrupt(other._interrupt),
    _statistics(std::move(other._statistics)),
    _timers(std::move(other._timers))
{
    other._prefix = nullptr;
    other._dqbf_prefix = nullptr;
    other._qbf_prefix = nullptr;
}

/**
 * \brief Creates a copy of a formula.
 */
inline Formula::Formula(const Formula& other) :
    _prefix(nullptr),
    _dqbf_prefix(nullptr),
    _qbf_prefix(nullptr),
    _clauses(other._clauses),
    _gates(other._gates),
    _occ_list(other._occ_list),
    _implications(other._implications),
    _implications_added(other._implications_added),
    _deleted_clause_numbers(other._deleted_clause_numbers),
    _deleted_var_numbers(other._deleted_var_numbers),
    _unit_stack(other._unit_stack),
    _assignment(other._assignment),
    _dont_touch(other._dont_touch),
    _variable_score(other._variable_score),
    _clause_sizes(other._clause_sizes),
    _candidates(_variable_score, other._candidates),
    _blocked_candidates(_clause_sizes, other._blocked_candidates),
    _removed_lits(other._removed_lits),
    _seen(other._seen),
    _seen2(other._seen2),
    _settings(other._settings),
    _elim_prob(other._elim_prob),
    _process_limit(other._process_limit),
    _interrupt(other._interrupt),
    _statistics(other._statistics),
    _timers(other._timers)
{
    if (other._dqbf_prefix) {
        _dqbf_prefix = new DQBFPrefix(*(other._dqbf_prefix));
        _prefix = _dqbf_prefix;
    } else if (other._qbf_prefix) {
        _qbf_prefix = new QBFPrefix(*(other._qbf_prefix));
        _prefix = _qbf_prefix;
    }
}

inline Formula& Formula::operator=(const Formula& other)
{
    if (&other == this) return *this;

    if (_prefix) {
        delete _prefix;
        _prefix = nullptr;
        _dqbf_prefix = nullptr;
        _qbf_prefix = nullptr;
    }

    if (other._dqbf_prefix) {
        _dqbf_prefix = new DQBFPrefix(*(other._dqbf_prefix));
        _prefix      = _dqbf_prefix;
    } else if (other._qbf_prefix) {
        _qbf_prefix = new QBFPrefix(*(other._qbf_prefix));
        _prefix     = _qbf_prefix;
    }

    _enforce_dqbf = other._enforce_dqbf;
    _enforce_stochastic = other._enforce_stochastic;
    _clauses = other._clauses;
    _gates = other._gates;
    _occ_list = other._occ_list;
    _implications = other._implications;
    _implications_added = other._implications_added;
    _deleted_clause_numbers = other._deleted_clause_numbers;
    _deleted_var_numbers = other._deleted_var_numbers;
    _unit_stack = other._unit_stack;
    _assignment = other._assignment;
    _dont_touch = other._dont_touch;
    _variable_score = other._variable_score;
    _clause_sizes = other._clause_sizes;
    _candidates.copy(_variable_score, other._candidates);
    _blocked_candidates.copy(_clause_sizes, other._blocked_candidates);
    _removed_lits = other._removed_lits;
    _seen = other._seen;
    _seen2 = other._seen2;
    _settings = other._settings;
    _elim_prob = other._elim_prob,
    _process_limit = other._process_limit;
    _interrupt = other._interrupt;
    _statistics = other._statistics;
    _timers = other._timers;

    belongIntoOneClause = other.belongIntoOneClause;
    max_resolveRmb_cost = other.max_resolveRmb_cost;
    nodes = other.nodes;

    return *this;
}


inline Formula& Formula::operator=(Formula&& other)
{
    if (&other == this) return *this;

    if (_prefix) {
        delete _prefix;
        _prefix = nullptr;
        _dqbf_prefix = nullptr;
        _qbf_prefix = nullptr;
    }
    _prefix = other._prefix;
    other._prefix = nullptr;

    _qbf_prefix = other._qbf_prefix;
    other._qbf_prefix = nullptr;

    _dqbf_prefix = other._dqbf_prefix;
    other._dqbf_prefix = nullptr;

    _enforce_dqbf = other._enforce_dqbf;
    _enforce_stochastic = other._enforce_stochastic;
    _clauses = std::move(other._clauses);
    _gates = std::move(other._gates);
    _occ_list = std::move(other._occ_list);
    _implications = std::move(other._implications);
    _implications_added = other._implications_added;
    _deleted_clause_numbers = std::move(other._deleted_clause_numbers);
    _deleted_var_numbers = std::move(other._deleted_var_numbers);
    _unit_stack = std::move(other._unit_stack);
    _assignment = std::move(other._assignment);
    _dont_touch = std::move(other._dont_touch);
    _variable_score = std::move(other._variable_score);
    _clause_sizes = std::move(other._clause_sizes);
    _candidates.move(_variable_score, std::move(other._candidates));
    _blocked_candidates.move(_clause_sizes, std::move(other._blocked_candidates));
    _removed_lits = std::move(other._removed_lits);
    _seen = std::move(other._seen);
    _seen2 = std::move(other._seen2);
    _settings = std::move(other._settings);
    _elim_prob = std::move(other._elim_prob),
    _process_limit = std::move(other._process_limit);
    _interrupt = other._interrupt;
    _statistics = std::move(other._statistics);
    _timers = std::move(other._timers);

    belongIntoOneClause = std::move(other.belongIntoOneClause);
    max_resolveRmb_cost = std::move(other.max_resolveRmb_cost);
    nodes = std::move(other.nodes);

    return *this;
}


/**
 * \brief Destroys the formula and frees the prefix
 */
inline Formula::~Formula() noexcept
{
//    val_assert(_prefix);
    delete _prefix;
}

/**
 * \brief Adds a clause to the formula.
 *
 * \note Calling this function is only allowed as long as none of
 *       the functions has been called which modify the formula
 *       (like unitPropagation(), findPure(), ...)
 * \param clause vector with the literals of the clause
 * \param needs_sorting if true, the clause is sorted and duplicate literals are
 * removed. \pre If needs_sorting is false, the caller has to guarantee that the
 * literals are ordered in increasing order and that no literal appears more
 * than once. \param check_subsumption if true, it is checked whether the clause
 * subsumes other clauses in database \param status Specifies whether the clause
 * is manadory or optional (deleted does not make much sense here!) \return the
 * clause ID if the clause was actually added; -2 if the clause was a tautology;
 * -1 if the clause was unit. \throw UNSATException if an empty clause is added.
 */
template <typename Container>
inline int
Formula::addClause(const Container& clause, bool needs_sorting, bool check_subsumption, ClauseStatus status)
{
    if (clause.empty()) {
        throw UNSATException("Trying to add an empty clause");
    }

    if (clause.size() == 1) {
        pushUnit(clause[0], PureStatus::UNIT);
        return -1;
    }

    ClauseID c_nr = 0;

    // Try to recycle clause numbers where possible.
    if (_deleted_clause_numbers.empty()) {
        _clauses.emplace_back(clause, needs_sorting, status);
        c_nr = _clauses.size() - 1;
    } else {
        c_nr = _deleted_clause_numbers.back();
        _deleted_clause_numbers.pop_back();
        _clauses[c_nr] = Clause(clause, needs_sorting, status);
    }

    return addClauseToLists(c_nr, check_subsumption);
}

/**
 * \brief Returns the number of active variables.
 *
 * A variable is active if it is contained in at least one clause.
 * \sa Formula::numUVars()
 * \sa Formula::numEVars()
 */
inline std::size_t
Formula::numVars() const noexcept
{
    val_assert(_prefix);
    return _prefix->numVars();
}

/**
 * \brief Returns the number of active universal variables.
 *
 * A variable is active if it has not been deleted.
 * \sa Formula::numVars()
 * \sa Formula::numEVars()
 */
inline std::size_t
Formula::numUVars() const noexcept
{
    val_assert(_prefix);
    return _prefix->numUVars();
}

/**
 * \brief Returns the number of active existential variables.
 *
 * A variable is active if it has not been deleted.
 * \sa Formula::numVars()
 * \sa Formula::numUVars()
 */
inline std::size_t
Formula::numEVars() const noexcept
{
    val_assert(_prefix);
    return _prefix->numEVars();
}

/**
 * \brief Returns the number of literals in the formula.
 *
 * \note The running time is linear in the size of the formula.
 */
inline std::size_t
Formula::numLiterals() const noexcept
{
    std::size_t result = 0u;
    const std::size_t clauses_size = _clauses.size();
    for (ClauseID c_nr = 0u; c_nr < clauses_size; ++c_nr) {
        if (!clauseDeleted(c_nr)) {
            result += _clauses[c_nr].size();
        }
    }
    return result;
}

/**
 * \brief Returns if a given variable is universal.
 *
 * A variable has one of three states: universal, existential, or deleted.
 * That means: deleted variables are neither universal nor existential.
 * \return true iff the variable is universal.
 * \sa Formula::isExistential(Variable)
 * \sa Formaul::varDeleted(Variable)
 */
inline bool
Formula::isUniversal(Variable var) const noexcept
{
    val_assert(_prefix);
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    return _prefix->isUniversal(var);
}

/**
 * \brief Returns if a given variable is existential.
 *
 * A variable has one of three states: universal, existential, or deleted.
 * That means deleted variables are neither universal nor existential.
 * \return true iff the variable is existential.
 * \sa Formula::isUniversal(Variable)
 * \sa Formaul::varDeleted(Variable)
 */
inline bool
Formula::isExistential(Variable var) const noexcept
{
    val_assert(_prefix);
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    return _prefix->isExistential(var);
}

/**
 * \brief Returns if a variable has been deleted.
 *
 * Variables are considered to be deleted if they do not
 * appear in the formula anymore.
 * \sa Formula::isUniversal(Variable)
 * \sa Formula::isExistential(Variable)
 */
inline bool
Formula::varDeleted(const Variable var) const noexcept
{
    val_assert(_prefix);
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    return _prefix->varDeleted(var);
}

/**
 * \brief Returns the maximal valid ID of a variable.
 *
 * Variable IDs always start at 1. Some of the IDs may be
 * inactive because the variables does not appear in the formula
 * anymore. Note that this function does not return the number
 * of active variables.
 * \sa Formula::numVars()
 * \sa Formula::maxLitIndex()
 */
inline Variable
Formula::maxVarIndex() const noexcept
{
    val_assert(_prefix);
    return _prefix->maxVarIndex();
}

/**
 * \brief Returns the minimal valid ID of a variable.
 *
 * Variable IDs always start at 1. Some of the IDs may be
 * inactive because the variables does not appear in the formula
 * anymore. Note that this function does not return the number
 * of active variables.
 * \sa Formula::numVars()
 * \sa Formula::maxVarIndex()
 * \sa Formula::maxLitIndex()
 * \sa Formula::minLitIndex()
 */
inline constexpr Variable
Formula::minVarIndex() noexcept
{
    return Prefix::minVarIndex();
}

/**
 * \brief Returns the minimal valid ID of a literal.
 *
 * Literal IDs always start at 2. Some of the IDs may be
 * inactive because the corresponding variables has been deleted.
 * Note that this function does not return the number
 * of active literals.
 * \sa Formula::numVars()
 * \sa Formula::maxVarIndex()
 * \sa Formula::minVarIndex()
 * \sa Formula::maxLitIndex()
 */
inline constexpr Literal
Formula::minLitIndex() noexcept
{
    return Prefix::minLitIndex();
}

/**
 * \brief Returns the maximal valid ID of a literal.
 *
 * Literal IDs always start at 2. Some of the IDs may be
 * inactive because the corresponding variables has been deleted.
 * Note that this function does not return the number
 * of active literals.
 * \sa Formula::numVars()
 * \sa Formula::maxVarIndex()
 */
inline Literal
Formula::maxLitIndex() const noexcept
{
    val_assert(_prefix);
    return _prefix->maxLitIndex();
}

/**
 * \brief Set the "don't touch" status of a variable to the given status
 *
 * Variables can be "don't touch" in context of incremental solving
 */
inline void
Formula::setDontTouch(const Variable var, const bool status) noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    _dont_touch[var] = status;
}

/**
 * \brief Returns if a variable is "don't touch".
 *
 * Variables can be "don't touch" in context of incremental solving
 */
inline bool
Formula::isDontTouch(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    return _dont_touch[var];
}

/**
 * \brief Removes a variable.
 *
 * \pre The variable may not occur in any clause anymore.
 */
inline void
Formula::removeVar(const Variable var)
{
    val_assert(_prefix);
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(!varDeleted(var));
    val_assert(_occ_list[var2lit(var, false)].empty());
    val_assert(_occ_list[var2lit(var, true)].empty());
    val_assert(_implications[var2lit(var, false)].empty());
    val_assert(_implications[var2lit(var, true)].empty());

    _prefix->removeVar(var);
    _deleted_var_numbers.push_back(var);
    _assignment[var] = TruthValue::UNKNOWN;
}

/**
 * \brief Removes a clause from the occurrence list of a literal
 */
inline void
Formula::removeFromOccList(Literal lit, ClauseID c_nr)
{
    val_assert(minLitIndex() <= lit && lit <= maxLitIndex());

    const auto position = std::find(_occ_list[lit].begin(), _occ_list[lit].end(), c_nr);
    if (position != _occ_list[lit].end()) {
        *position = _occ_list[lit].back();
        _occ_list[lit].pop_back();
    }
}

/**
 * \brief Returns the literal from the given clause with the shortest occurrence
 * list
 */
template <typename Container>
inline Literal
Formula::getMinOccLit(const Container& clause) const
{
    val_assert(!clause.empty());

    return *std::min_element(clause.cbegin(), clause.cend(),
        [this](const Literal a, const Literal b) -> bool { return _occ_list[a].size() < _occ_list[b].size(); });
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * The new variable is created in the right-most block, i.e., it
 * depends on all universal variables.
 */
inline Variable
Formula::addEVar()
{
    return setEVar(nextVar());
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * The new variable is created in the right-most block, i.e., it
 * depends on all universal variables.
 */
inline Variable
Formula::setEVar(Variable index)
{
    val_assert(_prefix);
    val_assert(_prefix->varDeleted(index));
    val_assert(index > 0);

    _prefix->addEVar(index);
    return index;
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * If the current prefix is a QBF prefix, it is automatically
 * converted to a DQBF prefix!
 * \param deps container with the dependency set of the created variable
 * \pre All variables in the dependency set must be universal (and be created
 * before).
 */
template <typename Container>
inline Variable
Formula::addEVar(const Container& deps)
{
    return this->setEVar(nextVar(), deps);
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * If the current prefix is a QBF prefix, it is automatically
 * converted to a DQBF prefix.
 * \param deps container with the dependency set of the created variable
 * \pre All variables in the dependency set must be universal (and be created
 * before).
 */
template <typename Container>
inline Variable
Formula::setEVar(Variable index, const Container& deps)
{
    val_assert(_prefix);
    val_assert(_prefix->varDeleted(index));
    val_assert(index > 0);

    if (_prefix->type() == PrefixType::QBF) {
        _dqbf_prefix = _qbf_prefix->convertToDQBF();
        delete _prefix;
        _prefix = _dqbf_prefix;
        _qbf_prefix = nullptr;
    }

    if (index > maxVarIndex()) setMaxVarIndex(index + 10);

    _dqbf_prefix->addEVar(index, deps);
    return index;
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * If the current prefix is a QBF prefix, it is automatically
 * converted to a DQBF prefix!
 * \param deps container with the dependency set of the created variable
 * \pre All variables in the dependency set must be universal (and be created
 * before).
 */
inline Literal
Formula::addEVar(std::set<Variable>&& deps)
{
    return setEVar(nextVar(), std::move(deps));
}

/**
 * \brief Creates a new existential variable and returns its ID.
 *
 * If the current prefix is a QBF prefix, it is automatically
 * converted to a DQBF prefix!
 * This function creates a variable with a predefined index. The caller
 * has to ensure that this variable does not exist yet. Additionally,
 * calling this function invalidates the list of unused variables
 * (Formula::_deleted_var_numbers); it needs to be re-created before
 * the next function (apart from Formula::addUVar(Variable) is called.
 *
 * \param index the index of the variable to be created.
 * \param deps container with the dependency set of the created variable
 * \pre All variables in the dependency set must be universal (and be created
 * before).
 */
inline Literal
Formula::setEVar(Variable index, std::set<Variable>&& deps)
{
    val_assert(_prefix);
    val_assert(_prefix->varDeleted(index));
    val_assert(index > 0);

    if (_prefix->type() == PrefixType::QBF) {
        _dqbf_prefix = _qbf_prefix->convertToDQBF();
        delete _prefix;
        _prefix = _dqbf_prefix;
        _qbf_prefix = nullptr;
    }

    if (index > maxVarIndex()) setMaxVarIndex(index + 10);

    _dqbf_prefix->addEVar(index, std::move(deps));
    return index;
}

/**
 * \brief Creates a new universal variable and returns its ID.
 */
inline Variable
Formula::addUVar()
{
    return setUVar(nextVar());
}

/**
 * \brief Creates a new universal variable and returns its ID.
 *
 * This function creates a variable with a predefined index. The caller
 * has to ensure that this variable does not exist yet. Additionally,
 * calling this function invalidates the list of unused variables
 * (Formula::_deleted_var_numbers); it needs to be re-created before
 * the next function (apart from Formula::addUVar(Variable) is called.
 *
 * \param index the index of the variable to be created.
 */
inline Variable
Formula::setUVar(Variable index)
{
    val_assert(_prefix);
    val_assert(_prefix->varDeleted(index));
    val_assert(index > 0);

    if (index > maxVarIndex()) setMaxVarIndex(index + 10);

    _prefix->addUVar(index);
    return index;
}

/**
 * \brief Returns the next available variable index.
 *
 * The returned variable index is marked as 'not available'.
 * \note This is the only function which may remove elements from
 *       `_deleted_var_numbers`.
 * \sa Formula::addEVar()
 * \sa Formula::addUVar()
 * \sa Formula::copyVar(Variable)
 * \sa Formula::removeVar(Variable)
 */
inline Variable
Formula::nextVar()
{
    if (_deleted_var_numbers.empty()) {
        setMaxVarIndex(maxVarIndex() + 10);
        val_assert(_deleted_var_numbers.size() >= 10);
    }

    const Variable result = _deleted_var_numbers.back();
    _deleted_var_numbers.pop_back();

    val_assert(minVarIndex() <= result && result <= maxVarIndex());
    val_assert(varDeleted(result));

    return result;
}

/**
 * \brief Updates data structures for new variables
 *
 * Resizes all variable-related data structures such that
 * variables up to the given index can be used. The newly
 * created variable numbers are inserted into `_deleted_var_numbers`.
 *
 * \param index the largest valid variable ID
 * \sa Formula::nextVar()
 */
inline void
Formula::setMaxVarIndex(const Variable index)
{
    if (index <= maxVarIndex()) return;

    const unsigned int old_max_index = maxVarIndex();

    _prefix->setMaxVarIndex(index);
    _dont_touch.resize(index + 1, false);
    _occ_list.resize(2 * index + 2);
    _implications.resize(2 * index + 2);
    _assignment.resize(index + 1, TruthValue::UNKNOWN);
    _seen.resize(2 * index + 2, false);
    _seen2.resize(2 * index + 2, false);
    _gates.resizeVars(index + 1);

    for (Variable var = index; var > old_max_index; --var) {
        val_assert(minVarIndex() <= var && var <= maxVarIndex());
        val_assert(varDeleted(var));

        _deleted_var_numbers.push_back(var);
    }

    val_assert(_deleted_var_numbers.size() >= maxVarIndex() - old_max_index);
    val_assert(maxVarIndex() >= index);
}

/**
 * \brief Returns the quantifier level of a variable in a QBF.
 *
 * \pre The formula must be a QBF.
 *
 * The quantifier levels are numbered from left to right,
 * starting with 0 being the left-most (outermost) block.
 */
inline std::size_t
Formula::getLevel(Variable var) const noexcept
{
    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::QBF);
    val_assert(_qbf_prefix);
    val_assert(minVarIndex() <= var && var <= maxVarIndex());

    return _qbf_prefix->getLevel(var);
}

/**
 * \brief Returns the number of clauses.
 */
inline std::size_t
Formula::numClauses() const noexcept
{
    return (_clauses.size() + _unit_stack.size()) - _deleted_clause_numbers.size();
}

/**
 * \brief Returns the maximum used clause ID.
 *
 * Clause IDs are numbered starting from 0. Some of the IDs may no
 * longer be used because clauses have been deleted. This can be checked
 * using Formula::clauseDeleted().
 */
inline ClauseID
Formula::maxClauseIndex() const noexcept
{
    val_assert(!_clauses.empty());

    return static_cast<ClauseID>(_clauses.size() - 1);
}

/**
 * Tries to find the given clause in the clause database.
 * \return the clause ID if the clause was found and -1 otherwise.
 * \note Unit clauses cannot be found this way as they are only
 *       contained in the unit_stack list.
 * \pre `clause` must be sorted and may not contain duplicate literals.
 */
template <typename Container>
inline int
Formula::findClause(const Container& clause)
{
    val_assert(clause.size() > 1);

    if (clause.size() == 2) {
        return getImplicationClause(negate(clause[0]), clause[1]);
    }

    const Literal lit = getMinOccLit(clause);

    for (const ClauseID c_nr : _occ_list[lit]) {
        if (_clauses[c_nr].size() != clause.size()) continue;
        if (std::equal(_clauses[c_nr].cbegin(), _clauses[c_nr].cend(), clause.cbegin(), clause.cend())) {
            return static_cast<int>(c_nr);
        }
    }

    return -1;
}

inline bool
Formula::isGateClause(const ClauseID c_nr) const noexcept
{
    val_assert(c_nr <= maxClauseIndex());

    return _gates.isGateClause(c_nr);
}

/**
 * \brief Returns a reference to the clause with the given ID.
 *
 * \pre The accessed clause may not be deleted.
 * \sa Formula::clauseDeleted()
 */
inline const Clause&
Formula::getClause(const ClauseID c_nr) const noexcept
{
    val_assert(c_nr <= maxClauseIndex());
    val_assert(!clauseDeleted(c_nr));

    return _clauses[c_nr];
}

/**
 * \brief Checks if a clause has been deleted.
 */
inline bool
Formula::clauseDeleted(const ClauseID c_nr) const noexcept
{
    val_assert(c_nr <= maxClauseIndex());

    return _clauses[c_nr].getStatus() == ClauseStatus::DELETED;
}

/**
 * \brief Checks if a clause is optional and can be deleted.
 */
inline bool
Formula::clauseOptional(const ClauseID c_nr) const noexcept
{
    val_assert(c_nr <= maxClauseIndex());

    return _clauses[c_nr].getStatus() == ClauseStatus::OPTIONAL;
}

/**
 * \brief Checks if a clause is mandatory and cannot be deleted.
 */
inline bool
Formula::clauseMandatory(const ClauseID c_nr) const noexcept
{
    val_assert(c_nr <= maxClauseIndex());

    return _clauses[c_nr].getStatus() == ClauseStatus::MANDATORY;
}

/**
 * \brief Returns the number of dependencies of a given variable.
 */
inline std::size_t
Formula::numDependencies(const Variable var) const noexcept
{
    val_assert(_prefix);

    return _prefix->numDependencies(var);
}

/**
 * \brief Checks if one variable depends on another one.
 * \pre One of the variables needs to be universal, the other one existential.
 * \return true if the dependency between var1 and var2 exists.
 */
inline bool
Formula::depends(const Variable var1, const Variable var2) const
{
    val_assert(_prefix);
    return _prefix->depends(var1, var2);
}

/**
 * \brief Removes a dependency between two variables.
 * \pre One of the variables needs to be universal, the other one existential.
 */
inline void
Formula::removeDependency(const Variable var1, const Variable var2)
{
    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::DQBF);

    _dqbf_prefix->removeDependency(var1, var2);
}

/**
 * \brief Adds a dependency between two variables.
 * \pre One of the variables needs to be universal, the other one existential.
 */
inline void
Formula::addDependency(const Variable var1, const Variable var2)
{
    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::DQBF);

    _dqbf_prefix->addDependency(var1, var2);
}

/**
 * \brief Returns a reference to the dependency set of a given variable.
 */
inline const std::set<Variable>&
Formula::getDependencies(const Variable var) const
{
    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::DQBF);

    return _dqbf_prefix->getDependencies(var);
}

/**
 * \brief Adds an implication lit1 => lit2 to the implication data structure
 * \note The equivalent implication ~lit2 => ~lit1 is not added automatically.
 */
inline void
Formula::addImplication(const Literal lit1, const Literal lit2, const ClauseID c_nr)
{
    val_assert(lit1 < _implications.size());
    val_assert(lit2 < _implications.size());

    _implications[lit1].insert(BinaryClause(lit2, c_nr));
    _implications_added = true;
}

/**
 * \brief Checks if the formula contains an implication lit1 => lit2.
 * \note This does not take transitivity into account.
 * \return -1 if implication is not found, and ID of existing binary clause
 * otherwise
 */
inline int
Formula::hasImplication(const Literal lit1, const Literal lit2) const noexcept
{
    val_assert(lit1 >= minLitIndex());
    val_assert(lit2 >= minLitIndex());
    val_assert(lit1 < _implications.size());
    val_assert(lit2 < _implications.size());

    auto found = _implications[lit1].find(BinaryClause(lit2));
    if (found != _implications[lit1].cend()) {
        return static_cast<int>(found->getClauseID());
    } else {
        return -1;
    }
}

/**
 * \brief Returns the clause ID of an implication.
 *
 * If there is the implication lit1 => lit2, the ID of the corresponding
 * binary clause is returned. Otherwise the return value is negative.
 */
inline int
Formula::getImplicationClause(const Literal lit1, const Literal lit2) const noexcept
{
    val_assert(lit1 >= minLitIndex());
    val_assert(lit2 >= minLitIndex());
    val_assert(lit1 < _implications.size());
    val_assert(lit2 < _implications.size());

    const auto found = _implications[lit1].find(BinaryClause(lit2));
    if (found == _implications[lit1].cend())
        return -1;
    else
        return static_cast<int>(found->getClauseID());
}

/**
 * \brief Removes the implication lit1 => lit2 from the implication data
 * structure. \note This does neither remove the converse implication ~lit2 =>
 * ~lit1 nor the corresponding binary clause (~lit1, lit2).
 */
inline void
Formula::removeImplication(const Literal lit1, const Literal lit2)
{
    val_assert(lit1 >= minLitIndex());
    val_assert(lit2 >= minLitIndex());
    val_assert(lit1 < _implications.size());
    val_assert(lit2 < _implications.size());

    _implications[lit1].erase(BinaryClause(lit2));
}

/**
 * \brief Returns the assignment of a variable.
 *
 * Value 0 corresponds to unassigned, +1 to TRUE, and
 * -1 to FALSE.
 */
inline std::int32_t
Formula::varAssignment(const Variable var) const noexcept
{
    if (_assignment[var] == TruthValue::UNKNOWN) {
        return 0;
    } else if (_assignment[var] == TruthValue::TRUE) {
        return 1;
    } else if (_assignment[var] == TruthValue::FALSE) {
        return -1;
    } else {
        val_assert(false);
        exit(0);
    }
    return 0;
}

inline void
Formula::enforceDQBF(const bool val) noexcept
{
    _enforce_dqbf = val;

    if (_enforce_dqbf && _prefix->type() == PrefixType::QBF) {
        _dqbf_prefix = _qbf_prefix->convertToDQBF();
        delete _qbf_prefix;
        _qbf_prefix = nullptr;
        _prefix = _dqbf_prefix;
    }
}

inline void
Formula::setInterrupt(const bool val) noexcept
{
    _interrupt = val;
}

inline void
Formula::clearSeen() const
{
    std::fill(_seen.begin(), _seen.end(), false);
}

inline std::size_t&
Formula::stat(const Statistics which)
{
    return _statistics[static_cast<unsigned int>(which)];
}

template <typename Container>
inline void
Formula::clearSeen(const Container& container) const
{
    for (const auto lit : container) {
        _seen[lit] = false;
    }
}

template <typename Container>
inline void
Formula::setSeen(const Container& container) const
{
    for (const auto lit : container) {
        val_assert(!_seen[lit]);
        _seen[lit] = true;
    }
}

template <typename Container>
inline void
Formula::clearSeen2(const Container& container) const
{
    for (const auto lit : container) {
        _seen2[lit] = false;
    }
}

template <typename Container>
inline void
Formula::setSeen2(const Container& container) const
{
    for (const auto lit : container) {
        val_assert(!_seen2[lit]);
        _seen2[lit] = true;
    }
}

}  // end namespace hqspre

#endif
