/**************************************************************
 *
 *       AIGPP Package // Manager_assignment.cc
 *
 *	Author:
 *         Florian Pigorsch
 *	  University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *	Last revision:
 *         $Revision: 649 $
 *	  $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "Manager.hh"

std::pair<bool, aigpp::VariableAssignment>
aigpp::Manager::satisfyingAssignment(const aigpp::Edge& e, bool minimize)
{
    /*
     * handle constant nodes
     */
    if (e.isConstant()) {
        /* const0 -> not satisfiable */
        if (structurallyFalse(e)) {
            return std::pair<bool, VariableAssignment>(false, VariableAssignment());
        }
        /* const1 -> satisfiable for all assignments */
        else {
            return std::pair<bool, VariableAssignment>(true, VariableAssignment());
        }
    }

    /*
     * clear the SAT solver instance
     */
    resetSolver(true);

    /*
     * create CNF for e's cone
     */
    createClauses(e.node());

    /*
     * create "miter" variable
     */
    Minisat::Lit f = Minisat::mkLit(_solver->newVar());
    if (e.isInverted()) {
        /*
         * create clauses for "!( m == e.node()->satVar )"
         */
        Minisat::vec<Minisat::Lit> clause1(2), clause2(2);
        clause1[0] = ~f;
        clause1[1] = ~Minisat::mkLit(e.node()->_satVar);
        clause2[0] = f;
        clause2[1] = Minisat::mkLit(e.node()->_satVar);
        _solver->addClause(clause1);
        _solver->addClause(clause2);
    } else {
        /*
         * create clauses for "f == e.node()->satVar"
         */
        Minisat::vec<Minisat::Lit> clause1(2), clause2(2);
        clause1[0] = ~f;
        clause1[1] = Minisat::mkLit(e.node()->_satVar);
        clause2[0] = f;
        clause2[1] = ~Minisat::mkLit(e.node()->_satVar);
        _solver->addClause(clause1);
        _solver->addClause(clause2);
    }

    /*
     * get a full satisfying assignment
     */
    bool                       solverState;
    Minisat::vec<Minisat::Lit> initialAssumptions;
    initialAssumptions.push(f);

    solverState = _solver->solve(initialAssumptions);

    /* not satisfiable */
    if (!solverState) {
        /* must not be the case when e.node() is reduced! */
        assert(!(e.node()->isReduced()));

        resetSolver(true);

        return std::pair<bool, VariableAssignment>(false, VariableAssignment());
    }

    /* satisfiable */
    VariableAssignment satAssignment = getCounterExample();

    /*
     * return the full assignment if minimize is false
     */
    if (!minimize) {
        resetSolver(true);
        return std::pair<bool, VariableAssignment>(true, satAssignment);
    }

    for (std::vector<Node*>::const_iterator v = _variables.begin(); v != _variables.end(); ++v) {
        if ((*v)->isSatVarValid()) {
            if (satAssignment.get((*v)->varIndex()) == VariableAssignment::Unassigned) {
                continue;
            }

            Minisat::vec<Minisat::Lit> assumptions;
            assumptions.push(~f);
            for (std::vector<Node*>::const_iterator v2 = _variables.begin(); v2 != _variables.end(); ++v2) {
                if (v2 == v) {
                    if (satAssignment.get((*v2)->varIndex()) == VariableAssignment::Positive) {
                        assumptions.push(~(Minisat::mkLit((*v2)->_satVar)));
                    } else if (satAssignment.get((*v2)->varIndex()) == VariableAssignment::Negative) {
                        assumptions.push(Minisat::mkLit((*v2)->_satVar));
                    } else {
                        abort();
                    }
                } else if (satAssignment.get((*v2)->varIndex()) == VariableAssignment::Positive) {
                    assumptions.push(Minisat::mkLit((*v2)->_satVar));
                } else if (satAssignment.get((*v2)->varIndex()) == VariableAssignment::Negative) {
                    assumptions.push(~(Minisat::mkLit((*v2)->_satVar)));
                } else {
                    // nothing to do. variable is nor assigned.
                }
            }

            bool satisfiable = _solver->solve(assumptions);
            if (!satisfiable) {
                satAssignment.setUnassigned((*v)->varIndex());
                assert(satAssignment.get((*v)->varIndex()) == VariableAssignment::Unassigned);
            }
        }
    }

    resetSolver(true);

    return std::pair<bool, VariableAssignment>(true, satAssignment);
}
