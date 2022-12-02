/**************************************************************
 *
 *       AIGPP Package // Manager_sat.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
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

#include "Manager.hh"

/* std */
#include <fstream>
#include <sstream>

/* LRABS utilities */
#include <lrabsutil/Resources.hh>

/* local */
#include "StatManager.hh"

bool
aigpp::Manager::satEqual(aigpp::Node* n1, aigpp::Node* n2)
{
    assert(n1 != nullptr && n2 != nullptr);

    if (verbosity() >= 4) {
        std::cout << "equivalent? " << n1->index() << " " << n2->index() << std::endl;
    }

    const double startTime = lrabs::cpuTime();

    ++_satChecksPresent;

    StatManager::instance().incSATChecks();

    createClauses(n1, n2);

    // miter m <=> n1 xor n2
    Minisat::Lit m  = Minisat::mkLit(_solver->newVar());
    Minisat::Lit f1 = Minisat::mkLit(n1->_satVar);
    Minisat::Lit f2 = Minisat::mkLit(n2->_satVar);

    Minisat::vec<Minisat::Lit> clause1(3), clause2(3), clause3(3), clause4(3);
    clause1[0] = ~m;
    clause1[1] = f1;
    clause1[2] = f2;
    clause2[0] = ~m;
    clause2[1] = ~f1;
    clause2[2] = ~f2;
    clause3[0] = m;
    clause3[1] = ~f1;
    clause3[2] = f2;
    clause4[0] = m;
    clause4[1] = f1;
    clause4[2] = ~f2;
    _solver->addClause(clause1);
    _solver->addClause(clause2);
    _solver->addClause(clause3);
    _solver->addClause(clause4);

    Minisat::vec<Minisat::Lit> assumptions(1, m);

    ++_mitersInCNF;

    StatManager::instance().updateSATMaxVars(_solver->nVars());
    StatManager::instance().incSATCreationTime(lrabs::cpuTime() - startTime);

    const bool satisfiable = _solver->solve(assumptions);
    StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

    if (satisfiable) {
        addCounterExample(getCounterExample());
        ++_openMitersInCNF;

        if (verbosity() >= 4) {
            std::cout << "no" << std::endl;
        }
    } else {
        Minisat::vec<Minisat::Lit> unit(1, ~m);
        _solver->addClause(unit);

        ++_closedMitersInCNF;

        StatManager::instance().incSATEquiv();

        if (verbosity() >= 4) {
            std::cout << "yes" << std::endl;
        }
    }

    return !satisfiable;
}

bool
aigpp::Manager::satEqualNegated(aigpp::Node* n1, aigpp::Node* n2)
{
    assert(n1 != nullptr && n2 != nullptr);

    // resetSolver();

    if (verbosity() >= 4) {
        std::cout << "antivalent? " << n1->index() << " " << n2->index() << std::endl;
    }

    const double startTime = lrabs::cpuTime();

    StatManager::instance().incSATChecks();
    ++_satChecksPresent;

    createClauses(n1, n2);

    // miter m <=> ~n1 xor n2
    Minisat::Lit m  = Minisat::mkLit(_solver->newVar());
    Minisat::Lit f1 = Minisat::mkLit(n1->_satVar);
    Minisat::Lit f2 = Minisat::mkLit(n2->_satVar);

    Minisat::vec<Minisat::Lit> clause1(3), clause2(3), clause3(3), clause4(3);
    clause1[0] = ~m;
    clause1[1] = ~f1;
    clause1[2] = f2;
    clause2[0] = ~m;
    clause2[1] = f1;
    clause2[2] = ~f2;
    clause3[0] = m;
    clause3[1] = f1;
    clause3[2] = f2;
    clause4[0] = m;
    clause4[1] = ~f1;
    clause4[2] = ~f2;
    _solver->addClause(clause1);
    _solver->addClause(clause2);
    _solver->addClause(clause3);
    _solver->addClause(clause4);

    Minisat::vec<Minisat::Lit> assumptions(1, m);

    ++_mitersInCNF;

    StatManager::instance().updateSATMaxVars(_solver->nVars());
    StatManager::instance().incSATCreationTime(lrabs::cpuTime() - startTime);

    if (_solver->solve(assumptions)) {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        addCounterExample(getCounterExample());

        ++_openMitersInCNF;

        if (verbosity() >= 4) {
            std::cout << "no" << std::endl;
        }

        return false;
    } else {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        Minisat::vec<Minisat::Lit> unit(1, ~m);
        _solver->addClause(unit);

        ++_closedMitersInCNF;

        StatManager::instance().incSATEquiv();

        if (verbosity() >= 4) {
            std::cout << "yes" << std::endl;
        }

        return true;
    }
}

bool
aigpp::Manager::satEqualZero(aigpp::Node* n1)
{
    assert(n1 != nullptr);

    if (verbosity() >= 4) {
        std::cout << "zero? " << n1->index() << std::endl;
    }

    const double startTime = lrabs::cpuTime();

    ++_satChecksPresent;
    StatManager::instance().incSATChecks();

    createClauses(n1);

    Minisat::Lit               f1 = Minisat::mkLit(n1->_satVar);
    Minisat::vec<Minisat::Lit> assumptions(1, f1);

    StatManager::instance().updateSATMaxVars(_solver->nVars());
    StatManager::instance().incSATCreationTime(lrabs::cpuTime() - startTime);

    if (_solver->solve(assumptions)) {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        addCounterExample(getCounterExample());

        if (verbosity() >= 4) {
            std::cout << "no" << std::endl;
        }

        return false;
    } else {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        Minisat::vec<Minisat::Lit> unit(1, ~f1);
        _solver->addClause(unit);

        StatManager::instance().incSATEquiv();

        if (verbosity() >= 4) {
            std::cout << "yes" << std::endl;
        }

        return true;
    }
}

bool
aigpp::Manager::satEqualOne(aigpp::Node* n1)
{
    assert(n1 != nullptr);

    if (verbosity() >= 4) {
        std::cout << "one? " << n1->index() << std::endl;
    }

    const double startTime = lrabs::cpuTime();

    ++_satChecksPresent;
    StatManager::instance().incSATChecks();

    createClauses(n1);

    Minisat::Lit               f1 = Minisat::mkLit(n1->_satVar);
    Minisat::vec<Minisat::Lit> assumptions(1, ~f1);

    StatManager::instance().updateSATMaxVars(_solver->nVars());
    StatManager::instance().incSATCreationTime(lrabs::cpuTime() - startTime);

    if (_solver->solve(assumptions)) {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        addCounterExample(getCounterExample());

        if (verbosity() >= 4) {
            std::cout << "no" << std::endl;
        }

        return false;
    } else {
        StatManager::instance().incSATTime(lrabs::cpuTime() - startTime);

        Minisat::vec<Minisat::Lit> unit(1, f1);
        _solver->addClause(unit);

        StatManager::instance().incSATEquiv();

        if (verbosity() >= 4) {
            std::cout << "yes" << std::endl;
        }

        return true;
    }
}

bool
aigpp::Manager::isBoolRedundant(aigpp::Node* node, aigpp::Node* var) const
{
    Minisat::Solver solver;

    std::stack<const Node*>             pending;
    std::map<const Node*, Minisat::Var> varmap1, varmap2;

    varmap1[var] = solver.newVar();
    varmap2[var] = solver.newVar();

    pending.push(nullptr);

    const Node* n = node;

    do {
        if (varmap1.find(n) != varmap1.end()) {
            n = pending.top();
            pending.pop();
        } else if (n->isVar()) {
            Minisat::Var v = solver.newVar();
            varmap1[n]     = v;
            varmap2[n]     = v;

            n = pending.top();
            pending.pop();
        } else if (varmap1.find(n->parent1().node()) == varmap1.end()) {
            pending.push(n);
            n = n->parent1().node();
        } else if (varmap1.find(n->parent2().node()) == varmap1.end()) {
            pending.push(n);
            n = n->parent2().node();
        } else {
            // n = (!)n1 & (!)n2  <=>  (!n | (!)n1) & (!n | (!)n2) & (n | !(!)n1 |
            // !(!)n2)

            Minisat::Var v1 = solver.newVar();
            Minisat::Var v2 = solver.newVar();
            varmap1[n]      = v1;
            varmap2[n]      = v2;

            {
                Minisat::Lit vn = Minisat::mkLit(v1, n->isNAND());

                Minisat::Var vn1 = varmap1[n->parent1().node()];
                Minisat::Var vn2 = varmap1[n->parent2().node()];

                Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

                clause1[0] = ~vn;
                clause1[1] = Minisat::mkLit(vn1, n->parent1().isInverted());
                clause2[0] = ~vn;
                clause2[1] = Minisat::mkLit(vn2, n->parent2().isInverted());
                clause3[0] = vn;
                clause3[1] = Minisat::mkLit(vn1, !(n->parent1().isInverted()));
                clause3[2] = Minisat::mkLit(vn2, !(n->parent2().isInverted()));

                solver.addClause(clause1);
                solver.addClause(clause2);
                solver.addClause(clause3);
            }

            {
                Minisat::Lit vn = Minisat::mkLit(v2, n->isNAND());

                Minisat::Var vn1 = varmap2[n->parent1().node()];
                Minisat::Var vn2 = varmap2[n->parent2().node()];

                Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

                clause1[0] = ~vn;
                clause1[1] = Minisat::mkLit(vn1, n->parent1().isInverted());
                clause2[0] = ~vn;
                clause2[1] = Minisat::mkLit(vn2, n->parent2().isInverted());
                clause3[0] = vn;
                clause3[1] = Minisat::mkLit(vn1, !(n->parent1().isInverted()));
                clause3[2] = Minisat::mkLit(vn2, !(n->parent2().isInverted()));

                solver.addClause(clause1);
                solver.addClause(clause2);
                solver.addClause(clause3);
            }

            n = pending.top();
            pending.pop();
        }
    } while (!pending.empty());

    Minisat::Lit               f1 = Minisat::mkLit(varmap1[node]);
    Minisat::Lit               f2 = Minisat::mkLit(varmap2[node]);
    Minisat::vec<Minisat::Lit> clause1(2), clause2(2);
    clause1[0] = f1;
    clause1[1] = f2;
    clause2[0] = ~f1;
    clause2[1] = ~f2;
    solver.addClause(clause1);
    solver.addClause(clause2);

    Minisat::vec<Minisat::Lit> clause3(1);
    clause3[0] = ~Minisat::mkLit(varmap2[var]);
    solver.addClause(clause3);

    const bool result = solver.solve();

    return !result;
}

std::vector<aigpp::Node*>
aigpp::Manager::boolRedundant(aigpp::Node* node, const std::vector<aigpp::Node*>& candidates) const
{
    Minisat::Solver solver;

    std::stack<const Node*> pending;
    using VarPair = std::pair<Minisat::Var, Minisat::Var>;
    std::map<const Node*, VarPair> varmap;

    pending.push(nullptr);

    const Node* n = node;

    do {
        if (varmap.find(n) != varmap.end()) {
            n = pending.top();
            pending.pop();
        } else if (n->isVar()) {
            varmap[n] = VarPair(solver.newVar(), solver.newVar());

            n = pending.top();
            pending.pop();
        } else if (varmap.find(n->parent1().node()) == varmap.end()) {
            pending.push(n);
            n = n->parent1().node();
        } else if (varmap.find(n->parent2().node()) == varmap.end()) {
            pending.push(n);
            n = n->parent2().node();
        } else {
            // n = (!)n1 & (!)n2  <=>  (!n | (!)n1) & (!n | (!)n2) & (n | !(!)n1 |
            // !(!)n2)

            Minisat::Var v1 = solver.newVar();
            Minisat::Var v2 = solver.newVar();
            varmap[n]       = VarPair(v1, v2);

            {
                Minisat::Lit vn = Minisat::mkLit(v1, n->isNAND());

                Minisat::Var vn1 = varmap[n->parent1().node()].first;
                Minisat::Var vn2 = varmap[n->parent2().node()].first;

                Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

                clause1[0] = ~vn;
                clause1[1] = Minisat::mkLit(vn1, n->parent1().isInverted());
                clause2[0] = ~vn;
                clause2[1] = Minisat::mkLit(vn2, n->parent2().isInverted());
                clause3[0] = vn;
                clause3[1] = Minisat::mkLit(vn1, !(n->parent1().isInverted()));
                clause3[2] = Minisat::mkLit(vn2, !(n->parent2().isInverted()));

                solver.addClause(clause1);
                solver.addClause(clause2);
                solver.addClause(clause3);
            }

            {
                Minisat::Lit vn = Minisat::mkLit(v2, n->isNAND());

                Minisat::Var vn1 = varmap[n->parent1().node()].second;
                Minisat::Var vn2 = varmap[n->parent2().node()].second;

                Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

                clause1[0] = ~vn;
                clause1[1] = Minisat::mkLit(vn1, n->parent1().isInverted());
                clause2[0] = ~vn;
                clause2[1] = Minisat::mkLit(vn2, n->parent2().isInverted());
                clause3[0] = vn;
                clause3[1] = Minisat::mkLit(vn1, !(n->parent1().isInverted()));
                clause3[2] = Minisat::mkLit(vn2, !(n->parent2().isInverted()));

                solver.addClause(clause1);
                solver.addClause(clause2);
                solver.addClause(clause3);
            }

            n = pending.top();
            pending.pop();
        }
    } while (!pending.empty());

    Minisat::Lit               f1 = Minisat::mkLit(varmap[node].first);
    Minisat::Lit               f2 = Minisat::mkLit(varmap[node].second);
    Minisat::vec<Minisat::Lit> clause1(2), clause2(2);
    clause1[0] = f1;
    clause1[1] = f2;
    clause2[0] = ~f1;
    clause2[1] = ~f2;
    solver.addClause(clause1);
    solver.addClause(clause2);

    const auto                 candidates_size_int = static_cast<int>(candidates.size());  // minisat works with ints!
    Minisat::vec<Minisat::Lit> miters(candidates_size_int);
    for (int i = 0; i < candidates_size_int; ++i) {
        Minisat::vec<Minisat::Lit> clause1(3), clause2(3), clause3(3), clause4(3);

        Minisat::Lit m  = Minisat::mkLit(solver.newVar());
        miters[i]       = m;
        Minisat::Lit v1 = Minisat::mkLit(varmap[candidates[i]].first);
        Minisat::Lit v2 = Minisat::mkLit(varmap[candidates[i]].second);

        clause1[0] = ~m;
        clause1[1] = v1;
        clause1[2] = v2;
        clause2[0] = ~m;
        clause2[1] = ~v1;
        clause2[2] = ~v2;
        clause3[0] = m;
        clause3[1] = ~v1;
        clause3[2] = v2;
        clause4[0] = m;
        clause4[1] = v1;
        clause4[2] = ~v2;

        solver.addClause(clause1);
        solver.addClause(clause2);
        solver.addClause(clause3);
        solver.addClause(clause4);
    }

    std::vector<Node*> redundant;

    Minisat::vec<Minisat::Lit> assumptions(candidates_size_int);

    for (int i = 0; i < candidates_size_int; ++i) {
        for (int j = 0; j < candidates_size_int; ++j) {
            assumptions[j] = ~miters[j];
        }
        assumptions[i] = miters[i];

        if (!(solver.solve(assumptions))) {
            redundant.push_back(candidates[i]);
        }
    }

    return redundant;
}

void
aigpp::Manager::resetSolver(bool force)
{
    assert(_solver != 0);

    /* no reset if solver is empty! */
    if ((_satChecksPresent == 0 || _noReset || _satChecksPresent < settings().getSATResetInterval()) && !force) {
        return;
    }

    StatManager::instance().incSATResets();

    delete _solver;
    _solver = new Minisat::Solver;

    for (Node* n : _nodesWithSATVars) {
        assert(n->isSatVarValid());
        n->invalidateSatVar();
    }
    _nodesWithSATVars.clear();

    _satChecksPresent     = 0;
    _nodesInCNF           = 0;
    _deletedNodesInCNF    = 0;
    _equivalentNodesInCNF = 0;
    _mitersInCNF          = 0;
    _openMitersInCNF      = 0;
    _closedMitersInCNF    = 0;
}

void
aigpp::Manager::createClauses(aigpp::Node* n1)
{
    createClauses(n1, n1);
}

#if 0
/* version that detects multi-ands and xors */
void
aigpp::Manager::
createClauses( aigpp::Node* n1, aigpp::Node* n2 )
{
    lockFlags<1>();

    std::stack<Node*>& pending = _universalNodeStack;
    assert( pending.empty() );

    pending.push( n1 );
    if( n1 != n2 )
    {
        pending.push( n2 );
    }

    Node* pt = 0;
    bool haveXOR = false;
    std::vector<Edge> inputs;
    Minisat::vec<Minisat::Lit> naryClause;
    Minisat::vec<Minisat::Lit> binClause;

    while( !pending.empty() )
    {
        pt = pending.top();

        if( pt->isSatVarValid() )
        {
            pending.pop();
            continue;
        }
        else if( pt->isVar() )
        {
            pt->setSatVar( _solver->newVar() );

            _nodesWithSATVars.insert( pt );

            pending.pop();
            continue;
        }

        if( ( haveXOR = collectXORInputs( pt, inputs ) ) )
        {
            assert( inputs.size() == 2 );
        }
        else
        {
            collectMultiAndInputs( pt, inputs );
        }


        for( const Edge& p : inputs)
        {
            if( !p.node()->isSatVarValid() )
            {
                pending.push( p.node() );
            }
        }

        if( pt != pending.top() )
        {
            continue;
        }

        pt->setSatVar( _solver->newVar() );
        _nodesWithSATVars.insert( pt );

        Minisat::Lit lit_pt = Minisat::mkLit( pt->_satVar, pt->isNAND() );
        ++_nodesInCNF;

        if( !haveXOR )
        {
            naryClause.clear();
            naryClause.push( lit_pt );

            for( std::size_t i = 0; i != inputs.size(); ++i )
            {
                binClause.clear();
                binClause.push( ~lit_pt );
                binClause.push( Minisat::mkLit( inputs[i].node()->_satVar, inputs[i].isInverted() ) );
                _solver->addClause( binClause );

                naryClause.push( Minisat::mkLit( inputs[i].node()->_satVar, !inputs[i].isInverted() ) );
            }

            _solver->addClause( naryClause );
        }
        else
        {
            /*
             * a = b ^ c
             * (!a+b+c)
             * (!a+!b+!c)
             * (a+b+!c)
             * (a+!b+c)
             */
            assert(inputs.size() == 2 );

            naryClause.clear();
            naryClause.push( ~lit_pt );
            naryClause.push( Minisat::mkLit( inputs[0].node()->_satVar, inputs[0].isInverted() ) );
            naryClause.push( Minisat::mkLit( inputs[1].node()->_satVar, inputs[1].isInverted() ) );
            _solver->addClause( naryClause );

            naryClause.clear();
            naryClause.push( ~lit_pt );
            naryClause.push( ~Minisat::mkLit( inputs[0].node()->_satVar, inputs[0].isInverted() ) );
            naryClause.push( ~Minisat::mkLit( inputs[1].node()->_satVar, inputs[1].isInverted() ) );
            _solver->addClause( naryClause );

            naryClause.clear();
            naryClause.push( lit_pt );
            naryClause.push( Minisat::mkLit( inputs[0].node()->_satVar, inputs[0].isInverted() ) );
            naryClause.push( ~Minisat::mkLit( inputs[1].node()->_satVar, inputs[1].isInverted() ) );
            _solver->addClause( naryClause );

            naryClause.clear();
            naryClause.push( lit_pt );
            naryClause.push( ~Minisat::mkLit( inputs[0].node()->_satVar, inputs[0].isInverted() ) );
            naryClause.push( Minisat::mkLit( inputs[1].node()->_satVar, inputs[1].isInverted() ) );
            _solver->addClause( naryClause );
        }

        pending.pop();
    }

    unlockFlags<1>();
}

#else
void
aigpp::Manager::createClauses(aigpp::Node* n1, aigpp::Node* n2)
{
    lockFlags<1>();

    std::stack<Node*>& pending = _universalNodeStack;
    assert(pending.empty());

    pending.push(nullptr);
    pending.push(n2);

    Node* n = n1;

    do {
        if (n->isSatVarValid()) {
            n = pending.top();
            pending.pop();
        } else if (n->isVar()) {
            n->setSatVar(_solver->newVar());

            _nodesWithSATVars.insert(n);

            n = pending.top();
            pending.pop();
        } else if (!(n->parent1().node()->isSatVarValid())) {
            pending.push(n);
            n = n->parent1().node();
        } else if (!(n->parent2().node()->isSatVarValid())) {
            pending.push(n);
            n = n->parent2().node();
        } else {
            // n = (!)n1 & (!)n2  <=>  (!n | (!)n1) & (!n | (!)n2) & (n | !(!)n1 |
            // !(!)n2)
            n->setSatVar(_solver->newVar());
            _nodesWithSATVars.insert(n);

            Minisat::Lit vn = Minisat::mkLit(n->_satVar, n->isNAND());

            Minisat::Var vn1 = n->parent1().node()->_satVar;
            Minisat::Var vn2 = n->parent2().node()->_satVar;

            Minisat::vec<Minisat::Lit> clause1(2), clause2(2), clause3(3);

            clause1[0] = ~vn;
            clause1[1] = Minisat::mkLit(vn1, n->parent1().isInverted());
            clause2[0] = ~vn;
            clause2[1] = Minisat::mkLit(vn2, n->parent2().isInverted());
            clause3[0] = vn;
            clause3[1] = Minisat::mkLit(vn1, !(n->parent1().isInverted()));
            clause3[2] = Minisat::mkLit(vn2, !(n->parent2().isInverted()));

            _solver->addClause(clause1);
            _solver->addClause(clause2);
            _solver->addClause(clause3);

            ++_nodesInCNF;

            n = pending.top();
            pending.pop();
        }
    } while (!pending.empty());

    unlockFlags<1>();
}
#endif

aigpp::VariableAssignment
aigpp::Manager::getCounterExample() const
{
    VariableAssignment c(variableCount());

    for (Node* v : _variables) {
        if (v->isSatVarValid()) {
            Minisat::lbool b = _solver->model[v->_satVar];

            if (b == Minisat::l_True) {
                c.setPositive(v->varIndex());
            } else if (b == Minisat::l_False) {
                c.setNegative(v->varIndex());
            }
        }
    }
    return c;
}
