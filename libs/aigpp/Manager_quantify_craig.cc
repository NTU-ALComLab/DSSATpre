#if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2) or defined(USE_INCMINICRAIG)

#    include "Manager.hh"

#    include <algorithm>

#    include <lrabsutil/Resources.hh>
#    include <lrabsutil/String.hh>

#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
#        include <minicraig/minisat_craig.hpp>
#    endif

#    if defined(USE_INCMINICRAIG)
#        include <incminicraig/incrementalminicraig.hpp>
#    endif

//#define DEBUG_CRAIG_QUANTIFY

#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
aigpp::InternalEdgeRef
simpleAIG2aigpp(const minicraig::SimpleAIG& aig, const minicraig::SimpleAIGEdge& root,
                const std::map<int, aigpp::InternalEdgeRef>& varMap)
{
    if (root == aig.getFalse()) {
        return aigpp::const0;
    } else if (root == aig.getTrue()) {
        return aigpp::const1;
    }

    std::stack<int>                       pending;
    std::map<int, aigpp::InternalEdgeRef> mapping;

    pending.push(minicraig::SimpleAIG::getIndex(root));

    while (!pending.empty()) {
        if (mapping.find(pending.top()) != mapping.end()) {
            pending.pop();
            continue;
        }

        const minicraig::SimpleAIGNode& n = aig.getNode(pending.top());

        /* current SimpleAIGNode is a variable -> lookup variable index in varMap */
        if (n.isVar()) {
            auto v = varMap.find(n.varIndex());
            assert(v != varMap.end());

            mapping[pending.top()] = v->second;

            pending.pop();
            continue;
        }

        std::map<int, aigpp::InternalEdgeRef>::const_iterator p1, p2;

        if ((p1 = mapping.find(minicraig::SimpleAIG::getIndex(n.p1()))) == mapping.end()) {
            pending.push(minicraig::SimpleAIG::getIndex(n.p1()));
        } else if ((p2 = mapping.find(minicraig::SimpleAIG::getIndex(n.p2()))) == mapping.end()) {
            pending.push(minicraig::SimpleAIG::getIndex(n.p2()));
        } else {
            aigpp::InternalEdgeRef result = p1->second.notIf(minicraig::SimpleAIG::isNegated(n.p1()))
                                            & p2->second.notIf(minicraig::SimpleAIG::isNegated(n.p2()));

            mapping[pending.top()] = result;
            pending.pop();
        }
    }

    return mapping[minicraig::SimpleAIG::getIndex(root)].notIf(minicraig::SimpleAIG::isNegated(root));
}
#    endif

#    if defined(USE_INCMINICRAIG)
aigpp::InternalEdgeRef
simpleAIG2aigpp(const simpaig::SimpleAIG& aig, const simpaig::SimpleAIGEdge& root,
                const std::map<int, aigpp::InternalEdgeRef>& varMap)
{
    if (root == aig.getFalse()) {
        return aigpp::const0;
    } else if (root == aig.getTrue()) {
        return aigpp::const1;
    }

    std::stack<int>                       pending;
    std::map<int, aigpp::InternalEdgeRef> mapping;

    pending.push(simpaig::SimpleAIG::getIndex(root));

    while (!pending.empty()) {
        if (mapping.find(pending.top()) != mapping.end()) {
            pending.pop();
            continue;
        }

        const simpaig::SimpleAIGNode& n = aig.getNode(pending.top());

        /* current SimpleAIGNode is a variable -> lookup variable index in varMap */
        if (n.isVar()) {
            std::map<int, aigpp::InternalEdgeRef>::const_iterator v = varMap.find(n.varIndex());
            assert(v != varMap.end());

            mapping[pending.top()] = v->second;

            pending.pop();
            continue;
        }

        std::map<int, aigpp::InternalEdgeRef>::const_iterator p1, p2;

        if ((p1 = mapping.find(simpaig::SimpleAIG::getIndex(n.p1()))) == mapping.end()) {
            pending.push(simpaig::SimpleAIG::getIndex(n.p1()));
        } else if ((p2 = mapping.find(simpaig::SimpleAIG::getIndex(n.p2()))) == mapping.end()) {
            pending.push(simpaig::SimpleAIG::getIndex(n.p2()));
        } else {
            aigpp::InternalEdgeRef result = p1->second.notIf(simpaig::SimpleAIG::isNegated(n.p1()))
                                            & p2->second.notIf(simpaig::SimpleAIG::isNegated(n.p2()));

            mapping[pending.top()] = result;
            pending.pop();
        }
    }

    return mapping[simpaig::SimpleAIG::getIndex(root)].notIf(simpaig::SimpleAIG::isNegated(root));
}
#    endif

void
aigpp::Manager::collectMultiAndInputs(aigpp::Node* root, std::vector<aigpp::Edge>& inputs) const
{
    assert(root != nullptr);
    assert(!root->isVar());

    static std::stack<Node*> pending;

    inputs.clear();

    pending.push(root);

    while (!pending.empty()) {
        Node* pt = pending.top();
        pending.pop();

        if ((pt->parent1().isInverted() != pt->parent1().node()->isNAND()) || pt->parent1().node()->isVar()
            || pt->parent1().node()->refCount() > 1) {
            inputs.push_back(pt->parent1());
        } else {
            pending.push(pt->parent1().node());
        }

        if ((pt->parent2().isInverted() != pt->parent2().node()->isNAND()) || pt->parent2().node()->isVar()
            || pt->parent2().node()->refCount() > 1) {
            inputs.push_back(pt->parent2());
        } else {
            pending.push(pt->parent2().node());
        }
    }

    assert(inputs.size() >= 2);
}

bool
aigpp::Manager::collectXORInputs(aigpp::Node* root, std::vector<aigpp::Edge>& inputs) const
{
    assert(root != nullptr);
    assert(!root->isVar());

    inputs.clear();

    const Edge& p1 = root->parent1();
    const Edge& p2 = root->parent2();

    if (p1.node()->isVar() || (p1.node()->refCount() > 1) || (p1.isInverted() == p1.node()->isNAND())) {
        return false;
    }
    if (p2.node()->isVar() || (p2.node()->refCount() > 1) || (p2.isInverted() == p2.node()->isNAND())) {
        return false;
    }

    if (p1.node()->parent1().structurallyAntivalent(p2.node()->parent1())
        && p1.node()->parent2().structurallyAntivalent(p2.node()->parent2())) {
        inputs.push_back(p1.node()->parent1());
        inputs.push_back(p1.node()->parent2());

        return true;
    } else {
        return false;
    }
}

#    ifdef USE_MINICRAIG
inline minicraig::Lit
createLit(minicraig::Var var, bool neg = false)
{
    return minicraig::Lit(var, neg);
}
#    endif
#    ifdef USE_MINICRAIG2
inline minicraig::Lit
createLit(minicraig::Var var, bool neg = false)
{
    return minicraig::mkLit(var, neg);
}
#    endif

#    ifdef USE_MINICRAIG
inline minicraig::Lit
negateLit(minicraig::Lit lit, bool neg)
{
    return minicraig::id(lit, neg);
}
#    endif
#    ifdef USE_MINICRAIG2
inline minicraig::Lit
negateLit(minicraig::Lit lit, bool neg)
{
    return lit ^ neg;
}
#    endif

/*
 *
 * Implementation of
 *
 * Jie-Hong R. Jiang
 * Quantifier Elimination via Functional Composition
 * CAV 2009
 *
 */
#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
std::pair<aigpp::InternalEdgeRef, aigpp::InternalEdgeRef>
aigpp::Manager::quantifierEliminationBySubstitution(const aigpp::InternalEdgeRef& root,
                                                    const aigpp::InternalEdgeRef& var, std::size_t wantedSize,
                                                    double timeout, bool* failed)
{
    std::pair<InternalEdgeRef, InternalEdgeRef> result;
    result.first = root;
    if (failed != nullptr) {
        *failed = false;
    }

    if (root.isConstant()) {
        result.second = root;
        return result;
    }

    assert(!(var.isConstant()));
    assert(var.node()->isVar());

    InternalEdgeRef cof0 = root.cofactor(!var);
    InternalEdgeRef cof1 = root.cofactor(var);
    // RW: Original version.
    //    if( cof0.isConstant() || cof1.isConstant() )
    //    {
    //        // at least one cofactor is constant
    //        result.first = cof0 + cof1;
    //        result.second = const0;
    //        *failed = false;
    //        return result;
    //    }
    // RW: new version BEGIN
    if (cof0.isConstant()) {
        result.second = !cof0;
        *failed       = false;
        return result;
    } else if (cof1.isConstant()) {
        result.second = cof1;
        *failed       = false;
        return result;
    }
    // RW: new version END

    /* compute union of the supports of both cofactors -> global variables for
     * interpolation! */
    std::set<Node*> support;
    {
        std::vector<Node*> tempsupport;

        tempsupport = cof0.node()->supportNodes();
        support.insert(tempsupport.cbegin(), tempsupport.cend());

        tempsupport = cof1.node()->supportNodes();
        support.insert(tempsupport.cbegin(), tempsupport.cend());
    }

    /* detection of multi-ands */

    const double startCraigTime = lrabs::cpuTime();

#        ifdef USE_MINICRAIG
    minicraig::MinisatCraig craigSolver;
#        endif

#        ifdef USE_MINICRAIG2
    const bool sym  = settings().quantifyUseDualInterpolation();
    const bool asym = true;
    // minicraig::MinisatCraig craigSolver( sym, asym );
    minicraig::MinisatCraig craigSolver(false, true);
#        endif

    /* allocate literals for global variables, setup mappings */
    int                                                        allocatedVariables = 0;
    std::map<Node*, std::pair<minicraig::Lit, minicraig::Lit>> mapping;
    std::map<int, InternalEdgeRef>                             backMapping;

    for (Node* const v : support) {
        assert(v != var.node());

        const minicraig::Lit vlit = createLit(craigSolver.newVar());

        ++allocatedVariables;

        mapping[v]                        = std::make_pair(vlit, vlit);
        backMapping[minicraig::var(vlit)] = InternalEdgeRef(v);
    }

    /* create clause sets */
    std::vector<Edge>                                                                       inputs;
    std::vector<std::map<Node*, std::pair<minicraig::Lit, minicraig::Lit>>::const_iterator> inputVariables;

    std::vector<std::vector<minicraig::Lit>> clausesA, clausesB;
    {
        std::size_t                 mergedNodes = 0;
        std::stack<Node*>           pending;
        std::vector<minicraig::Lit> binClause, naryClause;
        // RW:        std::vector<minicraig::Lit> binClause_b, naryClause_b;
        //        std::map<Node*, std::pair<minicraig::Lit, minicraig::Lit>
        //        >::const_iterator p1, p2;

        pending.push(cof0.node());
        pending.push(cof1.node());

        while (!pending.empty()) {
            Node* pt = pending.top();

            if (mapping.find(pt) != mapping.end()) {
                pending.pop();
                continue;
            }

            assert(!pt->isVar());
            {
                collectMultiAndInputs(pt, inputs);
                assert(inputs.size() >= 2);
                mergedNodes += inputs.size() - 2;
            }

            inputVariables.clear();

            for (const Edge& p : inputs) {
                std::map<Node*, std::pair<minicraig::Lit, minicraig::Lit>>::const_iterator p_it
                    = mapping.find(p.node());
                if (p_it == mapping.end()) {
                    pending.push(p.node());
                } else {
                    inputVariables.push_back(p_it);
                }
            }

            if (pt == pending.top()) {
                minicraig::Lit alit = createLit(craigSolver.newVar());
                minicraig::Lit blit = createLit(craigSolver.newVar());

                allocatedVariables += 2;

                mapping[pt] = std::pair<minicraig::Lit, minicraig::Lit>(alit, blit);

                if (pt->isNAND()) {
                    alit = ~alit;
                    blit = ~blit;
                }

                {
                    naryClause.clear();
                    naryClause.push_back(alit);

                    for (std::size_t i = 0; i != inputs.size(); ++i) {
                        binClause.clear();

                        binClause.push_back(~alit);
                        binClause.push_back(
                            negateLit(inputVariables[i]->second.first, (inputs[i].isInverted() ? 1 : 0)));
                        clausesA.push_back(binClause);

                        naryClause.push_back(
                            negateLit(inputVariables[i]->second.first, (inputs[i].isInverted() ? 0 : 1)));
                    }
                    clausesA.push_back(naryClause);

                    naryClause.clear();
                    naryClause.push_back(blit);

                    for (std::size_t i = 0; i != inputs.size(); ++i) {
                        binClause.clear();

                        binClause.push_back(~blit);
                        binClause.push_back(
                            negateLit(inputVariables[i]->second.second, (inputs[i].isInverted() ? 1 : 0)));

                        clausesB.push_back(binClause);

                        naryClause.push_back(
                            negateLit(inputVariables[i]->second.second, (inputs[i].isInverted() ? 0 : 1)));
                    }
                    clausesB.push_back(naryClause);
                }

                pending.pop();
            }
        }

        /* !cof0 => a */
        binClause.clear();
        binClause.push_back(~(negateLit(mapping[cof0.node()].first, cof0.isInverted())));
        clausesA.push_back(binClause);

        /*  cof1 => a */
        binClause.clear();
        binClause.push_back(negateLit(mapping[cof1.node()].first, cof1.isInverted()));
        clausesA.push_back(binClause);

        /*  cof0 => b */
        binClause.clear();
        binClause.push_back(negateLit(mapping[cof0.node()].second, cof0.isInverted()));
        clausesB.push_back(binClause);

        /* !cof1 => b */
        binClause.clear();
        binClause.push_back(~(negateLit(mapping[cof1.node()].second, cof1.isInverted())));
        clausesB.push_back(binClause);
    }

#        if 0
    {
        static int dump = -1;
        ++dump;

        std::cout << "\n\n\n\n\ndumping " << dump << std::endl;

        std::ofstream fileA( ( "a_" + lrabs::toString( dump ) + ".dimacs" ).c_str() );
        int maxA = -1;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesA.begin(); c != clausesA.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                if( minicraig::var(*p) > maxA ) maxA = minicraig::var(*p);
            }
        }
        fileA << "p cnf " << maxA << " " << clausesA.size() << std::endl;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesA.begin(); c != clausesA.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                fileA << ( minicraig::sign(*p) ? "-" : "" ) << minicraig::var(*p) << " ";
            }
            fileA << "0" << std::endl;
        }


        std::ofstream fileB( ( "b_" + lrabs::toString( dump ) + ".dimacs" ).c_str() );
        int maxB = -1;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesB.begin(); c != clausesB.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                if( minicraig::var(*p) > maxB ) maxB = minicraig::var(*p);
            }
        }
        fileB << "p cnf " << maxB << " " << clausesB.size() << std::endl;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesB.begin(); c != clausesB.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                fileB << ( minicraig::sign(*p) ? "-" : "" ) << minicraig::var(*p) << " ";
            }
            fileB << "0" << std::endl;
        }
    }
#        endif

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during initialization (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed       = true;
            result.second = const0;
            return result;
        }
    }

    std::size_t index                = 0;
    std::size_t timeoutCheckInterval = clausesA.size() / 10 + 1;

    /* add clause sets to the craig solver */
    index = 0;
    for (const std::vector<minicraig::Lit>& c : clausesB) {
        craigSolver.addBClause(c);

        ++index;
        if ((index % timeoutCheckInterval) == 0 && failed != nullptr) {
            if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
                std::cout << "timelimit hit during b clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
                *failed = true;

                result.second = const0;
                return result;
            }
        }
    }

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during b clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed = true;

            result.second = const0;
            return result;
        }
    }

    index = 0;
    for (const std::vector<minicraig::Lit>& c : clausesA) {
        craigSolver.addAClause(c);

        ++index;
        if ((index % timeoutCheckInterval) == 0 && failed != nullptr) {
            if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
                std::cout << "timelimit hit during a clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
                *failed = true;

                result.second = const0;
                return result;
            }
        }
    }

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during a clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed = true;

            result.second = const0;
            return result;
        }
    }

    if (failed != nullptr) {
        craigSolver.setResourceLimits(static_cast<int>(wantedSize), timeout, failed);
    }

#        ifdef DEBUG_CRAIG_QUANTIFY
    std::cout << "solving (line " << __LINE__ << ")" << std::endl;
#        endif

#        ifndef NDEBUG
    const bool satResult =
#        endif
        craigSolver.solve();
    if (failed != nullptr) {
        if (*failed) {
            result.second = const0;
            return result;
        }
    }

    assert(!satResult);

#        ifdef DEBUG_CRAIG_QUANTIFY
    const int size = craigSolver.getInterpolantAIG().getConeSize(craigSolver.getInterpolant());

    std::cout << "CRAIG:"
              << " learnts=" << craigSolver.internalSolver()->nLearnts()
              << " total_conflicts=" << craigSolver.internalSolver()->conflicts
              << " learnts_lits=" << craigSolver.internalSolver()->learnts_literals
              << " aig_nodes=" << craigSolver.getInterpolantAIG().size()
              << " and_nodes=" << craigSolver.getInterpolantAIG().andNodes() << " inter_simp=" << size << std::endl;
#        endif

#        ifdef DEBUG_CRAIG_QUANTIFY
    std::cout << "CRAIG: using original interpolant (7)" << std::endl;
#        endif

    result.second = simpleAIG2aigpp(craigSolver.getInterpolantAIG(), craigSolver.getInterpolant(), backMapping);

    return result;
}
#    endif

#    ifdef USE_INCMINICRAIG
std::pair<aigpp::InternalEdgeRef, aigpp::InternalEdgeRef>
aigpp::Manager::quantifierEliminationBySubstitution(const aigpp::InternalEdgeRef& root,
                                                    const aigpp::InternalEdgeRef& var, std::size_t wantedSize,
                                                    double timeout, bool* failed)
{
    std::pair<InternalEdgeRef, InternalEdgeRef> result;
    result.first = root;

    if (root.isConstant()) {
        result.second = root;
        return result;
    }

    assert(!(var.isConstant()));
    assert(var.node()->isVar());

    InternalEdgeRef cof0 = root.cofactor(!var);
    InternalEdgeRef cof1 = root.cofactor(var);

    if (cof0.isConstant() || cof1.isConstant()) {
#        ifdef DEBUG_CRAIG_QUANTIFY
        std::cout << "INFO: at least one cofactor is constant" << std::endl;
#        endif
        result.first  = cof0 + cof1;
        result.second = const0;
        return result;
    }

    /* compute union of the supports of both cofactors -> global variables for
     * interpolation! */
    std::set<Node*> support;
    {
        std::vector<Node*> tempsupport;

        tempsupport = cof0.node()->supportNodes();
        support.insert(tempsupport.begin(), tempsupport.end());

        tempsupport = cof1.node()->supportNodes();
        support.insert(tempsupport.begin(), tempsupport.end());
    }

    /* detection of multi-ands */

    double startCraigTime = lrabs::cpuTime();

    incminicraig::IncMiniCraig solver;
    unsigned int               nextVariable = 0;

    std::cout << "computing SYM and ASYM interpolants" << std::endl;
    // std::cout << "computing SYM interpolant" << std::endl;
    // std::cout << "computing ASYM interpolant" << std::endl;

    solver.setConstructionRules(incminicraig::BOTH, true);
    // solver.setConstructionRules( incminicraig::SYMMETRIC, false );
    // solver.setConstructionRules( incminicraig::ASYMMETRIC, false );

    if (failed != 0) {
        solver.setTimeout(timeout);
        solver.setSymNodeBound(static_cast<int>(wantedSize));
        solver.setAsymNodeBound(static_cast<int>(wantedSize));
    }

    /* allocate literals for global variables, setup mappings */
    int                                                    allocatedVariables = 0;
    std::map<Node*, std::pair<unsigned int, unsigned int>> mapping;
    std::map<int, InternalEdgeRef>                         backMapping;

    for (std::set<Node*>::const_iterator v = support.begin(); v != support.end(); ++v) {
        assert(*v != var.node());

        unsigned int vlit = (nextVariable++) << 1;

        ++allocatedVariables;

        mapping[*v]              = std::pair<unsigned int, unsigned int>(vlit, vlit);
        backMapping[(vlit >> 1)] = InternalEdgeRef(*v);
    }

    /* create clause sets */
    std::vector<Edge>                                                                   inputs;
    std::vector<std::map<Node*, std::pair<unsigned int, unsigned int>>::const_iterator> inputVariables;

    int                                    mergedNodes = 0;
    std::vector<std::vector<unsigned int>> clausesA, clausesB;
    {
        std::stack<Node*>                                                      pending;
        std::vector<unsigned int>                                              binClause, naryClause;
        std::vector<unsigned int>                                              binClause_b, naryClause_b;
        std::map<Node*, std::pair<unsigned int, unsigned int>>::const_iterator p1, p2;

        pending.push(cof0.node());
        pending.push(cof1.node());

        while (!pending.empty()) {
            Node* pt = pending.top();

            if (mapping.find(pt) != mapping.end()) {
                pending.pop();
                continue;
            }

            assert(!pt->isVar());
            {
                collectMultiAndInputs(pt, inputs);
                mergedNodes += inputs.size() - 2;
            }

            inputVariables.clear();

            for (std::vector<Edge>::const_iterator p = inputs.begin(); p != inputs.end(); ++p) {
                std::map<Node*, std::pair<unsigned int, unsigned int>>::const_iterator p_it = mapping.find(p->node());
                if (p_it == mapping.end()) {
                    pending.push(p->node());
                } else {
                    inputVariables.push_back(p_it);
                }
            }

            if (pt == pending.top()) {
                unsigned int alit = (nextVariable++) << 1;
                unsigned int blit = (nextVariable++) << 1;
                allocatedVariables += 2;

                mapping[pt] = std::pair<unsigned int, unsigned int>(alit, blit);

                if (pt->isNAND()) {
                    alit = alit ^ 1;
                    blit = blit ^ 1;
                }

                {
                    naryClause.clear();
                    naryClause.push_back(alit);

                    for (std::size_t i = 0; i != inputs.size(); ++i) {
                        binClause.clear();

                        binClause.push_back(alit ^ 1);
                        binClause.push_back((inputVariables[i]->second.first ^ (inputs[i].isInverted() ? 1 : 0)));
                        clausesA.push_back(binClause);

                        naryClause.push_back((inputVariables[i]->second.first ^ (inputs[i].isInverted() ? 0 : 1)));
                    }
                    clausesA.push_back(naryClause);

                    naryClause.clear();
                    naryClause.push_back(blit);

                    for (std::size_t i = 0; i != inputs.size(); ++i) {
                        binClause.clear();

                        binClause.push_back(blit ^ 1);
                        binClause.push_back((inputVariables[i]->second.second ^ (inputs[i].isInverted() ? 1 : 0)));
                        clausesB.push_back(binClause);

                        naryClause.push_back((inputVariables[i]->second.second ^ (inputs[i].isInverted() ? 0 : 1)));
                    }
                    clausesB.push_back(naryClause);
                }

                pending.pop();
            }
        }

        /* !cof0 => a */
        binClause.clear();
        binClause.push_back((mapping[cof0.node()].first ^ cof0.isInverted()) ^ 1);
        clausesA.push_back(binClause);

        /*  cof1 => a */
        binClause.clear();
        binClause.push_back((mapping[cof1.node()].first ^ cof1.isInverted()));
        clausesA.push_back(binClause);

        /*  cof0 => b */
        binClause.clear();
        binClause.push_back((mapping[cof0.node()].second ^ cof0.isInverted()));
        clausesB.push_back(binClause);

        /* !cof1 => b */
        binClause.clear();
        binClause.push_back((mapping[cof1.node()].second ^ cof1.isInverted()) ^ 1);
        clausesB.push_back(binClause);
    }

#        if 0
    {
        static int dump = -1;
        ++dump;

        std::cout << "\n\n\n\n\ndumping " << dump << std::endl;

        std::ofstream fileA( ( "a_" + lrabs::toString( dump ) + ".dimacs" ).c_str() );
        int maxA = -1;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesA.begin(); c != clausesA.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                if( minicraig::var(*p) > maxA ) maxA = minicraig::var(*p);
            }
        }
        fileA << "p cnf " << maxA << " " << clausesA.size() << std::endl;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesA.begin(); c != clausesA.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                fileA << ( minicraig::sign(*p) ? "-" : "" ) << minicraig::var(*p) << " ";
            }
            fileA << "0" << std::endl;
        }


        std::ofstream fileB( ( "b_" + lrabs::toString( dump ) + ".dimacs" ).c_str() );
        int maxB = -1;
        for( std::vector<std::vector<minicraig::Lit> >::const_iterator c = clausesB.begin(); c != clausesB.end(); ++c )
        {
            for( std::vector<minicraig::Lit>::const_iterator p = c->begin(); p != c->end(); ++p )
            {
                if( minicraig::var(*p) > maxB ) maxB = minicraig::var(*p);
            }
        }
        fileB << "p cnf " << maxB << " " << clausesB.size() << std::endl;
        for( const std::vector<minicraig::Lit>& c: clausesB )
        {
            for( minicraig::Lit p: c )
            {
                fileB << ( minicraig::sign(p) ? "-" : "" ) << minicraig::var(p) << " ";
            }
            fileB << "0" << std::endl;
        }
    }
#        endif
    for (const std::vector<unsigned int>& c : clausesB) {
        solver.addClause(c, incminicraig::B_CLAUSE);
    }
    for (const std::vector<unsigned int>& c : clausesA) {
        solver.addClause(c, incminicraig::A_CLAUSE);
    }

    std::vector<std::pair<unsigned int, incminicraig::clause_type>> assumptions;
    incminicraig::solver_result                                     satResult = solver.solve(assumptions);

    if (satResult != incminicraig::UNSAT) {
        assert(satResult != incminicraig::SAT);

        assert(failed != 0);
        *failed = true;
        return result;
    }

    int sizeSym  = -1;
    int sizeAsym = -1;

    if (solver.isInterpolantValid(true)) {
        sizeSym = solver.getAigManager(true)->getConeSize(solver.getAigInterpolant(true));
    }
    if (solver.isInterpolantValid(false)) {
        sizeAsym = solver.getAigManager(false)->getConeSize(solver.getAigInterpolant(false));
    }

#        ifdef DEBUG_CRAIG_QUANTIFY
    std::cout << "CRAIG:"
              << " inter_simp (sym/asym) = " << sizeSym << " " << sizeAsym << std::endl;
#        endif

    assert(sizeSym != -1 || sizeAsym != -1);
    if (sizeAsym == -1 || (sizeSym != -1 && sizeSym <= sizeAsym)) {
        result.second = simpleAIG2aigpp(*(solver.getAigManager(true)), solver.getAigInterpolant(true), backMapping);
    } else {
        result.second = simpleAIG2aigpp(*(solver.getAigManager(false)), solver.getAigInterpolant(false), backMapping);
    }

    if (failed != 0) {
        *failed = false;
    }

    return result;
}
#    endif

#endif /* USE_MINICRAIG or USE_MINICRAIG2 */
