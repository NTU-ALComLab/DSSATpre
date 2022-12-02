#if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2) or defined(USE_INCMINICRAIG)

#    include "Manager.hh"

#    include <algorithm>

#    include <lrabsutil/Resources.hh>
#    include <lrabsutil/String.hh>

#    define _unused(x) ((void)(x))

#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
#        include <minicraig/minisat_craig.hpp>
#    endif

#    if defined(USE_INCMINICRAIG)
#        include <incminicraig/incrementalminicraig.hpp>
#    endif

#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2) or defined(USE_INCMINICRAIG)
aigpp::InternalEdgeRef simpleAIG2aigpp(const minicraig::SimpleAIG& aig, const minicraig::SimpleAIGEdge& root,
                                       const std::map<int, aigpp::InternalEdgeRef>& varMap);
#    endif

namespace {
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
}  // end anonymous namespace

#    if defined(USE_MINICRAIG) or defined(USE_MINICRAIG2)
aigpp::InternalEdgeRef
aigpp::Manager::dcmin(const aigpp::InternalEdgeRef& function, const aigpp::InternalEdgeRef& care,
                      std::size_t wantedSize, double timeout, bool* failed)
{
    if (failed != nullptr) {
        *failed = false;
    }

    if (function.isConstant() || functionallyTrue(care)) {
        return function;
    }

    if (functionallyFalse(care)) {
        return const1;
    }

    /* compute union of the supports of both 'function' and 'care' -> global
     * variables for interpolation! */
    std::set<Node*> support;
    {
        std::vector<Node*> tempsupport;

        tempsupport = function.node()->supportNodes();
        support.insert(tempsupport.cbegin(), tempsupport.cend());

        tempsupport = care.node()->supportNodes();
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

        pending.push(function.node());  // TODO CHECK
        pending.push(care.node());      // TODO CHECK

        while (!pending.empty()) {
            Node* pt = pending.top();

            if (mapping.find(pt) != mapping.cend()) {
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
                const auto p_it = mapping.find(p.node());
                if (p_it == mapping.end()) {
                    pending.push(p.node());
                } else {
                    inputVariables.emplace_back(p_it);
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
        /* function => A */
        binClause.clear();
        //        binClause.push_back( ~( negateLit( mapping[cof0.node()].first,
        //        cof0.isInverted() ) ) );
        binClause.push_back(negateLit(mapping[function.node()].first, function.isInverted()));
        clausesA.push_back(binClause);

        /* cof1 => a */
        /* care => A */
        binClause.clear();
        //        binClause.push_back( negateLit( mapping[cof1.node()].first,
        //        cof1.isInverted() ) );
        binClause.push_back(negateLit(mapping[care.node()].first, care.isInverted()));
        clausesA.push_back(binClause);

        /*  cof0 => b */
        /* !function => B */
        binClause.clear();
        //        binClause.push_back( negateLit( mapping[cof0.node()].second,
        //        cof0.isInverted() ) );
        binClause.push_back(~(negateLit(mapping[function.node()].second, function.isInverted())));
        clausesB.push_back(binClause);

        /* !cof1 => b */
        /* care => B */
        binClause.clear();
        //        binClause.push_back( ~( negateLit( mapping[cof1.node()].second,
        //        cof1.isInverted() ) ) );
        binClause.push_back(negateLit(mapping[care.node()].second, care.isInverted()));
        clausesB.push_back(binClause);
    }

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during initialization (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed = true;
            return const0;
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
                return const0;
            }
        }
    }

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during b clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed = true;
            return const0;
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
                return const0;
            }
        }
    }

    if (failed != nullptr) {
        if ((timeout > 0) && (lrabs::cpuTime() > startCraigTime + timeout)) {
#        ifdef DEBUG_CRAIG_QUANTIFY
            std::cout << "timelimit hit during a clause insertion (line " << __LINE__ << ")" << std::endl;
#        endif
            *failed = true;
            return const0;
        }
    }

    if (failed != nullptr) {
        craigSolver.setResourceLimits(static_cast<int>(wantedSize), timeout, failed);
    }

#        ifdef DEBUG_CRAIG_QUANTIFY
    std::cout << "solving (line " << __LINE__ << ")" << std::endl;
#        endif

    const bool satResult = craigSolver.solve();
    _unused(satResult);

    if (failed != nullptr) {
        if (*failed) {
            return const0;
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

    return simpleAIG2aigpp(craigSolver.getInterpolantAIG(), craigSolver.getInterpolant(), backMapping);
}
#    endif

#endif /* USE_MINICRAIG or USE_MINICRAIG2 */
