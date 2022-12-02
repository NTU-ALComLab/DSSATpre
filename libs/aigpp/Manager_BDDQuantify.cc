/**************************************************************
 *
 *       aigpp // Manager_quantify.cc
 *
 *       Copyright (C) 2008 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 462 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include <fstream>
#include <sstream>

#include <cuddInt.h>
#include <lrabsutil/Resources.hh>

#include "CuddTools.hh"
#include "Manager.hh"

/* used for the competition version of aigqbf */
#define SKIP_FRAIGING_OF_BIG_AIGS
//#define DEBUG_BDD_QUANTIFY
#define FRAIGING_TIMEOUT (5)
#define LOWER_BDD_LIMIT (1000000.0)
#define UPPER_BDD_LIMIT (10000000.0)

#define _unused(x) ((void)(x))

aigpp::EdgeRef
aigpp::Manager::bddQuantify(const aigpp::EdgeRef& root, const std::vector<aigpp::EdgeRef>& variables, bool existential,
                            bool preferAIG, double timeout, bool* timeoutReached)
{
    aigpp::EdgeRef result = root;
    bddQuantifyInplace(result, variables, existential, preferAIG, true, timeout, timeoutReached);
    return result;
}

// prefer aig if 0<=size<=1.2*initial
void
aigpp::Manager::bddQuantifyInplace(aigpp::EdgeRef& root, const std::vector<aigpp::EdgeRef>& variables, bool existential,
                                   bool /* preferAIG */, bool optimizelast, double timeout, bool* timeoutReached)
{
    _unused(timeout);

    static std::size_t bddquantifications     = 0;
    static std::size_t bddquantifiedvariables = 0;
    static std::size_t aigquantifications     = 0;

    assert((timeout >= 0) == (timeoutReached != nullptr));
    if (timeoutReached != nullptr) {
        *timeoutReached = false;
    }

    /* setup bdd manager */
    cleanBDD();

    enum QuantificationMethod
    {
        QM_AIG,
        QM_BDD
    };

    static QuantificationMethod nextMethod = QM_BDD;

    /*
     * skip aig optimizations after 2 unsuccessful optimizations in sequence were
     * made. skip at most 20 optimizations, abort skipping if aig has grown too
     * much (20%)
     *
     * use a maximum of 20 seconds for initial bdd computation
     */
    const std::size_t maxFailedOptimizations = 2;
    const std::size_t maxSkipOptimizations   = 20;
    const double      skipRatio              = 1.2;
    const std::size_t maxSkipFRAIGOnTimeout  = 10;
    const double      maxBddComputationTime  = 20.0;
    const std::size_t maxNodesForFraiging    = 50000;
    // const int maxNodesForFraiging = 200000;

    static std::size_t failedOptimizations      = 0;
    static std::size_t skipOptimizations        = 0;
    static std::size_t maxSkipAigSize           = 0;
    static std::size_t skipFRAIGOnTimeout       = 0;
    static std::size_t consecutiveFRAIGTimeouts = 0;

    static double interpolationTryFactor = 1.5;

    /* get indices of variables that are to be quantified */
    std::set<int> variablesset;
    for (const auto& var : variables) {
        assert(var.isVariable());
        variablesset.insert(var.variableIndex());
    }

    /* get root's support */
    std::vector<int> supportvec = root.support();

    /* determine which variables (of root) have to be quantified */
    std::vector<int> remainingVariables;
    for (const int supp : supportvec) {
        if (variablesset.find(supp) != variablesset.end()) {
            remainingVariables.push_back(supp);
        }
    }

    /* configure scheduling */
    std::size_t nextSchedulingVarCount = remainingVariables.size();
    // int nextSchedulingVarCountInterval = remainingVariables.size() / 10;
    /* reschedule after each elimination */
    std::size_t nextSchedulingVarCountInterval = 1;

    std::size_t startAIGSize   = root.nodeCount();
    bool        lastStepWasAIG = false;

    while (!remainingVariables.empty()) {
        assert(!(root.isConstant()));

        /* simple case: root is a variable -> result is either const1 (existential)
         * or const0 (universal) */
        if (root.getInternal().node()->isVar()) {
#ifndef NDEBUG
            const int varIndex = root.getInternal().node()->varIndex();
            assert(std::find(remainingVariables.begin(), remainingVariables.end(), varIndex)
                   != remainingVariables.end());
#endif

            if (existential) {
                root = getConst1();
            } else {
                root = getConst0();
            }
            remainingVariables.clear();
            break;
        }

        // Current function is not only a variable
        std::size_t currentAIGSize = root.nodeCount();

        DdNode* bdd = nullptr;

        /* try to compute a bdd representation */
        if (nextMethod == QM_BDD) {
            /* try to build bdd */

            /* reset reordering limit */
            Cudd_SetNextReordering(_bddManager, (int)currentAIGSize);

            /* set bdd node limit */
            const double BDDLimit
                = std::max(LOWER_BDD_LIMIT, std::min(UPPER_BDD_LIMIT, 20.0 * static_cast<double>(currentAIGSize)));
            Cudd_SetMaxLive(_bddManager, static_cast<int>(BDDLimit));

            /*
             * build bdd.
             * abort if maxBddComputationTime is reached,
             * or the current number of bdd nodes is 200 times larger than the number
             * of processed aig nodes
             */
            Cudd_SetTimeLimit(_bddManager, static_cast<unsigned long>(maxBddComputationTime));
            bdd = buildBDDfromAIGlimited(root.getInternal().node(),
                                         200.0);  // remark: the returned BDD is already reference-counted.
            Cudd_UnsetTimeLimit(_bddManager);

            if (bdd != nullptr && root.isInverted()) {
                bdd = Cudd_Not(bdd);
            }

            Cudd_SetMaxLive(_bddManager, 1000000000);

            if (bdd == nullptr) {
                startAIGSize = currentAIGSize;
                nextMethod   = QM_AIG;
            }
        }

        /* check if expected bdd quantification result is larger that worst possible
         * aig quantification result */
        if (bdd != nullptr && remainingVariables.size() < 10) {
            const auto bddSize = static_cast<std::size_t>(Cudd_DagSize(bdd));

            if (bddSize > ((1 << remainingVariables.size()) * currentAIGSize)) {
                // BDD too large for the quantifications
                Cudd_RecursiveDeref(_bddManager, bdd);
                bdd        = nullptr;
                nextMethod = QM_AIG;
            }
        }

        /* building bdd failed */
        if (bdd == nullptr) {
            if (!lastStepWasAIG) {
                lastStepWasAIG = true;
                startAIGSize   = currentAIGSize;
            }

            //            if( _bddManager != nullptr )
            //            {
            //                cuddGarbageCollect( _bddManager, true );
            //            }

            // performing aig-based quantification

            /* reschedule */
            if (settings().quantifyUseScheduling() && remainingVariables.size() <= nextSchedulingVarCount) {
                nextSchedulingVarCount = remainingVariables.size() - nextSchedulingVarCountInterval;

                if (remainingVariables.size() > 1) {
                    std::vector<Node*> varNodes;
                    varNodes.reserve(remainingVariables.size());
                    std::transform(remainingVariables.cbegin(), remainingVariables.cend(), std::back_inserter(varNodes),
                                   [this](const int var) -> Node* { return variableInternal(var).node(); });
                    assert(varNodes.size() == remainingVariables.size());

                    std::vector<long> quality
                        = estimateRemainingNodesAfterQuantification(root.getInternal().node(), varNodes);

#ifndef NDEBUG
                    assert(quality.size() == remainingVariables.size());
                    for (long val : quality) {
                        assert(val >= 0);
                    }
#endif

                    /* TODO: better sorting */
                    long             best = -1, worst = -1;
                    std::vector<int> newSchedule(remainingVariables.size(), -1);

                    for (std::size_t i = 0; i != newSchedule.size(); ++i) {
                        long        min      = -1;
                        std::size_t minIndex = 0;

                        for (std::size_t j = 0; j != newSchedule.size(); ++j) {
                            if (quality[j] == -1) {
                                continue;
                            }

                            if (min == -1 || quality[j] < min) {
                                min      = quality[j];
                                minIndex = j;
                            }
                        }

                        assert(min != -1);

                        if (best == -1 || min < best) {
                            best = min;
                        }
                        if (worst == -1 || min > worst) {
                            worst = min;
                        }

                        newSchedule[i]    = remainingVariables[minIndex];
                        quality[minIndex] = -1;
                    }

                    std::swap(remainingVariables, newSchedule);
                }
            }

            InternalEdgeRef var = variableInternal(remainingVariables[0]);
            {
                if (remainingVariables.size() > 1) {
                    remainingVariables[0] = remainingVariables[1];
                } else {
                    remainingVariables.clear();
                }
            }

            ++aigquantifications;

            /* deep quantification */
            EdgeRef rootDeep = quantifyDeep(root, var.node()->varIndex(), existential);

            const std::size_t deepSize = rootDeep.nodeCount();

            /* try with interpolation */
            /* if big aig growth by deep quantification -> try interpolation */

            if (settings().quantifyUseInterpolation()
                && static_cast<double>(deepSize) >= interpolationTryFactor * static_cast<double>(currentAIGSize)) {
                bool            failed = false;
                InternalEdgeRef f      = root.getInternal().notIf(!existential);

                std::pair<InternalEdgeRef, InternalEdgeRef> interResult = quantifierEliminationBySubstitution(
                    f, var, 3 * currentAIGSize, settings().quantifyInterpolationTimeout(), &failed);
                // std::pair<InternalEdgeRef, InternalEdgeRef> interResult =
                // quantifierEliminationBySubstitution( f, var );
                if (!failed) {
                    // interpolation_succeeded
                    const std::size_t cumulatedSize = interResult.first.nodeCount() + interResult.second.nodeCount();

/* disabled for check equivalence! */
#if 1
                    if (cumulatedSize < deepSize) {
                        rootDeep.invalidate();
                    }
#endif

                    std::vector<InternalEdgeRef> varsVector;
                    varsVector.push_back(var);
                    std::vector<InternalEdgeRef> replVector;
                    replVector.push_back(interResult.second);

                    InternalEdgeRef interResult2 = interResult.first.compose(varsVector, replVector);
                    if (!existential) {
                        interResult2 = !interResult2;
                    }

                    const std::size_t interSize = interResult2.nodeCount();

                    if (interSize < deepSize || cumulatedSize < deepSize) {
                        interpolationTryFactor *= 0.9;
                        if (interpolationTryFactor < 1.0) {
                            interpolationTryFactor = 1.0;
                        }
                        rootDeep = DebugGetExternal(interResult2);

                    } else {
                        interpolationTryFactor *= 1.05;
                        if (interpolationTryFactor > 1.5) {
                            interpolationTryFactor = 1.5;
                        }
                    }
                } else {
                    interpolationTryFactor *= 1.1;
                    if (interpolationTryFactor > 1.5) {
                        interpolationTryFactor = 1.5;
                    }
                }
            }

            root = rootDeep;
            rootDeep.invalidate();

            var.clear();
            garbageCollect();

            /* optimization using rewriting and fraig */
            if ((optimizelast || !remainingVariables.empty()) && !(root.isConstant())) {
                if (skipOptimizations > 0 && root.nodeCount() < maxSkipAigSize) {
                    --skipOptimizations;
                } else {
                    std::size_t sizeBefore = root.nodeCount();
                    rewrite();

                    /* FRAIG */
                    if (nodeCount() < maxNodesForFraiging) {
                        if (skipFRAIGOnTimeout > 0) {
                            --skipFRAIGOnTimeout;
                        } else {
                            const std::size_t fraigSizeBefore = nodeCount();
                            makeFRAIG(FRAIGING_TIMEOUT);
                            const std::size_t fraigSizeAfter = nodeCount();

                            if (_timeoutDuringLastFRAIG) {
                                ++consecutiveFRAIGTimeouts;
                                skipFRAIGOnTimeout = consecutiveFRAIGTimeouts * maxSkipFRAIGOnTimeout;

                            } else if (fraigSizeAfter >= fraigSizeBefore) {
                                skipFRAIGOnTimeout = maxSkipFRAIGOnTimeout;
                            } else {
                                consecutiveFRAIGTimeouts = 0;
                            }
                        }
                    }

                    const std::size_t sizeAfter = root.nodeCount();

                    if (static_cast<double>(sizeAfter) >= 0.99 * static_cast<double>(sizeBefore)) {
                        ++failedOptimizations;
                        if (failedOptimizations == maxFailedOptimizations) {
                            failedOptimizations = 0;
                            skipOptimizations   = maxSkipOptimizations;
                            maxSkipAigSize      = static_cast<std::size_t>(skipRatio * static_cast<double>(sizeAfter));
                        }
                    } else {
                        failedOptimizations = 0;
                    }
                }
            }

            if (static_cast<double>(root.nodeCount()) <= 1.2 * static_cast<double>(startAIGSize)) {
                nextMethod = QM_AIG;
            } else {
                //                skipNextBdd = false;
                nextMethod = QM_BDD;
            }

            /* update remaining variables */
            if (root.isConstant()) {
                remainingVariables.clear();
            } else {
                std::vector<Node*> supportnodes = root.getInternal().node()->supportNodes();
                std::set<Node*>    current(supportnodes.begin(), supportnodes.end());

                std::vector<int> newRemaining;

                for (const int var : remainingVariables) {
                    Node* v = variableInternal(var).node();

                    if (current.find(v) != current.end()) {
                        newRemaining.push_back(var);
                        current.erase(v);
                    }
                }

                std::swap(newRemaining, remainingVariables);
            }
        }
        /* building bdd succeeded -> perform bdd quantification */
        else {
            lastStepWasAIG = false;

            EdgeRef     deepResult;
            std::size_t deepSize            = 0;
            bool        canPerformDeepQuant = false;
            if (remainingVariables.size() == 1) {
                canPerformDeepQuant = true;
                InternalEdgeRef var = variableInternal(remainingVariables[0]);
                deepResult          = quantifyDeep(root, var.node()->varIndex(), existential);
                deepSize            = deepResult.nodeCount();
            }

            /* TODO: disable bdd limits */
            Cudd_SetMaxLive(_bddManager, 1000000000);

            /* create cube of remaining variables */
            DdNode* cube = nullptr;
            {
                auto cubevars = new int[remainingVariables.size()];

                int* cubevarsit = cubevars;
                for (auto p = remainingVariables.begin(); p != remainingVariables.end(); ++p, ++cubevarsit) {
                    *cubevarsit = *p;
                }

                cube = Cudd_IndicesToCube(_bddManager, cubevars, static_cast<int>(remainingVariables.size()));
                delete[] cubevars;

                assert(cube != nullptr);

                Cudd_Ref(cube);
            }

            /* perform actual quantification */
            DdNode* quantifiedbdd = nullptr;
            if (existential) {
                quantifiedbdd = Cudd_bddExistAbstract(_bddManager, bdd, cube);
            } else {
                quantifiedbdd = Cudd_bddUnivAbstract(_bddManager, bdd, cube);
            }

            assert(quantifiedbdd != nullptr);

            Cudd_Ref(quantifiedbdd);
            Cudd_RecursiveDeref(_bddManager, bdd);
            Cudd_RecursiveDeref(_bddManager, cube);

            ++bddquantifications;
            bddquantifiedvariables += remainingVariables.size();

            /* translate resulting bdd to aig */
            {
                InternalEdgeRef quantifiedaig = createAIG(quantifiedbdd);
                Cudd_RecursiveDeref(_bddManager, quantifiedbdd);
                root = EdgeRef(_extRefTable, quantifiedaig);
            }

            const std::size_t resultSize = root.nodeCount();

            if (canPerformDeepQuant && deepSize < resultSize) {
                root       = deepResult;
                nextMethod = QM_AIG;
            } else if ((remainingVariables.size() < 10) && (resultSize > 1000)
                       && (resultSize >= (1 << remainingVariables.size()) * currentAIGSize)) {
                nextMethod = QM_AIG;
            } else {
                nextMethod = QM_BDD;
            }

            /* clear remaining variables */
            remainingVariables.clear();
        }
    }
}

void
aigpp::Manager::hybridQuantify(aigpp::Manager::HybridQuantifyResult& baseFunction, aigpp::EdgeRef& conjunctiveFunction,
                               const std::vector<aigpp::EdgeRef>& variables, bool existential, bool optimizelast)
{
    static std::size_t bddquantifications     = 0;
    static std::size_t bddquantifiedvariables = 0;
    static std::size_t aigquantifications     = 0;

    /* setup bdd manager */
    // cleanBDD();

    enum QuantificationMethod
    {
        QM_AIG,
        QM_BDD
    };

    static QuantificationMethod nextMethod = QM_BDD;

    /*
     * skip aig optimizations after 2 unsuccessful optimizations in sequence were
     * made. skip at most 20 optimizations, abort skipping if aig has grown too
     * much (20%)
     *
     * use a maximum of 20 seconds for initial bdd computation
     */
    const std::size_t maxFailedOptimizations = 2;
    const std::size_t maxSkipOptimizations   = 20;
    const double      skipRatio              = 1.2;
    const std::size_t maxSkipFRAIGOnTimeout  = 10;
    const double      maxBddComputationTime  = 20.0;
    const std::size_t maxNodesForFraiging    = 50000;

    static std::size_t failedOptimizations      = 0;
    static std::size_t skipOptimizations        = 0;
    static std::size_t maxSkipAigSize           = 0;
    static std::size_t skipFRAIGOnTimeout       = 0;
    static std::size_t consecutiveFRAIGTimeouts = 0;

    /* get indices of variables that are to be quantified */

    std::set<int> variablesset;
    for (const auto& edge : variables) {
        assert(edge.isVariable());
        variablesset.insert(edge.variableIndex());
    }

    /* compute complete root function */
    baseFunction.aig = baseFunction.aig & conjunctiveFunction;
    if (nextMethod != QM_BDD) {
        if (baseFunction.bdd != nullptr) {
            Cudd_RecursiveDeref(_bddManager, baseFunction.bdd);
            baseFunction.bdd = nullptr;
        }
    }
    if (baseFunction.bdd == nullptr) {
        conjunctiveFunction.invalidate();
    }

    /* get support */
    std::vector<int> supportvec = baseFunction.aig.support();

    /* determine which variables (of root) have to be quantified */
    std::set<int> remainingVariables;
    for (const int p : supportvec) {
        if (variablesset.find(p) != variablesset.end()) {
            remainingVariables.insert(p);
        }
    }

    std::size_t startAIGSize   = baseFunction.aig.nodeCount();
    bool        lastStepWasAIG = false;

    while (!remainingVariables.empty()) {
        assert(!(baseFunction.aig.isConstant()));

        /* simple case: root is a variable */
        if (baseFunction.aig.isVariable()) {
#ifndef NDEBUG
            const int varIndex = baseFunction.aig.variableIndex();
            assert(remainingVariables.find(varIndex) != remainingVariables.end());
#endif

            /* delete baseFunction.bdd */
            if (baseFunction.bdd != nullptr) {
                Cudd_RecursiveDeref(_bddManager, baseFunction.bdd);
                baseFunction.bdd = nullptr;
            }

            if (existential) {
                baseFunction.aig = getConst1();
            } else {
                baseFunction.aig = getConst0();
            }
            remainingVariables.clear();
            break;
        }

        std::size_t currentAIGSize = baseFunction.aig.nodeCount();

        DdNode* bdd      = nullptr;
        bool    triedBDD = false;

        /* try to compute a bdd representation */
        if (nextMethod == QM_BDD) {
            triedBDD = true;

            /* try to build bdd */

            /* reset reordering limit */
            if ((std::size_t)Cudd_ReadNextReordering(_bddManager) > 6 * currentAIGSize) {
                Cudd_SetNextReordering(_bddManager, static_cast<int>(6 * currentAIGSize));
            }

            double BDDLimit
                = std::max(LOWER_BDD_LIMIT, std::min(UPPER_BDD_LIMIT, 20.0 * static_cast<double>(currentAIGSize)));
            Cudd_SetMaxLive(_bddManager, (int)BDDLimit);

            /* build bdd */
            Cudd_SetTimeLimit(_bddManager, static_cast<unsigned long>(maxBddComputationTime));

            if (baseFunction.bdd == nullptr) {
                bdd = buildBDDfromAIG(baseFunction.aig.getInternal().node());
                if (bdd != nullptr && baseFunction.aig.isInverted()) {
                    bdd = Cudd_Not(bdd);
                }
            } else {
                DdNode* conj_bdd = buildBDDfromAIG(conjunctiveFunction.getInternal().node());
                if (conj_bdd != nullptr && conjunctiveFunction.isInverted()) {
                    conj_bdd = Cudd_Not(conj_bdd);
                }
                conjunctiveFunction.invalidate();

                if (conj_bdd == nullptr) {
                    startAIGSize = currentAIGSize;
                    nextMethod   = QM_AIG;

                } else {
                    bdd = Cudd_bddAnd(_bddManager, baseFunction.bdd, conj_bdd);

                    if (bdd == nullptr) {
                        startAIGSize = currentAIGSize;
                        nextMethod   = QM_AIG;
                    } else {
                        Cudd_Ref(bdd);
                    }

                    Cudd_RecursiveDeref(_bddManager, conj_bdd);
                }

                Cudd_RecursiveDeref(_bddManager, baseFunction.bdd);
                baseFunction.bdd = nullptr;
            }

            Cudd_UnsetTimeLimit(_bddManager);
            Cudd_SetMaxLive(_bddManager, 1000000000);
        }

        /* check if expected bdd quantification result is larger that worst possible
         * aig quantification result */
        if (bdd != nullptr && remainingVariables.size() < 10) {
            int bddSize = Cudd_DagSize(bdd);

            if (static_cast<std::size_t>(bddSize) > ((1 << remainingVariables.size()) * currentAIGSize)) {
                Cudd_RecursiveDeref(_bddManager, bdd);
                bdd        = nullptr;
                nextMethod = QM_AIG;
            }
        }

        /* building bdd failed */
        if (bdd == nullptr) {
            if (!lastStepWasAIG) {
                lastStepWasAIG = true;
                startAIGSize   = currentAIGSize;

            } else if (triedBDD) {
                startAIGSize = static_cast<std::size_t>(1.1 * static_cast<double>(startAIGSize));
            }

            if (baseFunction.bdd != nullptr) {
                Cudd_RecursiveDeref(_bddManager, baseFunction.bdd);
                baseFunction.bdd = nullptr;
            }

            if (_bddManager != nullptr) {
                cuddGarbageCollect(_bddManager, true);
            }

            InternalEdgeRef var = variableInternal(*(remainingVariables.begin()));
            remainingVariables.erase(remainingVariables.begin());

            ++aigquantifications;

            /* deep quantification */
            baseFunction.aig = quantifyDeep(baseFunction.aig, var.node()->varIndex(), existential);

            var.clear();
            garbageCollect();

#ifdef DEBUG_BDD_QUANTIFY
            std::cout << "resulting_deep_aig_size: " << baseFunction.aig.nodeCount() << std::endl;
#endif

            /* optimization using rewriting and fraig */
            if ((optimizelast || !remainingVariables.empty()) && !(baseFunction.aig.isConstant())) {
                if (skipOptimizations > 0 && baseFunction.aig.nodeCount() < maxSkipAigSize) {
#ifdef DEBUG_BDD_QUANTIFY
                    std::cout << "skipping rewriting and fraig" << std::endl;
#endif
                    --skipOptimizations;
                } else {
                    std::size_t sizeBefore = baseFunction.aig.nodeCount();

#ifdef DEBUG_BDD_QUANTIFY
                    std::cout << "rewriting..." << std::flush;
#endif
                    rewrite();

#ifdef DEBUG_BDD_QUANTIFY
                    std::cout << " -> " << (baseFunction.aig.nodeCount() - sizeBefore) << std::endl;
#endif

                    if (nodeCount() < maxNodesForFraiging) {
                        if (skipFRAIGOnTimeout > 0) {
#ifdef DEBUG_BDD_QUANTIFY
                            std::cout << "skipping fraiging due to previous timeouts " << skipFRAIGOnTimeout
                                      << std::endl;
#endif
                            --skipFRAIGOnTimeout;
                        } else {
#ifdef DEBUG_BDD_QUANTIFY
                            std::cout << "fraiging... " << FRAIGING_TIMEOUT << std::flush;
#endif
                            std::size_t fraigSizeBefore = nodeCount();
                            makeFRAIG(FRAIGING_TIMEOUT);
                            std::size_t fraigSizeAfter = nodeCount();
#ifdef DEBUG_BDD_QUANTIFY
                            std::cout << " -> " << (fraigSizeAfter - fraigSizeBefore) << std::endl;
#endif

                            if (_timeoutDuringLastFRAIG) {
                                ++consecutiveFRAIGTimeouts;
                                skipFRAIGOnTimeout = consecutiveFRAIGTimeouts * maxSkipFRAIGOnTimeout;

#ifdef DEBUG_BDD_QUANTIFY
                                std::cout << "consecutive timeouts " << consecutiveFRAIGTimeouts << std::endl;
                                std::cout << "skipping next fraigings " << skipFRAIGOnTimeout << std::endl;
#endif
                            } else if (fraigSizeAfter >= fraigSizeBefore) {
                                skipFRAIGOnTimeout = maxSkipFRAIGOnTimeout;
#ifdef DEBUG_BDD_QUANTIFY
                                std::cout << "fraig fail" << std::endl;
                                std::cout << "skipping next fraigings " << skipFRAIGOnTimeout << std::endl;
#endif
                            } else {
                                consecutiveFRAIGTimeouts = 0;
                            }
                        }
                    }

#ifdef DEBUG_BDD_QUANTIFY
                    std::cout << "optimized_aig_size: " << baseFunction.aig.nodeCount() << std::endl;
#endif

                    std::size_t sizeAfter = baseFunction.aig.nodeCount();

                    if (sizeAfter >= static_cast<std::size_t>(0.99 * static_cast<double>(sizeBefore))) {
                        ++failedOptimizations;
                        if (failedOptimizations == maxFailedOptimizations) {
                            failedOptimizations = 0;
                            skipOptimizations   = maxSkipOptimizations;
                            maxSkipAigSize      = static_cast<std::size_t>(skipRatio * static_cast<double>(sizeAfter));
#ifdef DEBUG_BDD_QUANTIFY
                            std::cout << "skipping optimization steps" << std::endl;
#endif
                        }
                    } else {
                        failedOptimizations = 0;
                    }
                }
            }

#ifdef DEBUG_BDD_QUANTIFY
            std::cout << "currentsize " << baseFunction.aig.nodeCount() << " " << currentAIGSize << std::endl;
#endif

            if (static_cast<double>(baseFunction.aig.nodeCount()) <= 1.2 * static_cast<double>(startAIGSize)) {
#ifdef DEBUG_BDD_QUANTIFY
                std::cout << "\nnext bdd should be skipped" << std::endl;
#endif
                nextMethod = QM_AIG;
            } else {
                nextMethod = QM_BDD;
            }

            /* update remaining variables */
            if (baseFunction.aig.isConstant()) {
                remainingVariables.clear();
            } else {
                std::vector<Node*> supportnodes = baseFunction.aig.getInternal().node()->supportNodes();

                std::set<int> newremaining;

                for (Node* n : supportnodes) {
                    int varindex = n->varIndex();

                    if (remainingVariables.find(varindex) != remainingVariables.cend()) {
                        newremaining.insert(varindex);
                    }
                }

                std::swap(newremaining, remainingVariables);
            }

        }
        /* building bdd succeeded -> perform bdd quantification */
        else {
            lastStepWasAIG = false;

#ifdef DEBUG_BDD_QUANTIFY
            std::cout << "bdd succeeded" << std::endl;
            std::cout << "current_bdd_size: " << Cudd_DagSize(bdd) << std::endl;
#endif
            /* TODO: disable bdd limits */
            Cudd_SetMaxLive(_bddManager, 1000000000);

#ifdef DEBUG_BDD_QUANTIFY
            std::cout << "performing bdd based quantification" << std::endl;
#endif

            /* create cube of remaining variables */
            DdNode* cube = nullptr;
            {
                auto cubevars = new int[remainingVariables.size()];

                int* cubevarsit = cubevars;
                for (auto p = remainingVariables.begin(); p != remainingVariables.end(); ++p, ++cubevarsit) {
                    *cubevarsit = *p;
                }

                cube = Cudd_IndicesToCube(_bddManager, cubevars, static_cast<int>(remainingVariables.size()));
                delete[] cubevars;

                assert(cube != 0);

                Cudd_Ref(cube);
            }

            /* perform actual quantification */
            DdNode* quantifiedbdd = nullptr;
            if (existential) {
                quantifiedbdd = Cudd_bddExistAbstract(_bddManager, bdd, cube);
            } else {
                quantifiedbdd = Cudd_bddUnivAbstract(_bddManager, bdd, cube);
            }

            assert(quantifiedbdd != 0);

            Cudd_Ref(quantifiedbdd);
            Cudd_RecursiveDeref(_bddManager, bdd);
            Cudd_RecursiveDeref(_bddManager, cube);

            ++bddquantifications;
            bddquantifiedvariables += remainingVariables.size();

            /* copy bdd to baseFunction.bdd */
            baseFunction.bdd = quantifiedbdd;
            Cudd_Ref(baseFunction.bdd);

            /* translate resulting bdd to aig */
            {
                InternalEdgeRef quantifiedaig = createAIG(quantifiedbdd);
                Cudd_RecursiveDeref(_bddManager, quantifiedbdd);
                baseFunction.aig = EdgeRef(_extRefTable, quantifiedaig);
            }

            const std::size_t resultSize = baseFunction.aig.nodeCount();

            if ((remainingVariables.size() < 10) && (resultSize > 1000)
                && (resultSize >= (1 << remainingVariables.size()) * currentAIGSize)) {
                nextMethod = QM_AIG;
            } else {
                nextMethod = QM_BDD;
            }

            /* clear remaining variables */
            remainingVariables.clear();
        }
    }

#ifdef DEBUG_BDD_QUANTIFY
    std::cout << "QUANT " << aigquantifications << " " << bddquantifications << " " << bddquantifiedvariables
              << std::endl;
#endif

    if (nextMethod != QM_BDD && baseFunction.bdd != nullptr) {
        Cudd_RecursiveDeref(_bddManager, baseFunction.bdd);
        baseFunction.bdd = nullptr;
    }
}

void
aigpp::Manager::hybridQuantify(aigpp::Manager::HybridQuantifyResult& baseFunction,
                               const std::vector<aigpp::EdgeRef>& variables, bool existential, bool optimizelast)
{
    EdgeRef conjunctiveFunction = getConst1();
    hybridQuantify(baseFunction, conjunctiveFunction, variables, existential, optimizelast);
}

int
aigpp::Manager::selectBestVariableToQuantify(const std::set<int>& variables, const aigpp::EdgeRef& root)
{
    if (variables.size() == 1) {
        return *(variables.begin());
    }

    int         bestVar   = -1;
    std::size_t bestSize  = 0;
    std::size_t worstSize = 0;

    for (int v : variables) {
        EdgeRef     f = quantifyDeep(root, v, true);
        std::size_t n = f.nodeCount();

        if (bestVar == -1 || n < bestSize) {
            bestVar  = v;
            bestSize = n;
        }

        if (bestVar == -1 || n > worstSize) {
            worstSize = n;
        }
    }

    return bestVar;
}

template <template <typename> class P = std::less>
struct compare_pair_second
{
    template <class T1, class T2>
    bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right)
    {
        return P<T2>()(left.second, right.second);
    }
};

std::vector<int>
aigpp::Manager::computeVariableOrdering(const std::set<int>& variables, const aigpp::EdgeRef& root)
{
    std::vector<int> result;
    result.reserve(variables.size());

    if (variables.size() == 1) {
        result.push_back(*(variables.begin()));
        return result;
    }

    std::vector<std::pair<int, std::size_t>> temp;

    for (const int v : variables) {
        EdgeRef           f = quantifyDeep(root, v, true);
        const std::size_t n = f.nodeCount();

        temp.emplace_back(v, n);
    }

    std::sort(temp.begin(), temp.end(), compare_pair_second<>());

    for (const auto& p : temp) {
        result.push_back(p.first);
    }

    return result;
}

aigpp::EdgeRef
aigpp::Manager::hybridQuantify2(aigpp::EdgeRef root, const std::vector<aigpp::EdgeRef>& variables)
{
    const double maxBDDComputationTime = 20.0;

    static std::size_t bddquantifications     = 0;
    static std::size_t bddquantifiedvariables = 0;
    static std::size_t aigquantifications     = 0;

    std::size_t startAIGSize   = 1;
    bool        lastStepWasAIG = false;

    /* determine which variables (of root) have to be quantified */
    std::set<int> remainingVariables;
    {
        std::set<int> variablesset;
        for (const auto& p : variables) {
            assert(p.isVariable());
            variablesset.insert(p.variableIndex());
        }

        /* get root's support */
        std::vector<int> supportvec = root.support();
        for (const int p : supportvec) {
            if (variablesset.find(p) != variablesset.cend()) {
                remainingVariables.insert(p);
            }
        }
    }

    /* setup bdd manager */
    cleanBDD();

    enum QuantificationMethod
    {
        QM_AIG,
        QM_BDD
    };

    std::vector<int> aigVariableOrder;

    static QuantificationMethod nextMethod = QM_BDD;

    while (!remainingVariables.empty()) {
        assert(!(root.isConstant()));

        std::cout << "remaining variables to be quantified: " << remainingVariables.size() << std::endl;
        std::cout << "current_aig_size: " << root.nodeCount() << std::endl;
        std::cout << "global_aig_size:  " << nodeCount() << std::endl;

        /* simple case: root is a variable */
        if (root.getInternal().node()->isVar()) {
#ifndef NDEBUG
            const int varIndex = root.getInternal().node()->varIndex();
            assert(remainingVariables.find(varIndex) != remainingVariables.end());
#endif
            root = getConst1();

            remainingVariables.clear();
            break;
        }

        std::size_t currentAIGSize = root.nodeCount();

        DdNode* bdd = nullptr;

        /* try to compute a bdd representation */
        if (nextMethod == QM_BDD) {
            if (static_cast<std::size_t>(Cudd_ReadNextReordering(_bddManager)) > 6 * currentAIGSize) {
                std::cout << "next reordering " << Cudd_ReadNextReordering(_bddManager) << " -> " << 6 * currentAIGSize
                          << std::endl;
                Cudd_SetNextReordering(_bddManager, static_cast<int>(6 * currentAIGSize));
            }

            double BDDLimit
                = std::max(LOWER_BDD_LIMIT, std::min(UPPER_BDD_LIMIT, 20.0 * static_cast<double>(currentAIGSize)));
            std::cout << "using a node limit of " << BDDLimit << std::endl;
            Cudd_SetMaxLive(_bddManager, static_cast<int>(BDDLimit));

            Cudd_SetTimeLimit(_bddManager, maxBDDComputationTime);
            bdd = buildBDDfromAIG(root.getInternal().node());
            Cudd_UnsetTimeLimit(_bddManager);

            if (bdd != nullptr && root.isInverted()) {
                bdd = Cudd_Not(bdd);
            }

            Cudd_SetMaxLive(_bddManager, 1000000000);

            if (bdd == nullptr) {
                std::cout << "bdd failed due to limits\n";

                startAIGSize = currentAIGSize;

                nextMethod = QM_AIG;
            }
        } else {
            std::cout << "bdd skipped (1)\n";
        }

        if (bdd == nullptr) {
            if (!lastStepWasAIG) {
                lastStepWasAIG = true;
                startAIGSize   = currentAIGSize;

#ifdef DEBUG_BDD_QUANTIFY
                std::cout << "last step was bdd -> resetting start aig size to " << startAIGSize << std::endl;
#endif
            }

            if (_bddManager != nullptr) {
                cuddGarbageCollect(_bddManager, true);
            }

            /* select var to quantify */
            // int bestVar = selectBestVariableToQuantify( remainingVariables, root );
            int bestVar = -1;
            {
                if (aigVariableOrder.empty()) {
                    aigVariableOrder = computeVariableOrdering(remainingVariables, root);
                }

                for (const int v : aigVariableOrder) {
                    if (remainingVariables.find(v) != remainingVariables.cend()) {
                        bestVar = v;
                        break;
                    }
                }
            }
            assert(bestVar != -1);

            ++aigquantifications;

            /* deep quantification */
            root = quantifyDeep(root, bestVar, true);

            root = _rwManager.localRewriteCone(root);

            // makeFRAIG( maxFraigingTime );

            std::size_t sizeAfter = root.nodeCount();

            if (static_cast<double>(sizeAfter) <= 1.2 * static_cast<double>(startAIGSize)) {
                nextMethod = QM_AIG;
            } else {
                nextMethod = QM_BDD;
            }

            /* update remaining variables */
            if (root.isConstant()) {
                remainingVariables.clear();
            } else {
                std::vector<int> supportvec = root.support();
                std::set<int>    newremaining;
                for (const int p : supportvec) {
                    if (remainingVariables.find(p) != remainingVariables.end()) {
                        newremaining.insert(p);
                    }
                }

                std::swap(newremaining, remainingVariables);
            }
        }
        /* building bdd succeeded -> perform bdd quantification */
        else {
            lastStepWasAIG = false;

            Cudd_ReduceHeap(_bddManager, CUDD_REORDER_SIFT, 0);
            Cudd_SetMaxLive(_bddManager, 1000000000);

            // performing bdd based quantification

            /* create cube of remaining variables */
            DdNode* cube = nullptr;
            {
                auto cubevars = new int[remainingVariables.size()];

                int* cubevarsit = cubevars;
                for (auto p = remainingVariables.begin(); p != remainingVariables.end(); ++p, ++cubevarsit) {
                    *cubevarsit = *p;
                }

                cube = Cudd_IndicesToCube(_bddManager, cubevars, static_cast<int>(remainingVariables.size()));
                delete[] cubevars;

                assert(cube != 0);

                Cudd_Ref(cube);
            }

            /* perform actual quantification */
            DdNode* quantifiedbdd = Cudd_bddExistAbstract(_bddManager, bdd, cube);
            assert(quantifiedbdd != 0);

            Cudd_Ref(quantifiedbdd);
            Cudd_RecursiveDeref(_bddManager, bdd);
            Cudd_RecursiveDeref(_bddManager, cube);

            ++bddquantifications;
            bddquantifiedvariables += remainingVariables.size();

            Cudd_ReduceHeap(_bddManager, CUDD_REORDER_SIFT, 0);

            /* translate resulting bdd to aig */
            {
                InternalEdgeRef quantifiedaig = createAIG(quantifiedbdd);
                Cudd_RecursiveDeref(_bddManager, quantifiedbdd);
                root = EdgeRef(_extRefTable, quantifiedaig);
            }

            const std::size_t resultSize = root.nodeCount();
            if ((remainingVariables.size() < 10) && (resultSize > 1000)
                && (resultSize >= (1 << remainingVariables.size()) * currentAIGSize)) {
                nextMethod = QM_AIG;
            } else {
                nextMethod = QM_BDD;
            }

            /* clear remaining variables */
            remainingVariables.clear();
        }

        printStats();
    }

    return root;
}
