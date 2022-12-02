/**************************************************************
 *
 *       AIGPP Package // Manager_sim.cc
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
#include <stack>

/* LRABS utilities */
#include <lrabsutil/Math.hh>
#include <lrabsutil/Resources.hh>
#include <lrabsutil/String.hh>

/* local */
#include "StatManager.hh"

void
aigpp::Manager::recomputeSimulation(bool warnings)
{
    assert(simCreation() == true);
    _simTable.clear();

    if (!warnings) {
        for (Node* n = _nodes; n != nullptr; n = n->next()) {
            if (n->isVar()) {
                continue;
            }

            n->updateSim();
            _simTable.insert(n);
        }
    } else {
        unsigned int diffCount = 0;

        SimVector s;

        for (Node* n = _nodes; n != nullptr; n = n->next()) {
            if (n->isVar()) {
                continue;
            }

            s = n->sim();

            n->updateSim();

            if (s != n->sim()) {
                ++diffCount;
            }

            _simTable.insert(n);
        }

        if (diffCount != 0) {
            std::cout << "WARNING: nodes with inconsistent simvectors " << diffCount << std::endl;
        }
    }
}

void
aigpp::Manager::updateSimTable(bool force)
{
    assert(simCreation() == true);

    /* update only if there are enaugh new sim vectors */
    if (!force && _newCounterExamples.size() < SimVector::BinSize) {
        return;
    }

    double startTime = lrabs::cpuTime();

    /* clear sim table */
    //_simTable.clear();
    _simTable.clearRememberClasses();

    while ((int)_newCounterExamples.size() >= (int)SimVector::BinSize) {
        auto addcounterex = (int)SimVector::BinSize;

        /* select counter examples to add */
        std::vector<VariableAssignment> counterex(addcounterex);
        for (int i = 0; i != addcounterex; ++i) {
            counterex[i] = _newCounterExamples.front();
            _counterExamples.push_back(_newCounterExamples.front());
            _newCounterExamples.pop_front();
        }

        /* erase counter examples if there are too many */
        while (_counterExamples.size() > SimVector::BinSize * SimVector::BinCount) {
            _counterExamples.pop_front();
        }

        /* update sim vectors of variables */
        for (Node* n : _variables) {
            auto               newcex = counterex.begin();
            SimVector::BinType bin    = 0ul;

            for (int cex = 0; cex != addcounterex; ++cex) {
                /* is this variable present in the counter example? */
                VariableAssignment::Assignment a = (*newcex).get(n->varIndex());

                /* occurs positive */
                if (a == VariableAssignment::Positive) {
                    bin |= (1ul << cex);
                }
                /* occurs negative */
                else if (a == VariableAssignment::Negative) {
                    bin &= ~(1ul << cex);
                }
                /* does not occur -> select random value */
                else {
                    if (_random.getBool()) {
                        bin |= (1ul << cex);
                    } else {
                        bin &= ~(1ul << cex);
                    }
                }

                ++newcex;
            }

            /* insert new bin into the sim vector */
            n->sim().setBin(_nextUpdateSimVectorBin, bin);
        }

        /* update sim vectors of the remaining nodes using the topological order of
         * the nodes list */
        for (Node* n = _nodes; n != nullptr; n = n->next()) {
            if (n->isVar()) {
                continue;
            }
            n->updateSim(_nextUpdateSimVectorBin);
        }

        ++_nextUpdateSimVectorBin;
        if (_nextUpdateSimVectorBin >= SimVector::BinCount) {
            _nextUpdateSimVectorBin = 0;
        }
    }

    /*
    for( Node* n = _nodes; n != 0; n = n->next() )
    {
        _simTable.insert( n );
    }
    */
    _simTable.rebuild(_nodes, true);

    StatManager::instance().incSimUpdates();
    StatManager::instance().incSimUpdateTime(lrabs::cpuTime() - startTime);
}

void
aigpp::Manager::updateSimTable(const aigpp::VariableAssignment& assignment)
{
    assert(simCreation() == true);

    double startTime = lrabs::cpuTime();

    /*  clear sim table */
    //_simTable.clear();
    _simTable.clearRememberClasses();

    const int binIndex = 0;

    /* update sim vectors of variables */
    for (Node* n : _variables) {
        SimVector::BinType bin = n->sim()[_nextUpdateSimVectorBin];

        // is this variable present in the counter example?
        VariableAssignment::Assignment a = assignment.get(n->varIndex());

        if (a == VariableAssignment::Positive) {
            bin |= (1ul << binIndex);
        } else if (a == VariableAssignment::Negative) {
            bin &= ~(1ul << binIndex);
        } else {
            if (_random.getBool()) {
                bin |= (1ul << binIndex);
            } else {
                bin &= ~(1ul << binIndex);
            }
        }

        n->sim().setBin(_nextUpdateSimVectorBin, bin);
    }

    /* update sim vectors of the remaining nodes using the topological order */
    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->isVar()) {
            continue;
        }
        n->updateSim(_nextUpdateSimVectorBin);
    }

    //    for( Node* n = _nodes; n != 0; n = n->next() )
    //    {
    //        _simTable.insert( n );
    //    }
    _simTable.rebuild(_nodes, true);

    StatManager::instance().incSimZeroUpdateTime(lrabs::cpuTime() - startTime);
}

void
aigpp::Manager::addCounterExample(const aigpp::VariableAssignment& counterex)
{
    StatManager::instance().incSimCounterEx();
    _newCounterExamples.push_back(counterex);

    if (settings().getSimVectorMerging()) {
        mergeCounterExamples();
    }
}

#ifdef USE_TEMPSIM
void
aigpp::Manager::updateTempSim(const aigpp::VariableAssignment& assignment, aigpp::Node* additionalNode)
{
    double startTime = lrabs::cpuTime();

    /* update sim vectors of variables */
    for (Node* n : _variables) {
        // is this variable present in the counter example?
        const VariableAssignment::Assignment a = assignment.get(n->varIndex());

        if (a == VariableAssignment::Positive) {
            n->setTempSim(~0ul);
        } else if (a == VariableAssignment::Negative) {
            n->setTempSim(0ul);
        } else {
            n->setTempSim(_randomTempSim.getULong());
        }
    }

    StatManager::instance().incTempSimUpdateTime(lrabs::cpuTime() - startTime);

    _tempSimUpToDate = false;
    propagateTempSim(additionalNode);
}

void
aigpp::Manager::propagateTempSim(aigpp::Node* additionalNode)
{
    if (!_tempSimUpToDate) {
        double startTime = lrabs::cpuTime();

        /* update sim vectors of the remaining nodes using the topological order */
        for (Node* n = _nodes; n != nullptr; n = n->next()) {
            if (n->isVar()) {
                continue;
            }
            n->updateTempSim();
        }

        if (additionalNode != nullptr) {
            additionalNode->updateTempSim();
        }

        _tempSimUpToDate = true;

        StatManager::instance().incTempSimUpdateTime(lrabs::cpuTime() - startTime);
    }
}

#endif

void
aigpp::Manager::addRandomCounterExample()
{
    addCounterExample(VariableAssignment());
}

void
aigpp::Manager::mergeCounterExamples()
{
    if (_newCounterExamples.size() <= 1) {
        return;
    }

    bool merged = true;
    while (merged) {
        merged = false;

        for (auto c1 = _newCounterExamples.begin(); c1 != _newCounterExamples.end(); ++c1) {
            for (auto c2 = c1 + 1; c2 != _newCounterExamples.end(); ++c2) {
                if (c1->compatible(*c2)) {
                    StatManager::instance().incSimCounterExMerges();

                    c1->mergeWith(*c2);

                    _newCounterExamples.erase(c2);
                    merged = true;
                    break;
                }
            }
            if (merged) {
                break;
            }
        }
    }
}
