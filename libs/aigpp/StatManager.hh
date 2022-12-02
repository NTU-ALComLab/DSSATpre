/**************************************************************
 *
 *       AIGPP Package // StatManager.hh
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
 *
 ***************************************************************/

#ifndef AIGPP_STATMANAGER_HH
#define AIGPP_STATMANAGER_HH

#include <iostream>

namespace aigpp {
class StatManager
{
   public:
    StatManager(const StatManager&) = delete;
    StatManager(StatManager&&)      = delete;
    StatManager& operator=(const StatManager&) = delete;
    StatManager& operator=(StatManager&&) = delete;

    static StatManager& instance();

    void updateMaxNodes(std::size_t n)
    {
        if (n > _maxNodes) _maxNodes = n;
    };
    void updateMaxDead(std::size_t n)
    {
        if (n > _maxDead) _maxDead = n;
    };
    void incCreated() { ++_nodesCreated; };
    void incDeleted() { ++_nodesDeleted; };
    void incRefs() { ++_refs; };
    void incDerefs() { ++_derefs; };
    void incReclaimed() { ++_reclaimed; };

    void incGarbageCollections() { ++_gc; };
    void incGarbageCollectionTime(double t) { _gcTime += t; };

    void incUniqueHits() { ++_uniqueHits; };
    void incComputedHitsAnd() { ++_computedHitsAnd; };
    void incComputedHitsCofactor() { ++_computedHitsCofactor; };
    void incComputedMissAnd() { ++_computedMissAnd; };
    void incComputedMissCofactor() { ++_computedMissCofactor; };
    void incComputedInsertions() { ++_computedInsertions; };
    void incComputedCollisions() { ++_computedCollisions; };

    void incAnd() { ++_and; };
    void incAndConstant() { ++_andConstant; };
    void incAndTrivial1() { ++_andTrivial1; };
    void incAndTrivial15() { ++_andTrivial15; };
    void decAndTrivial15() { --_andTrivial15; };
    void incAndTrivial2a() { ++_andTrivial2a; };
    void decAndTrivial2a() { --_andTrivial2a; };
    void incAndTrivial2b() { ++_andTrivial2b; };
    void decAndTrivial2b() { --_andTrivial2b; };
    void incAndStrashing() { ++_andStrashing; };
    void incAndKnown() { ++_andKnown; };
    void incAndSAT() { ++_andSAT; };
    void incAndSATZero() { ++_andSATZero; };
    void incAndNew() { ++_andNew; };

    void incReplacements() { ++_replacements; };
    void incReplacementsInverted() { ++_replacementsInverted; };
    void incReplacementTime(double t) { _replacementTime += t; };

    void incCofactor() { ++_cofactor; };
    void incCompose() { ++_compose; };
    void incQuantification() { ++_quantification; };

    void incSATChecks() { ++_satChecks; };
    void incSATEquiv() { ++_satEquiv; };
    void incSATResets() { ++_satResets; };
    void updateSATMaxVars(std::size_t v)
    {
        if (v > _satMaxVars) _satMaxVars = v;
    };
    void incSATTime(double t)
    {
        _satTime += t;
        if (t > _satMaxTime) _satMaxTime = t;
    };
    void incSATCreationTime(double t) { _satCreationTime += t; };

    void incSimZero() { ++_simZero; };
    void incSimRandomAdded() { ++_simRandom; };

    void incSimCounterEx() { ++_simCounterEx; };
    void incSimCounterExMerges() { ++_simCounterExMerges; };
    void incSimUpdates() { ++_simUpdates; };
    void incSimUpdateTime(double t) { _simUpdateTime += t; };
    void incTempSimUpdateTime(double t) { _simTempUpdateTime += t; };
    void incSimZeroUpdateTime(double t) { _simZeroUpdateTime += t; };

    void incBDDCalls() { ++_bddCalls; };
    void incBDDSuccess() { ++_bddSuccess; };
    void incBDDSkipped() { ++_bddSkipped; };
    void incBDDLimit() { ++_bddLimit; };
    void incBDDFailed() { ++_bddFailed; };
    void incBDDTime(double t) { _bddTime += t; };
    void incBDDAIGTime(double t) { _bddAIGTime += t; };

    void incRedundancyChecksCombined() { ++_redChecksCombined; };
    void incRedundancyChecksSingle(std::size_t n) { _redChecksSingle += n; };
    void incRedundanciesDetected(std::size_t n) { _redChecksDetected += n; };
    void incRedundancyCheckTime(double t) { _redCheckTime += t; };

   protected:
    /* methods */
    StatManager();

    /* variables */
    std::size_t _maxNodes, _maxDead;
    std::size_t _nodesCreated, _nodesDeleted;
    std::size_t _refs, _derefs;
    std::size_t _reclaimed;
    std::size_t _gc;
    double      _gcTime;

    std::size_t _uniqueHits;
    std::size_t _computedHitsAnd, _computedHitsCofactor;
    std::size_t _computedMissAnd, _computedMissCofactor;
    std::size_t _computedInsertions, _computedCollisions;

    std::size_t _and, _andConstant, _andTrivial1, _andTrivial15, _andTrivial2a, _andTrivial2b, _andStrashing, _andKnown,
        _andSAT, _andSATZero, _andNew;

    std::size_t _replacements, _replacementsInverted;
    double      _replacementTime;

    std::size_t _cofactor, _compose, _quantification;

    std::size_t _satChecks, _satEquiv, _satResets, _satMaxVars;
    double      _satTime, _satMaxTime, _satCreationTime;

    std::size_t _simZero, _simRandom;
    std::size_t _simCounterEx, _simCounterExMerges;
    std::size_t _simUpdates;
    double      _simUpdateTime;
    double      _simZeroUpdateTime;
    double      _simTempUpdateTime;

    std::size_t _bddCalls, _bddSuccess, _bddSkipped, _bddLimit, _bddFailed;
    double      _bddTime, _bddAIGTime;

    std::size_t _redChecksCombined, _redChecksSingle, _redChecksDetected;
    double      _redCheckTime;

    /* friends */
    friend class Manager;
    friend std::ostream& operator<<(std::ostream& os, const StatManager& m);
};

std::ostream& operator<<(std::ostream& os, const StatManager& m);
}  // namespace aigpp

#endif /* AIGPP_STATMANAGER_HH */
