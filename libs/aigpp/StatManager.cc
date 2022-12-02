/**************************************************************
 *
 *       AIGPP Package // StatManager.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *	Author:
 *         Florian Pigorsch
 *	  University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *	Last revision:
 *         $Revision: 377 $
 *	  $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "StatManager.hh"

aigpp::StatManager&
aigpp::StatManager::instance()
{
    static StatManager m;
    return m;
}

aigpp::StatManager::StatManager() :
    _maxNodes(0),
    _maxDead(0),
    _nodesCreated(0),
    _nodesDeleted(0),
    _refs(0),
    _derefs(0),
    _reclaimed(0),
    _gc(0),
    _gcTime(0.0),

    _uniqueHits(0),
    _computedHitsAnd(0),
    _computedHitsCofactor(0),
    _computedMissAnd(0),
    _computedMissCofactor(0),
    _computedInsertions(0),
    _computedCollisions(0),

    _and(0),
    _andConstant(0),
    _andTrivial1(0),
    _andTrivial15(0),
    _andTrivial2a(0),
    _andTrivial2b(0),
    _andStrashing(0),
    _andKnown(0),
    _andSAT(0),
    _andSATZero(0),
    _andNew(0),

    _replacements(0),
    _replacementsInverted(0),
    _replacementTime(0.0),

    _cofactor(0),
    _compose(0),
    _quantification(0),

    _satChecks(0),
    _satEquiv(0),
    _satResets(0),
    _satMaxVars(0),
    _satTime(0.0),
    _satMaxTime(0.0),
    _satCreationTime(0.0),

    _simZero(0),
    _simRandom(0),

    _simCounterEx(0),
    _simCounterExMerges(0),
    _simUpdates(0),
    _simUpdateTime(0.0),
    _simZeroUpdateTime(0.0),
    _simTempUpdateTime(0.0),

    _bddCalls(0),
    _bddSuccess(0),
    _bddSkipped(0),
    _bddLimit(0),
    _bddFailed(0),
    _bddTime(0.0),
    _bddAIGTime(0.0),

    _redChecksCombined(0),
    _redChecksSingle(0),
    _redChecksDetected(0),
    _redCheckTime(0.0)
{}

std::ostream&
aigpp::operator<<(std::ostream& os, const aigpp::StatManager& m)
{
    os << "aigpp package statistics\n"
       << "========================\n\n"

       << "node creation\n"
       << "-------------\n"
       << "  max nodes: " << m._maxNodes << '\n'
       << "  max dead: " << m._maxDead << '\n'
       << "  created: " << m._nodesCreated << '\n'
       << "  deleted: " << m._nodesDeleted << '\n'
       << "  referenced: " << m._refs << '\n'
       << "  dereferenced: " << m._derefs << '\n'
       << "  reclaimed: " << m._reclaimed << '\n'
       << "  garbage collections: " << m._gc << '\n'
       << "  gc time: " << m._gcTime << '\n'
       << '\n'

       << "caches\n"
       << "------\n"
       << "  unique table hits: " << m._uniqueHits << '\n'
       << "  computed table hits (and): " << m._computedHitsAnd << '\n'
       << "  computed table hits (cofactor): " << m._computedHitsCofactor << '\n'
       << "  computed table misses (and): " << m._computedMissAnd << '\n'
       << "  computed table misses (cofactor): " << m._computedMissCofactor << '\n'
       << "  computed table insertions: " << m._computedInsertions << '\n'
       << "  computed table collisions: " << m._computedCollisions << '\n'
       << '\n'

       << "and operations\n"
       << "--------------\n"
       << "  and calls: " << m._and << '\n'
       << "  constant operand: " << m._andConstant << '\n'
       << "  trivial operands 1: " << m._andTrivial1 << '\n'
       << "  trivial operands 15: " << m._andTrivial15 << '\n'
       << "  trivial operands 2a: " << m._andTrivial2a << '\n'
       << "  trivial operands 2b: " << m._andTrivial2b << '\n'
       << "  result by strashing: " << m._andStrashing << '\n'
       << "  result known: " << m._andKnown << '\n'
       << "  result by SAT: " << m._andSAT << '\n'
       << "  result by SAT (zero):  " << m._andSATZero << '\n'
       << "  new node: " << m._andNew << '\n'
       << '\n'

       << "sat\n"
       << "---\n"
       << "  checks: " << m._satChecks << '\n'
       << "  equivalence: " << m._satEquiv << '\n'
       << "  resets: " << m._satResets << '\n'
       << "  max vars: " << m._satMaxVars << '\n'
       << "  time: " << m._satTime << '\n'
       << "  max time: " << m._satMaxTime << '\n'
       << "  creation time: " << m._satCreationTime << '\n'
       << '\n'

       << "simulation\n"
       << "----------\n"
       << "  zero: " << m._simZero << '\n'
       << "  random added: " << m._simRandom << '\n'
       << "  counterex: " << m._simCounterEx << '\n'
       << "  counterex merges: " << m._simCounterExMerges << '\n'
       << "  updates: " << m._simUpdates << '\n'
       << "  update time: " << m._simUpdateTime << '\n'
       << "  zero update time: " << m._simZeroUpdateTime << '\n'
       << "  temp sim update time: " << m._simTempUpdateTime << '\n'
       << '\n'

       << "node replacements\n"
       << "-----------------\n"
       << "  positive: " << m._replacements << '\n'
       << "  inverted: " << m._replacementsInverted << '\n'
       << "  time: " << m._replacementTime << '\n'
       << '\n'

       << "bdd sweeping\n"
       << "------------\n"
       << "  bdd calls: " << m._bddCalls << '\n'
       << "  success: " << m._bddSuccess << '\n'
       << "  skipped: " << m._bddSkipped << '\n'
       << "  limit: " << m._bddLimit << '\n'
       << "  bdd time: " << m._bddTime << '\n'
       << "  aig time: " << m._bddAIGTime << '\n'
       << '\n'

       << "redundancy checks\n"
       << "-----------------\n"
       << "  combined: " << m._redChecksCombined << '\n'
       << "  single: " << m._redChecksSingle << '\n'
       << "  detected: " << m._redChecksDetected << '\n'
       << "  time: " << m._redCheckTime << '\n'
       << '\n'

       << "highlevel operations\n"
       << "--------------------\n"
       << "  cofactor: " << m._cofactor << '\n'
       << "  compose: " << m._compose << '\n'
       << "  quantification: " << m._quantification << '\n'
       << '\n';

    return os;
}
