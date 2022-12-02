/****************************************************************************************[Solver.C]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include "Solver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

#include <lrabsutil/Resources.hh>

#include "Sort.h"

// #define DEBUG_ANALYZE
// #define ONIONING

//=================================================================================================
// Constructor/Destructor:

namespace minicraig {

    Solver::Solver() :
#ifdef CRAIG_INTERPOLATION
        // Craig Interpolation Stuff
        no_more_b_clauses_(false),
        var_types_(),
        m_(),
        global_vars_(),
        max_craig_var_id_(0),
        partial_inter_on_dl_0_(),
        failed_( 0 ),
        maxAIGNodeCount_( 0 ),
        maxInterpolationTime_( 0 ),
        timeoutTime_( 0 ),
        rev_var_map_(NULL),
        var_map_(NULL),
        mapping_(),
        final_interpolants_(),
        approx_number_(0),
        status_(NULL),
#endif
#ifdef MY_DECISION
        pseudo_occurence_(),
#endif

        // Parameters: (formerly in 'SearchParams')
        var_decay(1 / 0.95), clause_decay(1 / 0.999), random_var_freq(0.02)
        , restart_first(100), restart_inc(1.5), learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

        // More parameters:
        //
        , expensive_ccmin  (false)
#if !defined(MY_DECISION) && !defined(MINISAT_Z09)
        , polarity_mode    (polarity_false)
#endif
#ifdef MY_DECISION
        , polarity_mode  (polarity_occ)
#endif
#ifdef MINISAT_Z09
        , polarity_mode  (polarity_stored)
#endif
        //, polarity_mode    (polarity_true)
        //, polarity_mode    (polarity_rnd)
        , verbosity        (0)

        // Statistics: (formerly in 'SolverStats')
        //
        , starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
        , clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

        , ok               (true)
        , cla_inc          (1)
        , var_inc          (1)
#ifdef CRAIG_INTERPOLATION
        , a_               (-1)
#endif
        , qhead            (0)
        , simpDB_assigns   (-1)
        , simpDB_props     (0)
        , order_heap       (VarOrderLt(activity))
        , random_seed      (91648253)
        , progress_estimate(0)
        , remove_satisfied (false)
    {
    }


    Solver::~Solver()
    {
        for (int i = 0; i < learnts.size(); i++) free(learnts[i]);
        for (int i = 0; i < clauses.size(); i++) free(clauses[i]);
    }


    //=================================================================================================
    // Minor methods:


    // Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
    // used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
    //
    Var Solver::newVar(bool sign, bool dvar)
    {
        int v = nVars();

        watches   .push();          // (list for positive literal)
        watches   .push();          // (list for negative literal)
        reason    .push(NULL);
        assigns   .push(toInt(l_Undef));
        level     .push(-1);
#ifdef CRAIG_INTERPOLATION
        assignment.push(-1);
#endif
        activity  .push(0);
        seen      .push(0);

#ifndef MINISAT_Z09
        polarity    .push((char)sign);
#endif
#ifdef MINISAT_Z09
        polarity.push((char)sign);
#endif
        decision_var.push((char)dvar);

        insertVarOrder(v);
        return v;
    }
#ifndef CRAIG_INTERPOLATION
    bool Solver::addClause(vec<Lit>& ps)
#endif
#ifdef CRAIG_INTERPOLATION
        bool Solver::addClause(vec<Lit>& ps,
                               ClauseType cl_type,
                               ClauseInfo info,
                               bool pure,
                               bool take_inter,
                               const SimpleAIGEdge inter )
#endif
    {
        assert(decisionLevel() == 0);
#ifdef CRAIG_INTERPOLATION
        // here we update the variable types e.g. we
        // set for all variables ocurring in any clauses
        // the corresponding variable type (A_LOCAL or B_LOCAL)
        // if one variable is both we know that the variable
        // is GLOBAL
        for( int t = 0 ; t < ps.size() ; ++t){
            if(cl_type == A_CLAUSE)
                setVarType( var(ps[t]) , A_LOCAL);
            else if(cl_type == B_CLAUSE)
                setVarType( var(ps[t]) , B_LOCAL);
            else{
                assert( take_inter );
            }
        }
#endif



#ifdef CRAIG_INTERPOLATION
        if( cl_type == A_CLAUSE && (!pure) ){
            no_more_b_clauses_ = true;
        }
        assert( (!no_more_b_clauses_) || cl_type == A_CLAUSE || pure || take_inter);

        SimpleAIGEdge partial_interpolant = m_.getFalse();

        if( take_inter ){
            partial_interpolant = inter;
        }
        else{
            if( cl_type == B_CLAUSE ){
                partial_interpolant = m_.getTrue();
            }
            else if( cl_type == A_CLAUSE && pure ){
                partial_interpolant = m_.getFalse();
            }
            else{
                std::vector<SimpleAIGEdge> lits_;
                for( int t = 0 ; t < ps.size() ; ++t){
                    if ( getVarType(ps[t]) == GLOBAL ){
                        lits_.push_back(getCraigVar(ps[t]));
                    }
                }
                if( lits_.size() > 0){
                    partial_interpolant = m_.MultiOr( lits_ );
                }
            }
        }
#endif

        if (!ok){
            // nothiing to do (we have already a valid interpolant)
            return false;
        }
        else{
            // Check if clause is satisfied and remove false/duplicate literals:
            sort(ps);
            Lit p; int i, j;
            for (i = j = 0, p = lit_Undef; i < ps.size(); i++){
                if (value(ps[i]) == l_True || ps[i] == ~p){
                    return true;
                }
                else if ( value(ps[i]) != l_False && ps[i] != p){
                    ps[j++] = p = ps[i];
                }
#ifdef CRAIG_INTERPOLATION
                else if(value(ps[i]) == l_False){

                    assert( var_types_[var(ps[i])] == A_LOCAL ||
                            var_types_[var(ps[i])] == B_LOCAL ||
                            var_types_[var(ps[i])] == GLOBAL );

                    if( isAlocal(ps[i]) ){
                        // or - case
                        partial_interpolant = m_.Or(partial_interpolant,getInterUnit(~ps[i]));
                    }
                    else{
                        // and - case
                        partial_interpolant = m_.And(partial_interpolant,getInterUnit(~ps[i]));
                    }
                }
#endif
            }
            ps.shrink(i - j);
        }

        if (ps.size() == 0){
#ifdef CRAIG_INTERPOLATION
            std::cerr << "Clause is empty!!!\n";
            // the computed interpolant is the final craig interpolant
            craig_interpolant_ = partial_interpolant;
#endif
            return ok = false;
        }
        else if (ps.size() == 1){
#ifdef CRAIG_INTERPOLATION
            // as minisat does not store unit clauses we have to store
            // out partial interpolant for this unit clause
            setPartialInter(ps[0],partial_interpolant);
#endif
            assert(value(ps[0]) == l_Undef);
            uncheckedEnqueue(ps[0]);
            return ok = (propagate() == NULL);
        }
        else{
            Clause* c = Clause_new(ps, false);
#ifdef CRAIG_INTERPOLATION
            c->setPartialInterpolant( partial_interpolant );
            c->setClauseType(cl_type);
            c->setClauseInfo(info);
#endif
            clauses.push(c);
            attachClause(*c);
        }
        return true;
    }


    void Solver::attachClause(Clause& c) {
        assert(c.size() > 1);
        assert(c.getPartialInterpolant().isValid());

#ifdef MY_DECISION
        for(int i = 0; i < c.size(); ++i){
            if((int)(pseudo_occurence_.size()) <= var(c[i])){
                pseudo_occurence_.resize(var(c[i]) + 1,0);
            }
            if(sign(c[i])){
                // negative literal
                --pseudo_occurence_[var(c[i])];
            }
            else{
                ++pseudo_occurence_[var(c[i])];
            }
        }
#endif
        watches[toInt(~c[0])].push(&c);
        watches[toInt(~c[1])].push(&c);
        if (c.learnt()) learnts_literals += c.size();
        else            clauses_literals += c.size();
    }


    void Solver::detachClause(Clause& c) {
        assert(c.size() > 1);
#ifdef MY_DECISION
        for(int i = 0; i < c.size(); ++i){
            assert((int)pseudo_occurence_.size() > var(c[i]));
            if(sign(c[i])){
                // negative literal
                ++pseudo_occurence_[var(c[i])];
            }
            else{
                --pseudo_occurence_[var(c[i])];
            }
        }
#endif
        assert(find(watches[toInt(~c[0])], &c));
        assert(find(watches[toInt(~c[1])], &c));
        remove(watches[toInt(~c[0])], &c);
        remove(watches[toInt(~c[1])], &c);
        if (c.learnt()){
#ifdef CRAIG_INTERPOLATION
            c.invalidatePartialInterpolant();
#endif
            learnts_literals -= c.size();
        }
        else
            clauses_literals -= c.size();
    }


    void Solver::removeClause(Clause& c) {
        detachClause(c);
        free(&c);
    }


    bool Solver::satisfied(const Clause& c) const {
        for (int i = 0; i < c.size(); i++)
            if (value(c[i]) == l_True)
                return true;
        return false; }


    // Revert to the state at given level (keeping all assignment at 'level' but not beyond).
    //
    void Solver::cancelUntil(int level) {
        if (decisionLevel() > level){
            for (int c = trail.size()-1; c >= trail_lim[level]; c--){
                Var     x  = var(trail[c]);
#ifdef CRAIG_INTERPOLATION
                assignment[x] = -1;
                a_--;
#endif
                assigns[x] = toInt(l_Undef);
                insertVarOrder(x); }
            qhead = trail_lim[level];
            trail.shrink(trail.size() - trail_lim[level]);
            trail_lim.shrink(trail_lim.size() - level);
        } }


    //=================================================================================================
    // Major methods:


    Lit Solver::pickBranchLit(int polarity_mode, double random_var_freq)
    {
        Var next = var_Undef;

        // Random decision:
        if (drand(random_seed) < random_var_freq && !order_heap.empty()){
            next = order_heap[irand(random_seed,order_heap.size())];
            if (toLbool(assigns[next]) == l_Undef && decision_var[next])
                rnd_decisions++; }

        // Activity based decision:
        while (next == var_Undef || toLbool(assigns[next]) != l_Undef || !decision_var[next])
            if (order_heap.empty()){
                next = var_Undef;
                break;
            }else
                next = order_heap.removeMin();

        bool sign = false;
        switch (polarity_mode){
        case polarity_true:  sign = false; break;
        case polarity_false: sign = true;  break;
        case polarity_user:  sign = polarity[next]; break;
#ifdef MINISAT_Z09
        case polarity_stored:  sign = polarity[next]; break; /* MODIFIED 2009 */
#endif
        case polarity_rnd:   sign = irand(random_seed, 2); break;
#ifdef MY_DECISION
        case polarity_occ:   sign = (pseudo_occurence_[next] < 0)?true:false; break;
#endif
        default: assert(false); }

        return next == var_Undef ? lit_Undef : Lit(next, sign);
    }


    /*_________________________________________________________________________________________________
      |
      |  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
      |
      |  Description:
      |    Analyze conflict and produce a reason clause.
      |
      |    Pre-conditions:
      |      * 'out_learnt' is assumed to be cleared.
      |      * Current decision level must be greater than root level.
      |
      |    Post-conditions:
      |      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
      |
      |  Effect:
      |    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
      |________________________________________________________________________________________________@*/
#ifndef CRAIG_INTERPOLATION
    void Solver::analyze(Clause* confl, vec<Lit>& out_learnt, int& out_btlevel)
#endif
#ifdef CRAIG_INTERPOLATION
    SimpleAIGEdge Solver::analyze(Clause* confl, vec<Lit>& out_learnt, int& out_btlevel )
#endif
    {
        int pathC = 0;
        Lit p     = lit_Undef;

        // Generate conflict clause:
        //
        out_learnt.push();      // (leave room for the asserting literal)
        int index   = trail.size() - 1;
        out_btlevel = 0;

#ifdef CRAIG_INTERPOLATION
        SimpleAIGEdge inter = (*confl).getPartialInterpolant();
        SimpleAIGEdge other;
        bool do_resolution = false;
        bool pivo_is_a_local = false;
#endif

#ifdef IMPROVED_CCMIN
        /// improved conflict clause minimization
        trace_reasons.clear();
#endif

        do{
            assert(confl != NULL);          // (otherwise should be UIP)
#ifdef IMPROVED_CCMIN
            /// improved conflict clause minimization
            trace_reasons.push(confl);
#endif
            Clause& c = *confl;

#ifdef CRAIG_INTERPOLATION
            if( do_resolution ){
                other = (c).getPartialInterpolant();
                if( pivo_is_a_local ){
                    inter = m_.Or(inter,other);
                }
                else{
                    inter = m_.And(inter,other);
                }
            }
#endif

#ifdef DEBUG_ANALYZE
            std::cout << "(";
            for ( int z = 0; z < c.size() ; ++z ){
                std::cout << toInt(c[z])  << " ";
            }
            std::cout << ") \n";
#endif
            if (c.learnt())
                claBumpActivity(c);

            for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
                Lit q = c[j];

                if (!seen[var(q)] && level[var(q)] > 0){
                    varBumpActivity(var(q));
                    seen[var(q)] = 1;
                    if (level[var(q)] >= decisionLevel()){
                        pathC++;
                    }
                    else{
                        out_learnt.push(q);
                        if (level[var(q)] > out_btlevel)
                            out_btlevel = level[var(q)];
                    }
                }
#ifdef CRAIG_INTERPOLATION
                else{
                    if( !seen[var(q)] ){
                        assert(level[var(q)] == 0);
                        if( isAlocal(q) ){
                            inter = m_.Or(inter,getInterUnit(~q));
                        }
                        else{
                            inter = m_.And(inter,getInterUnit(~q));
                        }
                    }
                }
#endif
            }

            // Select next clause to look at:
            while (!seen[var(trail[index--])]);

            p     = trail[index+1];

#ifdef CRAIG_INTERPOLATION
            pivo_is_a_local = isAlocal(p);
            do_resolution = true;
#endif

#ifdef DEBUG_ANALYZE
            std::cout << "   resolve pivot " << toInt(p) << "\n";
#endif
            confl = reason[var(p)];
            seen[var(p)] = 0;
            pathC--;

        }while (pathC > 0);

        out_learnt[0] = ~p;

        // Simplify conflict clause:
#ifndef IMPROVED_CCMIN
        int i, j;
#endif
#ifdef IMPROVED_CCMIN
        // improved conflict clause minimization
        int i, j, minim_res_ctr;
#endif

        if (expensive_ccmin){
            uint32_t abstract_level = 0;
            for (i = 1; i < out_learnt.size(); i++)
                abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels
                                                                     // involved in conflict)

            out_learnt.copyTo(analyze_toclear);
#ifndef IMPROVED_CCMIN
            for (i = j = 1; i < out_learnt.size(); i++){
#ifndef CRAIG_INTERPOLATION
                if (reason[var(out_learnt[i])] == NULL || !litRedundant(out_learnt[i], abstract_level)){
                    out_learnt[j++] = out_learnt[i];
                }
#endif
#ifndef CRAIG_INTERPOLATION
                // this is not working in craig modus
                assert(false);

                if (reason[var(out_learnt[i])] == NULL ){
                    // no reason no resolution no interpolant
                    out_learnt[j++] = out_learnt[i];
                }
                else if( litRedundant(out_learnt[i], abstract_level, result) ){
                    assert(reason[var(out_learnt[i])] != NULL);

                    if( isAlocal(out_learnt[i]) ){
                        inter = m_.Or(inter,result);
                    }
                    else{
                        inter = m_.And(inter,result);
                    }
                }
                else{
                    assert(reason[var(out_learnt[i])] != NULL );
                    out_learnt[j++] = out_learnt[i];
                }
#endif /* CRAIG_INTERPOLATION */
            }
#endif /* IMPROVED_CCMIN */

#ifdef IMPROVED_CCMIN

#ifndef CRAIG_INTERPOLATION
            minim_res_ctr = find_removable(out_learnt, abstract_level );
#endif /* CRAIG_INTERPOLATION */
#ifdef CRAIG_INTERPOLATION
            minim_res_ctr = find_removable(out_learnt, abstract_level , inter );
#endif /* CRAIG_INTERPOLATION */

            j = prune_removable(out_learnt);

#endif /* IMPROVED_CCMIN */
        }
        else{
            out_learnt.copyTo(analyze_toclear);
            //if(false){
#ifdef CRAIG_INTERPOLATION
            bool only_a   = true;
            bool only_b_g = true;
            std::vector<std::pair<int,Lit> > maybe_to_sort;
#endif
           
            for (i = j = 1; i < out_learnt.size(); i++){
                Clause& c = *reason[var(out_learnt[i])];
                if( reason[var(out_learnt[i])] == NULL ){
                    out_learnt[j++] = out_learnt[i];
                }
                else{

#ifdef CRAIG_INTERPOLATION
                    bool do_resolution = true;
#endif /* CRAIG_INTERPOLATION */

                    for (int k = 1; k < c.size(); k++){
                        if (!seen[var(c[k])] && level[var(c[k])] > 0){
                            out_learnt[j++] = out_learnt[i];

#ifdef CRAIG_INTERPOLATION
                            do_resolution = false;
#endif /* CRAIG_INTERPOLATION */

                            break;
                        }
                    }

#ifdef CRAIG_INTERPOLATION
                    if ( do_resolution ){
                        // update only_a and only_b_g
                        if( isAlocal(out_learnt[i]) ){
                            only_b_g = false;
                            maybe_to_sort.push_back(std::make_pair(assignment[var(out_learnt[i])],out_learnt[i]));
                        }
                        else{
                            only_a = false;
                            maybe_to_sort.push_back(std::make_pair(assignment[var(out_learnt[i])],out_learnt[i]) );
                        }
                        assert( assignment[var(out_learnt[i])] != -1 );
                    }
#endif /* CRAIG_INTERPOLATION */

                }
            }
#ifdef CRAIG_INTERPOLATION
            if( only_a ){
                // Or case
                for(unsigned h = 0 ; h < maybe_to_sort.size() ; ++h){
                    inter = m_.Or(inter,reason[var(maybe_to_sort[h].second)]->getPartialInterpolant());
                }
            }
            else if(only_b_g){
                // And case
                for(unsigned h = 0 ; h < maybe_to_sort.size() ; ++h){
                    inter = m_.And(inter,reason[var(maybe_to_sort[h].second)]->getPartialInterpolant());
                }
            }
            else{
                // we have to sort
                std::sort(maybe_to_sort.begin(), maybe_to_sort.end());
                //std::cout << "After sorting:\n";
                int h = (int)(maybe_to_sort.size()) - 1;
                //int h = 0;
                while ( h >= 0 ){
                    if ( isAlocal( maybe_to_sort[h].second ) ){
                        inter = m_.Or(inter,reason[var(maybe_to_sort[h].second)]->getPartialInterpolant());
                    }
                    else{
                        inter = m_.And(inter,reason[var(maybe_to_sort[h].second)]->getPartialInterpolant());
                    }
                    h--;
                }
            }
#endif
        }
        max_literals += out_learnt.size();
#ifndef IMPROVED_CCMIN
        out_learnt.shrink(i - j);
#endif
#ifdef IMPROVED_CCMIN
        out_learnt.shrink(out_learnt.size() - j);
#endif

#ifdef PARTIAL_INTERPOLANT_TEST
        createPartialTest(inter,out_learnt);
#endif

        tot_literals += out_learnt.size();

        // Find correct backtrack level:
        if (out_learnt.size() == 1){
            out_btlevel = 0;
        }
        else{
            int max_i = 1;
            for (int i = 2; i < out_learnt.size(); i++)
                if (level[var(out_learnt[i])] > level[var(out_learnt[max_i])])
                    max_i = i;
            Lit p             = out_learnt[max_i];
            out_learnt[max_i] = out_learnt[1];
            out_learnt[1]     = p;
            out_btlevel       = level[var(p)];
        }
        for (int j = 0; j < analyze_toclear.size(); j++)
            seen[var(analyze_toclear[j])] = 0;            // ('seen[]' is now cleared)

#ifdef DEBUG_ANALYZE
        std::cout << "Finall clause is: ";
        std::cout << "(";
        for ( int z = 0; z < out_learnt.size() ; ++z ){
            std::cout << toInt(out_learnt[z]) << " ";
        }
        std::cout << ") \n";
        //exit(0);
#endif
#ifdef CRAIG_INTERPOLATION
        return inter;
#endif


    }

    // Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
    // visiting literals at levels that cannot be removed later.
#ifndef CRAIG_INTERPOLATION
    bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
#endif
#ifdef CRAIG_INTERPOLATION
        bool Solver::litRedundant(Lit p, uint32_t abstract_levels, SimpleAIGEdge& result)
#endif
    {
        analyze_stack.clear();
        analyze_stack.push(p);
        int top = analyze_toclear.size();

        while (analyze_stack.size() > 0){
            assert(reason[var(analyze_stack.last())] != NULL);
            Clause& c = *reason[var(analyze_stack.last())];

            analyze_stack.pop();

            for (int i = 1; i < c.size(); i++){
                Lit p  = c[i];
                if (!seen[var(p)] && level[var(p)] > 0){
                    if (reason[var(p)] != NULL && (abstractLevel(var(p)) & abstract_levels) != 0){
                        seen[var(p)] = 1;
                        analyze_stack.push(p);
                        analyze_toclear.push(p);
#ifdef CRAIG_INTERPOLATION
                        if( result.isValid() ){
                            if( isAlocal(p) ){
                                result = m_.Or(result,reason[var(p)]->getPartialInterpolant());
                            }
                            else{
                                result = m_.And(result,reason[var(p)]->getPartialInterpolant());
                            }
                        }
                        else{
                            result = reason[var(p)]->getPartialInterpolant();
                        }
#endif
                    }
                    else{
                        for (int j = top; j < analyze_toclear.size(); j++){
                            seen[var(analyze_toclear[j])] = 0;
                        }
                        analyze_toclear.shrink(analyze_toclear.size() - top);
                        return false;
                    }
                }
            }
        }
        return true;
    }


    /*_________________________________________________________________________________________________
      |
      |  analyzeFinal : (p : Lit)  ->  [void]
      |
      |  Description:
      |    Specialized analysis procedure to express the final conflict in terms of assumptions.
      |    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
      |    stores the result in 'out_conflict'.
      |________________________________________________________________________________________________@*/
    void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
    {
        // std::cerr << "analyzeFinal!!\n";

        out_conflict.clear();
        out_conflict.push(p);

        if (decisionLevel() == 0){
            //std::cerr << "analyzeFinal decisionLevel() == 0 [ Solver.C]\n";
#ifdef CRAIG_INTERPOLATION
            SimpleAIGEdge final_craig = (reason[var(p)] != NULL )? reason[var(p)]->getPartialInterpolant() : getInterUnit(p);
            assert( final_craig.isValid() );
            craig_interpolant_ = final_craig;
#endif
            return;
        }

        seen[var(p)] = 1;

#ifdef CRAIG_INTERPOLATION
        SimpleAIGEdge final_craig = (reason[var(p)] != NULL)?reason[var(p)]->getPartialInterpolant() : getInterUnit(p);
#endif

        for (int i = trail.size()-1; i >= trail_lim[0]; i--){
            Var x = var(trail[i]);


#ifdef DEBUG_ANALYZE
            std::cout << "   Pivot: " << (x << 1) << "\n";
#endif

            if (seen[x]){
                if (reason[x] == NULL){
#ifdef DEBUG_ANALYZE
                    std::cout << "Reason for " << (x << 1) << " is Unit Clause\n";
                    std::cerr << "isAlocal( Lit(" << x << ")): " << isAlocal( Lit(x) ) << "\n";
                    std::cerr << "getInterUnit(trail[i]): " << getInterUnit(trail[i]) << "\n";
#endif

#ifdef CRAIG_INTERPOLATION
                    final_craig = (isAlocal( Lit(x) ))?
                        m_.Or(final_craig,getInterUnit(trail[i])) :
                        m_.And(final_craig,getInterUnit(trail[i]));
#endif
                    assert(level[x] > 0);
                    out_conflict.push(~trail[i]);
                }
                else{
                    Clause& c = *reason[x];
#ifdef CRAIG_INTERPOLATION
                    final_craig = (isAlocal( Lit(x) ))?
                        m_.Or(final_craig,c.getPartialInterpolant()) :
                            m_.And(final_craig,c.getPartialInterpolant());
#ifdef CLAUSE_ORIGN
                    const orign& other_orign = c.getClauseOrign();
                    orign::iterator iter = other_orign.begin();
                    while(  iter != other_orign.end() ){
                        final_orign.insert(*iter);
                        ++iter;
                    }
#endif

#endif

#ifdef DEBUG_ANALYZE
                    std::cout << "(";
                    for ( int z = 0; z < c.size() ; ++z ){
                        std::cout << toInt(c[z]) << " ";
                    }
                    std::cout << ") \n";
#endif
                    for (int j = 1; j < c.size(); j++)
                        if (level[var(c[j])] > 0)
                            seen[var(c[j])] = 1;
                }
                seen[x] = 0;
            }
        }
        seen[var(p)] = 0;
#ifdef CRAIG_INTERPOLATION
        craig_interpolant_ = final_craig;

#ifdef DEBUG_ANALYZE
        std::cerr << "craig_interpolant_: " << craig_interpolant_ << "\n";
#endif /* DEBUG_ANALYZE */

#endif /* CRAIG_INTERPOLATION */

    }


    void Solver::uncheckedEnqueue(Lit p, Clause* from)
    {
        assert(value(p) == l_Undef);
        assigns [var(p)] = toInt(lbool(!sign(p)));  // <<== abstract but not uttermost effecient
        level   [var(p)] = decisionLevel();
#ifdef CRAIG_INTERPOLATION
        assignment[var(p)] = a_++;
#endif
        reason  [var(p)] = from;
#ifdef MINISAT_Z09
        polarity[var(p)] = sign(p); /* Modified 2009 */
#endif
        trail.push(p);

#ifdef CRAIG_INTERPOLATION
        if(decisionLevel() == 0 && from != NULL){
            SimpleAIGEdge inter = from->getPartialInterpolant();
            SimpleAIGEdge other;

            for( int t = 1 ; t < (*from).size() ; ++t){
                other = getInterUnit(~((*from)[t]));
                if( isAlocal( (*from)[t] ) ){
                    inter = m_.Or(inter,other);
                }
                else{
                    inter = m_.And(inter,other);
                }
            }
            setPartialInter( p , inter);
        }
#endif
    }


    /*_________________________________________________________________________________________________
      |
      |  propagate : [void]  ->  [Clause*]
      |
      |  Description:
      |    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
      |    otherwise NULL.
      |
      |    Post-conditions:
      |      * the propagation queue is empty, even if there was a conflict.
      |________________________________________________________________________________________________@*/
    Clause* Solver::propagate()
    {
        Clause* confl     = NULL;
        int     num_props = 0;

        while (qhead < trail.size()){
            Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
            vec<Clause*>&  ws  = watches[toInt(p)];
            Clause         **i, **j, **end;
            num_props++;

            for (i = j = (Clause**)ws, end = i + ws.size();  i != end;){
                Clause& c = **i++;

                // Make sure the false literal is data[1]:
                Lit false_lit = ~p;
                if (c[0] == false_lit)
                    c[0] = c[1], c[1] = false_lit;

                assert(c[1] == false_lit);

                // If 0th watch is true, then clause is already satisfied.
                Lit first = c[0];
                if (value(first) == l_True){
                    *j++ = &c;
                }
                else{
                    // Look for new watch:
                    for (int k = 2; k < c.size(); k++)
                        if (value(c[k]) != l_False){
                            c[1] = c[k]; c[k] = false_lit;
                            watches[toInt(~c[1])].push(&c);
                            goto FoundWatch;
                        }

                    // Did not find watch -- clause is unit under assignment:
                    *j++ = &c;
                    if (value(first) == l_False){
                        confl = &c;
                        qhead = trail.size();
                        // Copy the remaining watches:
                        while (i < end)
                            *j++ = *i++;
                    }
                    else{
                        uncheckedEnqueue(first, &c);
                    }
                }
            FoundWatch:;
            }
            ws.shrink(i - j);
        }
        propagations += num_props;
        simpDB_props -= num_props;

        return confl;
    }

    /*_________________________________________________________________________________________________
      |
      |  reduceDB : ()  ->  [void]
      |
      |  Description:
      |    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
      |    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
      |________________________________________________________________________________________________@*/
    struct reduceDB_lt { bool operator () (Clause* x, Clause* y) { return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity()); } };
    void Solver::reduceDB()
    {
        int     i, j;
        double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

        sort(learnts, reduceDB_lt());
        for (i = j = 0; i < learnts.size() / 2; i++){
            if (learnts[i]->size() > 2 && !locked(*learnts[i]))
                removeClause(*learnts[i]);
            else
                learnts[j++] = learnts[i];
        }
        for (; i < learnts.size(); i++){
            if (learnts[i]->size() > 2 && !locked(*learnts[i]) && learnts[i]->activity() < extra_lim)
                removeClause(*learnts[i]);
            else
                learnts[j++] = learnts[i];
        }
        learnts.shrink(i - j);
#ifdef MINISAT_Z09
         nof_learnts   *= learntsize_inc; /* Modified 2009 */
#endif
    }

#ifdef CRAIG_INTERPOLATION
    void Solver::removeAllLearnts()
    {
        int i;
        int j;
        for (i = j = 0 ; i < learnts.size() ; i++){
            if ( !locked(*learnts[i]) ){
                (*learnts[i]).invalidatePartialInterpolant();
                removeClause(*learnts[i]);
            }
            else
                learnts[j++] = learnts[i];
        }
        learnts.shrink(i - j);
    }
#endif



    void Solver::removeSatisfied(vec<Clause*>& cs)
    {
        int i,j;
        for (i = j = 0; i < cs.size(); i++){
            if (satisfied(*cs[i]))
                removeClause(*cs[i]);
            else
                cs[j++] = cs[i];
        }
        cs.shrink(i - j);
    }


    /*_________________________________________________________________________________________________
      |
      |  simplify : [void]  ->  [bool]
      |
      |  Description:
      |    Simplify the clause database according to the current top-level assigment. Currently, the only
      |    thing done here is the removal of satisfied clauses, but more things can be put here.
      |________________________________________________________________________________________________@*/
    bool Solver::simplify()
    {
        // std::cout << "simplify()!!!" << std::endl;

        assert(decisionLevel() == 0);

        if (!ok || propagate() != NULL)
            return ok = false;

        if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
            return true;

        // Remove satisfied clauses:
        removeSatisfied(learnts);
        if (remove_satisfied)        // Can be turned off.
            removeSatisfied(clauses);

        // Remove fixed variables from the variable heap:
        order_heap.filter(VarFilter(*this));

        simpDB_assigns = nAssigns();
        simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

        return true;
    }


    /*_________________________________________________________________________________________________
      |
      |  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]
      |
      |  Description:
      |    Search for a model the specified number of conflicts, keeping the number of learnt clauses
      |    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
      |    indicate infinity.
      |
      |  Output:
      |    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
      |    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
      |    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
      |________________________________________________________________________________________________@*/
#ifndef MINISAT_Z09
    lbool Solver::search(int nof_conflicts, int nof_learnts)
#endif
#ifdef MINISAT_Z09
    lbool Solver::search()
#endif
    {
        assert(ok);
        int         backtrack_level;
        int         conflictC = 0;
        vec<Lit>    learnt_clause;

        starts++;

//        bool first = true;
        // std::cout << "There are " << nClauses() << " many clauses!!\n";

        // main DPLL loop
        // std::cout << "Entering main dpll loop \n";
        for (;;){

            Clause* confl = propagate();
            if (confl != NULL){
                // CONFLICT
                conflicts++; conflictC++;
#ifdef PARALLEL_MODE
                if(conflicts % 100 == 0 ){
                    if(*status_){
                        // in this case the correct outcome is not needed
                        // e.g. this result is false but a nother parallel
                        // running algorithm has already the corrcet solution
                        return l_Undef;
                    }
                }
#endif
                if (decisionLevel() == 0){
                    std::cout << "decisionLevel() == 0 (after propagation) [ Solver.C ]\n";
                    return l_False;
                }

//                first = false;

                learnt_clause.clear();
#ifndef CRAIG_INTERPOLATION
                analyze(confl, learnt_clause, backtrack_level);
#endif
#ifdef CRAIG_INTERPOLATION

                SimpleAIGEdge inter = analyze(confl,
                                               learnt_clause,
                                               backtrack_level
                );

                assert( inter.isValid() );

#endif
                cancelUntil(backtrack_level);

                assert(value(learnt_clause[0]) == l_Undef);

#ifdef MINISAT_Z09
                backtrackLevels[conflicts % restartMore]= backtrack_level;
#endif
                if (learnt_clause.size() == 1){

#ifdef CRAIG_INTERPOLATION
                    setPartialInter(learnt_clause[0], inter);
#endif
                    uncheckedEnqueue(learnt_clause[0]);
                }
                else{
                    Clause* c = Clause_new(learnt_clause, true);

#ifdef CRAIG_INTERPOLATION
                    c->setPartialInterpolant(inter);
                    c->setClauseType(L_CLAUSE);
                    c->setClauseInfo(LEARNT_CLAUSE);
#endif
                    learnts.push(c);

                    attachClause(*c);
                    claBumpActivity(*c);
                    uncheckedEnqueue(learnt_clause[0], c);
                }

                varDecayActivity();
                claDecayActivity();

            }
            else{
                // NO CONFLICT
#ifndef MINISAT_Z09
                if (nof_conflicts >= 0 && conflictC >= nof_conflicts){
                    // Reached bound on number of conflicts:
                    progress_estimate = progressEstimate();
                    cancelUntil(0);
                    return l_Undef;
                }
#endif /* MINISAT_Z09 */
#ifdef MINISAT_Z09
                if (conflictC >= restartMore) {  /* Modified 2009 */
                    // search and count local minimum
                    int LM= backtrackLevels[0];
                    int nofLM= 1;

                    for(int i=1; i< restartMore; i++) {
                        if(backtrackLevels[i]< LM) {
                            LM= backtrackLevels[i];
                            nofLM= 1;
                        } else if(backtrackLevels[i]== LM) {
                            nofLM++;
                        }
                    }
                    if(LM > restartTolerance && nofLM>= restartLess) { /* Modified 2009 */
                        // AVOIDANCE OF PLATEAUX
                        progress_estimate= progressEstimate();
                        cancelUntil(0);
                        return l_Undef;
                    }
                }

#endif /* MINISAT_Z09 */
                // Simplify the set of problem clauses:
                if (decisionLevel() == 0 && !simplify()){
                    std::cout << "decisionLevel() == 0 && !simplify() [ Solver.C ]\n";
                    return l_False;
                }

                if (nof_learnts >= 0 && learnts.size()-nAssigns() >= nof_learnts){
                    // Reduce the set of learnt clauses:
                    reduceDB();
                }

                Lit next = lit_Undef;
                while (decisionLevel() < assumptions.size()){
                    // std::cerr << "Last while loop!!!\n";
                    // Perform user provided assumption:
                    Lit p = assumptions[decisionLevel()];
                    if (value(p) == l_True){
                        // Dummy decision level:
                        newDecisionLevel();
                    }
                    else if (value(p) == l_False){

                        analyzeFinal(~p, conflict);

#ifdef DEBUG_ANALYZE
                        std::cerr << "Conflict is:\n( ";
                        for ( int z = 0; z < conflict.size() ; ++z ){
                            std::cerr << toInt(conflict[z])  << " ";
                            if( reason[var(conflict[z])] != NULL ){
                                std::cerr << "reason[var(conflict[z])]->getPartialInterpolant():\n";
                                std::cerr << reason[var(conflict[z])]->getPartialInterpolant() << "\n";
                            }
                            else{
                                std::cerr << "getInterUnit(~conflict[z]): \n"
                                          << getInterUnit(~conflict[z]) << "\n";
                            }

                        }
                        std::cerr << ") \n";
#endif
                        return l_False;
                    }
                    else{
                        next = p;
                        break;
                    }
                }

                if (next == lit_Undef){
                    // std::cerr << "new variable decision!!!\n";
                    // New variable decision:
                    decisions++;
                    next = pickBranchLit(polarity_mode, random_var_freq);

                    if (next == lit_Undef)
                        // Model found:
                        return l_True;
                }

                // Increase decision level and enqueue 'next'
                assert(value(next) == l_Undef);
                newDecisionLevel();
                uncheckedEnqueue(next);
            }
        }
    }


    double Solver::progressEstimate() const
    {
        double  progress = 0;
        double  F = 1.0 / nVars();

        for (int i = 0; i <= decisionLevel(); i++){
            int beg = i == 0 ? 0 : trail_lim[i - 1];
            int end = i == decisionLevel() ? trail.size() : trail_lim[i];
            progress += pow(F, i) * (end - beg);
        }

        return progress / nVars();
    }


    bool Solver::solve(const vec<Lit>& assumps){
        if( failed_ != 0 )
        {
            *failed_ = false;
            timeoutTime_ = lrabs::cpuTime() + maxInterpolationTime_;
        }

#ifdef CRAIG_INTERPOLATION
        for(int u = 0; u < assumps.size() ; ++u){
            // std::cerr << "Setting partial craig for assumption: " << toInt(assumps[u]) << "\n";
            if( isAlocal( assumps[u] )){
                // setVarType(int var, VarType type){
                setPartialInter(assumps[u],m_.getFalse());
            }
            else{
                setPartialInter(assumps[u],m_.getTrue());
            }
        }
#endif

        model.clear();
        conflict.clear();

        if (!ok) return false;

        assumps.copyTo(assumptions);

#ifndef MINISAT_Z09
        double  nof_conflicts = restart_first;
        double  nof_learnts   = nClauses() * learntsize_factor;
#endif
#ifdef MINISAT_Z09
        /* Modified 2009 */
        double cvr= (double)nClauses() / (double)nVars();
        nof_learnts= 300000 / cvr;
        restartLess= 5;
        restartMore= 42;
        restartTolerance= nVars() / 10000 +10;
        backtrackLevels= new int[restartMore];
#endif
        lbool   status        = l_Undef;

        if (verbosity >= 1){
            reportf("============================[ Search Statistics ]==============================\n");
            reportf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
            reportf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
            reportf("===============================================================================\n");
        }

        // Search:
        while (status == l_Undef){

            if ( failed_ != 0 )
            {
                if( (int)m_.andNodes() > maxAIGNodeCount_ )
                {
                    *failed_ = true;
                    //printf( "terminating interpolation due to excessive creation of aig nodes: %d/%d\n", m_.andNodes(), maxAIGNodeCount_ );
                    break;
                }

                if( lrabs::cpuTime() >= timeoutTime_ )
                {
                    *failed_ = true;
                    //printf( "terminating interpolation due to excessive cpu time usage\n" );
                    break;
                }
            }

#ifndef MINISAT_Z09
            if (verbosity >= 1){
                reportf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", (int)conflicts, order_heap.size(), nClauses(), (int)clauses_literals, (int)nof_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progress_estimate*100), fflush(stdout);
            }
            status = search((int)nof_conflicts, (int)nof_learnts);
            nof_conflicts *= restart_inc;
            nof_learnts   *= learntsize_inc;
#endif /* MINISAT_Z09 */
#ifdef MINISAT_Z09
            status = search();
#endif
        }

#ifndef MINISAT_Z09
        if (verbosity >= 1)
            reportf("===============================================================================\n");
#endif

        if (status == l_True){
            // Extend & copy model:
            model.growTo(nVars());
            for (int i = 0; i < nVars(); i++) model[i] = value(i);
            //#ifndef NDEBUG
            //verifyModel();
            //#endif
        }
        else if( status == l_False ){
#ifdef CRAIG_INTERPOLATION
            // here we know that the problem is unsat
            // now we put craig_interpolant_ to the vector
            // final_craig_interpolants_
            // this is only needed for the fixed point check
            // when you are trying to solve bmc problems
            //
            // NOTE: It is not allowed to delete minisat between
            //       the different unsat bmc instances. If the current
            //       considered bmc instance is SAT you can delete minisat
            //       and create a new one
            final_interpolants_.push_back(craig_interpolant_);
#endif
            if (conflict.size() == 0)
                ok = false;
        }
        else
        {
                /* failed */
        }

        cancelUntil(0);
        return status == l_True;
    }

    //=================================================================================================
    // Debug methods:


    void Solver::verifyModel()
    {
#ifndef NDEBUG
        bool failed = false;
#endif
        for (int i = 0; i < clauses.size(); i++){
            assert(clauses[i]->mark() == 0);
            Clause& c = *clauses[i];
            for (int j = 0; j < c.size(); j++)
                if (modelValue(c[j]) == l_True)
                    goto next;

            reportf("unsatisfied clause: ");
            printClause(*clauses[i]);
            reportf("\n");
#ifndef NDEBUG
            failed = true;
#endif
        next:;
        }
#ifndef NDEBUG
        assert(!failed);
#endif
        reportf("Verified %d original clauses.\n", clauses.size());
    }


    void Solver::checkLiteralCount()
    {
        // Check that sizes are calculated correctly:
        int cnt = 0;
        for (int i = 0; i < clauses.size(); i++)
            if (clauses[i]->mark() == 0)
                cnt += clauses[i]->size();

        if ((int)clauses_literals != cnt){
            fprintf(stderr, "literal count: %d, real value = %d\n", (int)clauses_literals, cnt);
            assert((int)clauses_literals == cnt);
        }
    }


#ifdef IMPROVED_CCMIN

    /// NEW CCMIN PROCEDURES FOLLOW
    int  Solver::prune_removable(vec<Lit>& out_learnt) {
        int i, j, sz = out_learnt.size();
        j = 1;
        for (i = 1; i < sz; i++) {
            if ((seen[var(out_learnt[i])] & (1|2)) == (1|2)) {
                assert((seen[var(out_learnt[i])] & (4|8)) == 0);
                out_learnt[j++] = out_learnt[i];
            }
        }
        //if (j != sz) {
        //  reportf("prune out_learnt[%d]=%d\n", j-1, toInt(out_learnt[j-1]));}
        return j;
    }

#ifndef CRAIG_INTERPOLATION
    int  Solver::find_removable(vec<Lit>& out_learnt, uint32_t abstract_level)
#endif
#ifdef CRAIG_INTERPOLATION
    int  Solver::find_removable(vec<Lit>& out_learnt, uint32_t abstract_level, SimpleAIGEdge& inter)
#endif
    {
        int found_some;
        found_some = 0;
        int sz = out_learnt.size();
        int i;
        trace_lits_minim.clear();
        for (i = 1; i < sz; i++) {
            Lit curLit = out_learnt[i];
            if (level[var(curLit)] <= 0)
                continue;

            if ((seen[var(curLit)] & (2|4|8)) == 0) {
                found_some |= dfs_removable(curLit, abstract_level);
            }
        }
        int minim_res_ctr;
        if (found_some){
#ifndef CRAIG_INTERPOLATION
            minim_res_ctr = res_removable();
#endif
#ifdef CRAIG_INTERPOLATION
            minim_res_ctr = res_removable(inter);
#endif
        }
        else{
            minim_res_ctr = 0;
        }
        return minim_res_ctr;
    }

    int  Solver::quick_keeper(Lit p, uint32_t abstract_level, int maykeep) {
        // See if I can kill myself right away.
        // maykeep == 1 if I am in the original conflict clause.
        if (reason[var(p)] == NULL) {
            return (maykeep ? 2 : 8);
        } else if ((abstractLevel(var(p)) & abstract_level) == 0) {
            assert(maykeep == 0);
            return 8;
        } else {
            return 0;
        }
    }

    int  Solver::dfs_removable(Lit p, uint32_t abstract_level) {
        int pseen = seen[var(p)];
        assert((pseen & (2|4|8)) == 0);
        int maykeep = pseen & (1);
        int pstatus;
        pstatus = quick_keeper(p, abstract_level, maykeep);
        if (pstatus) {
            seen[var(p)] |= (char) pstatus;
            if (pseen == 0) analyze_toclear.push(p);
            return 0;
        }

        int found_some;
        found_some = 0;
        pstatus = 4;
        Clause& rp = *reason[var(p)];
        int sz = rp.size();
        int i;
        // rp[0] is p.  The rest of rp are predecessors of p.
        for (i = 1; i < sz; i++) {
            Lit q = rp[i];
            if (level[var(q)] <= 0)
                continue;

            if ((seen[var(q)] & (2|4|8)) == 0) {
                found_some |= dfs_removable(q, abstract_level);
            }
            int qseen = seen[var(q)];
            if (qseen & (8)) {
                pstatus = (maykeep ? 2 : 8);
                break;
            }
            assert((qseen & (2|4)));
        }
        if (pstatus == 4) {
            // We might want to resolve p out.  See res_removable().
            trace_lits_minim.push(p);
        }
        seen[var(p)] |= (char) pstatus;
        if (pseen == 0) analyze_toclear.push(p);
        found_some |= maykeep;
        return found_some;
    }

    void  Solver::mark_needed_removable(Lit p) {
        Clause& rp = *reason[var(p)];
        for (int i = 1; i < rp.size(); i++){
            Lit q  = rp[i];
            if (level[var(q)] <= 0)
                continue;

            int qseen = seen[var(q)];
            if ((qseen & (1)) == 0 && reason[var(q) ] != NULL) {
                seen[var(q)] |= 1;
                if (qseen == 0) analyze_toclear.push(q);
            }
        }
        return;
    }

#ifndef CRAIG_INTERPOLATION
    int  Solver::res_removable(void)
#endif
#ifdef CRAIG_INTERPOLATION
    int Solver::res_removable(SimpleAIGEdge& inter)
#endif
    {
        int minim_res_ctr = 0;
        while (trace_lits_minim.size() > 0){
            Lit p = trace_lits_minim.last();
            assert(reason[var(p)] != NULL);
            trace_lits_minim.pop();
            int pseen = seen[var(p)];
            if (pseen & (1)) {
                minim_res_ctr ++;
                trace_reasons.push(reason[var(p)]);
                mark_needed_removable(p);
#ifdef CRAIG_INTERPOLATION
                if(isAlocal(p)){
                    inter = m_.Or(inter,reason[var(p)]->getPartialInterpolant());
                }
                else{
                    inter = m_.And(inter,reason[var(p)]->getPartialInterpolant());
                }
#endif
            }
        }
        return minim_res_ctr;
    }

#endif /* IMPROVED_CCMIN */



}
