/****************************************************************************************[Solver.h]
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

#ifndef MINICRAIG_Solver_h
#define MINICRAIG_Solver_h

#include <cstdio>
#include <vector>

#include "Vec.h"
#include "Heap.h"
#include "Alg.h"

#include "SolverTypes.h"
#include "craig_define.hpp"
#include "proof_define.hpp"
#include "parallel_define.hpp"

#ifdef CRAIG_INTERPOLATION
#include "types.hpp"
#include "simple_aig.hpp"
#include <map>

#ifdef PARTIAL_INTERPOLANT_TEST
#include <fstream>
#include <sstream>
#endif


// #define IMPROVED_CCMIN
// #define MINISAT_Z09
#endif /* CRAIG_INTERPOLATION */

//=================================================================================================
// Solver -- the main class:

namespace minicraig {

class Solver {
public:

    // Constructor/Destructor:
    //
    Solver();
    ~Solver();

#ifdef CRAIG_INTERPOLATION
    // after adding all B-Clauses it is not allowed
    // to add further problem clauses that are not pure
    // e.g. not lying totally in the A or B clause set
    bool no_more_b_clauses_;

    // here we store for each variable  the corresponding variable type
    std::vector< VarType > var_types_;

    // aig manager manages the partial craig interpolants
    SimpleAIG m_;
    const SimpleAIG& getInterpolantAIG() const
    {
        return m_;
    }
    

    // here we store the aig node for a global variable
    // we store always positive variables
    std::vector<SimpleAIGEdge> global_vars_;

    unsigned max_craig_var_id_;

    unsigned getMaxCraigVarId(){
        return max_craig_var_id_;
    }


    // in partial_inter_on_dl_0_ vector we store the partial interpolant
    // of a unit clause on decision level 0
    std::vector<SimpleAIGEdge> partial_inter_on_dl_0_;

    // final craig interpolant
    SimpleAIGEdge craig_interpolant_;
    const SimpleAIGEdge& getInterpolant() const
    {
        return craig_interpolant_;
    }
    
    bool* failed_;
    int maxAIGNodeCount_;
    double maxInterpolationTime_;
    double timeoutTime_;

    void setResourceLimits( int nodeLimit, double timeLimit, bool* failed )
    {
        failed_ = failed;
        maxAIGNodeCount_ = nodeLimit;
        maxInterpolationTime_ = timeLimit;
    }

    // reverse var map
    std::vector<int>* rev_var_map_;

    // var map
    std::vector<Var>* var_map_;

    // mapping is used to construct out of the aig
    // representation of the craig interpolant a
    // cnf formule e.g. vector< vector<int> >
    // this map  is filled in method addCraigVar

    std::map<unsigned, SimpleAIGEdge> mapping_;

    std::vector<SimpleAIGEdge> final_interpolants_;

    int approx_number_;

    bool* status_;

    void setStatus(bool* status){
        status_ = status;
    }


    SimpleAIGEdge& getCraigInterpolant(){
        return craig_interpolant_;
    }

    void setRevVarMap(std::vector<int>* rev_var_map){
        rev_var_map_ = rev_var_map;
    }

    std::vector<int>* getRevVarMap(){
        return rev_var_map_;
    }


    void setVarMap(std::vector<Var>* var_map){
        var_map_ = var_map;
    }

    std::vector< std::vector<unsigned> > getCnfCraig(unsigned& nextfree);

    // helper methods
    bool isAlocal(Lit lit)const{
        if( var_types_.size() <=  (unsigned)var(lit) ){
            std::cerr << "variabble: " << var(lit) << "has no type!!\n";
        }

        assert(var_types_.size() > (unsigned)var(lit) );
        return ( A_LOCAL == var_types_[var(lit)] );
    }

    void setVarType(int var, VarType type){

        if( var_types_.size() <= (unsigned)var ){
            // e.g. vector not big enough
            // create room and insert
            var_types_.resize(var + 1,UNDEF);
            var_types_[var] = type;
        }
        else{
            // e.g. vector has enough room
            if( var_types_[var] == UNDEF ){
                var_types_[var] = type;
            }
            else{
                if(  !(var_types_[var] == GLOBAL) && !( var_types_[var] == type) ){
                    // here a varibel became a "global" variable
                    // var is the minisat var e.g. just an integer
                    addCraigVar(var);
                    var_types_[var] = GLOBAL;
                    // std::cerr << "var: " << var << " is global variable!!\n";
                }
            }
        }
    }

    VarType getVarType(Lit lit){
        assert( var_types_.size() > (unsigned)(var(lit)) );
        assert( var_types_[var(lit)] != UNDEF );
        return var_types_[var(lit)];
    }

    void setPartialInter(Lit lit,const SimpleAIGEdge& inter){
        if( partial_inter_on_dl_0_.size() <= (unsigned)(toInt(lit)) ){
             partial_inter_on_dl_0_.resize(toInt(lit) + 1);
        }
        if( partial_inter_on_dl_0_[toInt(lit)].isValid() )
        {
            //std::cout << "WARNING: overwriting partial (variable) interpolant" << std::endl;
            //std::cout << "variable index = " << var(lit) << std::endl;
            //m_.print( partial_inter_on_dl_0_[toInt(lit)] );
            //m_.print( inter );
        }
        partial_inter_on_dl_0_[toInt(lit)] = inter;
        reason[var(lit)] = NULL;
    }

    SimpleAIGEdge& getInterUnit(Lit lit){
        assert( partial_inter_on_dl_0_.size() > (unsigned)(toInt(lit)) );
        assert( partial_inter_on_dl_0_[toInt(lit)].isValid() );
        return partial_inter_on_dl_0_[toInt(lit)];
    }

    SimpleAIGEdge getCraigVar(Lit lit){
        // assert that the global var is already inserted
        // and the corresponding craig variable has been
        // created
        assert( global_vars_.size() > (unsigned)(var(lit)) );
        assert( !m_.isConstant( global_vars_[ var( lit ) ] ) );
        
        if( sign(lit) ){
            // negative literal
            return !(global_vars_[var(lit)]);

        }
        else{
            // positive literal
            return (global_vars_[var(lit)]);
        }
    }
    void addCraigVar(int var){

        // std::cerr << "approx_number_ is : " << approx_number_ << "\n";

        if( global_vars_.size() <= (unsigned)(var) ){
            // e.g. vector not big enough
            // create room and insert
            global_vars_.resize(var + 1, m_.getFalse() );
        }
        assert( m_.isConstant( global_vars_[var] ) );

        global_vars_[var] = m_.addVar(var);
        mapping_[ var ] = global_vars_[var];

    }

    void clearInterpolants(){
        final_interpolants_.clear();
    }

    void setApproxNumber(int value){
        approx_number_ = value;
    }

    void removeAllLearnts();

#endif /* CRAIG_INTERPOLATION */

    const vec<lbool>& getModel(){
            return model;
    }

#ifdef MY_DECISION
    std::vector<int> pseudo_occurence_;
#endif
    // Problem specification:
    //
    Var     newVar    (bool polarity = true, bool dvar = true); // Add
                                                                // a
                                                                // new
                                                                // variable
                                                                // with
                                                                // parameters
                                                                // specifying variable mode.
#ifndef CRAIG_INTERPOLATION
    bool    addClause (vec<Lit>& ps);                           // Add
                                                                // a
                                                                // clause
                                                                // to
                                                                // the
                                                                // solver. NOTE!
                                                                // 'ps'
                                                                // may
                                                                // be
                                                                // shrunk
                                                                // by
                                                                // this
                                                                // method!

#endif
#ifdef CRAIG_INTERPOLATION

    bool    addClause (vec<Lit>& ps,
                       ClauseType type,
                       ClauseInfo info,
                       bool pure,
                       bool take_inter,
                       const SimpleAIGEdge inter);
    
#ifdef PARTIAL_INTERPOLANT_TEST
    unsigned int test_next_free;
    std::vector<Var>* test_var_map;
    std::vector<int>* test_rev_var_map;
    std::string a_clauses_;
    std::string b_clauses_;
    unsigned int a_num_;
    unsigned int b_num_;
    
    void set_var_map(std::vector<Var>* v_map){
        test_var_map = v_map;
    }

    void set_rev_var_map( std::vector<int>* r_v_map ){
        test_rev_var_map = r_v_map;
    }
    
    void set_a_clauses(const std::string& a_clauses,unsigned int a_num){
        a_clauses_ = a_clauses;
        a_num_ = a_num;
    }
    
    void set_b_clauses(const std::string& b_clauses,unsigned int b_num){
        b_clauses_ = b_clauses;
        b_num_ = b_num;
    }

    void set_next_free(unsigned n_free){
        test_next_free = n_free;
    }
    
    void createPartialTest(const SimpleAIGEdge& partial_inter,
                           const vec<Lit>& learnt_clause ,
                           const unsigned name = 0) const{
        
        std::vector<int> a_local_lits_negated;
        
        std::vector<int> b_global_or_b_local_lits_negated;
        
        for(int t = 0; t <  learnt_clause.size() ; ++t ){
            // we do not want assumptions
            if( var(learnt_clause[t]) != 0 && var(learnt_clause[t]) != 1){
                if(isAlocal(learnt_clause[t])){
                    if( sign(learnt_clause[t]) ){
                        // negated literal
                        assert( test_rev_var_map->at(var(learnt_clause[t])) != -1);
                        a_local_lits_negated.push_back( (test_rev_var_map->at(var(learnt_clause[t]))) );
                    }
                    else{
                        assert( test_rev_var_map->at(var(learnt_clause[t])) != -1);
                        a_local_lits_negated.push_back( (-1)* (test_rev_var_map->at(var(learnt_clause[t]))) );
                    }
                }
                else{
                    if( sign(learnt_clause[t]) ){
                        // negated literal
                        assert( test_rev_var_map->at(var(learnt_clause[t])) != -1);
                        b_global_or_b_local_lits_negated.push_back( (test_rev_var_map->at(var(learnt_clause[t]))) );
                    }
                    else{
                        assert( test_rev_var_map->at(var(learnt_clause[t])) != -1);
                        b_global_or_b_local_lits_negated.push_back( (-1) * (test_rev_var_map->at(var(learnt_clause[t]))) );
                    } 
                }
            }
        }
        
        std::vector< std::vector< unsigned > > craig;
        std::stringstream craig_stream;
        std::stringstream craig_a_trigger;
        std::stringstream craig_b_trigger;
        unsigned int next_tseitin = test_next_free + 1;
        
        for ( unsigned int i = 0; i < a_local_lits_negated.size() ; ++i){
            craig_a_trigger << a_local_lits_negated[i] << " 0\n";
        }

        for ( unsigned int i = 0; i < b_global_or_b_local_lits_negated.size() ; ++i){
            craig_b_trigger << b_global_or_b_local_lits_negated[i] << " 0\n";
        }
        
        if( aigpp::structurallyFalse( partial_inter ) ){
            
            // craig.push_back(std::vector<unsigned>());
            std::cerr << "partial interpolant is FALSE!!!\n";
            craig_b_trigger << "-1 0\n" << "1 0\n";
        }
        else if( aigpp::structurallyTrue( partial_inter ) ){
            std::cerr << "partial interpolant is TRUE!!!\n";
            craig_a_trigger << "-1 0\n" << "1 0\n";
        }
        else{
            
            unsigned outputLit = 0;
            unsigned nextfree = test_next_free;

            std::map<unsigned int, SimpleAIGEdge> map2 = mapping_;
            
            craig = m_.createCNF( partial_inter,
                                  map2,
                                  outputLit,
                                  nextfree,
                                  false);
            
            std::vector< unsigned > clause;
            clause.push_back(outputLit);
            craig.push_back(clause);

            
            std::vector<int> test_rev_var_map2 = *test_rev_var_map;
            
            for(unsigned int t = 0; t < craig.size()-1 ; ++t){
                for(unsigned int j = 0; j < craig[t].size() ; ++j){
                    
                    bool is_neg = craig[t][j]&1;
                    
                    if( test_rev_var_map2.size() <= (craig[t][j] >> 1) ){
                        test_rev_var_map2.resize((craig[t][j] >> 1)+1,-1);
                    }
                    if( test_rev_var_map2.at((craig[t][j] >> 1)) == -1){
                        test_rev_var_map2.at((craig[t][j] >> 1)) = next_tseitin;
                        next_tseitin++;
                    }
                    if( is_neg ){
                        craig_stream << "-" << test_rev_var_map2.at((craig[t][j] >> 1)) << " ";
                    }
                    else{
                        craig_stream << test_rev_var_map2.at((craig[t][j] >> 1)) << " ";
                    }
                }
                craig_stream << "0\n";
            }
            
            assert( craig.back().size() == 1);
            
            bool is_neg = craig.back()[0]&1;
            if( test_rev_var_map2.size() <= (craig.back()[0] >> 1) ){
                test_rev_var_map2.resize((craig.back()[0] >> 1) + 1,-1);
            }
            if( test_rev_var_map2.at((craig.back()[0] >> 1)) == -1){
                test_rev_var_map2.at((craig.back()[0] >> 1)) = next_tseitin;
                next_tseitin++;
            }
            if( is_neg ){
                craig_a_trigger << test_rev_var_map2.at((craig.back()[0] >> 1)) << " 0\n";
                craig_b_trigger << "-" << test_rev_var_map2.at((craig.back()[0] >> 1)) << " 0\n";
            }
            else{
                craig_b_trigger << test_rev_var_map2.at((craig.back()[0] >> 1)) << " 0\n";
                craig_a_trigger << "-" << test_rev_var_map2.at((craig.back()[0] >> 1)) << " 0\n";
            }
        }
        // now build two files
        std::stringstream sstream_a;
        std::stringstream sstream_b;

        if( name == 0){
            sstream_a << "cnfa-" << conflicts << ".dimacs";
            sstream_b << "cnfb-" << conflicts << ".dimacs";
        }
        else{
            sstream_a << "cnfa-" << conflicts << "-" << name << ".dimacs";
            sstream_b << "cnfb-" << conflicts << "-" << name << ".dimacs";
        }

        std::string file_a = sstream_a.str();;
        std::string file_b = sstream_b.str();;
       
        const char* filename_a = file_a.c_str();
        const char* filename_b = file_b.c_str();
        
        std::ofstream ofstream_a(filename_a, std::ios::out);
        std::ofstream ofstream_b(filename_b, std::ios::out);
        
        ofstream_a <<  "p cnf " << next_tseitin << " " << a_num_ + craig.size() << "\n";
        ofstream_b <<  "p cnf " << next_tseitin << " " << b_num_ + craig.size() << "\n";
        
        ofstream_a << a_clauses_ << craig_stream.str() << craig_a_trigger.str() << std::flush;
        
        ofstream_b << b_clauses_;
        ofstream_b << craig_stream.str();
        ofstream_b << craig_b_trigger.str() << std::flush;
    }
#endif /* PARTIAL_INTERPOLANT_TEST */

    
#endif /* CRAIG_INTERPOLATION */

    // Solving:
    //
    bool    simplify     ();                        // Removes already satisfied clauses.
    bool    solve        (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    bool    solve        ();                        // Search without assumptions.
    bool    okay         () const;                  // FALSE means solver is in a conflicting state

    // Variable mode:
    //
    void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

    // Read state:
    //
    lbool   value      (Var x) const;       // The current value of a variable.
    lbool   value      (Lit p) const;       // The current value of a literal.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model.
                                            // The last call to solve must have been satisfiable.
    int     nAssigns   ()      const;       // The current number of assigned literals.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nVars      ()      const;       // The current number of variables.

    // Extra results: (read-only member variable)
    //
    vec<lbool> model;             // If problem is satisfiable, this vector contains the model (if any).
    vec<Lit>   conflict;          // If problem is unsatisfiable (possibly under assumptions),
                                  // this vector represent the final conflict clause expressed in the assumptions.

    // Mode of operation:
    //
    double    var_decay;          // Inverse of the variable activity decay factor.                                            (default 1 / 0.95)
    double    clause_decay;       // Inverse of the clause activity decay factor.                                              (1 / 0.999)
    double    random_var_freq;    // The frequency with which the decision heuristic tries to choose a random variable.        (default 0.02)
    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)
    bool      expensive_ccmin;    // Controls conflict clause minimization.                                                    (default TRUE)
    int       polarity_mode;      // Controls which polarity the decision heuristic chooses. See enum below for allowed modes. (default polarity_false)
    int       verbosity;          // Verbosity level. 0=silent, 1=some progress report                                         (default 0)
#ifdef MINISAT_Z09
    /* Modified 2009 */
	int restartLess;
	int restartMore;
	float restartTolerance;
	double nof_learnts;
	int* backtrackLevels;

#endif /* MINISAT_Z09 */

    enum { polarity_true = 0,
           polarity_false = 1,
           polarity_user = 2,
           polarity_rnd = 3
#ifdef MY_DECISION
           ,polarity_occ = 4
#endif /* MY_DECISION */
#ifdef MINISAT_Z09
           ,polarity_stored = 5
#endif
    };

    // Statistics: (read-only member variable)
    //
    uint64_t starts, decisions, rnd_decisions, propagations, conflicts;
    uint64_t clauses_literals, learnts_literals, max_literals, tot_literals;

    void     varBumpActivities( const vec<Var>& vs );
    
protected:

    // Helper structures:
    //
    struct VarOrderLt {
        const vec<double>&  activity;
        bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
        explicit VarOrderLt(const vec<double>&  act) : activity(act) { }
    };

    friend class VarFilter;
    struct VarFilter {
        const Solver& s;
        explicit VarFilter(const Solver& _s) : s(_s) {}
        bool operator()(Var v) const { return toLbool(s.assigns[v]) == l_Undef && s.decision_var[v]; }
    };

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable.
                                          // No part of the solver state may be used!
    vec<Clause*>        clauses;          // List of problem clauses.
    vec<Clause*>        learnts;          // List of learnt clauses.
    double              cla_inc;          // Amount to bump next clause with.
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.
    vec<vec<Clause*> >  watches;          // 'watches[lit]' is a list of constraints watching 'lit'
                                          // (will go there if literal becomes true).
    vec<char>           assigns;          // The current assignments (lbool:s stored as char:s).
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision_var;     // Declares if a variable is eligible for selection in the decision heuristic.
    vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
    vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<Clause*>        reason;           // 'reason[var]' is the clause that implied the variables current value, or 'NULL' if none.
    vec<int>            level;            // 'level[var]' contains the level at which the assignment was made.
#ifdef CRAIG_INTERPOLATION
    vec<int>            assignment;       // assignment[var] contains information about when var has been assigned a value
    int                 a_;
#endif
    int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double              random_seed;      // Used by the random variable selection.
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            add_tmp;

    // Main internal methods:
    //
    void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
    Lit      pickBranchLit    (int polarity_mode, double random_var_freq);             // Return the next decision variable.
    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     uncheckedEnqueue (Lit p, Clause* from = NULL);                            // Enqueue a literal. Assumes value of literal is undefined.
    bool     enqueue          (Lit p, Clause* from = NULL);                            // Test if fact 'p' contradicts current state, enqueue otherwise.
    Clause*  propagate        ();                                                      // Perform unit propagation. Returns possibly conflicting clause.
    void     cancelUntil      (int level);                                             // Backtrack until a certain level.
#ifndef CRAIG_INTERPOLATION
    void     analyze          (Clause* confl, vec<Lit>& out_learnt, int& out_btlevel); // (bt = backtrack)
#endif
#ifdef CRAIG_INTERPOLATION
    SimpleAIGEdge analyze  (Clause* confl,
                             vec<Lit>& out_learnt,
                             int& out_btlevel  // (bt = backtrack)
    );
#endif
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY
                                                                                       // THE ORDINARIY "analyze" BY SOME
                                                                                       // REASONABLE
                                                                                       // GENERALIZATION?
#ifndef CRAIG_INTERPOLATION
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper
                                                                                       // method
                                                                                       // for
                                                                                       // 'analyze()')
#endif  /* CRAIG_INTERPOLATION */
#ifdef CRAIG_INTERPOLATION
    bool     litRedundant     (Lit p, uint32_t abstract_levels, SimpleAIGEdge& inter);
#endif /* CRAIG_INTERPOLATION */
#ifndef MINISAT_Z09
    lbool    search           (int nof_conflicts, int nof_learnts);                    // Search for a given number of conflicts.
#endif /* MINISAT_Z09 */
#ifdef MINISAT_Z09
    lbool    search           ();
#endif /* MINISAT_Z09 */
    void     reduceDB         ();                                                      // Reduce the set of learnt clauses.


    void     removeSatisfied  (vec<Clause*>& cs);                                      // Shrink 'cs' to contain only non-satisfied clauses.

    // Maintaining Variable/Clause activity:
    //
    void     varDecayActivity ();                      // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void     varBumpActivity  (Var v);                 // Increase a variable with the current 'bump' value.
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.

    // Operations on clauses:
    //
    void     attachClause     (Clause& c);             // Attach a clause to watcher lists.
    void     detachClause     (Clause& c);             // Detach a clause to watcher lists.
    void     removeClause     (Clause& c);             // Detach and free a clause.
    bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    // Misc:
    //
    int      decisionLevel    ()      const; // Gives the current decisionlevel.
    uint32_t abstractLevel    (Var x) const; // Used to represent an abstraction of sets of decision levels.
    double   progressEstimate ()      const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    // Debug:
    void     printLit         (Lit l);
    template<class C>
    void     printClause      (const C& c);
    void     verifyModel      ();
    void     checkLiteralCount();

    // Static helpers:
    //

    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size);
    }

#ifdef IMPROVED_CCMIN
    //// improved conflict clause minimizations
    //   AVG routines for minimization.
    vec<Clause*>        trace_reasons;                      // clauses to resolve to give CC
    vec<Lit>            trace_lits_minim;                   // lits maybe used in minimization
    int  prune_removable(vec<Lit>& out_learnt);
#ifndef CRAIG_INTERPOLATION
    int  find_removable(vec<Lit>& out_learnt, uint32_t abstract_level);
#endif
#ifdef CRAIG_INTERPOLATION
    int  find_removable(vec<Lit>& out_learnt, uint32_t abstract_level, SimpleAIGEdge& inter);
#endif
    int  quick_keeper(Lit p, uint32_t abstract_level, int maykeep);
    int  dfs_removable(Lit p, uint32_t abstract_level);
    void  mark_needed_removable(Lit p);
#ifndef CRAIG_INTERPOLATION
    int  res_removable(void);
#endif
#ifdef CRAIG_INTERPOLATION
    int  res_removable(SimpleAIGEdge& inter);
#endif
#endif
};


//=================================================================================================
// Implementation of inline methods:


inline void Solver::insertVarOrder(Var x) {
    if (!order_heap.inHeap(x) && decision_var[x]) { order_heap.insert(x); 
    }
}

inline void Solver::varDecayActivity() { var_inc *= var_decay; }
inline void Solver::varBumpActivity(Var v) {
    if ( (activity[v] += var_inc) > 1e100 ) {
        // Rescale:
        for (int i = 0; i < nVars(); i++) {
            activity[i] *= 1e-100;
        }
        var_inc *= 1e-100;
    }

    // Update order_heap with respect to new activity:
    if (order_heap.inHeap(v)) {
        order_heap.decrease(v); 
    }
}

inline void Solver::varBumpActivities( const vec<Var>& vs ) 
{
    bool rescale = false;
    
    for( int i = 0; i != vs.size(); ++i )
    {
        Var v = vs[i];
        if ( ( activity[v] += var_inc ) > 1e100 ) 
        {
            rescale = true;
        }
        
        if( order_heap.inHeap( v ) )
        {
            order_heap.decrease( v );
        }
    }
    
    if( rescale )
    {
        for( int i = 0; i < nVars(); i++ )
        {
            activity[i] *= 1e-100;
        }
        var_inc *= 1e-100; 
    }
}

inline void Solver::claDecayActivity() { cla_inc *= clause_decay; }
inline void Solver::claBumpActivity (Clause& c) {
        if ( (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++) {
                learnts[i]->activity() *= 1e-20;
            }
            cla_inc *= 1e-20; }
}

inline bool     Solver::enqueue         (Lit p, Clause* from)   { return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true); }
inline bool     Solver::locked          (const Clause& c) const { return reason[var(c[0])] == &c && value(c[0]) == l_True; }
inline void     Solver::newDecisionLevel()                      { trail_lim.push(trail.size()); }

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }
inline uint32_t Solver::abstractLevel (Var x) const   {
    return 1u << (level[x] & 31u);
}
inline lbool    Solver::value         (Var x) const   { return toLbool(assigns[x]); }
inline lbool    Solver::value         (Lit p) const   { return toLbool(assigns[var(p)]) ^ sign(p); }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return assigns.size(); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity    [v] = (char)b; }
inline void     Solver::setDecisionVar(Var v, bool b) { decision_var[v] = (char)b; if (b) { insertVarOrder(v); } }
inline bool     Solver::solve         ()              { vec<Lit> tmp; return solve(tmp); }
inline bool     Solver::okay          ()      const   { return ok; }



//=================================================================================================
// Debug + etc:


#define reportf(format, args...) ( fflush(stdout), fprintf(stderr, format, ## args), fflush(stderr) )

static inline void logLit(FILE* f, Lit l)
{
    fprintf(f, "%sx%d", sign(l) ? "~" : "", var(l)+1);
}

static inline void logLits(FILE* f, const vec<Lit>& ls)
{
    fprintf(f, "[ ");
    if (ls.size() > 0){
        logLit(f, ls[0]);
        for (int i = 1; i < ls.size(); i++){
            fprintf(f, ", ");
            logLit(f, ls[i]);
        }
    }
    fprintf(f, "] ");
}

static inline const char* showBool(bool b) { return b ? "true" : "false"; }


// Just like 'assert()' but expression will be evaluated in the release version as well.
static inline void check(bool expr) { (void)expr; assert(expr); }


inline void Solver::printLit(Lit l)
{
    reportf("%s%d:%c", sign(l) ? "-" : "", var(l)+1, value(l) == l_True ? '1' : (value(l) == l_False ? '0' : 'X'));
}


template<class C>
inline void Solver::printClause(const C& c)
{
    for (int i = 0; i < c.size(); i++){
        printLit(c[i]);
        fprintf(stderr, " ");
    }
}

}
//=================================================================================================
#endif
