/**********************************************************************
 * minisat_bmc.hpp
 *
 * Stefan Kupferschmid <skupfers@informatik.uni-freiburg.de>
 *
 **********************************************************************/
#ifndef MINICRAIG_MINISAT_BMC_HPP
#define MINICRAIG_MINISAT_BMC_HPP

#include <utility>
#include "bmc_var_map.hpp"
#include "bmc_problem.hpp"
#include "bmc_instance_generator.hpp"
#include "Vec.h"
#include "SimpSolver.h"
#include <vector>

namespace minicraig {

#ifdef COLLECT_LATCH_VARIABLES
    enum Result{UNS = 0,
                COUNTEREXAMPLE = 1,
                UNKNOWN = 2
    };
#endif

    class MinisatBmc{
    public:
        MinisatBmc();
        MinisatBmc(const MinisatBmc& other) = delete;
        MinisatBmc& operator=(const MinisatBmc& other) = delete;

        ~MinisatBmc();

        void setInstanceGenerator(BmcInstGenerator* gen);

        // solve until depth
        std::pair<bool,int> solve(int depth);

        // bmc solve procedure with unsat core
        std::pair<bool,int> solve_bmc_unsat(int depth);


        // some printing functions
        void printVarMap();

        void printRevVarMap();

        void resetMiniSatBmc(int depth);

#ifdef PARALLEL_MODE
        void setStatus(bool* status){
            status_ = status;
            minisat_->setStatus(status_);
        }
#endif

#ifdef COLLECT_LATCH_VARIABLES
        std::vector<std::set<unsigned> >& getLatchVarsInCraig(){
            return latch_vars_in_craig_;
        }
        Result getResult(){
            return result_;
        }
#endif

        bool startFirstCheck();

        void printTrace(unsigned depth);


    private:
        BmcInstGenerator* gen_;

        int current_bmc_depth_;

        vec<Lit> assumptions_;

        Var last_target_trigger_;

        Var last_craig_trigger_;

        Var init_trigger_;

        SimpSolver* minisat_;

        std::vector<Var> var_map_;

        std::vector<int> rev_var_map_;

        bool craig_added_;

        unsigned index_to_last_target_;

        unsigned index_to_last_craig_;

        unsigned craig_counter_;

        int depth_of_last_trans_;

        // approx_number is the number of
        // transition steps that will be overapproximated
        // by the craig interpolant
        int approx_number_;

        bool backward_;

        bool exp_back_;

        // status is needed if you are using a parallel mode
        // if not status_ is set to NULL;
        bool* status_;


#ifdef CLAUSE_ORIGN
        orign last_empty_clause_orign_;

        bool use_modified_trans_;

        std::vector<std::vector< BmcLit >* > mod_trans_;
#endif

#ifdef COLLECT_LATCH_VARIABLES
        // here we collect all global variables that occur
        // in the craig interpolants
        std::vector<std::set<unsigned> > latch_vars_in_craig_;
        uint64_t latch_set_counter_;
        Result result_;
#endif

        // addInitial adds the initial state to the problem
        // and returns the assumption minisat variable to trigger
        // these clauses
        Var addInitial(int depth=0);

        Var addInvariant(int depth);

        // add tarnsition realtione clauses with
        // e.g. addTrans(2,5) the following will be added
        // T(5,6) /\ T(4,5) /\ T(3,4)
        void addTrans(int before, int depth);

        Var addTarget(int depth);

        Var addExpTarget(int depth);

        Var addCraig();

        bool solveCurrentProblem(void);

        void deactivateLastTarget(void);

        void deactivateLastCraig(void);

        void addNextVariables(void);

    };
}
#endif /* MINISAT_BMC_HPP */
