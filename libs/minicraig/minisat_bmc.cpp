#include "minisat_bmc.hpp"
#include "craig_define.hpp"
#include "proof_define.hpp"
#include "parallel_define.hpp"
#include <cassert>

//#define DEBUG_MINISAT_BMC
//#define DEBUG_BMC_PROBLEM
//#define DECIDE_ONLY_INPUTS_AND_LATCHES

namespace minicraig {

    MinisatBmc::MinisatBmc():gen_(NULL),
                             current_bmc_depth_(-1),
                             assumptions_(),
                             last_target_trigger_(0),
                             last_craig_trigger_(0),
                             init_trigger_(0),
                             minisat_(NULL),
                             var_map_(),
                             rev_var_map_(),
                             craig_added_(false),
                             index_to_last_target_(0),
                             index_to_last_craig_(0),
                             craig_counter_(0),
                             depth_of_last_trans_(0),
                             approx_number_(1),
                             backward_(false),
                             exp_back_(false),
                             status_(NULL)
#ifdef CLAUSE_ORIGN
                            ,last_empty_clause_orign_(),
                             use_modified_trans_(false),
                             mod_trans_()
#endif
#ifdef COLLECT_LATCH_VARIABLES
                            ,latch_vars_in_craig_(1),
                             latch_set_counter_(0),
                             result_()
#endif
    {
        minisat_ = new SimpSolver();
#ifdef CRAIG_INTERPOLATION
        minisat_->setRevVarMap( &(rev_var_map_) );
        minisat_->setVarMap( &(var_map_) );
#endif
    }


    MinisatBmc::~MinisatBmc()
    {
        delete minisat_;
    }


    void MinisatBmc::setInstanceGenerator(BmcInstGenerator* gen)
    {
        assert( gen != NULL );
        gen_ = gen;
#ifdef CRAIG_INTERPOLATION
        minisat_->setBmcVarMap(gen_->getBmcVarMap() );
#endif
        approx_number_ = (-1)*(gen_->getOverApprox() - 1);

        backward_ = gen->isBackward();

        exp_back_ = gen->isExpBackward();

        // not both backward options at the same time!!
        assert(!(exp_back_ && backward_));


#ifdef CRAIG_INTERPOLATION
        minisat_->setApproxNumber( approx_number_ );
#endif
    }

    std::pair<bool,int> MinisatBmc::solve(int depth)
    {

        assert( gen_ != NULL );

        gen_->init(0);

        addNextVariables();

        last_craig_trigger_ = addInitial(approx_number_);

        init_trigger_ = last_craig_trigger_;

        last_target_trigger_ = addTarget(approx_number_);

        // important to add first the initial trigger
        assumptions_.push(Lit( init_trigger_, false));
        index_to_last_craig_ = assumptions_.size() - 1;

        // activating target through assumption
        // and store the index to the last target trigger
        // in assuimption vector
        assumptions_.push(Lit( last_target_trigger_, false));
        index_to_last_target_ = assumptions_.size()-1;

        current_bmc_depth_ = 0;

        // solve for depth 0
#ifndef PARALLEL_MODE
        std::cout << "Solving depth 0       ";
#endif

        bool result = solveCurrentProblem();

        if( result ){
#ifndef PARALLEL_MODE
            std::cout << "sat\n";
            std::cout << "Result:                 sat\n";
#endif
#ifdef COLLECT_LATCH_VARIABLES
            result_ = COUNTEREXAMPLE;
#endif
            return std::make_pair(result,current_bmc_depth_);
        }
#ifndef PARALLEL_MODE
        std::cout << "unsat\n";
#endif

        // reset bm instance generater and compute
        // the desired number of disjunctions
        // gen_->resetBmcInstGen();

        // now we start bmc + craig for the Problem
        // I_0 /\ T_(0,1) /\ not(P_1)

        // if(startFirstCheck()){
        //    // not correct
        //    return std::make_pair(result,current_bmc_depth_);
        // }

        gen_->resetBmcInstGen();

        gen_->setOverApprox(gen_->getOverApprox());

        // delete old minisat
        delete minisat_;

        minisat_ = new SimpSolver();

#ifdef PARALLEL_MODE
        minisat_->setStatus(status_);
#endif

        gen_->init((-1)*approx_number_);

        var_map_.clear();
        rev_var_map_.clear();

#ifdef CRAIG_INTERPOLATION
        minisat_->setApproxNumber( approx_number_ );
        minisat_->setRevVarMap( &(rev_var_map_) );
        minisat_->setVarMap( &(var_map_) );
        minisat_->setBmcVarMap(gen_->getBmcVarMap() );
#endif

        addNextVariables();
        assumptions_.clear();

        last_craig_trigger_ = addInitial(approx_number_);

        init_trigger_ = last_craig_trigger_;

        assumptions_.push(Lit( init_trigger_, false));
        index_to_last_craig_ = assumptions_.size() - 1;

        // main bmc loop

        for(int d = gen_->getOverApprox() ; d < depth ; d+=gen_->getOverApprox() )
        {

#ifdef PARALLEL_MODE
            if( *status_ ){
                return std::make_pair(result,current_bmc_depth_);
            }
#endif

            if( !craig_added_ ){
                gen_->init( d + 1 - approx_number_ );
                // std::cout << "d + 1 - approx_number_:" << d + 1 - approx_number_ << "\n";

                addNextVariables();

                // first adding B-Clauses
                if( backward_ ){
                    last_target_trigger_ = addTarget(d + 1);
                }
                else if(exp_back_){
                    last_target_trigger_ = addExpTarget(d + approx_number_ + 1);
                }
                else{
                    last_target_trigger_ = addTarget(d + approx_number_ + 1);
                }

                // std::cout << "d: " << d <<"\n";
                if( exp_back_ ){
                    // std::cout << "gen_->getOverApprox(): " <<
                    // gen_->getOverApprox() << "\n";
                    // std::cout << "addTrans("
                    //          << d - depth_of_last_trans_
                    //          << "," << d + approx_number_ << ")\n";

                    addTrans( d - depth_of_last_trans_ , d + approx_number_ );

                    depth_of_last_trans_ = d + 1;
                    // + gen_->getDepthOffset() + approx_number_;
                }
                else{
                    addTrans( d - depth_of_last_trans_ - approx_number_ , d );
                    depth_of_last_trans_ = d + 1 +  gen_->getDepthOffset();
                }

                // activating target through assumption
                assumptions_.push(Lit( last_target_trigger_, false));
                index_to_last_target_ = assumptions_.size() - 1;

                // here we add the initial or if backward is enabled the
                // target as first craig interpolant to the set of craig
                // interpolants used for fixed point computation
#ifdef CRAIG_INTERPOLATION
                minisat_->addFirstCraig( gen_->getBmcProblem() );
#endif
            }
#ifndef PARALLEL_MODE
            std::cout << "Solving depth " << d + 1 + gen_->getDepthOffset() << std::flush;
#endif
            result = solveCurrentProblem();

            current_bmc_depth_ = d + 1;

            if( result ){
#ifdef CRAIG_INTERPOLATION
                if( craig_added_
#ifdef CLAUSE_ORIGN
                    || use_modified_trans_
#endif
                ){
#ifndef PARALLEL_MODE
                    std::cout << "       sat (over-approx.)\n";
#endif
#ifdef COLLECT_LATCH_VARIABLES
                    if ( latch_set_counter_ < 4){
                        latch_vars_in_craig_.push_back(std::set<unsigned>());
                        latch_set_counter_++;
                    }
                    else{
                        result_ = UNKNOWN;
                        std::cout << "There are " << latch_vars_in_craig_.size() << " many latch sets\n";
                        std::set<unsigned>::iterator iter;
                        for( unsigned int i = 0 ; i  < latch_vars_in_craig_.size() ; ++i){
                            std::cout << "Set " << (i+1) << ":\n";
                            iter = latch_vars_in_craig_[i].begin();
                            while( iter !=  latch_vars_in_craig_[i].end() ){
                                std::cout << "     L" << (*iter) << "\n";
                                iter++;
                            }
                        }
                        return std::make_pair(result,current_bmc_depth_);
                    }
#endif
                    // minisat_->removeCraigClauses();
                    // at the moment it is not that fast if we reset
                    // minisatbmc after every fith generated interpolant
                    if( craig_counter_ > 10 && false){
                        // reset the old minisat instance
                        // std::cerr << "Doing solver reset!!!\n";
                        resetMiniSatBmc(d);
                        craig_counter_ = 0;
                    }
                    else{

                        // some clause removal

                        // minisat_->removeCraigClauses();
                        // minisat_->removeTargetClauses();

                        // delete the old craig interpolants
                        minisat_->clearInterpolants();
                        //minisat_->removeAllLearnts();

                        // folowing 3 lines are used if you do no reset
                        deactivateLastTarget();
                        deactivateLastCraig();

                        craig_added_ = false;

                        // trigger initial clauses
                        last_craig_trigger_ = init_trigger_;
                        assumptions_[0] = (Lit( init_trigger_, false));

                        // uppdate index_to_last_craig_
                        index_to_last_craig_ = 0;

                        // craig_counter_ = 0;

                    }
                    // decrement d (depth) as we have to check this
                    // current depth again
                    d-= gen_->getOverApprox();
                }
                else{
#endif
#ifndef PARALLEL_MODE
                    std::cout << "       sat\n";
                    std::cout << "Result:                 sat\n";
                    // printTrace(d + 1 + gen_->getDepthOffset() );
#ifdef COLLECT_LATCH_VARIABLES
                    result_ = COUNTEREXAMPLE;
#endif

#endif
                    return std::make_pair(result,current_bmc_depth_);
#ifdef CRAIG_INTERPOLATION
                }
#endif
            }
            else{
#ifdef CLAUSE_ORIGN
                std::cerr << "\n   Number of transition clauses needed for unsat core: "
                          << minisat_->getEmptyClauseOrign().size() << "\n";
#endif

#ifdef CRAIG_INTERPOLATION
                if( craig_added_ ){
#ifndef PARALLEL_MODE
                    std::cout << "       unsat (over-approx.)\n";
#endif
                }
                else{
#endif
#ifndef PARALLEL_MODE
                    std::cout << "       unsat\n";
#endif
                    // resetMiniSatBmc(d--);
                    minisat_->removeAllLearnts();
#ifndef CRAIG_INTERPOLATION
                    deactivateLastTarget();
#endif
#ifdef CRAIG_INTERPOLATION
                }
                //minisat_->removeAllLearnts();
                deactivateLastCraig();

                // perform fixed point check
                if( minisat_->performFixedPointCheck() ){
#ifndef PARALLEL_MODE
                    std::cout << "Result:                 uns\n";
#endif

#ifdef COLLECT_LATCH_VARIABLES
                    result_ = UNS;
                    unsigned max_craig = minisat_->getMaxCraigVarId();
                    unsigned next_free = max_craig + 1;
                    unsigned old_next_free = next_free;
                    std::vector< std::vector<unsigned> > craig_cnf = minisat_->getCnfCraig(next_free);

                    if( craig_cnf.size() == 0 || craig_cnf[0].size() == 0 ){
                        std::cerr << "Interpolant is const (true/false)\n";
                    }
                    else{

                        for(unsigned t = 0; t < craig_cnf.size() ; ++t){
                            for(unsigned int s = 0; s < craig_cnf[t].size() ; ++s){
                                if( ( craig_cnf[t][s] >> 1 ) < old_next_free ){
                                    // global variables
                                    latch_vars_in_craig_[latch_set_counter_].insert(
                                        (gen_->getBmcVarMap()->
                                         getBmcPair(minisat_->getRevVarMap()->
                                                    at(craig_cnf[t][s]>>1))).first << 1);
                                }
                            }
                        }
                    }
                    //  std::set<unsigned>::iterator it = latch_vars_in_craig_.begin();
                    //  std::cout << "The following global variables are needed\n"
                    //            << "for the fixed point check!!!\n";

                    //  while(it != latch_vars_in_craig_.end() ){
                    //     std::cout << "Litid " << *it << "\n";
                    //     ++it;
                    //  }
#endif /* COLLECT_LATCH_VARIABLES */

                    // std::cout << "       safe system\n";
                    return std::make_pair(false,depth);
                }

#ifdef PARALLEL_MODE
                if(*status_){
                    return std::make_pair(false,depth);
                }
#endif
                last_craig_trigger_ = addCraig();

                // increment craig counter
                ++craig_counter_;

                // activating craig through assumption
                // and store the index to the last craig trigger
                // in assuimption vector
                assumptions_.push(Lit( last_craig_trigger_, false));
                index_to_last_craig_ = assumptions_.size() - 1;
#endif
            }

        }
        return std::make_pair(false,depth);
    }

    Var MinisatBmc::addInitial(int depth)
    {
#ifdef DEBUG_BMC_PROBLEM
        std::cerr << "adding I(" << depth << ")\n";
        // assert that gen_ is not NULL
#endif
        assert(gen_ != NULL);

        const std::vector<std::vector< BmcLit>* >& cl = gen_->getBmcProblem()->getInit();

        //std::cerr << "dl.size(): " << cl.size() << "\n";

        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

        Var trigger = minisat_->newVar();
        Lit neg_trigger = Lit(trigger,true);

        vec<Lit> clause;

        for( unsigned int t = 0 ; t < cl.size() ; ++t){
            clause.clear();
            for(unsigned int s = 0; s < cl[t]->size() ; ++s ){
                // assert( cl[t]->at(s).getBmcInstance() == 0);
                if(cl[t]->at(s).getVarId() & 1 ){ // negative literal
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth)>>1],true));
                }
                else{
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth)>>1],false));
                }
            }

            clause.push(neg_trigger);

#ifndef CRAIG_INTERPOLATION
            minisat_->addClause(clause);
#endif
#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
            // we are not intersted in the orign of clauses bleonging to the
            // initial state. We want all the clauses that originates from the
            // transition relation and belonging to the unsat core.
            orign clo;
#endif
            minisat_->addClause(clause,
                                A_CLAUSE,
                                INITIAL_CLAUSE,
                                true,
                                false,
                                ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                ,clo
#endif
            );
#endif
        }
#ifdef DEBUG_MINISAT_BMC
        std::cerr << "initial trigger is: " << toInt(~(Lit(trigger))) << "\n";
#endif
        return trigger;
    }


    void MinisatBmc::addTrans(int before, int depth)
    {

        // assert that gen_ is not NULL
        assert(gen_ != NULL);

#ifdef CLAUSE_ORIGN

        //     if( use_modified_trans_ ){

        //         std::cout << "Modified trans clauses look like:" << std::endl;
        //         orign::iterator iter =  last_empty_clause_orign_.begin();

        //         while( iter != last_empty_clause_orign_.end() ){

        //             std::cout << "( " << std::flush;
        //             for ( unsigned s = 0; s < (**iter).size() ; ++s){
        //                 (**iter)[s].printBmcLit();
        //                 std::cout << " " << std::flush;
        //             }
        //              std::cout << ")" << std::endl;
        //             ++iter;
        //         }
        //     }

        if( use_modified_trans_ ){
            std::vector<std::vector< BmcLit>* > tmp;
            std::cerr << "Computing modified transition!!!\n";
            orign::iterator iter =  last_empty_clause_orign_.begin();
            while( iter != last_empty_clause_orign_.end() ){
                //std::vector< BmcLit> test = *(*iter);
                tmp.push_back(*iter);
                ++iter;
            }
            mod_trans_ = tmp;
        }
#endif

        const std::vector<std::vector< BmcLit>* >& cl =
#ifdef CLAUSE_ORIGN
            use_modified_trans_? mod_trans_ :
#endif
            gen_->getBmcProblem()->getTrans();



        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

//        int instance;

        while( before >= 0 ){

            int instance = depth;

#ifdef DEBUG_BMC_PROBLEM
            std::cerr << "adding T(" << depth
                      << "," << (depth + 1) << ") ";
            if( depth <= 0 ){
                std::cerr << "[ A-CLAUSE ]\n";
            }
            else if(depth == 1){
                std::cerr << "[ not pure B-CLAUSE]\n";
            }
            else{
                std::cerr << "[ pure B-CLAUSE]\n";
            }
#endif

            vec<Lit> clause;

            for( unsigned int t = 0 ; t < cl.size() ; ++t){
                clause.clear();
                for(unsigned int s = 0; s < cl[t]->size() ; ++s ){
                    assert(  cl[t]->at(s).getBmcInstance() == 0 ||  cl[t]->at(s).getBmcInstance() == 1 );

                    if(cl[t]->at(s).getVarId() & 1 ){ // negative literal
                        clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                         cl[t]->at(s).getBmcInstance() +
                                                                         instance )>>1],true));
                    }
                    else{
                        clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                         cl[t]->at(s).getBmcInstance() +
                                                                         instance )>>1],false));
                    }
                }
#ifndef CRAIG_INTERPOLATION
                minisat_->addClause(clause);
#endif
#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
                orign clo;
                clo.insert(const_cast<std::vector<BmcLit>*>(cl[t]));
#endif
                if( depth <= 0){
#ifdef CLAUSE_ORIGN
                    // clo.insert(const_cast<std::vector<BmcLit>*>(cl[t]));
#endif
                    minisat_->addClause(clause,
                                        A_CLAUSE,
                                        TRANS_CLAUSE,
                                        false,
                                        false,
                                        ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                        ,clo
#endif
                    );
                }
                else if(depth == 1 ){
                    minisat_->addClause(clause,
                                        B_CLAUSE,
                                        TRANS_CLAUSE,
                                        false,
                                        false,
                                        ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                        ,clo
#endif
                    );
                }
                else{
                    minisat_->addClause(clause,
                                        B_CLAUSE,
                                        TRANS_CLAUSE,
                                        true,
                                        false,
                                        ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                        ,clo
#endif
                    );
                }
#endif
            }
            --depth;
            --before;
        }
    }


    Var MinisatBmc::addTarget(int depth)
    {
#ifdef DEBUG_BMC_PROBLEM
        std::cerr << "adding Target(" << depth << ")\n";
#endif

        // assert that gen_ is not NULL
        assert(gen_ != NULL);

        const std::vector<std::vector< BmcLit>* >& cl = gen_->getBmcProblem()->getTarget();

        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

        vec<Lit> clause;

        Var trigger = minisat_->newVar();
        Lit neg_trigger = Lit(trigger,true);

        for( unsigned int t = 0 ; t < cl.size() ; ++t){
            clause.clear();
            for(unsigned int s = 0; s < cl[t]->size() ; ++s ){

                if(cl[t]->at(s).getVarId() & 1 ){ // negative literal
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth )>>1],true));
                }
                else{
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth )>>1],false));
                }
            }
            clause.push(neg_trigger);
#ifndef CRAIG_INTERPOLATION
            minisat_->addClause(clause);
#endif
#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
            orign clo;
#endif
            if( depth <= 0 ){
                minisat_->addClause(clause,
                                    A_CLAUSE,
                                    TARGET_CLAUSE,
                                    true,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
            else if( depth == 1 ){
                minisat_->addClause(clause,
                                    B_CLAUSE,
                                    TARGET_CLAUSE,
                                    false,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
            else{
                minisat_->addClause(clause,
                                    B_CLAUSE,
                                    TARGET_CLAUSE,
                                    true,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
#endif
        }
#ifdef DEBUG_MINISAT_BMC
        std::cerr << "target trigger is: " << toInt(~(Lit(trigger))) << "\n";
#endif
        return trigger;
    }

    Var MinisatBmc::addExpTarget(int depth)
    {

#ifdef DEBUG_BMC_PROBLEM
        std::cerr << "adding expensive Target(" << depth << ")\n";
#endif

        // assert that gen_ is not NULL
        assert(gen_ != NULL);

        const std::vector<std::vector< BmcLit>* >& cl = gen_->getBmcProblem()->getExpensiveBackwardTarget();

        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

        vec<Lit> clause;

        Var trigger = minisat_->newVar();
        Lit neg_trigger = Lit(trigger,true);

        for( unsigned int t = 0 ; t < cl.size() ; ++t){
            clause.clear();
            for(unsigned int s = 0; s < cl[t]->size() ; ++s ){

                if(cl[t]->at(s).getVarId() & 1 ){ // negative literal
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth )>>1],true));
                }
                else{
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(),
                                                                     cl[t]->at(s).getBmcInstance() + depth )>>1],false));
                }
            }
            clause.push(neg_trigger);
#ifndef CRAIG_INTERPOLATION
            minisat_->addClause(clause);
#endif
#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
            orign clo;
#endif
            if( depth <= 0 ){
                minisat_->addClause(clause,
                                    A_CLAUSE,
                                    TARGET_CLAUSE,
                                    true,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
            else if( depth == 1 ){
                minisat_->addClause(clause,
                                    B_CLAUSE,
                                    TARGET_CLAUSE,
                                    false,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
            else{
                minisat_->addClause(clause,
                                    B_CLAUSE,
                                    TARGET_CLAUSE,
                                    true,
                                    false,
                                    ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                    ,clo
#endif
                );
            }
#endif
        }
#ifdef DEBUG_MINISAT_BMC
        std::cerr << "expensive target trigger is: " << toInt(~(Lit(trigger))) << "\n";
#endif
        return trigger;



    }


    bool MinisatBmc::solveCurrentProblem(void)
    {
        // activate minisat with simplification and turn it off
        // std::cout << "Solving current subproblem!!\n";
#ifdef PARALLEL_MODE
        if(*status_){
            return true;
        }
#endif
        return minisat_->solve(assumptions_, true , true);
    }


    void MinisatBmc::deactivateLastTarget(void)
    {
#ifdef DEBUG_BMC_PROBLEM
        std::cerr << "Deactivating last target!!!\n";
#endif
        assumptions_[index_to_last_target_] = Lit(last_target_trigger_, true);
    }

    void MinisatBmc::deactivateLastCraig(void)
    {
#ifdef DEBUG_BMC_PROBLEM
        if(index_to_last_craig_ == 0 ){
            std::cerr << "Deactivating initial!!!\n";
        }
        else{
            std::cerr << "Deactivating last craig!!!\n";
        }
#endif
        assumptions_[index_to_last_craig_] = Lit(last_craig_trigger_, true);
    }

    void MinisatBmc::addNextVariables()
    {
        std::vector<int>& next_vars = gen_->getNewVariables();

#ifdef DECIDE_ONLY_INPUTS_AND_LATCHES
        std::vector<BMC_vartype>& next_var_types_ = gen_->getNewVariableTypes();
#endif

        for( unsigned int z = 0; z < next_vars.size() ; ++z ){

            // first create if needed a new entry for the
            // current concidered variable
            if( (int)var_map_.size() <= next_vars[z] )
                var_map_.resize( next_vars[z] + 1);

            assert( var_map_[next_vars[z]] == 0 );

#ifndef DECIDE_ONLY_INPUTS_AND_LATCHES
            // insert variable
            var_map_[ next_vars[z] ] = minisat_->newVar();
#endif
            // the follwing is comented out as at the momment it is
            // not clear whether all hwmcc benchmarks can be decided if
            // you only dercide INPUT_VARIABLES
#ifdef DECIDE_ONLY_INPUTS_AND_LATCHES
            if( next_var_types_[z] == INPUT_VAR || next_var_types_[z] == LATCH_VAR){
                var_map_[ next_vars[z] ] = minisat_->newVar();
            }
            else{
            // not a decision variable
                var_map_[ next_vars[z] ] = minisat_->newVar(true,false);
            }
#endif

            // now fill the reverse var map rev_var_map_

            if( (int)rev_var_map_.size() <=  var_map_[next_vars[z]] ){
                rev_var_map_.resize( var_map_[ next_vars[z] ] + 1 , -1);
            }

            rev_var_map_[ var_map_[next_vars[z]] ] = next_vars[z];
        }
    }
#ifdef CRAIG_INTERPOLATION
    //#define START_ONLY_FROM_NEW_STATES
    Var MinisatBmc::addCraig()
    {
#ifdef DEBUG_BMC_PROBLEM
        std::cerr << "adding Craig(" << 0 << ")\n";
#endif
        unsigned max_craig = minisat_->getMaxCraigVarId();

        unsigned next_free = max_craig + 1;
        unsigned old_next_free = next_free;

        std::vector< std::vector<unsigned> > craig_cnf = minisat_->getCnfCraig(next_free);

        // The folowing line simplifys a given craig iterpolant
        // it seems that only about ~3 clauses can be eleminated
        // that is why we do not do so
        //
        // gen_->simplifyCraig(craig_cnf,max_craig,next_free);

        // std::cerr << "There are " << craig_cnf.size() << " many craig clauses!\n";

#ifdef DEBUG_MINISAT_BMC
        std::cerr << "Craig clauses:\n{\n";
        for(unsigned t = 0; t < craig_cnf.size() ; ++t){
            std::cerr << "  ( ";
            for(unsigned h = 0 ; h < craig_cnf[t].size() ; ++h){
                std::cerr << craig_cnf[t][h] << " ";
            }
            std::cerr << ")\n";
        }
        std::cerr << "}\n";
#endif

        if( craig_cnf.size() == 0 || craig_cnf[0].size() == 0 ){
            std::cerr << "Fixed point reached --> system is safe\n";
            return 0;
        }

        Var map[next_free];

        for( unsigned z = 0; z < next_free ; ++z){
            map[z] = 0;
        }
#ifndef START_ONLY_FROM_NEW_STATES
        Var trigger = minisat_->newVar();
        Lit neg_trigger = Lit(trigger,true);
#endif
#ifdef START_ONLY_FROM_NEW_STATES
        Var trigger;
#endif
        // minisat clause that is filled and inserted below
        vec<Lit> clause;

        for(unsigned t = 0; t < craig_cnf.size() ; ++t){
            clause.clear();
            for(unsigned int s = 0; s < craig_cnf[t].size() ; ++s){
                if( ( craig_cnf[t][s] >> 1 ) >= old_next_free ){
                    // std::cerr << "Inserting craig tseitin!!\n";
                    // variable is tseitin variable
                    if( map[ (craig_cnf[t][s] >> 1) ] == 0 ){
                        map[ (craig_cnf[t][s] >> 1) ] = minisat_->newVar();
                    }

                    if(craig_cnf[t][s] & 1 ){ // negative literal
                        clause.push(Lit( (map[ craig_cnf[t][s] >> 1]) , true ));
                    }
                    else{
                        clause.push(Lit( (map[ craig_cnf[t][s] >> 1]) , false ));
                    }
                }
                else{
                    // global variable
                    if(craig_cnf[t][s] & 1 ){ // negative literal
                        clause.push(Lit(craig_cnf[t][s] >> 1,true));
                    }
                    else{
                        clause.push(Lit(craig_cnf[t][s] >> 1,false));
                    }
#ifdef COLLECT_LATCH_VARIABLES
                    latch_vars_in_craig_[latch_set_counter_].insert( (gen_->getBmcVarMap()->getBmcPair(minisat_->getRevVarMap()->at(craig_cnf[t][s]>>1))).first << 1);
#endif
                }

            }
#ifndef START_ONLY_FROM_NEW_STATES
            clause.push(neg_trigger);
#endif /* START_ONLY_FROM_NEW_STATES */
#ifdef START_ONLY_FROM_NEW_STATES
            if( t + 1 == craig_cnf.size() ){
                //assert last clause is unit clause e.g. the trigger
                //for this craig interpolant
                assert(clause.size() == 1);
                trigger = var(clause[0]);
                craig_added_ = true;
                return trigger;
            }

#endif /* START_ONLY_FROM_NEW_STATES */


#ifndef CRAIG_INTERPOLATION
            minisat_->addClause(clause);
#endif

#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
            orign clo;
#endif /* CLAUSE_ORIGN */
            minisat_->addClause(clause,
                                A_CLAUSE,
                                CRAIG_CLAUSE,
                                true,
                                false,
                                ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                ,clo
#endif /* CLAUSE_ORIGN */
            );
#endif /* CRAIG_INTERPOLATION */
        }

        //     std::cout << "Following latch variables are in the craig interpolant:\n";
        //     std::set<unsigned>::iterator iter = latch_vars_in_craig.begin();
        //     while(iter != latch_vars_in_craig.end() ) {
        //         std::cout << "Var: " << *iter << "\n";
        //         ++iter;
        //     }
        //     std::cout << "Number of different latch variables is: " << latch_vars_in_craig.size() << "\n";

        //     std::cout << "Latch variables contained in Target are:\n";
        //     for( unsigned i = 0 ; i < gen_->getGlobalVars().size() ; ++i){
        //         std::cout << "Var " << var_map_[gen_->getGlobalVars()[i]] << "\n";
        //     }

        craig_added_ = true;
#ifdef DEBUG_MINISAT_BMC
        std::cerr << "craig trigger is: " << toInt(~(Lit(trigger))) << "\n";
#endif
        return trigger;
    }
#endif

    void MinisatBmc::printVarMap()
    {
        std::cerr << "VarMap ( mapping from varibles to minisat variables)\n";
        for(unsigned int z = 0 ; z < var_map_.size() ; ++z){
            std::cerr << "( " << z << " ) --> ( " << var_map_[z] << " )\n";
        }
    }


    void MinisatBmc::printRevVarMap()
    {
        std::cerr << "RevVarMap ( mapping from minisat variables to variables)\n";
        for(unsigned int z = 0 ; z < rev_var_map_.size() ; ++z){
            std::cerr << "( " << z << " ) --> ( " << rev_var_map_[z] << " )\n";
        }
    }

    // we do a solver reset.
    // e.g.
#ifdef CRAIG_INTERPOLATION
    void MinisatBmc::resetMiniSatBmc(int depth)
    {
        //std::cerr << "Reset MinisatBMC!!!\n";
        // we do a solver reset
        gen_->resetBmcInstGen();

        current_bmc_depth_ = depth;

        assumptions_.clear();

        last_target_trigger_ = 0;

        last_craig_trigger_ = 0;

        init_trigger_ = 0;
#ifdef CLAUSE_ORIGN
        //last_empty_clause_orign_ = minisat_->getEmptyClauseOrign();

        //use_modified_trans_ = !use_modified_trans_;
        //std::cout << "use_modified_trans_ : " <<  use_modified_trans_ << "\n";
        // use_modified_trans_ = true;
#endif

        delete minisat_;

        minisat_ = new SimpSolver();

#ifdef PARALLEL_MODE
        minisat_->setStatus(status_);
#endif

        minisat_->setApproxNumber( approx_number_ );

        var_map_.clear();

        rev_var_map_.clear();

        craig_added_ = false;

        index_to_last_target_ = 0;

        index_to_last_craig_ = 0;

        craig_counter_ = 0;

        depth_of_last_trans_ = 0;

        minisat_->setRevVarMap( &(rev_var_map_) );
        minisat_->setVarMap( &(var_map_) );
        minisat_->setBmcVarMap(gen_->getBmcVarMap() );

        gen_->init( (-1)*approx_number_ );

        // now we have to add all initial state clauses
        addNextVariables();

        last_craig_trigger_ = addInitial(approx_number_);
        init_trigger_ = last_craig_trigger_;

        // important to add first the initial trigger
        assumptions_.push(Lit( init_trigger_, false));
        index_to_last_craig_ = assumptions_.size() - 1;

    }
#endif


#ifdef CLAUSE_ORIGN
    std::pair<bool,int> MinisatBmc::solve_bmc_unsat(int depth)
    {
        assert( gen_ != NULL );

        // declare needed variables
        addNextVariables();

        // add the initial state
        init_trigger_ = addInitial();

        // add the target state
        last_target_trigger_ = addTarget(0);

        // important to add first the initial trigger
        assumptions_.push(Lit( init_trigger_, false));

        assumptions_.push(Lit( init_trigger_, false));

        // activating target through assumption
        // and store the index to the last target trigger
        // in assuimption vector
        assumptions_.push(Lit( last_target_trigger_, false));
        index_to_last_target_ = assumptions_.size()-1;

        current_bmc_depth_ = 0;

        // solve for depth 0
        std::cout << "Solving depth 0       ";
        bool result = solveCurrentProblem();

        if( result ){
            std::cout << "sat\n";
            return std::make_pair(result,current_bmc_depth_);
        }


        std::cout << "unsat\n";

        deactivateLastTarget();

        for(int d = 0 ; d < depth ; ++d){

            // declare nest variables
            gen_->init( d + 1 );

            addNextVariables();

            last_target_trigger_ = addTarget(d + 1);

            addTrans( d - depth_of_last_trans_ , d );

            depth_of_last_trans_ = d + 1;

            // activating target through assumption
            assumptions_.push(Lit( last_target_trigger_, false));

            index_to_last_target_ = assumptions_.size()-1;

            std::cout << "Solving depth " << d + 1  << std::flush;

            // solve that problem
            result = solveCurrentProblem();

            current_bmc_depth_ = d + 1;

            if( result){
                if( use_modified_trans_ ){
                    std::cout << "\n sat by using mod. trans\n";
                    resetMiniSatBmc(d--);
                    use_modified_trans_ = false;
                }
                else{
                    std::cout << "       sat\n";
                    return std::make_pair(result,current_bmc_depth_);
                }
            }
            else{


                std::cout << "       unsat\n";

                deactivateLastTarget();

                use_modified_trans_ = true;

                if( !last_empty_clause_orign_.empty() ){
                    orign& new_empty = minisat_->getEmptyClauseOrign();
                    orign::iterator iter = new_empty.begin();

                    while( iter !=  new_empty.end() ){
                        last_empty_clause_orign_.insert(*iter);
                        iter++;
                    }



                }
                else{
                    last_empty_clause_orign_ = minisat_->getEmptyClauseOrign();
                }
#ifdef CLAUSE_ORIGN
                std::cerr << "\n   Number of transition clauses needed for unsat core: "
                          << last_empty_clause_orign_.size() << "\n";
#endif
            }
        }
        return std::make_pair(false,depth);
    }
#endif


    Var MinisatBmc::addInvariant(int depth)
    {
        std::cerr << "adding Invariant(" << depth << ")\n";


        const std::vector<std::vector< BmcLit>* >& cl = gen_->getBmcProblem()->getOrigTarget();

        std::cerr << "There are " << cl.size() << " clauses for the assumed invariant!!\n";

        assert(cl[0]->size() == 1 );

        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

        int instance = 0 ;

        Var trigger = minisat_->newVar();
        Lit neg_trigger = Lit(trigger,true);

        vec<Lit> clause;

        if( cl[0]->at(0).getVarId() & 1 ){
            clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[0]->at(0).getVarId(), instance + depth)>>1], false));
        }
        else{
            clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[0]->at(0).getVarId(), instance + depth)>>1], true ));
        }
        clause.push(neg_trigger);

#ifndef CRAIG_INTERPOLATION
        minisat_->addClause(clause);
#endif

#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
        // we are not intersted in the orign of clauses bleonging to the
        // initial state. We want all the clauses that originates from the
        // transition relation and belonging to the unsat core.
        orign clo;
#endif
        minisat_->addClause(clause,
                            A_CLAUSE,
                            INITIAL_CLAUSE,
                            true,
                            false,
                            ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                            ,clo
#endif
        );
#endif

        for( unsigned int t = 1 ; t < cl.size() ; ++t){
            clause.clear();
            for(unsigned int s = 0; s < cl[t]->size() ; ++s ){
                if(cl[t]->at(s).getVarId() & 1 ){ // negative literal
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(), instance + depth)>>1],true));
                }
                else{
                    clause.push(Lit( var_map_[bmc_var_map->getSatLit(cl[t]->at(s).getVarId(), instance + depth )>>1],false));
                }
            }
            clause.push(neg_trigger);

#ifndef CRAIG_INTERPOLATION
            minisat_->addClause(clause);
#endif
#ifdef CRAIG_INTERPOLATION
#ifdef CLAUSE_ORIGN
            // we are not intersted in the orign of clauses bleonging to the
            // initial state. We want all the clauses that originates from the
            // transition relation and belonging to the unsat core.
            orign clo;
#endif
            minisat_->addClause(clause,
                                A_CLAUSE,
                                INITIAL_CLAUSE,
                                true,
                                false,
                                ((*minisat_).m_).getConst1()
#ifdef CLAUSE_ORIGN
                                ,clo
#endif
            );
#endif
        }
#ifdef DEBUG_MINISAT_BMC
        std::cerr << "initial trigger is: " << toInt(~(Lit(trigger))) << "\n";
#endif
        return trigger;
    }

#ifdef CRAIG_INTERPOLATION
    bool MinisatBmc::startFirstCheck()
    {
        // delete old solver and instantiate new solver
        if( minisat_ != NULL ){
            delete minisat_;
            minisat_ = new SimpSolver();
        }

        // setup bmc problem for depth 1
        assert( gen_ != NULL );
        gen_->init(1);

        var_map_.clear();
        rev_var_map_.clear();

        minisat_->setApproxNumber( 0 );
        minisat_->setRevVarMap( &(rev_var_map_) );
        minisat_->setVarMap( &(var_map_) );
        minisat_->setBmcVarMap(gen_->getBmcVarMap() );

        // add desired variables
        addNextVariables();

        // clear the assumpptions
        assumptions_.clear();

        // add pure A-Clauses (e.g. the initial state)
        last_craig_trigger_ = addInitial(0);
        init_trigger_ = last_craig_trigger_;
        assumptions_.push(Lit( init_trigger_, false));
        index_to_last_craig_ = assumptions_.size() - 1;

        // now add the B-Clauses in this case only one Target
        // e.g. target Tar_1
        last_target_trigger_ = addTarget(1);

        assumptions_.push(Lit( last_target_trigger_, false));
        index_to_last_target_ = assumptions_.size() - 1;

        // now add Trans relation T_(0,1)
        addTrans(0,0);

        // add the first craig interpolant
        minisat_->addFirstCraig( gen_->getBmcProblem() );

        std::cout << "Solving depth 1       ";
        if( solveCurrentProblem() ){
            std::cout << "sat\n";
            std::cout << "Result:                 sat\n";
            return true;
        }
        else{
            std::cout << "unsat\n";

            bool go_on = true;
            unsigned int depth_counter = 1;

            while( go_on ){

                ++depth_counter;

                std::cerr << "Starting fixed point check " << std::endl;
                if( minisat_->performFixedPointCheck() ){
                    std::cout << "Result:                 uns\n";
                    return true;
                }

                deactivateLastCraig();
                last_craig_trigger_ = addCraig();

                // activating craig through assumption
                // and store the index to the last craig trigger
                // in assuimption vector
                assumptions_.push(Lit( last_craig_trigger_, false));
                index_to_last_craig_ = assumptions_.size() - 1;

                std::cout << "Solving depth " << depth_counter << std::flush;

                go_on = !(solveCurrentProblem());

                if( !go_on ){

                    minisat_->clearInterpolants();
                    // folowing 3 lines are used if you do no reset
                    deactivateLastTarget();
                    deactivateLastCraig();

                    craig_added_ = false;

                    // trigger initial clauses
                    last_craig_trigger_ = init_trigger_;
                    assumptions_[0] = (Lit( init_trigger_, false));

                    // uppdate index_to_last_craig_
                    index_to_last_craig_ = 0;

                    std::cout << "       sat (over-approx.)\n";
                    return false;
                }
                else{
                    std::cout << "       unsat (over-approx.)\n";

                }

            }
        }
    }
#endif

    void MinisatBmc::printTrace(unsigned depth)
    {
        std::cout << "There exists a trace by depth " << depth << "\n";
        minisat_->extendModel();
        // minisat_->verifyModel();

        const vec<lbool>& model = minisat_->getModel();

        std::vector<BMC_vartype>& variable_type_map = gen_->getBmcProblem()->getVarTypeMap();
        BmcVarMap* bmc_var_map = gen_->getBmcVarMap();

        //bmc_var_map->printBmcVarMap();
        // first the inputs
        int start_depth = -(gen_->getDepthOffset());
        //std::cout << "start_depth: " << start_depth << "\n";

        for(unsigned int i = 0; i <  variable_type_map.size() ; ++i){
            if( variable_type_map[i] == INPUT_VAR ){
                std::cout << "    Input " << i << ":\n";
                for(unsigned int d = 0; d <= depth ; d++){
                    std::cout << "        " << "@" << d << ":  ";
                    unsigned minisat_var = var_map_[ (bmc_var_map->getSatLit( i << 1,start_depth + d)) >> 1 ];
                    if( (int)minisat_var >= model.size() ){
                        std::cout << "1";
                    }
                    else{
                        assert( l_Undef != model[minisat_var] );
                        if( model[minisat_var] == l_True ){
                            std::cout << "1";
                        }
                        else{
                            std::cout << "0";
                        }
                    }
                    std::cout << "\n";
                }
            }
        }
        for(unsigned int i = 0; i <  variable_type_map.size() ; ++i){
            if( variable_type_map[i] == LATCH_VAR ){
                std::cout << "    Latch " << i << ":\n";
                for(unsigned int d = 0; d <= depth ; d++){
                    std::cout << "        " << "@" << d << ":  ";
                    unsigned minisat_var = var_map_[ (bmc_var_map->getSatLit( i << 1,start_depth + d)) >> 1 ];
                    if( (int)minisat_var >= model.size() ){
                        std::cout << "1";
                    }
                    else{
                        assert( l_Undef != model[minisat_var] );
                        if( model[minisat_var] == l_True ){
                            std::cout << "1";
                        }
                        else{
                            std::cout << "0";
                        }
                    }
                    std::cout << "\n";
                }
            }
        }

    }

}

