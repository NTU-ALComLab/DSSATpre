/**********************************************************************
 * minisat_craig.hpp
 *
 * Stefan Kupferschmid <skupfers@informatik.uni-freiburg.de>
 *
 **********************************************************************/
#ifndef MINICRAIG_MINISAT_CRAIG_HPP
#define MINICRAIG_MINISAT_CRAIG_HPP

#include "craig_define.hpp"
#include "SimpSolver.h"
#include "Solver.h"
#include "Vec.h"
#include <cassert>
#include <fstream>
#include <vector>
#include "SolverTypes.h"


namespace minicraig {

    class MinisatCraig{
    public:
        MinisatCraig();
        MinisatCraig(const MinisatCraig& other) = delete;
        ~MinisatCraig();
        MinisatCraig& operator=(const MinisatCraig& other) = delete;

        // functions for generating variables, literals and clauses
        Var newVar();
        static Lit genLit(Var var, bool neg);

        void addAClause( const std::vector<Lit>& clause );
        void addAClause( Lit lit1 );
        void addAClause( Lit lit1, Lit lit2 );
        void addAClause( Lit lit1, Lit lit2, Lit lit3 );

        void addBClause( const std::vector<Lit>& clause );
        void addBClause( Lit lit1 );
        void addBClause( Lit lit1, Lit lit2 );
        void addBClause( Lit lit1, Lit lit2, Lit lit3 );

        // dumping to file
        void startDump( int fileIndex = 0 );
        void stopDump();
        void dumpClause( bool a, const std::vector<Lit>& clause );
        void dumpClause( bool a, Lit lit1 );
        void dumpClause( bool a, Lit lit1, Lit lit2 );
        void dumpClause( bool a, Lit lit1, Lit lit2, Lit lit3 );

        void bumpActivitiesA( int count = 1 );
        void bumpActivitiesB( int count = 1 );

        // solve method
        bool solve();

        void setResourceLimits( int nodeLimit, double timeLimit, bool* failed )
        {
            minisat_->setResourceLimits( nodeLimit, timeLimit, failed );
        }

        const SimpleAIG& getInterpolantAIG() const
        {
            return minisat_->getInterpolantAIG();
        }

        const SimpleAIGEdge& getInterpolant() const
        {
            return minisat_->getInterpolant();
        }
        
#ifdef PARTIAL_INTERPOLANT_TEST    
        void setVarMap(std::vector<Var>* var_map){
            minisat_->set_var_map(var_map);
        }
        void setRevVarMap(std::vector<int>* rev_var_map){
            minisat_->set_rev_var_map(rev_var_map);
        }

        void setAclauses(const std::string& a_clauses,unsigned int a_num){
            minisat_->set_a_clauses(a_clauses,a_num);
        }
        void setBclauses(const std::string& b_clauses,unsigned int b_num){
            minisat_->set_b_clauses(b_clauses,b_num);
        }
        void setNextFree(unsigned n_free){
            minisat_->set_next_free(n_free);
        }
#endif

        // method for receiving craig interpolant in cnf
        std::vector< std::vector<unsigned> > getCnfCraig(unsigned& next_free);

        VarType getVarType( Lit x ) const
        {
            return minisat_->getVarType( x );
        }

        SimpSolver* internalSolver()
        {
            return minisat_;
        }

    private:
        Lit pos_a_trigger_;
        Lit neg_a_trigger_;

        Lit pos_b_trigger_;
        Lit neg_b_trigger_;

        // a pointer to minisat simp solver 
        SimpSolver* minisat_;

        //Solver* minisat_;
        bool no_more_b_clauses_;

        bool dumpToFile_;
        std::ofstream* dumpFile_;
        int dumpMaxIndex_;
        int dumpClauses_;
    };
}
#endif /* MINISAT_CRAIG_HPP */
