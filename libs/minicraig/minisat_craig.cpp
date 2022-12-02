#include "minisat_craig.hpp"

#include <sstream>

namespace minicraig {

    MinisatCraig::MinisatCraig()
    : pos_a_trigger_(0),
      neg_a_trigger_(0),
      pos_b_trigger_(0),
      neg_b_trigger_(0),
      minisat_(NULL),
      no_more_b_clauses_(false),
      dumpToFile_( false ),
      dumpFile_( 0 ),
      dumpMaxIndex_(0),
      dumpClauses_(0)
      {

        minisat_ = new SimpSolver();
        Var trigger_var_a = minisat_->newVar();
        Var trigger_var_b = minisat_->newVar();

        pos_a_trigger_ = Lit(trigger_var_a,false);
        neg_a_trigger_ = Lit(trigger_var_a,true);

        pos_b_trigger_ = Lit(trigger_var_b,false);
        neg_b_trigger_ = Lit(trigger_var_b,true);

    }

    MinisatCraig::~MinisatCraig(){
        delete minisat_;
        
        if( dumpToFile_ )
        {
            stopDump();
        }
    }

    Var MinisatCraig::newVar(){
        return minisat_->newVar();
    }

    Lit MinisatCraig::genLit(Var var, bool neg){
        return Lit(var,neg);
    }

    void 
    MinisatCraig::
    addAClause( const std::vector<Lit>& clause )
    {
        no_more_b_clauses_ = true;
        
        /* 
         * create a minisat clause with the right size
         */
        vec<Lit> minisat_clause( clause.size() + 1 );

        minisat_clause[0] = pos_a_trigger_;
        for( unsigned int i = 0; i != clause.size(); ++i)
        {
            minisat_clause[i+1] = clause[i];
        }

        minisat_->addClause(
            minisat_clause,
            A_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( true, clause );
        }
    }
    
    void 
    MinisatCraig::
    addAClause( Lit lit1 )
    {
        no_more_b_clauses_ = true;
        
        vec<Lit> minisat_clause( 2 );
        minisat_clause[0] = pos_a_trigger_;
        minisat_clause[1] = lit1;
        
        minisat_->addClause(
            minisat_clause,
            A_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( true, lit1 );
        }
    }
    
    void 
    MinisatCraig::
    addAClause( Lit lit1, Lit lit2 )
    {
        no_more_b_clauses_ = true;
        
        vec<Lit> minisat_clause( 3 );
        minisat_clause[0] = pos_a_trigger_;
        minisat_clause[1] = lit1;
        minisat_clause[2] = lit2;
        
        minisat_->addClause(
            minisat_clause,
            A_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( true, lit1, lit2 );
        }
    }
    
    void 
    MinisatCraig::
    addAClause( Lit lit1, Lit lit2, Lit lit3 )
    {
        no_more_b_clauses_ = true;
        
        vec<Lit> minisat_clause( 4 );
        minisat_clause[0] = pos_a_trigger_;
        minisat_clause[1] = lit1;
        minisat_clause[2] = lit2;
        minisat_clause[3] = lit3;
        
        minisat_->addClause(
            minisat_clause,
            A_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( true, lit1, lit2, lit3 );
        }
    }
    
    
    void 
    MinisatCraig::
    addBClause( const std::vector<Lit>& clause )
    {
        assert( !no_more_b_clauses_ );
        
        /* 
         * create a minisat clause with the right size
         */
        vec<Lit> minisat_clause( clause.size() + 1 );
        
        minisat_clause[0] = pos_b_trigger_;
        for( unsigned int i = 0; i != clause.size(); ++i )
        {
            minisat_clause[i+1] = clause[i];
        }

        minisat_->addClause(
            minisat_clause,
            B_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( false, clause );
        }
    }
    
    void 
    MinisatCraig::
    addBClause( Lit lit1 )
    {
        assert( !no_more_b_clauses_ );
        
        vec<Lit> minisat_clause( 2 );
        minisat_clause[0] = pos_b_trigger_;
        minisat_clause[1] = lit1;
        
        minisat_->addClause(
            minisat_clause,
            B_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( false, lit1 );
        }
    }
    
    void 
    MinisatCraig::
    addBClause( Lit lit1, Lit lit2 )
    {
        assert( !no_more_b_clauses_ );
        
        vec<Lit> minisat_clause( 3 );
        minisat_clause[0] = pos_b_trigger_;
        minisat_clause[1] = lit1;
        minisat_clause[2] = lit2;
        
        minisat_->addClause(
            minisat_clause,
            B_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( false, lit1, lit2 );
        }
    }
    
    void 
    MinisatCraig::
    addBClause( Lit lit1, Lit lit2, Lit lit3 )
    {
        assert( !no_more_b_clauses_ );
        
        vec<Lit> minisat_clause( 4 );
        minisat_clause[0] = pos_b_trigger_;
        minisat_clause[1] = lit1;
        minisat_clause[2] = lit2;
        minisat_clause[3] = lit3;
        
        minisat_->addClause(
            minisat_clause,
            B_CLAUSE,
            INITIAL_CLAUSE,
            false,
            false,
            ((*minisat_).m_).getTrue() );
        
        if( dumpToFile_ )
        {
            dumpClause( false, lit1, lit2, lit3 );
        }
    }
    
    
    void MinisatCraig::startDump( int fileIndex )
    {
        assert( !dumpToFile_ );
        dumpToFile_ = true;
        std::ostringstream oss; 
        oss << "interpolation_" << fileIndex << ".cnf";
        dumpFile_ = new std::ofstream( oss.str().c_str() );
        *dumpFile_ << "p cnf MAXINDEX CLAUSES\n";
        
        dumpMaxIndex_ = 0;
        dumpClauses_ = 0;
    }
    
    void MinisatCraig::stopDump()
    {
        assert( dumpToFile_ );
        
        *dumpFile_ << "c STATISTICS " << dumpMaxIndex_ << " " << dumpClauses_ << "\n" << std::flush;
        
        delete dumpFile_;
        dumpFile_ = 0;
    }

    void MinisatCraig::dumpClause( bool a, const std::vector<Lit>& clause )
    {
        assert( dumpToFile_ );
        *dumpFile_ << "c" << ( a ? "A" : "B" ) << "\n";
        for( std::vector<Lit>::const_iterator p = clause.begin(); p != clause.end(); ++p )
        {
            if( var(*p) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(*p) + 1;
            *dumpFile_ << ( sign(*p) ? "-" : "" ) << 1+var(*p) << " ";
        }
        *dumpFile_ << "0\n";
        
        ++dumpClauses_;
    }
    
    void MinisatCraig::dumpClause( bool a, Lit lit1 )
    {
        assert( dumpToFile_ );
        *dumpFile_ << "c" << ( a ? "A" : "B" ) << "\n";
        *dumpFile_ << ( sign(lit1) ? "-" : "" ) << 1+var(lit1) << " ";
        *dumpFile_ << "0\n";
        
        if( var(lit1) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit1) + 1;
        ++dumpClauses_;
    }
    
    void MinisatCraig::dumpClause( bool a, Lit lit1, Lit lit2 )
    {
        assert( dumpToFile_ );
        *dumpFile_ << "c" << ( a ? "A" : "B" ) << "\n";
        *dumpFile_ << ( sign(lit1) ? "-" : "" ) << 1+var(lit1) << " ";
        *dumpFile_ << ( sign(lit2) ? "-" : "" ) << 1+var(lit2) << " ";
        *dumpFile_ << "0\n";
        
        if( var(lit1) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit1) + 1;
        if( var(lit2) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit2) + 1;
        ++dumpClauses_;
    }
    
    void MinisatCraig::dumpClause( bool a, Lit lit1, Lit lit2, Lit lit3 )
    {
        assert( dumpToFile_ );
        *dumpFile_ << "c" << ( a ? "A" : "B" ) << "\n";
        *dumpFile_ << ( sign(lit1) ? "-" : "" ) << 1+var(lit1) << " ";
        *dumpFile_ << ( sign(lit2) ? "-" : "" ) << 1+var(lit2) << " ";
        *dumpFile_ << ( sign(lit3) ? "-" : "" ) << 1+var(lit3) << " ";
        *dumpFile_ << "0\n";
        
        if( var(lit1) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit1) + 1;
        if( var(lit2) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit2) + 1;
        if( var(lit3) + 1 > dumpMaxIndex_ ) dumpMaxIndex_ = var(lit3) + 1;
        ++dumpClauses_;
    }
    
        
    void MinisatCraig::bumpActivitiesA( int count )
    {
        vec<Var> vs;
        for( int i = 0; i != (int)minisat_->var_types_.size(); ++i )
        {
            if( minisat_->var_types_[i] == A_LOCAL )
            {
                vs.push( Var( i ) );
            }
        }
        
        for( int i = 0; i < count; ++i )
        {
            minisat_->varBumpActivities( vs );
        }
    }
    
    void MinisatCraig::bumpActivitiesB( int count )
    {
        vec<Var> vs;
        for( int i = 0; i != (int)minisat_->var_types_.size(); ++i )
        {
            if( minisat_->var_types_[i] == B_LOCAL )
            {
                vs.push( Var( i ) );
            }
        }
        
        for( int i = 0; i < count; ++i )
        {
            minisat_->varBumpActivities( vs );
        }
    }
    
    bool MinisatCraig::solve(){
        // activate minisat with simplification and turn it off
        vec<Lit> assumptions;
        assumptions.push(neg_a_trigger_);
        assumptions.push(neg_b_trigger_);
        //return minisat_->solve(assumptions, true , true);
        return minisat_->solve(assumptions, false, false);
        //return minisat_->solve(assumptions);
    }
#if 0
    std::vector< std::vector<unsigned> > MinisatCraig::getCnfCraig(unsigned& next_free){
        return minisat_->getCnfCraig(next_free);
    }
#endif
}
